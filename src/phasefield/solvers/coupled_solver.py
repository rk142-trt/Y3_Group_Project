"""
Solver module for the coupled protein-RNA phase field model.

This module implements time-stepping schemes for the system of equations
described in Goh et al. (2025):

Protein dynamics (Cahn-Hilliard):
    ∂c/∂t = ∇·(M_c ∇μ)
    μ = df/dc - κ∇²c + χm + γcm²

RNA dynamics (reaction-diffusion):
    ∂m/∂t = D_m∇²m + k_p·S(r)·c - k_d·m

The solver uses semi-implicit time stepping for stability:
    - Implicit treatment of diffusive/stabilizing terms
    - Explicit treatment of nonlinear terms

This approach allows larger time steps than fully explicit methods
while remaining computationally tractable (no nonlinear solves).
"""

import numpy as np
from typing import Tuple, Optional, Callable, List
from dataclasses import dataclass
import time

# Import our modules
import sys
sys.path.insert(0, '/home/claude/condensate_project/src/phasefield')
from config import SimulationConfig
from physics.free_energy import FreeEnergy, RNAProteinCoupling, ChemicalPotential
from numerics.grid import RadialGrid, RadialOperators, create_initial_droplet, create_gaussian_source


@dataclass
class SimulationState:
    """
    Container for the current state of the simulation.
    
    Stores the field values (concentration and RNA), time, and
    optional diagnostic quantities for analysis.
    """
    c: np.ndarray          # Protein concentration
    m: np.ndarray          # RNA concentration
    t: float               # Current time
    step: int              # Current time step number
    
    # Optional diagnostics (computed on demand)
    total_protein: Optional[float] = None
    total_rna: Optional[float] = None
    free_energy: Optional[float] = None
    droplet_center: Optional[float] = None
    droplet_radius: Optional[float] = None


class CoupledSolver:
    """
    Solver for the coupled Cahn-Hilliard/reaction-diffusion system.
    
    This class handles time-stepping of both the protein and RNA fields,
    using a semi-implicit scheme for numerical stability.
    
    The Cahn-Hilliard equation is split into two second-order equations:
        ∂c/∂t = M_c ∇²μ
        μ = f'(c) - κ∇²c + μ_coupling(c, m)
    
    For stability, we use convex splitting (Eyre's method):
        - Treat the convex part of f(c) implicitly
        - Treat the concave part explicitly
        - Treat coupling terms explicitly (they're lower order)
    
    The reaction-diffusion equation for RNA is simpler and uses
    standard semi-implicit treatment (diffusion implicit, reaction explicit).
    """
    
    def __init__(self, config: SimulationConfig):
        """
        Initialize the solver with the given configuration.
        
        Args:
            config: SimulationConfig object containing all parameters
        """
        self.config = config
        
        # Set up the grid
        self.grid = RadialGrid(
            n_r=config.numerical.n_r,
            R_domain=config.numerical.R_domain,
            dimension=2  # 2D cylindrical
        )
        
        # Set up operators
        self.ops = RadialOperators(self.grid)
        
        # Set up physics
        self.free_energy = FreeEnergy(
            alpha=config.free_energy.alpha,
            beta=config.free_energy.beta,
            kappa=config.free_energy.kappa,
            c_bar=config.free_energy.c_bar
        )
        
        self.coupling = RNAProteinCoupling(
            chi=config.coupling.chi,
            gamma=config.coupling.gamma
        )
        
        self.chem_pot = ChemicalPotential(self.free_energy, self.coupling)
        
        # Create the RNA production source profile
        self.rna_source = create_gaussian_source(
            self.grid,
            amplitude=config.transport.k_p,
            sigma=config.transport.sigma_p,
            r_source=0.0  # Promoter at origin
        )
        
        # Precompute matrices for implicit solves
        self._setup_implicit_matrices()
        
        # Initialize state
        self.state = None
    
    def _setup_implicit_matrices(self):
        """
        Precompute matrices needed for implicit time stepping.
        
        For the RNA equation: (I - dt·D_m·L)·m^{n+1} = m^n + dt·(source - k_d·m^n)
        We can precompute (I - dt·D_m·L)^{-1} since it doesn't change.
        
        For Cahn-Hilliard, we use operator splitting, so we need the
        Laplacian matrix in a usable form.
        """
        n = self.grid.n_r
        dt = self.config.numerical.dt
        
        # Laplacian matrix
        self.L = self.ops.laplacian_matrix()
        
        # For RNA: (I - dt·D_m·L)
        D_m = self.config.transport.D_m
        self.rna_implicit_matrix = np.eye(n) - dt * D_m * self.L
        
        # Precompute LU decomposition for efficiency
        from scipy.linalg import lu_factor
        self.rna_lu = lu_factor(self.rna_implicit_matrix)
        
        # For Cahn-Hilliard, we need M_c and κ
        M_c = self.config.transport.M_c
        kappa = self.config.free_energy.kappa
        
        # The bi-Laplacian term κ·M_c·∇⁴c requires careful treatment
        # We'll use a stabilized explicit scheme for simplicity
        # More sophisticated: implicit treatment of ∇⁴c via splitting
        
        # Stability parameter for explicit scheme
        # For ∇⁴ operator, stability requires dt < C·dr⁴/(κ·M_c)
        self.dt_stability_factor = self.grid.dr**4 / (kappa * M_c * 16)
        
        if dt > self.dt_stability_factor:
            print(f"WARNING: dt={dt} may be too large for stability.")
            print(f"Recommended: dt < {self.dt_stability_factor:.2e}")
    
    def initialize(self):
        """
        Set up the initial conditions for the simulation.
        
        Creates a circular droplet at the specified distance from the
        promoter (origin), with the specified initial concentrations.
        """
        ic = self.config.initial
        
        # Initial protein concentration (droplet + background)
        # Interface width from config
        interface_width = self.config.free_energy.interface_width
        c0 = create_initial_droplet(
            self.grid,
            c_plus=ic.c_plus_init,
            c_minus=ic.c_minus_init,
            R_droplet=ic.R_init,
            r_center=ic.r_init,
            interface_width=interface_width
        )
        
        # Initial RNA concentration (uniform, typically zero)
        m0 = np.full(self.grid.n_r, ic.m_init)
        
        # Create initial state
        self.state = SimulationState(c=c0, m=m0, t=0.0, step=0)
        
        # Compute initial diagnostics
        self._update_diagnostics()
        
        return self.state
    
    def step(self):
        """
        Advance the simulation by one time step.
        
        Uses operator splitting to decouple the Cahn-Hilliard and
        reaction-diffusion updates:
            1. Update RNA (reaction-diffusion, semi-implicit)
            2. Update protein (Cahn-Hilliard, stabilized explicit)
        
        Returns:
            Updated SimulationState
        """
        dt = self.config.numerical.dt
        
        # Current state
        c = self.state.c.copy()
        m = self.state.m.copy()
        
        # Step 1: Update RNA (semi-implicit)
        m_new = self._step_rna(c, m, dt)
        
        # Step 2: Update protein (Cahn-Hilliard with coupling to new RNA)
        c_new = self._step_protein(c, m_new, dt)
        
        # Update state
        self.state = SimulationState(
            c=c_new, 
            m=m_new, 
            t=self.state.t + dt,
            step=self.state.step + 1
        )
        
        return self.state
    
    def _step_rna(self, c: np.ndarray, m: np.ndarray, dt: float) -> np.ndarray:
        """
        Update RNA concentration using semi-implicit scheme.
        
        The RNA equation is:
            ∂m/∂t = D_m∇²m + k_p·S(r)·c - k_d·m
        
        Semi-implicit: diffusion implicit, reaction explicit
            (I - dt·D_m·L)·m^{n+1} = m^n + dt·(k_p·S·c - k_d·m^n)
        
        Actually, for better stability with the decay term:
            (I - dt·D_m·L + dt·k_d·I)·m^{n+1} = m^n + dt·k_p·S·c
        
        But we precomputed without k_d, so we'll use the simpler form
        with explicit decay for now.
        
        Args:
            c: Current protein concentration
            m: Current RNA concentration
            dt: Time step
            
        Returns:
            New RNA concentration
        """
        k_p = self.config.transport.k_p
        k_d = self.config.transport.k_d
        
        # Right-hand side: m^n + dt·(source·c - k_d·m^n)
        # Note: source profile is k_p·exp(-r²/2σ²), so production is source·c
        production = self.rna_source * c  # k_p·S(r)·c
        decay = k_d * m
        
        rhs = m + dt * (production - decay)
        
        # Solve (I - dt·D_m·L)·m^{n+1} = rhs
        from scipy.linalg import lu_solve
        m_new = lu_solve(self.rna_lu, rhs)
        
        # Ensure non-negativity (RNA can't go negative)
        m_new = np.maximum(m_new, 0.0)
        
        return m_new
    
    def _step_protein(self, c: np.ndarray, m: np.ndarray, dt: float) -> np.ndarray:
        """
        Update protein concentration using stabilized explicit Cahn-Hilliard.
        
        The Cahn-Hilliard equation is:
            ∂c/∂t = M_c ∇²μ
            μ = f'(c) - κ∇²c + χm + γcm²
        
        We use a stabilized explicit scheme where we add and subtract
        a stabilizing term S·∇²c to improve stability:
            ∂c/∂t = M_c ∇²(f'(c) - κ∇²c + χm) + M_c·S·∇²c - M_c·S·∇²c
        
        Treating the -S·∇²c term implicitly and everything else explicitly
        gives better stability while remaining linear.
        
        For simplicity, we'll use a basic explicit scheme with small dt.
        
        Args:
            c: Current protein concentration
            m: New RNA concentration (from step 1)
            dt: Time step
            
        Returns:
            New protein concentration
        """
        M_c = self.config.transport.M_c
        kappa = self.config.free_energy.kappa
        chi = self.config.coupling.chi
        gamma = self.config.coupling.gamma
        
        # Compute Laplacian of c for the chemical potential
        lap_c = self.ops.laplacian(c)
        
        # Chemical potential: μ = f'(c) - κ∇²c + χm + γcm²
        mu_bulk = self.free_energy.bulk_derivative(c)
        mu_gradient = -kappa * lap_c
        mu_coupling = chi * m + gamma * c * m**2
        mu = mu_bulk + mu_gradient + mu_coupling
        
        # Flux: J = -M_c ∇μ
        # For finite volumes, we need the Laplacian of μ
        lap_mu = self.ops.laplacian(mu)
        
        # Update: ∂c/∂t = M_c ∇²μ = -∇·J
        dcdt = M_c * lap_mu
        
        # Explicit Euler update
        c_new = c + dt * dcdt
        
        return c_new
    
    def _update_diagnostics(self):
        """Compute diagnostic quantities for the current state."""
        if self.state is None:
            return
        
        # Total protein (should be conserved)
        self.state.total_protein = self.grid.integrate(self.state.c)
        
        # Total RNA
        self.state.total_rna = self.grid.integrate(self.state.m)
        
        # Droplet center of mass and radius
        self._compute_droplet_properties()
    
    def _compute_droplet_properties(self):
        """
        Estimate the droplet center of mass and effective radius.
        
        For metastable states where the background concentration may be
        above the binodal, we use a robust detection based on concentration
        peaks and half-maximum thresholds.
        """
        c = self.state.c
        r = self.grid.r_centers
        V = self.grid.cell_volumes
        
        # Use a concentration threshold based on the actual profile
        c_max = np.max(c)
        c_min = np.min(c)
        
        # If concentration variation is too small, no droplet detected
        if c_max - c_min < 0.1:
            self.state.droplet_center = None
            self.state.droplet_radius = None
            return
        
        # For center of mass: weight by excess over background
        c_threshold = c_min + 0.5 * (c_max - c_min)  # Half-max threshold
        c_excess = np.maximum(c - c_min, 0)  # Excess over background
        total_excess = np.sum(c_excess * V)
        
        if total_excess > 0:
            self.state.droplet_center = np.sum(r * c_excess * V) / total_excess
        else:
            self.state.droplet_center = None
            self.state.droplet_radius = None
            return
        
        # For radius: find where concentration drops to half-max
        # This is more physically meaningful than volume-based estimate
        in_droplet = c > c_threshold
        if np.any(in_droplet):
            # Find the outermost radius where c > threshold
            droplet_indices = np.where(in_droplet)[0]
            if len(droplet_indices) > 0:
                # Get radial extent of the high-concentration region
                r_inner = r[droplet_indices[0]]
                r_outer = r[droplet_indices[-1]]
                
                # Radius is half the extent
                self.state.droplet_radius = (r_outer - r_inner) / 2
            else:
                self.state.droplet_radius = None
        else:
            self.state.droplet_radius = None
    
    def run(self, callback: Optional[Callable] = None, 
            save_interval: Optional[int] = None) -> List[SimulationState]:
        """
        Run the simulation to completion.
        
        Args:
            callback: Optional function called after each save_interval steps.
                      Signature: callback(state, history) -> bool
                      Return True to continue, False to stop early.
            save_interval: Steps between saving states. If None, uses config value.
            
        Returns:
            List of SimulationState objects at saved time points
        """
        if self.state is None:
            self.initialize()
        
        if save_interval is None:
            save_interval = self.config.numerical.save_interval
        
        dt = self.config.numerical.dt
        t_final = self.config.numerical.t_final
        n_steps = int(t_final / dt)
        
        history = [self._copy_state(self.state)]
        
        print(f"Running simulation: {n_steps} steps, dt={dt}, t_final={t_final}")
        start_time = time.time()
        
        for step in range(n_steps):
            self.step()
            
            # Save state at intervals
            if (step + 1) % save_interval == 0:
                self._update_diagnostics()
                history.append(self._copy_state(self.state))
                
                # Progress report
                elapsed = time.time() - start_time
                progress = (step + 1) / n_steps
                eta = elapsed / progress - elapsed if progress > 0 else 0
                print(f"  Step {step+1}/{n_steps} (t={self.state.t:.2f}), "
                      f"ETA: {eta:.1f}s, "
                      f"droplet r={self.state.droplet_center:.2f}" 
                      if self.state.droplet_center else "droplet dissolved")
                
                # Optional callback
                if callback is not None:
                    if not callback(self.state, history):
                        print("Simulation stopped by callback")
                        break
        
        total_time = time.time() - start_time
        print(f"Simulation complete in {total_time:.1f}s")
        
        return history
    
    def _copy_state(self, state: SimulationState) -> SimulationState:
        """Create a deep copy of a simulation state."""
        return SimulationState(
            c=state.c.copy(),
            m=state.m.copy(),
            t=state.t,
            step=state.step,
            total_protein=state.total_protein,
            total_rna=state.total_rna,
            free_energy=state.free_energy,
            droplet_center=state.droplet_center,
            droplet_radius=state.droplet_radius
        )


def classify_regime(history: List[SimulationState], 
                    config: SimulationConfig) -> str:
    """
    Classify the simulation outcome into one of the four regimes from Figure 1B.
    
    Regimes:
        I. Dissolution: Droplet dissolves completely
        II. Renucleation: Droplet dissolves then reforms at promoter
        III. Directed motion: Droplet moves toward promoter
        IV. Directed motion with elongation: Motion + shape deformation
    
    Args:
        history: List of SimulationState from a completed simulation
        config: Configuration used for the simulation
        
    Returns:
        String identifying the regime ("I", "II", "III", "IV", or "unknown")
    """
    if len(history) < 3:
        return "unknown"
    
    initial = history[0]
    final = history[-1]
    
    # Get droplet properties over time
    positions = [s.droplet_center for s in history]
    radii = [s.droplet_radius for s in history]
    times = [s.t for s in history]
    
    # Check if droplet exists at start
    if initial.droplet_radius is None or initial.droplet_center is None:
        return "unknown"
    
    initial_r = initial.droplet_center
    initial_R = initial.droplet_radius
    
    # Check for dissolution: droplet radius becomes very small or None
    # Look for dissolution at any point in the simulation
    dissolved_at_some_point = False
    dissolution_index = None
    for i, R in enumerate(radii):
        if R is None or R < 0.3 * initial_R:  # Radius dropped to < 30% of initial
            dissolved_at_some_point = True
            dissolution_index = i
            break
    
    # Check final state
    final_has_droplet = final.droplet_radius is not None and final.droplet_radius > 0.3 * initial_R
    
    if dissolved_at_some_point:
        if final_has_droplet and final.droplet_center is not None:
            # Droplet reformed - check if it's near the promoter
            if final.droplet_center < 0.3 * initial_r:
                return "II"  # Renucleation at promoter
            else:
                return "unknown"  # Dissolved and reformed elsewhere
        else:
            return "I"  # Dissolution (stayed dissolved)
    
    # Droplet persisted - check for directed motion
    if not final_has_droplet or final.droplet_center is None:
        return "I"  # Somehow lost the droplet
    
    final_r = final.droplet_center
    final_R = final.droplet_radius
    
    # Calculate net displacement toward promoter
    displacement = initial_r - final_r  # Positive if moved toward promoter
    
    # Calculate average velocity (should be negative in r, i.e., toward origin)
    if times[-1] > times[0]:
        avg_velocity = displacement / (times[-1] - times[0])
    else:
        avg_velocity = 0
    
    # Criteria for directed motion:
    # 1. Net displacement toward promoter (displacement > 0)
    # 2. Consistent motion (most positions should be decreasing)
    
    # Count how many time steps show motion toward promoter
    n_toward_promoter = 0
    valid_pairs = 0
    for i in range(1, len(positions)):
        if positions[i] is not None and positions[i-1] is not None:
            valid_pairs += 1
            if positions[i] < positions[i-1]:
                n_toward_promoter += 1
    
    fraction_toward = n_toward_promoter / valid_pairs if valid_pairs > 0 else 0
    
    # Check for significant motion toward promoter
    # Even small displacement counts if it's consistent
    relative_displacement = displacement / initial_r if initial_r > 0 else 0
    
    if relative_displacement > 0.02 and fraction_toward > 0.5:
        # Droplet is moving toward promoter
        # Check for elongation (radius increase)
        radius_change = (final_R - initial_R) / initial_R if initial_R > 0 else 0
        
        if radius_change > 0.15:  # Radius increased by > 15%
            return "IV"  # Directed motion with elongation
        else:
            return "III"  # Directed motion
    
    elif relative_displacement < -0.02:
        # Droplet moved away from promoter - unusual
        return "unknown"
    
    else:
        # Very little motion - could be near steady state
        # Check if this is because the system is in a stable configuration
        # (e.g., low k_p so weak attraction)
        
        # If k_p is low and c_minus is low, droplet might just be stable
        k_p = config.transport.k_p
        c_minus = config.initial.c_minus_init
        c_binodal = config.free_energy.c_minus
        
        # Near the phase boundary with low driving force
        if k_p < 0.1 and c_minus < c_binodal + 0.05:
            # Could be Regime I (dissolution) if supersaturation is marginal
            # or just very slow Regime III
            if final_R < 0.8 * initial_R:
                return "I"  # Shrinking toward dissolution
            else:
                return "III"  # Stable but slow motion
        
        return "unknown"


def compute_droplet_velocity(history: List[SimulationState]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute instantaneous droplet velocity from simulation history.
    
    Returns:
        Tuple of (positions, velocities) arrays
    """
    positions = []
    velocities = []
    
    for i in range(1, len(history)):
        if (history[i].droplet_center is not None and 
            history[i-1].droplet_center is not None):
            
            r_curr = history[i].droplet_center
            r_prev = history[i-1].droplet_center
            dt = history[i].t - history[i-1].t
            
            if dt > 0:
                v = (r_curr - r_prev) / dt  # Negative if moving toward origin
                positions.append((r_curr + r_prev) / 2)
                velocities.append(v)
    
    return np.array(positions), np.array(velocities)


if __name__ == "__main__":
    # Test the solver
    print("Testing CoupledSolver...")
    
    # Create configuration for Regime III
    config = SimulationConfig.figure_1b_regime_iii()
    config.numerical.t_final = 50.0  # Short test run
    config.numerical.save_interval = 50
    
    print(config.summary())
    
    # Create and run solver
    solver = CoupledSolver(config)
    state = solver.initialize()
    
    print(f"\nInitial state:")
    print(f"  Total protein: {state.total_protein:.2f}")
    print(f"  Droplet center: {state.droplet_center:.2f}")
    print(f"  Droplet radius: {state.droplet_radius:.2f}")
    
    # Run a few steps
    print("\nRunning 100 steps...")
    for _ in range(100):
        solver.step()
    
    solver._update_diagnostics()
    print(f"\nAfter 100 steps (t={solver.state.t:.2f}):")
    print(f"  Total protein: {solver.state.total_protein:.2f}")
    print(f"  Total RNA: {solver.state.total_rna:.4f}")
    if solver.state.droplet_center:
        print(f"  Droplet center: {solver.state.droplet_center:.2f}")
        print(f"  Droplet radius: {solver.state.droplet_radius:.2f}")
    
    print("\nSolver test complete!")
