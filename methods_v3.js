const docxLib = require("docx");
const fss = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType, ImageRun
} = docxLib;

const t = (text, opts = {}) => new TextRun({ text, ...opts });
const tb = (text, opts = {}) => new TextRun({ text, bold: true, ...opts });
const ti = (text, opts = {}) => new TextRun({ text, italics: true, ...opts });
const p = (children, opts = {}) => new Paragraph({ children: Array.isArray(children) ? children : [children], ...opts });

const spc = { before: 100, after: 100, line: 276 };
const eqSpc = { before: 140, after: 140, line: 276 };
const capSpc = { before: 40, after: 180, line: 240 };

// Load images
const imgExprVsEP = fss.readFileSync("/sessions/wonderful-lucid-euler/mnt/condensate_project/images/Expression_vs_EP.png");
const imgDistDep = fss.readFileSync("/sessions/wonderful-lucid-euler/mnt/condensate_project/images/Distance_Dependent_Expression.png");
const imgSuperTyp = fss.readFileSync("/sessions/wonderful-lucid-euler/mnt/condensate_project/images/Super_vs_Typical.png");
const imgAnalysis = fss.readFileSync("/sessions/wonderful-lucid-euler/mnt/condensate_project/data/experimental/analysis_summary.png");

function placeholder(label) {
  return [
    new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 180, after: 40 },
      border: {
        top: { style: BorderStyle.DASHED, size: 1, color: "AAAAAA", space: 6 },
        bottom: { style: BorderStyle.DASHED, size: 1, color: "AAAAAA", space: 6 },
        left: { style: BorderStyle.DASHED, size: 1, color: "AAAAAA", space: 6 },
        right: { style: BorderStyle.DASHED, size: 1, color: "AAAAAA", space: 6 },
      },
      children: [t(`[ INSERT ${label} ]`, { bold: true, color: "999999", size: 22 })]
    })
  ];
}

function img(buf, w, h, maxW) {
  maxW = maxW || 5486400;
  const aspect = h / w;
  const rw = Math.round(maxW / 9525);
  const rh = Math.round(maxW * aspect / 9525);
  return new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 180, after: 40 },
    children: [new ImageRun({ data: buf, transformation: { width: rw, height: rh }, type: "png" })]
  });
}

function cap(text) {
  return p([ti(text, { size: 20 })], { spacing: capSpc, alignment: AlignmentType.CENTER });
}

const doc = new Document({
  styles: {
    default: { document: { run: { font: "Times New Roman", size: 24 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, italics: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 160, after: 80 }, outlineLevel: 2 } },
    ]
  },
  numbering: { config: [] },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }
      }
    },
    children: [

      // ========= TITLE =========
      new Paragraph({ heading: HeadingLevel.HEADING_1, children: [t("3. Methods")] }),

      // ========= 3.1 PHASE FIELD MODEL =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.1 Phase field model")] }),

      p([
        t("We model the transcriptional condensate as a binary protein-solvent mixture coupled to a non-conserved RNA field. The protein concentration "),
        ti("c"),
        t("("),
        tb("r"),
        t(", "),
        ti("t"),
        t(") evolves by conserved Cahn-Hilliard dynamics; the RNA concentration "),
        ti("m"),
        t("("),
        tb("r"),
        t(", "),
        ti("t"),
        t(") evolves by reaction-diffusion with localised production and uniform degradation. This follows Goh "),
        ti("et al."),
        t(" [1], where electrostatic attraction between positively charged IDR proteins and negatively charged nascent RNA provides the coupling."),
      ], { spacing: spc }),

      // --- FIGURE: Schematic of the model ---
      ...placeholder("FIGURE 1: Schematic \u2014 condensate droplet, promoter at origin, RNA gradient, and the feedback loop (RNA attracts condensate \u2192 condensate lands on promoter \u2192 more RNA produced)"),
      cap("Figure 1. Schematic of the coupled model. A protein condensate (blue) sits at distance r from the promoter (green, at origin). The promoter produces RNA (red gradient) proportional to local protein concentration. The RNA gradient creates an asymmetric chemical potential that pulls the condensate toward the source."),

      p([t("The total free energy of the system is")], { spacing: { before: 100, after: 30 } }),

      // Eq 1
      p([
        ti("F"),
        t("["),
        ti("c"),
        t(", "),
        ti("m"),
        t("] = \u222B [ (\u03B1/4)("),
        ti("c"),
        t(" \u2013 "),
        ti("c\u0304"),
        t(")\u2074 + (\u03B2/2)("),
        ti("c"),
        t(" \u2013 "),
        ti("c\u0304"),
        t(")\u00B2 + (\u03BA/2)|\u2207"),
        ti("c"),
        t("|\u00B2 + \u03C7"),
        ti("cm"),
        t(" + (\u03B3/2)"),
        ti("c"),
        t("\u00B2"),
        ti("m"),
        t("\u00B2 ] d"),
        tb("r"),
        t("          (1)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("The first two terms form a double-well potential centred at "),
        ti("c\u0304"),
        t(" = 4.0. With \u03B1 = 1.0 and \u03B2 = \u20130.25, the binodal concentrations are "),
        ti("c"),
        t("\u207B = 3.5 (dilute) and "),
        ti("c"),
        t("\u207A = 4.5 (dense). The gradient term (\u03BA = 0.05) penalises sharp interfaces, giving an equilibrium interface width of \u221A(2\u03BA/|\u03B2|) \u2248 0.63. The \u03C7 term (\u03C7 = \u20130.1) is an attractive protein-RNA coupling that lowers the chemical potential in RNA-rich regions. We set \u03B3 = 0 (no reentrant repulsion) throughout, consistent with Figure 1B of [1]."),
      ], { spacing: spc }),

      // --- FIGURE: Double well potential plot ---
      ...placeholder("FIGURE 2: Double-well free energy f(c) showing the two minima at c\u207B = 3.5 and c\u207A = 4.5, with the critical point at c\u0304 = 4.0"),
      cap("Figure 2. The bulk free energy density f(c) = (\u03B1/4)(c \u2013 c\u0304)\u2074 + (\u03B2/2)(c \u2013 c\u0304)\u00B2, showing the two stable phases at c\u207B = 3.5 and c\u207A = 4.5. The condensate (dense phase) sits at c\u207A; the surrounding solvent at c\u207B."),

      p([t("Protein obeys the Cahn-Hilliard equation, which conserves total mass:")], { spacing: { before: 100, after: 30 } }),

      // Eq 2
      p([
        t("\u2202"),
        ti("c"),
        t("/\u2202"),
        ti("t"),
        t(" = \u2207 \u00B7 ("),
        ti("M"),
        t("\u2080 \u2207\u03BC),     \u03BC = \u03B1("),
        ti("c"),
        t(" \u2013 "),
        ti("c\u0304"),
        t(")\u00B3 + \u03B2("),
        ti("c"),
        t(" \u2013 "),
        ti("c\u0304"),
        t(") \u2013 \u03BA\u2207\u00B2"),
        ti("c"),
        t(" + \u03C7"),
        ti("m"),
        t("          (2)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("The divergence form guarantees mass conservation. The mobility "),
        ti("M"),
        t("\u2080 = 1.0 is constant. The cubic term drives phase separation; the \u2013\u03BA\u2207\u00B2"),
        ti("c"),
        t(" term opposes sharp gradients; and the \u03C7"),
        ti("m"),
        t(" term biases protein flux toward high-RNA regions. RNA obeys:"),
      ], { spacing: spc }),

      // Eq 3
      p([
        t("\u2202"),
        ti("m"),
        t("/\u2202"),
        ti("t"),
        t(" = "),
        ti("D"),
        t("\u2098\u2207\u00B2"),
        ti("m"),
        t(" + "),
        ti("k"),
        t("\u209A"),
        ti("S"),
        t("("),
        tb("r"),
        t(")"),
        ti("c"),
        t(" \u2013 "),
        ti("k"),
        t("\u2091"),
        ti("m"),
        t("          (3)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("where "),
        ti("S"),
        t("("),
        tb("r"),
        t(") = exp(\u2013|"),
        tb("r"),
        t("|\u00B2/2\u03C3\u209A\u00B2) is a Gaussian source at the promoter (\u03C3\u209A = 2.5). Production requires both proximity to the promoter and protein presence. The RNA diffusion length \u221A("),
        ti("D"),
        t("\u2098/"),
        ti("k"),
        t("\u2091) = 1.0 sets how far the gradient reaches. This is the key feedback loop: RNA attracts the condensate (\u03C7 < 0), and the condensate amplifies RNA production when it reaches the promoter."),
      ], { spacing: spc }),


      // ========= 3.2 FOUR REGIMES =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.2 Dynamical regimes")] }),

      p([
        t("The interplay of two control parameters \u2014 the RNA production rate "),
        ti("k"),
        t("\u209A and the dilute-phase protein concentration "),
        ti("c"),
        t("\u207B \u2014 gives rise to four distinct dynamical regimes, illustrated in Figure 3."),
      ], { spacing: spc }),

      // --- FIGURE: 2D regime snapshots (4-panel) ---
      ...placeholder("FIGURE 3: Four-panel showing 2D simulation snapshots for Regimes I\u2013IV (protein top row, RNA bottom row) at selected time points \u2014 dissolution, renucleation, directed motion, elongation"),
      cap("Figure 3. Snapshots from 2D simulations of the four dynamical regimes. Each panel shows the protein field (top) and RNA field (bottom) at selected time points. I: dissolution \u2014 the droplet shrinks and vanishes. II: renucleation \u2014 the droplet dissolves but reforms at the promoter. III: directed motion \u2014 the droplet migrates intact toward the promoter. IV: elongation \u2014 the droplet stretches as it moves inward."),

      // Parameter table
      (() => {
        const border = { style: BorderStyle.SINGLE, size: 1, color: "888888" };
        const borders = { top: border, bottom: border, left: border, right: border };
        const cm = { top: 50, bottom: 50, left: 80, right: 80 };
        const hdr = { fill: "D9E2F3", type: ShadingType.CLEAR };
        const hc = (txt, w) => new TableCell({ borders, width: { size: w, type: WidthType.DXA }, shading: hdr, margins: cm, children: [p([tb(txt, { size: 20 })], { alignment: AlignmentType.CENTER })] });
        const c = (txt, w) => new TableCell({ borders, width: { size: w, type: WidthType.DXA }, margins: cm, children: [p([t(txt, { size: 20 })], { alignment: AlignmentType.CENTER })] });
        return new Table({
          width: { size: 9360, type: WidthType.DXA },
          columnWidths: [1800, 2520, 2520, 2520],
          rows: [
            new TableRow({ children: [hc("Regime", 1800), hc("k\u209A", 2520), hc("c\u207B", 2520), hc("Behaviour", 2520)] }),
            new TableRow({ children: [c("I", 1800), c("0.05", 2520), c("3.51", 2520), c("Dissolution", 2520)] }),
            new TableRow({ children: [c("II", 1800), c("0.40", 2520), c("3.51", 2520), c("Renucleation", 2520)] }),
            new TableRow({ children: [c("III", 1800), c("0.08", 2520), c("3.53", 2520), c("Directed motion", 2520)] }),
            new TableRow({ children: [c("IV", 1800), c("0.25", 2520), c("3.60", 2520), c("Elongation", 2520)] }),
          ]
        });
      })(),

      p([
        ti("Table 1. Regime-specific parameters. Fixed across all runs: \u03B1 = 1, \u03B2 = \u20130.25, \u03BA = 0.05, c\u0304 = 4.0, \u03C7 = \u20130.1, \u03B3 = 0, M\u2080 = 1, D\u2098 = 1, k\u2091 = 1, \u03C3\u209A = 2.5, c\u207A(0) = 5.5, R\u2080 = 2, r\u2080 = 10."),
        t(" 2D grid: 256\u00B2 cells on [\u221225, 25]\u00B2, \u0394t = 5\u00D710\u207B\u2074."),
      ], { spacing: { before: 40, after: 100 }, size: 20 }),

      p([
        t("When "),
        ti("c"),
        t("\u207B is close to the lower binodal (3.5), there is barely enough protein to sustain a condensate; even weak RNA-induced redistribution dissolves it (Regime I). At the same low "),
        ti("c"),
        t("\u207B but high "),
        ti("k"),
        t("\u209A, the RNA gradient is strong enough to nucleate a new condensate directly at the promoter (Regime II). At intermediate "),
        ti("c"),
        t("\u207B the droplet is stable and the RNA gradient pulls it inward (Regime III), and at higher "),
        ti("c"),
        t("\u207B the droplet elongates as it moves (Regime IV) because protein accumulates at the leading edge faster than the trailing edge dissolves."),
      ], { spacing: spc }),

      // --- FIGURE: Phase diagram ---
      ...placeholder("FIGURE 4: Phase diagram in (k\u209A, c\u207B) plane, coloured by regime, reproducing Figure 1B of Goh et al."),
      cap("Figure 4. Phase diagram in the (k\u209A, c\u207B) parameter plane. Each point is coloured by its classified regime. The four regimes occupy distinct regions consistent with Figure 1B of [1]."),


      // ========= 3.3 NUMERICAL IMPLEMENTATION =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.3 Numerical implementation")] }),

      p([
        t("We solve the coupled system on a 2D Cartesian grid (256\u00D7256 cells, domain [\u221225, 25]\u00B2) using a split-step scheme. The promoter sits at the origin and the initial droplet at (10, 0), initialised with a smooth tanh interface profile:"),
      ], { spacing: spc }),

      // Eq 4
      p([
        ti("c"),
        t("("),
        tb("r"),
        t(", 0) = ("),
        ti("c"),
        t("\u207A + "),
        ti("c"),
        t("\u207B)/2 + ("),
        ti("c"),
        t("\u207A \u2013 "),
        ti("c"),
        t("\u207B)/2 \u00B7 tanh[("),
        ti("R"),
        t("\u2080 \u2013 |"),
        tb("r"),
        t(" \u2013 "),
        tb("r"),
        t("\u2080|) / "),
        ti("w"),
        t("]          (4)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("with "),
        ti("R"),
        t("\u2080 = 2.0, "),
        tb("r"),
        t("\u2080 = (10, 0), and interface width "),
        ti("w"),
        t(" = \u221A(2\u03BA/|\u03B2|) \u2248 0.63. The Laplacian uses a five-point stencil assembled as a sparse matrix (SciPy CSR, 65,536\u00D765,536 with ~325k nonzeros). No-flux boundary conditions are applied by zeroing boundary face fluxes. Figure 5 illustrates the grid setup and initial condition."),
      ], { spacing: spc }),

      // --- FIGURE: Grid/IC illustration ---
      ...placeholder("FIGURE 5: Initial condition \u2014 2D heatmap of protein field at t = 0 showing the circular droplet at (10, 0) with the promoter marked at the origin"),
      cap("Figure 5. Initial protein concentration field. The circular condensate droplet (c\u207A = 5.5) is placed at (10, 0) in a dilute background (c\u207B). The promoter at the origin is marked by a cross."),

      p([
        t("At each time step, we first update RNA semi-implicitly:"),
      ], { spacing: { before: 100, after: 30 } }),

      // Eq 5
      p([
        t("("),
        ti("I"),
        t(" \u2013 \u0394"),
        ti("t"),
        t(" "),
        ti("D"),
        t("\u2098"),
        tb("L"),
        t(") "),
        ti("m"),
        t("\u207F\u207A\u00B9 = "),
        ti("m"),
        t("\u207F + \u0394"),
        ti("t"),
        t("("),
        ti("k"),
        t("\u209A"),
        ti("S"),
        t("\u2299"),
        ti("c"),
        t("\u207F \u2013 "),
        ti("k"),
        t("\u2091"),
        ti("m"),
        t("\u207F)          (5)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("The system matrix is constant, so we precompute its sparse LU factorisation once and reuse it. We then update protein explicitly:"),
      ], { spacing: spc }),

      // Eq 6
      p([
        ti("c"),
        t("\u207F\u207A\u00B9 = "),
        ti("c"),
        t("\u207F + \u0394"),
        ti("t"),
        t(" "),
        ti("M"),
        t("\u2080\u2207\u00B2\u03BC          (6)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("The explicit Cahn-Hilliard step requires \u0394"),
        ti("t"),
        t(" < \u0394"),
        ti("x"),
        t("\u2074/(32\u03BA"),
        ti("M"),
        t("\u2080) \u2248 9\u00D710\u207B\u2074; we use \u0394"),
        ti("t"),
        t(" = 5\u00D710\u207B\u2074. Mass conservation holds to ~10\u207B\u00B9\u2076 relative error."),
      ], { spacing: spc }),

      p([
        t("We also built a 1D radial solver (200 cell-centred finite volumes, cylindrical Laplacian) for rapid parameter sweeps. However, the 1D geometry treats the condensate as a ring, producing an artificial outward drift from curvature asymmetry that competes with the RNA attraction. All reported results therefore use the 2D solver."),
      ], { spacing: spc }),


      // ========= 3.4 DROPLET TRACKING =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.4 Droplet tracking")] }),

      p([
        t("We define the droplet as cells where "),
        ti("c"),
        t(" > "),
        ti("c\u0304"),
        t(" and compute its centre of mass as:"),
      ], { spacing: spc }),

      // Eq 7
      p([
        tb("r"),
        t("\u209c\u2098 = \u03A3 "),
        tb("r"),
        t("\u1D62("),
        ti("c"),
        t("\u1D62 \u2013 "),
        ti("c\u0304"),
        t(")"),
        ti("V"),
        t("\u1D62 / \u03A3 ("),
        ti("c"),
        t("\u1D62 \u2013 "),
        ti("c\u0304"),
        t(")"),
        ti("V"),
        t("\u1D62          (7)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("summing only over dense-phase cells. The droplet aspect ratio is extracted from the inertia tensor eigenvalues. Figure 6 shows representative tracking results."),
      ], { spacing: spc }),

      // --- FIGURE: Distance vs time ---
      ...placeholder("FIGURE 6: Droplet distance from promoter vs. time for all four regimes, showing dissolution (I), renucleation (II), steady inward drift (III), and elongated inward motion (IV)"),
      cap("Figure 6. Droplet centre-of-mass distance from the promoter over time for each regime. Regime III (directed motion) shows steady inward drift; Regime IV shows similar drift accompanied by increasing aspect ratio."),


      // ========= 3.5 ANALYTICAL VELOCITY =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.5 Analytical droplet velocity")] }),

      p([
        t("In the sharp-interface limit (thin interface, constant interior concentration), the droplet velocity can be derived analytically [1, Eq. 11]. The idea is that the RNA gradient creates an asymmetric chemical potential across the droplet surface; the side closer to the promoter sees more RNA, so \u03BC is lower there, driving a net protein flux toward the source. The velocity in 3D spherical geometry is:"),
      ], { spacing: spc }),

      // Eq 8
      p([
        ti("v"),
        t(" = \u2013("),
        ti("M"),
        t("\u2080\u03C7"),
        ti("d"),
        t(" / \u0394"),
        ti("cV"),
        t(") \u222B\u2091 \u2207"),
        ti("m"),
        t("("),
        tb("r"),
        t(") \u00B7 "),
        tb("e\u0302"),
        t("\u1D65 d\u1D48"),
        tb("r"),
        t("          (8)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("where "),
        ti("d"),
        t(" is the spatial dimension, \u0394"),
        ti("c"),
        t(" = "),
        ti("c"),
        t("\u207A \u2013 "),
        ti("c"),
        t("\u207B, "),
        ti("V"),
        t(" is droplet volume, and the integral is over the droplet domain. The velocity is proportional to the RNA gradient asymmetry across the condensate. It peaks at "),
        ti("r"),
        t(" \u2248 "),
        ti("R"),
        t(" (when the leading edge just touches the source) and decays at larger distances. Figure 7 compares the analytical prediction with simulation data."),
      ], { spacing: spc }),

      // --- FIGURE: Velocity vs distance ---
      ...placeholder("FIGURE 7: Velocity vs. distance from promoter \u2014 analytical curve (Eq. 8) overlaid with 2D simulation measurements, showing the non-monotonic profile peaking near r \u2248 R"),
      cap("Figure 7. Droplet velocity as a function of distance from the promoter. The analytical prediction (solid line) peaks near r \u2248 R and decays at large distances. Simulation data (points) follow the same trend. The non-monotonic profile arises because the RNA gradient is most asymmetric when the droplet just contacts the source."),


      // ========= 3.6 ROUSE POLYMER =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.6 Rouse polymer model")] }),

      p([
        t("To connect condensate motion to gene regulation, we model the chromatin fibre as a Rouse chain: "),
        ti("N"),
        t(" beads connected by harmonic springs, each undergoing overdamped Brownian dynamics [1, Section V]. The enhancer is pinned at one end, the promoter at the other. A condensate exerts a harmonic trap on the promoter bead, pulling it toward the condensate position predicted by the phase field model. Figure 8 illustrates this setup."),
      ], { spacing: spc }),

      // --- FIGURE: Rouse polymer schematic ---
      ...placeholder("FIGURE 8: Schematic of Rouse polymer model \u2014 chain of beads with enhancer at one end, promoter at other, condensate pulling the promoter bead via a harmonic trap (similar to Figure 3 from your notes)"),
      cap("Figure 8. Rouse polymer model for enhancer-promoter contacts. The chromatin is modelled as N beads (grey) connected by springs. The enhancer (red, bead 1) is fixed; the condensate (blue) pulls the promoter bead (green, bead N) via a harmonic potential. Contact occurs when the enhancer-promoter distance falls below a threshold d\u209C."),

      p([t("The equation of motion for each bead is:")], { spacing: { before: 80, after: 30 } }),

      // Eq 9
      p([
        t("\u03B6 d"),
        ti("x"),
        t("\u1D62/d"),
        ti("t"),
        t(" = "),
        ti("k"),
        t("\u209B("),
        ti("x"),
        t("\u1D62\u208A\u2081 + "),
        ti("x"),
        t("\u1D62\u208B\u2081 \u2013 2"),
        ti("x"),
        t("\u1D62) + "),
        ti("F"),
        t("\u1D62 + \u03B6\u221A(2"),
        ti("D"),
        t("/\u0394"),
        ti("t"),
        t(")\u03B7\u1D62          (9)"),
      ], { spacing: eqSpc, alignment: AlignmentType.CENTER }),

      p([
        t("where the spring force term "),
        ti("k"),
        t("\u209B("),
        ti("x"),
        t("\u1D62\u208A\u2081 + "),
        ti("x"),
        t("\u1D62\u208B\u2081 \u2013 2"),
        ti("x"),
        t("\u1D62) is a discrete Laplacian reflecting Hookean interactions between neighbours, and \u03B7\u1D62 is Gaussian thermal noise satisfying the fluctuation-dissipation theorem. The external force "),
        ti("F"),
        t("\u1D62 is zero except at the promoter bead, where it pulls toward the condensate. We run 10\u2075 steps, discard the first half for equilibration, and count contacts (enhancer-promoter distance < "),
        ti("d"),
        t("\u209C = 0.5) in the remainder. By sweeping genomic distance (10\u2013150 beads) and trap strength, we produce the contact probability heatmap in Figure 9."),
      ], { spacing: spc }),

      // --- FIGURE: Contact probability heatmap ---
      ...placeholder("FIGURE 9: Contact probability heatmap (YlOrRd colourmap) with genomic distance on y-axis and trap strength on x-axis"),
      cap("Figure 9. Contact probability heatmap from Rouse polymer simulations. Contact probability increases with trap strength (proportional to k\u209A) and decreases with genomic distance, consistent with the model\u2019s prediction that condensate-mediated looping is strongest for highly transcribed genes at intermediate separations."),


      // ========= 3.7 EXPERIMENTAL VALIDATION =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.7 Experimental validation")] }),

      p([
        t("We tested the model predictions against mESC Micro-C data from Hsieh "),
        ti("et al."),
        t(" [2] and matched RNA-seq expression data. The model predicts that condensate-mediated looping should increase contact frequency preferentially for highly transcribed genes (high "),
        ti("k"),
        t("\u209A) and that this effect should be strongest at intermediate genomic distances. Figures 10\u201313 show the experimental results."),
      ], { spacing: spc }),

      // --- FIGURE 10: Expression vs E-P contact ---
      img(imgExprVsEP, 640, 480),
      cap("Figure 10. Gene expression (log\u2082(TPM + 1)) versus enhancer-promoter contact frequency (log O/E) for mESC enhancer-promoter pairs. The positive correlation (Spearman \u03C1 = 0.372) supports the model prediction that more highly transcribed genes (higher k\u209A) show elevated contact."),

      // --- FIGURE 11: Distance-dependent correlation ---
      img(imgDistDep, 1000, 500),
      cap("Figure 11. Spearman correlation between expression and contact frequency, binned by genomic distance. The correlation peaks at 75\u2013150 kb, consistent with the model prediction that condensate-mediated looping is most effective at intermediate distances where the RNA gradient can bridge the gap."),

      // --- FIGURE 12: Super vs Typical ---
      img(imgSuperTyp, 1000, 400),
      cap("Figure 12. Expression vs. contact for super-enhancers (left, \u03C1 = 0.389) and typical enhancers (right, \u03C1 = 0.308). The stronger correlation at super-enhancers is consistent with larger condensates and steeper RNA gradients at multi-gene transcriptional hubs."),

      // --- FIGURE 13: Analysis summary ---
      img(imgAnalysis, 800, 800),
      cap("Figure 13. Summary of experimental analysis. (A) Expression vs. E-P contacts. (B) Distribution of E-P distances. (C) Distance-stratified correlation (orange: significant). (D) Expression distribution."),


      // ========= 3.8 SOFTWARE =========
      new Paragraph({ heading: HeadingLevel.HEADING_2, children: [t("3.8 Software")] }),

      p([
        t("All code was written in Python 3 (NumPy, SciPy, Matplotlib). The codebase is modular: separate packages for physics, numerics, solvers, and analysis. All parameters are stored as serialisable dataclasses with factory methods for each regime preset. A 2D simulation of 500 time units runs in approximately 30 minutes on a single CPU core."),
      ], { spacing: spc }),

    ]
  }]
});

Packer.toBuffer(doc).then(buffer => {
  fss.writeFileSync("/sessions/wonderful-lucid-euler/mnt/condensate_project/methods_section.docx", buffer);
  console.log("Done.");
});
