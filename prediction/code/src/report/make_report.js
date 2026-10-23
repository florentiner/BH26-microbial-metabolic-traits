const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer,
} = require('docx');
const fs = require('fs');

const ACCENT = '2A78D6';
const MUTED = '52514E';
const HEAD_BG = 'EDF2FA';
const ALT_BG = 'F7F8FA';
const TOTAL = 9360; // DXA, US Letter with 1" margins

const P = (text, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 120, line: 276 },
  alignment: opts.align,
  children: [new TextRun({ text, bold: opts.bold, italics: opts.italics, size: opts.size ?? 21,
    color: opts.color, font: 'Calibri' })],
});

const Rich = (runs, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 120, line: 276 },
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i, size: 21, color: r.c, font: 'Calibri' })),
});

const H = (text, level) => new Paragraph({
  heading: level,
  spacing: { before: level === HeadingLevel.HEADING_1 ? 360 : 260, after: 140 },
  children: [new TextRun({ text, font: 'Calibri', bold: true,
    size: level === HeadingLevel.HEADING_1 ? 30 : 24,
    color: level === HeadingLevel.HEADING_1 ? ACCENT : '1A1A1A' })],
});

const bullet = (text, opts = {}) => new Paragraph({
  numbering: { reference: 'bul', level: 0 },
  spacing: { after: 80, line: 276 },
  children: [new TextRun({ text, size: 21, font: 'Calibri', bold: opts.bold })],
});

// text scheme (Figure 1): bordered text boxes joined by arrows, built from an invisible 5-column layout table
const SC = { shared: ['F1F1F1', '8C8C8C'], ko: ['E3EEFA', '2A78D6'], lm: ['FDEBD6', 'D98200'], out: ['E2F2E6', '2E7D32'] };
const SW = [1300, 2990, 780, 2990, 1300];                     // side | left box | gap | right box | side = 9360
const NOB = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
const spanW = (a, n) => SW.slice(a, a + n).reduce((x, y) => x + y, 0);
const schemeCell = (a, n, lines, kind) => {
  const [fill, edge] = kind ? SC[kind] : [undefined, undefined];
  const b = kind ? { style: BorderStyle.SINGLE, size: 8, color: edge } : NOB;
  return new TableCell({ columnSpan: n, width: { size: spanW(a, n), type: WidthType.DXA },
    shading: fill ? { type: ShadingType.CLEAR, fill, color: 'auto' } : undefined,
    borders: { top: b, bottom: b, left: b, right: b },
    margins: { top: kind ? 60 : 0, bottom: kind ? 60 : 0, left: 100, right: 100 },
    children: (lines.length ? lines : ['']).map((l, k) => new Paragraph({ alignment: AlignmentType.CENTER,
      spacing: { after: 0, line: 240 },
      children: [new TextRun({ text: l, bold: kind ? k === 0 : l.length <= 3, size: kind ? (k === 0 ? 19 : 17) : (l.length <= 3 ? 22 : 17),
        color: kind ? (k === 0 ? '1A1A1A' : '404040') : '6B6B6B', font: 'Calibri' })] })) });
};
// row types: centred box, two boxes side by side, arrows
const rowCentre = (lines, kind) => new TableRow({ cantSplit: true, children: [schemeCell(0, 1, [], null),
  schemeCell(1, 3, lines, kind), schemeCell(4, 1, [], null)] });
const rowPair = (left, right) => new TableRow({ cantSplit: true, children: [
  schemeCell(0, 2, left ? left[0] : [], left ? left[1] : null), schemeCell(2, 1, [], null),
  schemeCell(3, 2, right ? right[0] : [], right ? right[1] : null)] });
const rowArrow = (t = '↓') => new TableRow({ cantSplit: true, children: [schemeCell(0, 1, [], null),
  schemeCell(1, 3, [t], null), schemeCell(4, 1, [], null)] });
const rowArrows = (l = '↓', r = '↓') => new TableRow({ cantSplit: true, children: [schemeCell(0, 2, [l], null),
  schemeCell(2, 1, [], null), schemeCell(3, 2, [r], null)] });
const scheme = rows => new Table({ columnWidths: SW, width: { size: 9360, type: WidthType.DXA },
  borders: { top: NOB, bottom: NOB, left: NOB, right: NOB, insideHorizontal: NOB, insideVertical: NOB }, rows });

// LLM baseline comparison, read from results/baseline2_llm_summary.csv (src/21_baseline2_llm.py score)
const LLM = (() => {
  const [head, ...lines] = fs.readFileSync('./results/baseline2_llm_summary.csv', 'utf8').trim().split('\n');
  const cols = head.split(',');
  return lines.map(l => Object.fromEntries(l.split(',').map((v, k) => [cols[k], v === '' ? null : (isNaN(+v) ? v : +v)])));
})();
const llmRow = (arm, model) => LLM.find(r => r.llm_arm === arm && r.model === model) || {};
const sgn = x => (x >= 0 ? '+' : '−') + Math.abs(x).toFixed(3);
const llmDelta = (arm, model) => { const r = llmRow(arm, model);
  if (r.llm_minus_model === null || r.llm_minus_model === undefined) return '—';
  return `${sgn(r.llm_minus_model)}${(r.ci_low > 0 || r.ci_high < 0) ? ' *' : ''}`; };
function llmTable() {
  const W5 = [3000, 1150, 1500, 1650, 2060];
  const ref = 'H_genome';                                 // the 96-genome reference set for the model rows
  const rows = [
    ['Taxonomic null', llmRow(ref, 'null_taxNN').macro_auprc, 'null'],
    ['Claude — gene evidence, anonymised', llmRow('H_genome', 'H_genome').macro_auprc, 'H_genome'],
    ['GPT-6 Astra — gene evidence, anonymised', llmRow('H_gpt_genome', 'H_gpt_genome').macro_auprc, 'H_gpt_genome'],
    ['B2 (released model)', llmRow(ref, 'B2').macro_auprc, 'B2'],
    ['B2XL (best KO-only model)', llmRow(ref, 'B2XL').macro_auprc, 'B2XL'],
    ['B2XL + Bacformer-large (best model)', llmRow(ref, 'B2XLE_large').macro_auprc, 'B2XLE_large'],
    ['GPT-6 Astra + organism name (n = 84)', llmRow('H_gpt_named', 'H_gpt_named').macro_auprc, 'H_gpt_named'],
    ['Claude + organism name', llmRow('H_named', 'H_named').macro_auprc, 'H_named'],
  ].map(([name, v, key]) => {
    let d;
    if (key.startsWith('H_')) d = ['B2', 'B2XL', 'B2XLE_large'].map(m => llmDelta(key, m));
    else if (key === 'null') d = ['B2', 'B2XL', 'B2XLE_large'].map(m => sgn(v - llmRow(ref, m).macro_auprc));
    else d = ['—', '—', '—'];
    return [name, v.toFixed(3), ...d];
  });
  return table(['Baseline', 'Macro AUPRC', 'vs released B2', 'vs B2XL (KO only)', 'vs B2XL + Bacformer-large (best)'], rows, W5,
    { boldRows: ['B2 (released model)', 'B2XL (best KO-only model)', 'B2XL + Bacformer-large (best model)'] });
}

function table(headers, rows, widths, opts = {}) {
  const cell = (text, { bold, bg, align } = {}) => new TableCell({
    width: { size: widths[arguments.callee.i ?? 0], type: WidthType.DXA },
    shading: bg ? { type: ShadingType.CLEAR, fill: bg, color: 'auto' } : undefined,
    margins: { top: 60, bottom: 60, left: 90, right: 90 },
    children: [new Paragraph({
      alignment: align,
      spacing: { after: 0, line: 240 },
      children: [new TextRun({ text: String(text), bold, size: 19, font: 'Calibri' })],
    })],
  });
  const mkRow = (cells, { bg, bold } = {}) => new TableRow({
    children: cells.map((t, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: bg ? { type: ShadingType.CLEAR, fill: bg, color: 'auto' } : undefined,
      margins: { top: 60, bottom: 60, left: 90, right: 90 },
      children: [new Paragraph({
        alignment: i === 0 ? AlignmentType.LEFT : AlignmentType.RIGHT,
        spacing: { after: 0, line: 240 },
        children: [new TextRun({ text: String(t), bold: bold || (opts.boldRows || []).includes(t),
          size: 19, font: 'Calibri' })],
      })],
    })),
  });
  return new Table({
    columnWidths: widths,
    width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: {
      top: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
      bottom: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: 'E4E7EC' },
      insideVertical: { style: BorderStyle.NONE },
    },
    rows: [mkRow(headers, { bg: HEAD_BG, bold: true }),
      ...rows.map((r, k) => mkRow(r, { bg: k % 2 ? ALT_BG : undefined }))],
  });
}

const W4 = [3600, 1920, 1920, 1920];
const W3 = [4680, 2340, 2340];
const W2 = [6000, 3360];

const doc = new Document({
  creator: 'bio_pred',
  title: 'Electron donor/acceptor prediction from genomes — results report',
  numbering: {
    config: [{
      reference: 'bul',
      levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 360, hanging: 200 } } } }],
    }],
  },
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: {
      default: new Footer({ children: [new Paragraph({
        alignment: AlignmentType.RIGHT,
        children: [new TextRun({ children: ['Page ', PageNumber.CURRENT, ' of ', PageNumber.TOTAL_PAGES],
          size: 18, color: MUTED, font: 'Calibri' })],
      })] }),
    },
    children: [
      new Paragraph({
        spacing: { after: 60 },
        children: [new TextRun({ text: 'Predicting electron donors and acceptors from prokaryotic genomes',
          bold: true, size: 34, color: ACCENT, font: 'Calibri' })],
      }),
      new Paragraph({
        spacing: { after: 260 },
        border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: ACCENT } },
        children: [new TextRun({ text: 'Status report on results, baselines and publication readiness — 18 September 2026',
          size: 22, color: MUTED, font: 'Calibri' })],
      }),

      H('1. Summary', HeadingLevel.HEADING_1),
      P('586 prokaryotic genomes (358 with curated electron donor/acceptor annotations) were used to train and evaluate models that predict 14 metabolic categories (7 electron donors, 7 electron acceptors) from the genome alone. Evaluation holds out whole genera and, separately, whole families, so reported numbers describe performance on lineages the model has never seen.'),
      P('The released model uses KEGG orthologue (KO) annotation of an 852-gene energy-metabolism panel and reaches macro AUPRC 0.576 on held-out genera against 0.472 for a baseline that copies labels from the closest relatives. On 562 genomes from genera absent from training, labelled by an independent curation, it reaches 0.511 against 0.380 for that baseline (difference +0.131, 95% CI +0.100 to +0.167).'),
      Rich([{ t: 'Headline finding: ', b: true }, { t: 'six of fourteen categories are predicted well out of clade (AUPRC 0.82–0.89), three moderately, and five poorly. The poorly predicted ones are precisely those whose database labels disagree between close relatives, and an independent expert curation of the same organisms agrees with ours only at κ = 0.61 — barely more than our model agrees with that curation (κ = 0.52). The limit is the label definition, not the model.' }]),

      Rich([{ t: 'Input is the genome and nothing else. ', b: true }, { t: 'Features come only from the genome\u2019s proteins \u2014 KofamScan KO annotation and, for the fusion models, protein-language-model embeddings of the proteins in coordinate order. No taxonomy, organism name or literature enters the model (the taxonomic baseline uses taxonomy, but it is a baseline the model must beat). An accession uses the assembly\u2019s NCBI annotation; an unannotated nucleotide FASTA is gene-called first, which on a test genome gave 3 242 genes against NCBI\u2019s 3 163, a mean score difference of 0.015 and identical binary calls.' }], { after: 200 }),

      H('2. Baselines and models', HeadingLevel.HEADING_1),
      P('Macro AUPRC over 14 traits, out-of-fold, phylogeny-blocked cross-validation. Higher is better; the prevalence null is the floor.', { color: MUTED, size: 19 }),
      table(['Model / baseline', 'Held-out genus', 'Held-out family'], [
        ['Prevalence null (floor)', '0.185', '0.178'],
        ['Expert rules on marker gene complexes', '0.317', '0.317'],
        ['Taxonomic null (copy from relatives)', '0.472', '0.380'],
        ['Pooled ESM-2 35M embeddings', '0.544', '0.497'],
        ['Attention MIL over ESM-2 proteins', '0.544', '0.486'],
        ['Attention MIL, clean contrast (D_nopf)', '0.548', '0.495'],
        ['Bacformer-base (genome language model)', '0.551', '0.522'],
        ['Pooled ESM-C 300M embeddings', '0.574', '0.509'],
        ['Bacformer-large', '0.577', '0.531'],
        ['B2: 852-KO energy panel (RELEASED)', '0.576', '0.560'],
        ['B2 + second learner (B2XL)', '0.602', '0.573'],
        ['B2 + Bacformer-base', '0.604', '0.583'],
        ['B2XL + ESM-2', '0.616', '0.584'],
        ['B2XL + Bacformer-large (best internal)', '0.625', '0.585'],
      ], W3, { boldRows: ['B2: 852-KO energy panel (RELEASED)', 'B2XL + Bacformer-large (best)'] }),
      P('Protein language models alone do not beat gene content: the best of them (Bacformer-large, 0.577) ties the KO panel (0.576). Combined with gene content they add +0.03 to +0.05, and after family-wise correction across all 28 variants tried, two combinations remain significant (B2XL + Bacformer-large, +0.049; B2XL + ESM-2, +0.041). Their advantage over a control that merely averages two learners on the same KO features is not statistically separable at this sample size.', { after: 200 }),


      Rich([{ t: 'Figure 1. ', b: true }, { t: 'Scheme of the released model (a, B2) and the best internal model (b, B2XL + Bacformer-large). Each box names the step and the model or tool used for it. Grey: shared steps; blue: gene-content (KO) branch; orange: genome-language-model branch; green: output. Hyper-parameters and cut-offs are chosen inside the training folds only.' }], { after: 100 }),
      Rich([{ t: 'a  B2 — released model (CPU only, about 2 minutes per genome)', b: true }], { after: 60 }),
      scheme([
        rowCentre(['Genome', 'input: RefSeq accession or nucleotide FASTA', 'tool: NCBI Datasets API v2 (download)'], 'shared'),
        rowArrow(),
        rowCentre(['Proteins', 'model: NCBI gene-calling pipeline PGAP (genes already called in RefSeq);', 'for raw DNA: Pyrodigal gene finder (Prodigal model)'], 'shared'),
        rowArrow(),
        rowCentre(['Protein table', 'tool: own script — one protein per gene, in genome order, MD5 per sequence'], 'shared'),
        rowArrow(),
        rowCentre(['KEGG orthologue annotation', 'model: KEGG KOfam profile HMMs (hidden Markov models)', 'searched with pyhmmer (HMMER3), KofamScan per-KO thresholds'], 'ko'),
        rowArrow(),
        rowCentre(['Feature vector', '852-KO energy panel scanned; the 585 KOs that vary across the 358 training', 'genomes are used as features (present = 1 / absent = 0)'], 'ko'),
        rowArrow(),
        rowCentre(['Classifier, one per trait', 'model: L2 logistic regression (scikit-learn) or XGBoost', '(300 trees, depth 2 or 4, 3 seeds); best of 6 chosen in inner CV'], 'ko'),
        rowArrow(),
        rowCentre(['Prediction', 'probability → yes / no at the cut-off maximising MCC; 14 traits', '(a Platt-calibrated probability is also reported)'], 'out'),
      ]),
      P('', { after: 120 }),
      Rich([{ t: 'b  B2XL + Bacformer-large — best internal model (language models run on a GPU)', b: true }], { after: 20 }),
      P('Evaluated in cross-validation and on the external genomes; not packaged in the predictor (src/predict.py computes only the smaller ESM-2 / Bacformer-base embeddings).', { size: 17, color: MUTED, after: 60 }),
      scheme([
        rowCentre(['Genome → Proteins → Protein table', 'RefSeq accessions only: NCBI Datasets download, PGAP genes, genome order', '(no raw-FASTA route: the language-model vectors were computed for these genomes)'], 'shared'),
        rowArrows('↓   gene-content branch (CPU)', '↓   language-model branch (GPU: Kaggle, 2 × T4)'),
        rowPair([['KEGG orthologue annotation', 'model: KOfam profile HMMs', '(pyhmmer, KofamScan thresholds)'], 'ko'],
                [['Protein language model', 'model: ESM-C 300M (ESM++ small)', 'each protein → 960 numbers'], 'lm']),
        rowArrows(),
        rowPair([['Feature vector', '585 varying KOs, 1 / 0'], 'ko'],
                [['Genome language model', 'model: Bacformer-large — 30-layer transformer', 'over proteins in genome order (6 000 per window)'], 'lm']),
        rowArrows(),
        rowPair([['KO classifier', 'model: L2 logistic regression + XGBoost,', 'probabilities averaged'], 'ko'],
                [['Genome vector', 'mean over all proteins → 960 numbers'], 'lm']),
        rowArrows('│', '↓'),
        new TableRow({ cantSplit: true, children: [schemeCell(0, 2, ['│'], null), schemeCell(2, 1, [], null),
          schemeCell(3, 2, ['Language-model classifier', 'model: L2 logistic regression'], 'lm')] }),
        rowArrows('↘', '↙'),
        rowCentre(['Combine and predict', 'probability = ½ KO classifier + ½ language-model classifier', 'yes / no at the MCC cut-off; 14 traits'], 'out'),
      ]),
      P('', { after: 160 }),
      H('2b. Energy-metabolism pathways', HeadingLevel.HEADING_1),
      P('The project\u2019s third deliverable. The supplied labels are organism-level categories, so no pathway-level model can be trained from them; pathway presence is instead read out deterministically from the same KO annotation. All 71 KEGG modules in the energy-metabolism, carbon-fixation, methane, nitrogen and sulfur classes are parsed and scored for completeness per genome (median 6 modules complete per genome), and predict.py emits them alongside the trait predictions. These are annotations rather than predictions, but they track the curated labels closely \u2014 mean completeness in labelled positives versus negatives: dissimilatory sulfate reduction 0.81 vs 0.15, SOX thiosulfate oxidation 0.60 vs 0.09, nitrification 0.52 vs 0.07, methanogenesis (CO2) 0.49 vs 0.06, denitrification 0.46 vs 0.11. They also give resolution the labels lack: denitrification and DNRA are one label category but separate modules, as are the three methanogenesis/acetogenesis routes.', { after: 200 }),

      H('3. How to read these numbers', HeadingLevel.HEADING_1),
      Rich([{ t: 'AUPRC (area under the precision\u2013recall curve) ', b: true }, { t: 'is the primary metric. For each trait the model scores every genome; AUPRC summarises how well positives are ranked above negatives across all thresholds. It runs from the trait\u2019s prevalence (random) to 1.0 (perfect). It is used instead of accuracy or AUROC because the traits are rare: a trait with 3% positives is 97% accurate if the model always says \u201cno\u201d, and AUROC stays optimistically high under that imbalance.' }]),
      bullet('Read AUPRC against the prevalence, not against 1.0. Nitrite oxidation at 0.34 with 3.4% prevalence is a 10\u00d7 lift over chance; H2 oxidation at 0.87 with 58% prevalence is a 1.5\u00d7 lift.'),
      bullet('Macro AUPRC is the unweighted mean over the 14 traits, so rare traits count as much as common ones.'),
      Rich([{ t: 'MCC (Matthews correlation coefficient) ', b: true }, { t: 'describes the actual yes/no calls at the chosen threshold: \u22121 perfectly wrong, 0 chance, 1 perfect. Rough reading: above 0.7 strong, 0.4\u20130.7 useful, 0.2\u20130.4 weak, below 0.2 not usable. Thresholds are always chosen on inner cross-validation folds, never on test data.' }]),
      Rich([{ t: 'The two baselines every claim is measured against. ', b: true }, { t: 'The prevalence null always predicts the trait\u2019s frequency (macro 0.185 \u2014 the floor). The taxonomic null copies the average label of the closest relatives in training; because physiology is often conserved within a lineage, a model is only useful if it beats what taxonomy alone already tells you. At 0.472 (genus) and 0.380 (family) it is a demanding baseline, and the gap to it is this project\u2019s headline quantity.' }]),
      Rich([{ t: 'Confidence intervals ', b: true }, { t: 'are 95% cluster bootstraps: genomes are resampled in whole genera or families, never individually, because genomes of one genus are not independent observations. If an interval for a difference includes zero, the comparison is a tie at this sample size, whatever the point estimate says.' }]),
      Rich([{ t: 'Multiplicity. ', b: true }, { t: '28 model variants were tried, so the report also gives max-T simultaneous intervals that hold jointly across every variant, and Benjamini\u2013Hochberg q-values for per-trait tests. A result with q below 0.05 survives that correction; an unadjusted interval excluding zero with q above 0.05 is suggestive only.' }]),
      Rich([{ t: 'Hold-out levels. ', b: true }, { t: '\u201cHeld-out genus\u201d means no genome of that genus was in training; \u201cheld-out family\u201d means no genome of that family was. Family-level numbers describe a genuinely novel organism and are always lower; the gap measures how much performance is lineage memorisation.' }]),
      Rich([{ t: 'What counts as good here. ', b: true }, { t: 'Reproducing a curated database from genome sequence is bounded by how consistently that database is curated: two independent expert curations of these organisms agree at \u03ba 0.61, and our model agrees with the independent one at \u03ba 0.52. Per-trait AUPRC around 0.85 on marker-defined traits is therefore close to the practical limit, while traits below 0.3 are categories the databases themselves do not define consistently.' }, { t: '', b: false }], { after: 200 }),
      H('4. LLM baselines (arm H)', HeadingLevel.HEADING_1),
      P('An LLM was asked to infer donors and acceptors from the same genome evidence the models use: significant KofamScan hits retrieved from the RDF Portal BH26 endpoint that TogoMCP queries. Two variants were run: anonymised (gene evidence only) and named (evidence plus the organism name). Scored on the same 96 genomes as the models, paired.', { color: MUTED, size: 19 }),
      llmTable(),
      P('Differences are LLM minus model on the same genomes (paired). * 95% genus-cluster bootstrap confidence interval excludes zero. The GPT-6 Astra run with the organism name returned complete answers for 84 of the 96 genomes; its comparisons use those 84.', { color: MUTED, size: 18, after: 140 }),
      bullet(`Reading the same gene evidence, both LLMs clearly beat the taxonomic null (${llmDelta('H_genome', 'null_taxNN')} and ${llmDelta('H_gpt_genome', 'null_taxNN')}) but lose to the trained models — to the released B2 (${llmDelta('H_genome', 'B2')}, ${llmDelta('H_gpt_genome', 'B2')}) and more clearly to the best model, B2XL + Bacformer-large (${llmDelta('H_genome', 'B2XLE_large')}, ${llmDelta('H_gpt_genome', 'B2XLE_large')}; Claude, GPT-6 Astra).`),
      bullet('They struggle where marker genes are ambiguous: nitrite oxidation 0.22–0.28 versus 0.73 for the model, because nitrite oxidase and nitrate reductase share KO identifiers.'),
      bullet(`Given the organism name, Claude overtakes even the best model (${llmDelta('H_named', 'B2XLE_large')}); GPT-6 Astra is ahead by ${llmDelta('H_gpt_named', 'B2XLE_large')}, within noise. That measures recall of published physiology — the same literature the labels come from — not prediction from a genome, and it is unavailable for the uncharacterised organisms this project targets.`),
      P('The published ChatGPT/TogoMCP example also searched BacDive and species descriptions, so it operates in this "named" regime.', { after: 200 }),

      H('5. Per-trait results (released model, held-out genus)', HeadingLevel.HEADING_1),
      table(['Trait', 'Positives', 'AUPRC', 'MCC'], [
        ['Sulfate reduction', '63', '0.887', '0.76'],
        ['Sulfur oxidation', '131', '0.876', '0.70'],
        ['H2 oxidation', '208', '0.873', '0.50'],
        ['Methanogens / acetogens', '75', '0.863', '0.85'],
        ['O2 (aerobic respiration)', '177', '0.828', '0.68'],
        ['Ammonia oxidation', '32', '0.817', '0.74'],
        ['Nitrate reduction', '102', '0.678', '0.48'],
        ['Sulfur reduction', '50', '0.530', '0.44'],
        ['Fe(II) oxidation', '43', '0.482', '0.35'],
        ['CO oxidation', '45', '0.442', '0.35'],
        ['Nitrite oxidation', '12', '0.337', '0.43'],
        ['Fe(III) reduction', '32', '0.214', '0.07'],
        ['"Other" acceptors', '23', '0.136', '0.19'],
        ['"Other" donors', '23', '0.095', '0.08'],
      ], W4),
      P('Predictions for the 228 unannotated genomes are released only for the traits whose calls show out-of-clade skill; the rest are reported as scores without a call.', { after: 200 }),

      H('6. Evidence that the results are real', HeadingLevel.HEADING_1),
      bullet('External validation: 562 genomes from genera never seen in training, independently curated (Madin et al. 2020). Macro AUPRC 0.511 versus 0.380 for the taxonomic null (+0.131, CI +0.100 to +0.167). Largest margins: sulfate reduction 0.76 vs 0.38, sulfur oxidation 0.40 vs 0.13, methanogenesis 0.83 vs 0.63.'),
      bullet('Genuine novelty: maximum protein-sequence overlap between any external and any training genome is 0.009.'),
      bullet('Leakage quantified: a random (non-phylogenetic) split would have reported 0.647 instead of 0.576 — the inflation typical of papers that do not block by lineage.'),
      bullet('Robustness: conclusions hold under five re-drawn fold allocations, GTDB taxonomy, restriction to complete genomes, and treating empty label blocks as unknown.'),
      bullet('Pre-registration written before any model was scored, with its six predictions scored afterwards: two supported, two partly, one mostly, and one (“gene content beats every embedding arm”) not supported.'),
      bullet('All 28 model variants tried are logged with multiplicity-controlled inference; nothing was silently dropped.'),

      H('7. The ceiling: what limits performance', HeadingLevel.HEADING_1),
      P('Two independent expert curations of the same 179 species were compared directly, like for like:', { after: 100 }),
      table(['Comparison', 'Agreement (κ)'], [
        ['The two curations against each other', '0.61'],
        ['Our model against the independent curation', '0.52'],
        ['Same, excluding the H2 definitional mismatch', '0.68 (curations)'],
      ], W2),
      P('The model reproduces an independent curation almost as well as two expert curations reproduce each other. Trait by trait, curation agreement tracks model performance: methanogens κ 0.82 (model AUPRC 0.90), nitrate 0.73 (0.79), sulfate 0.61, Fe(III) reduction 0.47 (0.19). For H2 oxidation the two databases disagree structurally — ours calls 125 of 179 species positive, the other 33 — a difference of definition rather than noise.'),
      bullet('Within-genus label heterogeneity mirrors this: every genus containing an iron reducer or an “other” category also contains members labelled negative, versus 4 of 18 for methanogens.'),
      bullet('Direct label contradictions: 5 genomes carry the complete methanogenesis gene set but are labelled non-methanogen; 9 obligate aerobes are labelled oxygen-negative; 248 confident out-of-fold "false positives" are shared by two independent models, 131 carrying a complete marker complex.'),
      bullet('More data helps slowly: doubling the number of labelled genera adds about +0.02 macro AUPRC, so reaching 0.65 would need on the order of a thousand labelled genera.'),

      H('8. Publication assessment', HeadingLevel.HEADING_1),
      Rich([{ t: 'Recommendation: ', b: true }, { t: 'publishable, framed as a benchmark and honest-evaluation paper rather than as a new predictor. ' },
        { t: 'Suitable venues: Bioinformatics, mSystems, Microbiome, Nucleic Acids Research (Methods), ISME Communications. A top-tier machine-learning venue is not a good fit: the models are standard and the contribution lies in the evaluation, the baselines and the ceiling analysis.' }]),
      P('What carries the paper:', { bold: true, after: 80 }),
      bullet('Phylogeny-blocked evaluation with frozen, hashed splits and two nulls — most genome-to-phenotype papers report random-split numbers, which we show are inflated by ~0.07.'),
      bullet('External validation on 562 novel-genus genomes with independent labels.'),
      bullet('An LLM/TogoMCP baseline, including the leakage demonstration when the organism is named.'),
      bullet('A quantitative ceiling from cross-curation agreement, which reframes "can this be predicted" as "which curation are you predicting".'),
      bullet('Pre-registration with scored predictions, a complete variant ledger, and four independent audit rounds.'),
      P('What must not be claimed:', { bold: true, after: 80 }),
      bullet('That protein language models solve this task. Alone they tie gene content; in combination their contribution cannot be separated from simple learner ensembling at N = 358.'),
      bullet('That the model predicts metabolism. It predicts database annotations of metabolism, which is what the ceiling analysis quantifies.'),
      bullet('Anything about iron metabolism or the "other" categories beyond "not predictable from these labels".'),
      Rich([{ t: 'Deployment as part of TogoMCP. ', b: true }, { t: 'The released model needs only KO annotations \u2014 exactly what the RDF Portal BH26 endpoint already serves, and our KofamScan reproduction matches its significant hits to 1 cell in 43 344. Any of the ~23 400 genomes in that endpoint can therefore be scored with no download and no GPU, and any other genome from an accession or a raw FASTA. Exposed as predict / explain / list_traits, each answer would carry the KO evidence behind the call, that trait\u2019s measured accuracy on unseen genera and families, and flags for whether the genome, genus or family was in training, with calls returned only for the traits with demonstrated out-of-clade skill. This is the deliverable that helps where the literature-based route fails \u2014 organisms nobody has characterised \u2014 and it is roughly one to two days of work, since the models, thresholds and calibration are already serialised.' }]),
      P('Open items before submission:', { bold: true, after: 80 }),
      bullet('Integrity audit (in progress) — an adversarial check for leakage, cherry-picking and metric gaming.'),
      bullet('Attention-MIL rerun without the annotation-based prefilter, for a clean attention-versus-pooling comparison.'),
      bullet('External validation of the two best combined models, awaiting GPU features for the large encoder.'),


      H('9. Methods and models: completed checks, and what to try next', HeadingLevel.HEADING_1),
      P('Completed since the previous revision \u2014 both settled, nothing is pending', { bold: true, after: 80 }),
      table(['Experiment', 'What it settled', 'Outcome'], [
        ['Arm D_nopf \u2014 attention MIL with the KO bag-prefilter forced off', 'protocol Q3: attention vs mean pooling, with no curated annotation entering the arm', 'answered "no": 0.548 vs 0.544, \u0394 +0.004 [\u22120.024, +0.045]'],
        ['External validation of the two best internal models (ESM-C / Bacformer-large features for the 562 novel-genus genomes)', 'whether the best internal model (0.625) transfers to independently curated labels', 'it does not: \u22120.012 on all traits, \u22120.080 on the no-FAPROTAX subset. B2 stays released and "internal gains do not transfer" now covers every variant'],
      ], [3200, 3400, 2760]),
      P('Methods and models worth trying next', { bold: true, after: 60 }),
      P('Probability = chance of a measurable, honestly evaluated macro-AUPRC gain on the held-out-genus benchmark, grounded in this project\u2019s own measurements: 28 variants plateau within \u00b10.05, the learning curve gives \u2248+0.02 per doubling of training genera, and two expert curations agree only at \u03ba 0.61.', { color: MUTED, size: 19 }),
      table(['Method / model', 'Expected \u0394', 'Probability', 'Difficulty'], [
        ['Protein-level retrieval (arm G): ESM-C similarity to clade-resolved seed sets + rules', '+0.01\u20130.03 macro; +0.1\u20130.3 on ambiguous traits', '60%', 'Medium-high'],
        ['Trait-specific curated HMM panels (As/Se/DMSO/TMAO, Cyc2/Mtr/Omc)', '+0.00\u20130.02 macro; large per trait', '60% per trait', 'Medium'],
        ['Gene-neighbourhood / operon features (nar+narK, dsr+apr+qmo, mtr cluster)', '+0.01 to +0.03', '45%', 'Medium'],
        ['Semi-supervised use of unlabelled genomes (pseudo-labels, contrastive pretraining)', '+0.01 to +0.03', '40%', 'Medium'],
        ['Phylogeny-aware models (random effects, clade-equalising weights)', 'up to +0.03 (family)', '40%', 'Medium'],
        ['Rare-trait architectures (prototype / few-shot, shared-trunk multi-task)', '+0.00\u20130.02 macro', '35%', 'Medium'],
        ['Hierarchical / constrained multi-label decoding', '+0.00 to +0.02', '30%', 'Medium'],
        ['Full 7 900-KO profile instead of the 852-KO panel (tested as B3)', '+0.00 to +0.02', '30%', 'Low'],
        ['Label-noise-robust / PU learning', '~0 measured', '25%', 'Medium'],
        ['Fine-tuning Bacformer-large (LP-FT / LoRA)', '+0.00 to +0.02', '20%', 'High (GPU days)'],
        ['Stacking more encoders together', '+0.00 to +0.01', '20%', 'Low'],
        ['Conformal / selective prediction (deployment, not accuracy)', 'usable precision, not AUPRC', '80%', 'Low'],
      ], W4),
      Rich([{ t: 'Data-side lanes, for comparison. ', b: true }, { t: 'Two non-modelling lanes dominate the table above: more labelled genera (+0.04 to +0.07 for ~1 000\u20131 500 genera, 85% likely, months of curation) and harmonising labels across databases (+0.03 to +0.08 measured, 70% likely). The model search is in its flat region, so effort is better spent on retrieval, curated HMMs, selective prediction and the label lanes than on bigger encoders or more ensembling.' }], { after: 200 }),
      H('10. Deliverables', HeadingLevel.HEADING_1),
      bullet('Trained models per arm, with calibration, in results/models/; the released model is documented in RELEASED.json.'),
      bullet('Metric tables: per trait, macro, per fold, paired comparisons, simultaneous inference, sensitivity analyses, research ledger, external validation, LLM baselines, calibration, learning curve.'),
      bullet('Predictions for the 228 unannotated genomes, with reliability flags and both raw and calibrated scores.'),
      bullet('A command-line tool that takes a genome accession and returns predictions, tested end to end.'),
      bullet('Documentation: pre-registration, deviation log, research-loop record, external-validation design, LLM-baseline protocol, and a results summary in which every number cites its source file.'),
    ],
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync('./docs/inorganic_prediction_report.docx', b);
  console.log('written');
});
