const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat,
} = require('docx');
const fs = require('fs');

const ACCENT = '2A78D6';
const MUTED = '52514E';
const HEAD_BG = 'EDF2FA';
const ALT_BG = 'F7F8FA';

const S = 17;           // 8.5 pt body, keeps the summary on one page
const P = (text, o = {}) => new Paragraph({
  spacing: { after: o.after ?? 60, line: 230 },
  children: [new TextRun({ text, size: o.size ?? S, bold: o.bold, color: o.color, font: 'Calibri' })],
});
const Rich = (runs, o = {}) => new Paragraph({
  spacing: { after: o.after ?? 60, line: 230 },
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, size: o.size ?? S, color: r.c, font: 'Calibri' })),
});
const H = (text) => new Paragraph({
  spacing: { before: 120, after: 60 },
  children: [new TextRun({ text, bold: true, size: 21, color: ACCENT, font: 'Calibri' })],
});
const bullet = (text) => new Paragraph({
  numbering: { reference: 'b', level: 0 },
  spacing: { after: 30, line: 230 },
  children: [new TextRun({ text, size: S, font: 'Calibri' })],
});

const W = [3960, 1500, 1500, 2400];
const row = (cells, { bg, bold } = {}) => new TableRow({
  children: cells.map((t, i) => new TableCell({
    width: { size: W[i], type: WidthType.DXA },
    shading: bg ? { type: ShadingType.CLEAR, fill: bg, color: 'auto' } : undefined,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    children: [new Paragraph({
      alignment: i === 0 || i === 3 ? AlignmentType.LEFT : AlignmentType.RIGHT,
      spacing: { after: 0, line: 220 },
      children: [new TextRun({ text: String(t), bold, size: 17, font: 'Calibri' })],
    })],
  })),
});

const doc = new Document({
  creator: 'bio_pred',
  title: 'Electron donor/acceptor prediction — one-page summary',
  numbering: { config: [{ reference: 'b', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•',
    alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 300, hanging: 180 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: S } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 },
      margin: { top: 900, bottom: 900, left: 1080, right: 1080 } } },
    children: [
      new Paragraph({
        spacing: { after: 40 },
        children: [new TextRun({ text: 'Predicting electron donors and acceptors from genomes — one-page summary',
          bold: true, size: 26, color: ACCENT, font: 'Calibri' })],
      }),
      new Paragraph({
        spacing: { after: 140 },
        border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: ACCENT } },
        children: [new TextRun({ text: '586 genomes (358 labelled) · 14 traits · evaluation holds out whole genera and whole families · 18 September 2026',
          size: 16, color: MUTED, font: 'Calibri' })],
      }),

      H('Methods in brief'),
      Rich([{ t: 'Best model (B2XL + Bacformer-large). ', b: true }, { t: 'Two views of each genome, averaged. View 1: presence/absence of 852 energy-metabolism KEGG orthologues (KofamScan) through a logistic regression and a gradient-boosted tree. View 2: a logistic probe on a 960-dimensional genome vector from Bacformer-large, a genome language model that reads all proteins in genome order via ESM-C 300M. One model per trait; hyper-parameters and thresholds chosen inside training folds only.' }]),
      Rich([{ t: 'Released model (B2). ', b: true }, { t: 'The same KO view alone — no GPU, no language model. It is released because no variant beat it on an independent external test set (see below).' }]),
      Rich([{ t: 'Baselines. ', b: true }, { t: '(i) prevalence \u2014 always predict the trait frequency; (ii) expert rules \u2014 call a trait if its marker complex is complete; (iii) taxonomic null \u2014 copy the average label of the closest relatives in training, the demanding baseline since physiology is lineage-conserved.' }]),
      Rich([{ t: 'LLM baseline (arm H). ', b: true }, { t: 'Claude and GPT-6 Astra got exactly the evidence the model gets \u2014 a genome\u2019s significant KofamScan hits from the RDF Portal BH26 endpoint that TogoMCP queries \u2014 and returned a probability per trait, anonymised or with the organism name. No tools, no web; 96 genomes, paired with the models.' }], { after: 120 }),

      H('Results'),
      new Table({
        columnWidths: W,
        width: { size: W.reduce((a, b) => a + b, 0), type: WidthType.DXA },
        borders: {
          top: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
          bottom: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
          left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
          insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: 'E4E7EC' },
          insideVertical: { style: BorderStyle.NONE },
        },
        rows: [
          row(['Model / baseline', 'New genus', 'New family', 'Note'], { bg: HEAD_BG, bold: true }),
          row(['Prevalence (floor)', '0.185', '0.178', 'guessing']),
          row(['Expert marker-gene rules', '0.317', '0.317', 'no learning']),
          row(['Taxonomic null (copy relatives)', '0.472', '0.380', 'the baseline to beat'], { bg: ALT_BG }),
          row(['LLM on gene evidence (Claude / GPT)', '0.578 / 0.588', '—', '96 genomes, anonymised']),
          row(['B2 — released KO model', '0.576', '0.560', 'ships as the predictor'], { bg: ALT_BG }),
          row(['B2XL + Bacformer-large — best', '0.625', '0.585', 'not confirmed externally']),
          row(['LLM given the organism name', '0.731 / 0.777', '—', 'literature recall, not prediction'], { bg: ALT_BG }),
        ],
      }),
      P('Macro AUPRC over 14 traits (0.185 = guessing, 1.0 = perfect). LLM rows are on a 96-genome subset where the released model scores 0.666.', { size: 15, color: MUTED, after: 90 }),
      bullet('Six traits work out of clade (sulfate reduction 0.89, sulfur oxidation 0.88, H2 0.87, methanogens 0.86, O2 0.83, ammonia oxidation 0.82); three are moderate (0.44\u20130.68) and five poor (iron traits and the two "other" categories).'),
      bullet('It generalises: on 562 genomes from genera never seen in training, independently labelled, 0.511 versus 0.380 for the taxonomic null (+0.131, CI +0.100 to +0.167).'),
      bullet('Versus the LLM baselines: reading the same gene evidence, both LLMs beat the taxonomic null (+0.09 to +0.10) but lose to the trained model by 0.08–0.11 (CIs exclude zero for three of the four comparisons). They fail hardest where genes are ambiguous, e.g. nitrite oxidation 0.22–0.28 versus 0.73, because nitrite oxidase and nitrate reductase share KO identifiers. Given the organism name they overtake the model — but that is recalling published physiology, which is where the labels come from, and is unavailable for uncharacterised organisms.'),
      bullet('Protein language models do not replace gene content: internally they tie it (0.577 vs 0.576) and fused they add +0.03 to +0.05, but on the external set they lose (Bacformer-large \u22120.071, ESM-C \u22120.116 vs the KO model) and no fusion transferred \u2014 which is why the simple KO model ships.'),

      H('Is this good?'),
      Rich([{ t: 'Yes for the well-defined half of the traits, and the limit is the labels rather than the model. ', b: true }, { t: 'Two independent expert curations of the same organisms agree at κ 0.61; our model agrees with the independent curation at κ 0.52. So on marker-defined traits the model is already close to the level at which experts reproduce each other, and the traits it fails on are those the databases themselves define inconsistently — for iron reduction, every genus containing a positive also contains members labelled negative.' }]),
      Rich([{ t: 'What this is worth publishing for: ', b: true }, { t: 'an honest out-of-clade benchmark (a random split would have reported 0.647 instead of 0.576), external validation on novel genera, an LLM/TogoMCP baseline including its leakage when the organism is named, and a quantified ceiling. Not supported: that genome language models solve this task, or that the model predicts metabolism rather than database annotations of it.' }]),
      Rich([{ t: 'To go further, ', b: true }, { t: 'more labelled genera help slowly (+0.02 per doubling) and better label definitions help more; curated gene models fix specific failures, as one did for nitrite oxidation (0.34 \u2192 0.64).' }]),
      Rich([{ t: 'Good enough to serve inside TogoMCP. ', b: true }, { t: 'The released model needs only KO annotations \u2014 what the RDF Portal BH26 endpoint already holds, and our KofamScan reproduction matches its significant hits to 1 cell in 43 344 \u2014 so any of the ~23 400 genomes there can be scored with no download and no GPU, and anything else from an accession or a raw FASTA. Exposed as predict / explain / list_traits, each answer would carry the KO evidence behind the call, the measured accuracy of that trait on unseen genera and families, and flags for whether the genome, genus or family was in training; calls are returned only for the 8 traits with demonstrated out-of-clade skill. That helps exactly where the literature route fails \u2014 organisms nobody has characterised yet.' }]),
    ],
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync('./docs/inorganic_prediction_one_page_report.docx', b);
  console.log('written');
});
