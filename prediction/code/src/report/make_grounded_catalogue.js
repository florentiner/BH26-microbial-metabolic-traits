// Grounded rule catalogue: every gene rule judged on (1) how often it holds and (2) whether it makes biological sense.
// Every number is read from results/rulebook/grounded_catalogue.json (src/55_grounded_catalogue_data.py).
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer,
} = require('docx');
const fs = require('fs');

const ROOT = '.';
const D = JSON.parse(fs.readFileSync(`${ROOT}/results/rulebook/grounded_catalogue.json`, 'utf8'));
const ACCENT = '2A78D6', MUTED = '52514E', HEAD_BG = 'EDF2FA', ALT_BG = 'F7F8FA';
const VCOL = { 'Use': 'CFEBD6', 'Use with care': 'FFF1CC', 'Sound, not reliable enough': 'E3EEF8', 'Do not rely on it': 'FBE0D6', 'Reject': 'EEEEEE' };
const VTXT = { 'Use': '1E6B34', 'Use with care': '8A5A00', 'Sound, not reliable enough': '1F5B9A', 'Do not rely on it': 'A33A1F', 'Reject': '555555' };
const GCOL = { 'grounded': 'CFEBD6', 'partly grounded': 'FFF1CC', 'not grounded': 'FBE0D6' };
const ORDER = ['Use', 'Use with care', 'Sound, not reliable enough', 'Do not rely on it', 'Reject'];

const pct = x => (x === null || x === undefined || Number.isNaN(x)) ? '–' : `${Math.round(100 * x)} %`;
const nice = t => t.replace(/^donor_/, 'donor: ').replace(/^acceptor_/, 'acceptor: ').replace('uses_O2', 'uses O₂').replace(/_/g, ' ');
const genesOf = rule => rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, '').replace(/ AND /g, ' + ');
const P = (text, o = {}) => new Paragraph({ spacing: { after: o.after ?? 110, line: 270 }, keepNext: o.keepNext,
  children: [new TextRun({ text, bold: o.bold, italics: o.italics, size: o.size ?? 21, color: o.color, font: 'Calibri' })] });
const Rich = (runs, o = {}) => new Paragraph({ spacing: { after: o.after ?? 110, line: 270 }, keepNext: o.keepNext,
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i, size: r.s ?? 20, color: r.c, font: 'Calibri' })) });
const H = (text, level) => new Paragraph({ heading: level, keepNext: true,
  spacing: { before: level === HeadingLevel.HEADING_1 ? 340 : 220, after: 120 },
  children: [new TextRun({ text, font: 'Calibri', bold: true, size: level === HeadingLevel.HEADING_1 ? 30 : 24,
    color: level === HeadingLevel.HEADING_1 ? ACCENT : '1A1A1A' })] });
const bullet = (text, o = {}) => new Paragraph({ numbering: { reference: 'bul', level: 0 }, spacing: { after: 50, line: 255 },
  children: [new TextRun({ text, size: o.size ?? 19, font: 'Calibri' })] });
const spacer = () => new Paragraph({ spacing: { after: 80 }, children: [] });

function table(headers, rows, widths, o = {}) {
  const mk = (cells, { bg, bold, fills, colors } = {}) => new TableRow({ cantSplit: true, children: cells.map((t, i) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    shading: (fills && fills[i]) || bg ? { type: ShadingType.CLEAR, fill: (fills && fills[i]) || bg, color: 'auto' } : undefined,
    margins: { top: 45, bottom: 45, left: 75, right: 75 },
    children: [new Paragraph({ alignment: i === 0 || (o.left || []).includes(i) ? AlignmentType.LEFT : AlignmentType.RIGHT,
      spacing: { after: 0, line: 235 }, children: [new TextRun({ text: String(t), bold: bold || (colors && colors[i] ? true : false), color: colors && colors[i], size: 17, font: 'Calibri' })] })],
  })) });
  return new Table({ columnWidths: widths, width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: { top: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' }, bottom: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: 'E4E7EC' }, insideVertical: { style: BorderStyle.NONE } },
    rows: [mk(headers, { bg: HEAD_BG, bold: true }), ...rows.map((r, k) => mk(r.cells || r, { bg: k % 2 ? ALT_BG : undefined, fills: r.fills, colors: r.colors }))] });
}
// a rule card header: verdict badge + the rule, in one shaded box
const card = (r) => new Table({ columnWidths: [1900, 7460], width: { size: 9360, type: WidthType.DXA },
  borders: { top: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' }, bottom: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' },
    left: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' }, right: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' },
    insideHorizontal: { style: BorderStyle.NONE }, insideVertical: { style: BorderStyle.SINGLE, size: 2, color: 'B7C7DE' } },
  rows: [new TableRow({ cantSplit: true, children: [
    new TableCell({ width: { size: 1900, type: WidthType.DXA }, shading: { type: ShadingType.CLEAR, fill: VCOL[r.verdict], color: 'auto' },
      margins: { top: 80, bottom: 80, left: 100, right: 100 },
      children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: r.verdict, bold: true, size: 20, color: VTXT[r.verdict], font: 'Calibri' })] })] }),
    new TableCell({ width: { size: 7460, type: WidthType.DXA }, shading: { type: ShadingType.CLEAR, fill: 'FAFBFD', color: 'auto' },
      margins: { top: 80, bottom: 80, left: 140, right: 120 },
      children: r.rule.replace(/^IF /, '').split(/ THEN /).flatMap((part, i) => i === 0
        ? part.split(' AND ').map((c, j) => new Paragraph({ spacing: { after: 10 }, children: [
          new TextRun({ text: j === 0 ? 'IF   ' : 'AND  ', bold: true, size: 19, font: 'Consolas', color: ACCENT }),
          new TextRun({ text: c, size: 19, font: 'Consolas' })] }))
        : [new Paragraph({ spacing: { after: 0 }, children: [new TextRun({ text: 'THEN ', bold: true, size: 19, font: 'Consolas', color: ACCENT }),
          new TextRun({ text: nice(part), bold: true, size: 19, font: 'Consolas' })] })]) }),
  ] })] });

const tasks = Object.entries(D.tasks).filter(([, b]) => b.status === 'done');
const NEXT = D.n_external_novel ?? 562;
const pending = Object.entries(D.tasks).filter(([, b]) => b.status !== 'done');
const C = [];
C.push(new Paragraph({ spacing: { after: 60 }, children: [new TextRun({ text: 'Grounded rule catalogue', bold: true, size: 40, font: 'Calibri', color: '1A1A1A' })] }));
C.push(P('Gene rules for microbial physiology, each judged on two questions: how often does it hold, and does it make biological sense?', { size: 23, color: MUTED, after: 60 }));
C.push(P(`Generated from results/rulebook/ by src/report/make_grounded_catalogue.js. Included: ${tasks.map(([, b]) => b.name).join(', ')}${pending.length ? `. Still training: ${pending.map(([, b]) => b.name).join(', ')}` : ''}.`, { size: 18, color: MUTED, after: 80 }));
C.push(P('Exploratory, defined after the rules were built: the grounding axis, the five-verdict matrix, the breadth comparison and the “decisive” criterion for forest rules were not pre-registered (docs/preregistration_trees_xai.md, Addendum L). The rules and their probabilities come from the pre-registered procedure (Addendum G); the verdicts are a reading aid, not a test.', { size: 18, color: MUTED, after: 200 }));

// ------------------------------------------------------------------ 1. summary
C.push(H('1. Summary', HeadingLevel.HEADING_1));
const ov = D.overall || {};
const total = Object.values(ov).reduce((a, b) => a + b, 0);
C.push(P(`${total} rules were judged. ${ov['Use'] || 0} can be used as they are — they hold often and every gene in them has a known role that connects it to the prediction. ${ov['Use with care'] || 0} hold often but rest partly on genes with no known role. ${ov['Do not rely on it'] || 0} are accurate on these data but have no biological reason behind them, so they most likely track lineage and are expected to fail on new organisms. ${ov['Sound, not reliable enough'] || 0} make biological sense but did not hold often enough, and ${ov['Reject'] || 0} are neither.`));
C.push(table(['Verdict', ...tasks.map(([, b]) => b.name), 'All'],
  ORDER.map(v => ({ cells: [v, ...tasks.map(([, b]) => b.counts[v] || 0), ov[v] || 0], fills: [VCOL[v]], colors: [VTXT[v]] })),
  [2800, ...tasks.map(() => Math.floor(5360 / tasks.length)), 1200].map((w, i, a) => i === a.length - 2 ? w + (5360 - Math.floor(5360 / tasks.length) * tasks.length) : w)));
C.push(spacer());
const tg = D.transfer_by_grounding || [];
if (tg.length) {
  const g = n => tg.find(r => r.grounding === n);
  const gr = g('grounded'), pg = g('partly grounded'), ng = g('not grounded');
  C.push(P(`Does grounding predict transfer? On the external genomes (exploratory; small groups): ${[['grounded', gr], ['partly grounded', pg], ['not grounded', ng]].filter(([, x]) => x).map(([n, x]) => `${n} ${pct(x.heldout)} held-out → ${pct(x.external)} external (${x.n} rules)`).join('; ')}. Every group lost precision on new organisms, and the groups are too small to show that grounded rules transfer better. Grounding is a reason to prefer a rule, not evidence that it holds on new lineages.`));
}
C.push(P('The rules to use:', { bold: true, after: 60 }));
const useRows = tasks.flatMap(([, b]) => b.rules.filter(r => r.verdict === 'Use'));
C.push(table(['Trait', 'Rule (genes)', 'Holds', 'Base rate', 'Fires on', 'Catches', 'Why it makes sense'],
  useRows.map(r => [nice(r.trait) + (r.cls === 'lacks trait' ? ' (rule-out)' : ''), genesOf(r.rule), pct(r.probability), pct(r.base_rate), pct(r.fires_on_share), pct(r.coverage),
    r.grounding_reasons.split('; ').slice(0, 2).join('; ')]),
  [1650, 2150, 650, 700, 700, 700, 2810], { left: [1, 6] }));
C.push(P('Holds = share of genomes the rule fires on where it is right; compare it with the base rate of the same genomes — the share you would get right by always guessing the same answer (high for rule-out rules, because most traits are rare). Fires on = share of the training genomes the rule applies to. Catches = share of training genomes with (or, for rule-out rules, without) the trait that the rule picks up.', { size: 17, color: MUTED }));

// ------------------------------------------------------------------ 2. how rules are judged
C.push(H('2. How each rule is judged', HeadingLevel.HEADING_1));
C.push(P('Question 1 — how often does it hold? (the numbers)', { bold: true, after: 50 }));
C.push(bullet(`Holds: of the genomes the rule fires on, the share where its prediction is right — measured on ${NEXT} external genomes from genera absent from training where available, otherwise on genera held out while the rule was learned. Held-out genera are closer to the training data than new lineages, so held-out values are usually higher than external ones (rulebook, Fig. 2 and Supplementary Note 5). Shown with its 95 % interval.`));
C.push(bullet('Fires on: the share of the training genomes the rule applies to — how often the pattern exists at all.'));
C.push(bullet('Catches: the share of training genomes with the trait (or without it, for rules that rule it out) that the rule picks up.'));
C.push(bullet('Trusted: holds at least 80 % of the time, the interval’s lower end is at least 60 %, and it clearly beats guessing from the base rate.'));
C.push(P('Question 2 — does it make biological sense? (the mechanism)', { bold: true, after: 50 }));
C.push(bullet('Each gene condition is checked against known biology, in the direction the rule uses it: a gene of the trait’s own pathway, present in a rule that predicts the trait or absent in one that rules it out; an absence that argues against this particular trait, in a rule that rules it out; or a marker of aerobic or anaerobic life whose signal fits the trait’s oxygen lifestyle. A pathway gene required in the opposite state, or a lifestyle signal at odds with the trait, counts against grounding.'));
C.push(bullet('Grounded: every condition connects to the trait this way. Partly grounded: at least one does, and at least one has no known role, contradicts the biology, or is uninformative. Not grounded: none does — the rule is a pure statistical pattern.'));
C.push(bullet('Neutral (counts for neither side): an absent aerobic marker, because most organisms lack it anyway; and any oxygen-lifestyle marker for a trait found in both aerobes and anaerobes (e.g. H₂ oxidation, fermentation).'));
C.push(P('The verdict combines both:', { bold: true, after: 50 }));
C.push(table(['', 'Grounded', 'Partly grounded', 'Not grounded'], [
  { cells: ['Trusted', 'Use', 'Use with care', 'Do not rely on it'], fills: [undefined, VCOL['Use'], VCOL['Use with care'], VCOL['Do not rely on it']], colors: [undefined, VTXT['Use'], VTXT['Use with care'], VTXT['Do not rely on it']] },
  { cells: ['Not trusted', 'Sound, not reliable enough', 'Reject', 'Reject'], fills: [undefined, VCOL['Sound, not reliable enough'], VCOL['Reject'], VCOL['Reject']], colors: [undefined, VTXT['Sound, not reliable enough'], VTXT['Reject'], VTXT['Reject']] },
], [1800, 2520, 2520, 2520]));
C.push(spacer());

// ------------------------------------------------------------------ 3. catalogue
C.push(H('3. The catalogue, trait by trait', HeadingLevel.HEADING_1));
C.push(P('Rules that can be used (Use, Use with care) and rules that make sense but are not reliable yet are shown in full; rules that should not be relied on or are rejected are listed briefly at the end of each trait.'));
for (const [k, b] of tasks) {
  C.push(H(b.name, HeadingLevel.HEADING_2));
  const traits = [...new Set(b.rules.map(r => r.trait))];
  for (const t of traits) {
    const rs = b.rules.filter(r => r.trait === t);
    C.push(P(nice(t), { bold: true, size: 22, color: ACCENT, after: 60, keepNext: true }));
    const full = rs.filter(r => ['Use', 'Use with care', 'Sound, not reliable enough'].includes(r.verdict));
    const brief = rs.filter(r => !['Use', 'Use with care', 'Sound, not reliable enough'].includes(r.verdict));
    if (!full.length) C.push(P('No rule for this trait both holds and makes sense.', { italics: true, color: MUTED, size: 19 }));
    for (const r of full) {
      C.push(card(r));
      C.push(Rich([{ t: 'How often it holds: ', b: true }, { t: `${pct(r.probability)} [${pct(r.probability_ci_low)}–${pct(r.probability_ci_high)}] on ${r.probability_source} (base rate ${pct(r.base_rate)}${r.probability_source === 'external genomes' && r.external_n_pos !== null && r.external_n_pos !== undefined ? `; ${r.external_n_pos} external genomes with the trait${r.external_n_pos < 5 ? ' — too few for a reliable external test' : ''}` : ''}). Fires on ${pct(r.fires_on_share)} of training genomes (${r.fires_on_n}); catches ${pct(r.coverage)} of training genomes ${r.cls === 'has trait' ? 'with' : 'without'} the trait.` }], { after: 40 }));
      C.push(Rich([{ t: `Does it make biological sense? ${r.grounding.charAt(0).toUpperCase() + r.grounding.slice(1)}. `, b: true, c: r.grounding === 'grounded' ? '1E6B34' : (r.grounding === 'partly grounded' ? '8A5A00' : 'A33A1F') }, { t: r.logic }], { after: 30 }));
      r.conditions.split(' | ').forEach(c => C.push(bullet(c)));
      C.push(Rich([{ t: 'Organisms: ', b: true }, { t: `${r.captures || '–'}.${r.misfires_training || r.misfires_external ? ` Misfires: ${[r.misfires_training, r.misfires_external ? `external — ${r.misfires_external}` : ''].filter(x => x).join('; ')}.` : ''}` }], { after: 150 }));
    }
    if (brief.length) {
      C.push(table(['Verdict', 'Rule (genes)', 'Holds', 'Why not'],
        brief.map(r => ({ cells: [r.verdict, genesOf(r.rule), pct(r.probability),
          r.verdict === 'Do not rely on it' ? `accurate here, but ${r.grounding_reasons.split('; ').filter(x => x.includes('no known role') || x.includes('required') || x.includes('at odds') || x.includes('says little') || x.includes('weak evidence')).join('; ') || 'no gene has a known role'} — likely lineage` :
          (r.trusted ? r.grounding_reasons : `${r.probability >= 0.8 && r.probability_ci_low >= 0.6 ? `holds ${pct(r.probability)} but does not clearly beat the base rate (${pct(r.base_rate)})` : r.probability >= 0.8 ? `holds ${pct(r.probability)} but the lower confidence bound (${pct(r.probability_ci_low)}) is below 60 %` : `holds only ${pct(r.probability)}`}${r.grounding !== 'grounded' ? '; ' + (r.grounding_reasons.split('; ').filter(x => x.includes('no known role') || x.includes('required') || x.includes('at odds') || x.includes('says little') || x.includes('weak evidence'))[0] || '') : ''}`)],
          fills: [VCOL[r.verdict]], colors: [VTXT[r.verdict]] })),
        [1600, 2700, 700, 4360], { left: [1, 3] }));
      C.push(spacer());
    }
  }
}

// ------------------------------------------------------------------ 4. forest questions
C.push(H('4. The random forest’s own questions, judged the same way', HeadingLevel.HEADING_1));
C.push(P('The single-gene questions the random forest asks first most often (rulebook, Supplementary Table S3). A question that raises the probability of the trait is judged as a rule that predicts it; one that lowers it as a rule that rules it out. “Trusted” here (a post-hoc criterion) means it is decisive — at least 80 % one way — and clearly beats the base rate of the same genomes; P(trait) is external where available, otherwise held-out.'));
for (const [k, b] of tasks) {
  const fr = (b.forest || []).filter(r => r.verdict !== 'Reject');
  if (!fr.length) continue;
  C.push(P(b.name, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Question', 'P(trait)', 'Base rate', 'Sense', 'Verdict'],
    fr.map(r => ({ cells: [nice(r.trait), r.rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, ''), pct(r.p_trait), pct(r.base_rate), r.grounding, r.verdict],
      fills: [undefined, undefined, undefined, undefined, GCOL[r.grounding], VCOL[r.verdict]], colors: [undefined, undefined, undefined, undefined, undefined, VTXT[r.verdict]] })),
    [1900, 2500, 900, 900, 1400, 1760], { left: [1, 4, 5] }));
  C.push(spacer());
}

// ------------------------------------------------------------------ 5. limits
C.push(H('5. Limits', HeadingLevel.HEADING_1));
C.push(bullet('“Makes biological sense” is judged against a curated table of gene functions (src/gene_mechanisms.py). A gene missing from the table counts as having no known role, so a rule can be more grounded than the table shows — never less.'));
C.push(bullet('The external test exists for donor/acceptor rules only; oxygen-use rules are measured on held-out genera.'));
C.push(bullet('Most reliable rules draw their hits mainly from one phylum and describe that lineage’s version of the physiology; use the set of rules for a trait, not one rule alone (rulebook, Supplementary Note 3).'));
C.push(bullet('Labels are closed-world: a trait not recorded for an organism counts as absent, so some “wrong” calls are missing labels.'));

const doc = new Document({ creator: 'bio_pred', title: 'Grounded rule catalogue',
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 360, hanging: 200 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [new TextRun({
      children: ['Page ', PageNumber.CURRENT, ' of ', PageNumber.TOTAL_PAGES], size: 18, color: MUTED, font: 'Calibri' })] })] }) },
    children: C }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(`${ROOT}/docs/grounded_rule_catalogue.docx`, b); console.log('written docs/grounded_rule_catalogue.docx'); });
