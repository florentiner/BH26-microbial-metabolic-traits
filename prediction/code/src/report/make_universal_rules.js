// Universal rules: gene rules that hold across lineages, with biological interpretation.
// Every number is read from results/rulebook/universal_data.json (src/57_universal_data.py).
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer,
} = require('docx');
const fs = require('fs');

const ROOT = '.';
const D = JSON.parse(fs.readFileSync(`${ROOT}/results/rulebook/universal_data.json`, 'utf8'));
const ACCENT = '2A78D6', MUTED = '52514E', HEAD_BG = 'EDF2FA', ALT_BG = 'F7F8FA';
const GOOD = 'CFEBD6', CARE = 'FFF1CC', GREY = 'EEEEEE';
const f = (x, n = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–' : Number(x).toFixed(n);
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
  const mk = (cells, { bg, bold, fills } = {}) => new TableRow({ cantSplit: true, children: cells.map((t, i) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    shading: (fills && fills[i]) || bg ? { type: ShadingType.CLEAR, fill: (fills && fills[i]) || bg, color: 'auto' } : undefined,
    margins: { top: 45, bottom: 45, left: 75, right: 75 },
    children: [new Paragraph({ alignment: i === 0 || (o.left || []).includes(i) ? AlignmentType.LEFT : AlignmentType.RIGHT,
      spacing: { after: 0, line: 235 }, children: [new TextRun({ text: String(t), bold, size: 17, font: 'Calibri' })] })],
  })) });
  return new Table({ columnWidths: widths, width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: { top: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' }, bottom: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: 'E4E7EC' }, insideVertical: { style: BorderStyle.NONE } },
    rows: [mk(headers, { bg: HEAD_BG, bold: true }), ...rows.map((r, k) => mk(r.cells || r, { bg: k % 2 ? ALT_BG : undefined, fills: r.fills }))] });
}
const card = (badge, fill, color, rule) => new Table({ columnWidths: [1900, 7460], width: { size: 9360, type: WidthType.DXA },
  borders: { top: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' }, bottom: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' },
    left: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' }, right: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' },
    insideHorizontal: { style: BorderStyle.NONE }, insideVertical: { style: BorderStyle.SINGLE, size: 2, color: 'B7C7DE' } },
  rows: [new TableRow({ cantSplit: true, children: [
    new TableCell({ width: { size: 1900, type: WidthType.DXA }, shading: { type: ShadingType.CLEAR, fill, color: 'auto' },
      margins: { top: 80, bottom: 80, left: 100, right: 100 },
      children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ text: badge, bold: true, size: 20, color, font: 'Calibri' })] })] }),
    new TableCell({ width: { size: 7460, type: WidthType.DXA }, shading: { type: ShadingType.CLEAR, fill: 'FAFBFD', color: 'auto' },
      margins: { top: 80, bottom: 80, left: 140, right: 120 },
      children: rule.replace(/^IF /, '').split(/ THEN /).flatMap((part, i) => i === 0
        ? part.split(' AND ').map((c, j) => new Paragraph({ spacing: { after: 10 }, children: [
          new TextRun({ text: j === 0 ? 'IF   ' : 'AND  ', bold: true, size: 19, font: 'Consolas', color: ACCENT }),
          new TextRun({ text: c, size: 19, font: 'Consolas' })] }))
        : [new Paragraph({ spacing: { after: 0 }, children: [new TextRun({ text: 'THEN ', bold: true, size: 19, font: 'Consolas', color: ACCENT }),
          new TextRun({ text: nice(part), bold: true, size: 19, font: 'Consolas' })] })]) }),
  ] })] });

const tasks = Object.entries(D.tasks).filter(([, b]) => b.status === 'done');
const allRules = tasks.flatMap(([k, b]) => b.rules.map(r => ({ ...r, task: k, taskName: b.name })));
const uniPos = allRules.filter(r => r.universal && r.cls === 'has trait');
const direct = uniPos.filter(r => !r.indirect), indirect = uniPos.filter(r => r.indirect);
const uniNeg = allRules.filter(r => r.universal && r.cls === 'lacks trait');
const LC = tasks.flatMap(([k, b]) => (b.label_completeness || []).map(r => ({ ...r, task: k })));

const C = [];
const NEXT = D.n_external_novel ?? 562;
C.push(new Paragraph({ spacing: { after: 60 }, children: [new TextRun({ text: 'Universal gene rules', bold: true, size: 40, font: 'Calibri', color: '1A1A1A' })] }));
C.push(P('Rules that are not specific to one organism, family or phylum: built only from genes whose function connects them to the trait, and tested on lineages they never saw', { size: 23, color: MUTED, after: 60 }));
C.push(P('Generated from results/rulebook/ by src/report/make_universal_rules.js. Method: docs/preregistration_trees_xai.md, Addendum J (pre-registered), K (post-hoc revision) and L (audit corrections), disclosed in section 6.', { size: 18, color: MUTED, after: 200 }));

// ------------------------------------------------------------------ 1. results
C.push(H('1. Results', HeadingLevel.HEADING_1));
C.push(P(`Rules were searched for ${tasks.map(([, b]) => `${b.n_traits} ${b.name.toLowerCase()} trait${b.n_traits > 1 ? 's' : ''}`).join(', ')}, using only genes with a known function for each trait. A rule counts as universal when (i) it is right at least 80 % of the time, pooled over the lineages left out while it was learned (each test genome counted once), with a lower bound of at least 60 % that beats the base rate; (ii) in all data it is at least 80 % correct in at least three lineages (two when the trait occurs in only two); (iii) the same enzyme-level rule is learned again in at least half of the leave-one-out runs; (iv) its genes connect to the trait in the direction the rule uses them; and (v) where an independent external test exists, it beats the external base rate.`));
C.push(P(`${direct.length} direct universal rule${direct.length === 1 ? '' : 's'} ${direct.length === 1 ? 'was' : 'were'} found, each built around pathway genes of the trait (some with an added lifestyle condition): ${direct.map(r => `${nice(r.trait)} (${genesOf(r.rule)}: ${pct(r.unseen_precision)} correct in unseen lineages, ${r.unseen_hits} distinct test genomes in the ${r.unseen_units_fired ?? '–'} held-out lineages where it fired${r.large_units_never_tested ? `; never tested unseen in ${r.large_units_never_tested}` : ''})`).join('; ')}.${indirect.length ? ` ${indirect.length} further rule${indirect.length > 1 ? 's are' : ' is'} statistically universal but indirect (section 3).` : ''}`));
const nif = LC.find(r => r.trait === 'nitrogen_fixation'), den = LC.find(r => r.trait === 'denitrification');
if (nif) C.push(P(`For several traits the biology gives a universal rule but the labels cannot confirm it. The complete nitrogenase (NifHDK) is present in ${nif.carriers} genomes across ${nif.n_phyla_with_marker} phyla, yet only ${pct(nif.share_of_carriers_labelled)} of them are labelled nitrogen fixers${den ? `; for NO + N₂O reductase and denitrification it is ${pct(den.share_of_carriers_labelled)}` : ''}. Traits that are rarely tested are recorded as absent, so the measured precision of such rules is capped by the labels, not by the biology (section 4).`));
const mcr = direct.find(r => r.trait === 'methanogenesis'), mcrL = LC.find(r => r.trait === 'methanogenesis');
if (mcr && mcrL) C.push(P(`Methyl-CoM reductase is the cleanest case: ${pct(mcrL.share_of_carriers_labelled)} of its ${mcrL.carriers} carriers are labelled methanogens and ${pct(mcrL.share_of_trait_with_marker)} of methanogens carry it; the rule was relearned in ${pct(mcr.lopo_stability_enzyme)} of leave-one-lineage-out runs and was right ${pct(mcr.unseen_precision)} of the time in the lineages it had not seen.`));

// ------------------------------------------------------------------ 2. method
C.push(H('2. What “universal” means here, and how it was tested', HeadingLevel.HEADING_1));
C.push(bullet('Only meaningful genes: each trait’s own pathway enzymes, plus markers of aerobic or anaerobic life, from a curated table (src/gene_mechanisms.py). Caveat: some lifestyle markers in that table (korA/C/D, hcp, the PFL-activating enzyme, aceE, coxL, cooS) were added after they appeared in earlier analyses of overlapping genomes, so the candidate set is not fully independent of the data (Addendum L).'));
C.push(bullet('Leave one lineage out: for every lineage with at least five organisms carrying the trait, rules were learned on all other lineages and tested on the one left out. The lineage is the phylum when the trait occurs in three or more phyla, otherwise the class or order.'));
C.push(bullet('The rule is judged by its worst lineage, not its average: among rules that are at least 80 % correct in every lineage where they fire on at least three genomes (measured in the training data), the one that catches the most organisms is chosen; if no eligible rule reaches 80 % everywhere, the one with the best worst-lineage precision is taken.'));
C.push(bullet('Rules are compared at the enzyme level: subunits of one enzyme (for example narG and narH) count as the same rule.'));
C.push(spacer());

// ------------------------------------------------------------------ 3. rules
C.push(H('3. The universal rules', HeadingLevel.HEADING_1));
const showRule = (r, badge, fill, color) => {
  C.push(P(`${nice(r.trait)} — ${r.taskName}`, { bold: true, size: 22, color: ACCENT, after: 60, keepNext: true }));
  C.push(card(badge, fill, color, r.rule));
  C.push(Rich([{ t: 'In lineages it never saw: ', b: true }, { t: `${pct(r.unseen_precision)} correct [${pct(r.unseen_ci_low)}–${pct(r.unseen_ci_high)}], ${r.unseen_hits} distinct test genomes. It was relearned in ${r.unseen_units_tested ?? '–'} leave-one-out runs but fired in only ${r.unseen_units_fired ?? '–'} held-out ${r.lineage_unit === 'class' ? 'classes' : r.lineage_unit === 'order' ? 'orders' : 'phyla'}: ${r.unseen_by_unit || '–'} (at least 80 % correct in ${r.unseen_units_precision_ge_80 ?? '–'} with ≥ 3 hits); the same enzyme rule was learned again in ${pct(r.lopo_stability_enzyme)} of the runs.${r.large_units_never_tested ? ` Never tested on unseen genomes of the large lineages ${r.large_units_never_tested}, so its generality there rests on in-sample data only.` : ''}` }], { after: 30 }));
  C.push(Rich([{ t: 'Across lineages, in all data: ', b: true }, { t: `at least ${pct(r.worst_unit_precision)} correct in every ${r.lineage_unit} where it fires on at least three genomes (${r.n_units} ${r.lineage_unit === 'class' ? 'classes' : r.lineage_unit === 'order' ? 'orders' : 'phyla'}): ${r.per_unit}.` }], { after: 30 }));
  C.push(Rich([{ t: 'Reach: ', b: true }, { t: `catches ${pct(r.coverage)} of all organisms with the trait (base rate ${pct(r.base_rate)}).${r.external_precision !== undefined && r.external_precision !== null ? ` On ${NEXT} independent external genomes from new genera: ${pct(r.external_precision)} correct (n = ${r.external_support}; external base rate ${pct(r.external_base_rate)}).` : ''}` }], { after: 40 }));
  C.push(P('Why it makes sense:', { bold: true, size: 20, after: 20, keepNext: true }));
  r.conditions.split(' | ').forEach(c => C.push(bullet(c)));
  if (r.why_general) C.push(Rich([{ t: r.indirect ? 'Caution: ' : 'Why it could hold in every lineage (interpretation, from textbook physiology): ', b: true, c: r.indirect ? '8A5A00' : '1E6B34' }, { t: r.why_general }], { after: 160 }));
  else C.push(spacer());
};
for (const r of direct) showRule(r, 'Universal', GOOD, '1E6B34');
if (indirect.length) {
  C.push(H('Statistically universal, but indirect', HeadingLevel.HEADING_2));
  C.push(P('These rules hold across lineages, but they predict the trait from the absence of its canonical enzyme together with lifestyle markers — they recognise the kind of organism that does the job with other enzymes, rather than the enzyme doing it.'));
  for (const r of indirect) showRule(r, 'Universal (indirect)', CARE, '8A5A00');
}

// ------------------------------------------------------------------ 4. limited by labels
C.push(H('4. Universal in biology, but not confirmable with these labels', HeadingLevel.HEADING_1));
C.push(P('For each trait’s core enzyme: how many genomes carry the complete enzyme, how many of those carriers are labelled with the trait, and how many labelled organisms carry it. When almost every carrier is labelled (methyl-CoM reductase), a universal rule can be confirmed. When most carriers are unlabelled, a rule built on the enzyme looks imprecise even if the biology is sound — the labels are incomplete, not the rule.'));
C.push(table(['Trait', 'Core enzyme', 'Carriers', 'Carriers labelled', 'Labelled organisms with it', 'Phyla'],
  LC.map(r => ({ cells: [nice(r.trait), r.core_marker, r.carriers, `${r.carriers_labelled} (${pct(r.share_of_carriers_labelled)})`, pct(r.share_of_trait_with_marker), r.n_phyla_with_marker],
    fills: [undefined, undefined, undefined, r.share_of_carriers_labelled >= 0.8 ? GOOD : CARE, undefined, undefined] })),
  [2000, 2600, 900, 1500, 1500, 860], { left: [1] }));
C.push(spacer());
for (const r of LC.filter(x => x.why)) C.push(Rich([{ t: `${nice(r.trait)}: `, b: true }, { t: r.why }], { after: 60 }));

// ------------------------------------------------------------------ 5. rule-out rules
C.push(H('5. Universal rules that rule a trait out', HeadingLevel.HEADING_1));
C.push(P('From the post-hoc run (Addendum K). Rules that say an organism lacks a trait and hold across lineages. Most traits are rare, so check each against its base rate — a rule is only useful when it clearly beats always answering “no”.'));
C.push(table(['Trait', 'Rule (genes)', 'Unseen lineages', 'Base rate', 'Lift', 'Lineages fired', 'Catches'],
  uniNeg.map(r => ({ cells: [nice(r.trait), genesOf(r.rule), pct(r.unseen_precision), pct(r.base_rate), r.base_rate ? f(r.unseen_precision / r.base_rate) : '–', `${r.unseen_units_fired ?? '–'}${(r.unseen_units_fired ?? 0) < 2 ? ' ⚠' : ''}`, pct(r.coverage)] })),
  [2000, 3000, 1000, 800, 700, 900, 960], { left: [1] }));
C.push(P('Lift = unseen precision ÷ base rate; for common “lacks” classes it is close to 1 even for a correct rule. ⚠ = tested unseen in fewer than two lineages.', { size: 17, color: MUTED }));
C.push(spacer());

// ------------------------------------------------------------------ 6. transparency
C.push(H('6. How the analysis changed — disclosed', HeadingLevel.HEADING_1));
const preRec = tasks.map(([, b]) => b.preregistered_recorded_n_neg || 0).reduce((a, b) => a + b, 0);
const preNow = tasks.map(([, b]) => b.preregistered_rerun_n_neg || 0).reduce((a, b) => a + b, 0);
C.push(P(`As first specified (pre-registered, Addendum J), the search found no universal rule that predicts a trait — only ${preRec} rules that rule one out, as recorded in Addendum K (re-running that specification with the corrected grounding of Addendum L leaves ${preNow}). Pre-registered prediction J-P1 — ${D.j_p1 ? D.j_p1.prediction : ''} — is therefore not supported as specified; after the post-hoc revision, ${D.j_p1 ? D.j_p1.posthoc : ''}. The output showed three structural reasons, and the method was revised after seeing them (Addendum K):`));
C.push(bullet('Some traits occur in fewer than three phyla (every methanogen is in one NCBI phylum), so no rule could span three phyla. The lineage unit now adapts: phylum, or class, or order.'));
C.push(bullet('Subunits of one enzyme counted as different rules; rules are now compared at the enzyme level.'));
C.push(bullet('A trait found in exactly three phyla left only two after holding one out, so its rule could never be relearned; the requirement now adapts to the lineages that remain. And ranking by worst-lineage precision alone picked tiny rules; among rules at least 80 % correct in every lineage, the broadest is now chosen.'));
C.push(bullet('Two thresholds were relaxed in the same revision and are disclosed here: the lineage requirement can fall to two lineages when the trait occurs in only two, and “learned again in most runs” became “in at least half of the runs”.'));
C.push(P('Later audit rounds (Addendum L) also changed which genes are candidates: sulfur oxidation, iron oxidation and methylotrophy no longer carry an oxygen-lifestyle expectation, so both aerobic and anaerobic markers became candidates for them. A universal claim now also needs unseen evidence in at least two lineages where the rule fired and at least 20 distinct unseen genomes, and, where an informative external test exists, an external lower confidence bound of at least 60 %.', { after: 80 }));
C.push(P('After an independent audit (Addendum L), three corrections were made: hits of rules for the same enzyme learned in the same run are merged, so each test genome counts once (earlier versions double-counted); a rule the independent external test contradicts is no longer called universal (the external test counts only when it has at least five genomes of each class; otherwise it is reported as uninformative); and the external genomes are restricted to genera absent from training.', { after: 80 }));
C.push(P('These revisions were made because of a result, so the rules in section 3 are post hoc. They are still tested only on lineages left out while learning, which is what makes them credible.', { size: 19, color: MUTED }));

// ------------------------------------------------------------------ 7. limits
C.push(H('7. Limits', HeadingLevel.HEADING_1));
C.push(bullet('“All bacteria” here means all prokaryotes in the data. Methanogenesis is an archaeal trait, so its universal rule spans archaeal classes, not bacterial phyla.'));
C.push(bullet('Lineages follow NCBI taxonomy; GTDB would split some phyla (e.g. Proteobacteria) and could change which rules count as spanning three lineages.'));
C.push(bullet('The donor/acceptor data (358 genomes) have few lineages per trait, so few of their rules could be tested on unseen lineages at all.'));
C.push(bullet('Labels are closed-world: an untested trait counts as absent, which lowers the measured precision of rules for under-reported traits (section 4).'));

const doc = new Document({ creator: 'bio_pred', title: 'Universal gene rules',
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 360, hanging: 200 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [new TextRun({
      children: ['Page ', PageNumber.CURRENT, ' of ', PageNumber.TOTAL_PAGES], size: 18, color: MUTED, font: 'Calibri' })] })] }) },
    children: C }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(`${ROOT}/docs/universal_gene_rules.docx`, b); console.log('written docs/universal_gene_rules.docx'); });
