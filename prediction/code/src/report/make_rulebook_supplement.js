// Rulebook: high-probability gene rules with explanations. Every number is read from
// results/rulebook/rulebook_data.json (built by src/52_rulebook_report_data.py).
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer,
} = require('docx');
const fs = require('fs');

const ROOT = '.';
const D = JSON.parse(fs.readFileSync(`${ROOT}/results/rulebook/rulebook_data.json`, 'utf8'));
const ACCENT = '2A78D6', MUTED = '52514E', HEAD_BG = 'EDF2FA', ALT_BG = 'F7F8FA';
const GOOD = 'DCEFE0', MID = 'FFF3D6', BAD = 'FBE3D6', GREY = 'EEEEEE';
const LEVEL_BG = { mechanistic: GOOD, physiological: GOOD, 'contains unexplained conditions': MID };
const VBG = { supported: GOOD, 'not supported': BAD, pending: GREY, 'not run': GREY };

const f = (x, n = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–' : Number(x).toFixed(n);
const pct = x => (x === null || x === undefined) ? '–' : `${Math.round(100 * x)} %`;
const nice = t => t.replace(/^donor_/, 'donor: ').replace(/^acceptor_/, 'acceptor: ').replace(/_/g, ' ');
const P = (text, o = {}) => new Paragraph({ spacing: { after: o.after ?? 110, line: 270 }, keepNext: o.keepNext,
  children: [new TextRun({ text, bold: o.bold, italics: o.italics, size: o.size ?? 21, color: o.color, font: o.font ?? 'Calibri' })] });
const Rich = (runs, o = {}) => new Paragraph({ spacing: { after: o.after ?? 110, line: 270 }, keepNext: o.keepNext,
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i, size: r.s ?? 21, color: r.c, font: r.f ?? 'Calibri' })) });
const H = (text, level) => new Paragraph({ heading: level, keepNext: true,
  spacing: { before: level === HeadingLevel.HEADING_1 ? 340 : 220, after: 120 },
  children: [new TextRun({ text, font: 'Calibri', bold: true, size: level === HeadingLevel.HEADING_1 ? 30 : 24,
    color: level === HeadingLevel.HEADING_1 ? ACCENT : '1A1A1A' })] });
const bullet = (text, o = {}) => new Paragraph({ numbering: { reference: 'bul', level: 0 }, spacing: { after: 60, line: 260 },
  children: [new TextRun({ text, size: o.size ?? 20, font: 'Calibri', color: o.color })] });
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
// the rule itself, shown as IF / AND / THEN lines in a shaded box
const ruleBox = (rule, fill) => new Table({ columnWidths: [9360], width: { size: 9360, type: WidthType.DXA },
  borders: { top: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' }, bottom: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' },
    left: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' }, right: { style: BorderStyle.SINGLE, size: 4, color: 'B7C7DE' },
    insideHorizontal: { style: BorderStyle.NONE }, insideVertical: { style: BorderStyle.NONE } },
  rows: [new TableRow({ cantSplit: true, children: [new TableCell({ width: { size: 9360, type: WidthType.DXA },
    shading: { type: ShadingType.CLEAR, fill, color: 'auto' }, margins: { top: 80, bottom: 80, left: 140, right: 140 },
    children: rule.replace(/^IF /, '').split(/ THEN /).flatMap((part, i) => i === 0
      ? part.split(' AND ').map((c, j) => new Paragraph({ spacing: { after: 20 }, children: [
        new TextRun({ text: j === 0 ? 'IF   ' : 'AND  ', bold: true, size: 20, font: 'Consolas', color: ACCENT }),
        new TextRun({ text: c, size: 20, font: 'Consolas' })] }))
      : [new Paragraph({ spacing: { after: 0 }, children: [new TextRun({ text: 'THEN ', bold: true, size: 20, font: 'Consolas', color: ACCENT }),
        new TextRun({ text: nice(part.replace(/^(has|lacks) /, m => m)), bold: true, size: 20, font: 'Consolas' })] })]) })] })] });

const TITLE = require('./rulebook_title.js');
const C = [];
C.push(new Paragraph({ spacing: { after: 60 }, children: [new TextRun({ text: 'Supplementary Information',
  bold: true, size: 38, font: 'Calibri', color: '1A1A1A' })] }));
C.push(P(`for: ${TITLE}`, { size: 23, color: MUTED, after: 60 }));
const doneTasks = Object.entries(D.tasks).filter(([, b]) => b.status === 'done');
const pendTasks = Object.entries(D.tasks).filter(([, b]) => b.status !== 'done');
C.push(P('Contents: Supplementary Notes 1–5 and Tables S1–S5. Every number is read from results/rulebook/ (src/report/make_rulebook_supplement.js).', { size: 18, color: MUTED, after: 200 }));

// ------------------------------------------------------------------ helpers
const hasK = (r, ko) => !!r && r.rule.includes(`(${ko})`);
const NEXT = D.n_external_novel ?? 562;
const baseFor = r => r.probability_source === 'external genomes' ? r.external_base_rate : r.base_rate;
const task = k => D.tasks[k] && D.tasks[k].status === 'done' ? D.tasks[k] : null;
const DA = task('da'), OX = task('ox');
const reliable = b => b ? b.rules.filter(r => r.cls === 'has trait' && r.high_probability_and_informative) : [];
const bestRel = (b, trait) => reliable(b).filter(r => r.trait === trait).sort((a, c) => c.coverage - a.coverage)[0];
const markerMcc = (b, trait) => { const r = (b && b.rules_vs_best || []).find(x => x.trait === trait && x.rule_set.startsWith('marker')); return r ? r.rule_mcc : null; };
const modelMcc = (b, trait) => { const r = (b && b.rules_vs_best || []).find(x => x.trait === trait); return r ? r.best_mcc : null; };
const fRule = (b, trait, gene, present) => (b && b.forest_rules || []).find(r => r.trait === trait && r.genes === gene && r.rule.includes(`) ${present ? 'present' : 'absent'}`));
const pExt = r => r ? (r.p_external ?? r.p_heldout) : null;
const genesOf = r => r ? r.rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, '').replace(/ AND /g, ' + ') : '–';
const relTraits = b => [...new Set(reliable(b).map(r => r.trait))];
const physRel = (b, trait, cls = 'has trait') => (b ? b.rules : []).filter(r => r.trait === trait && r.cls === cls && r.high_probability_and_informative && r.interpretation_level !== 'contains unexplained conditions').sort((a, c) => c.coverage - a.coverage)[0];


// ------------------------------------------------------------------ 4. how to read
C.push(H('Supplementary Note 1. How to read a rule', HeadingLevel.HEADING_1));
C.push(P('Each rule names up to three genes that must be present or absent. The numbers under it say how far to trust it:'));
C.push(table(['Number', 'Meaning'], [
  ['Probability', `Of the genomes the rule fires on, the share that truly have (or lack) the trait. Measured on the ${NEXT} lineage-novel external genomes where available, otherwise on genera held out while the rule was learned. Held-out genera are closer to the training data than truly new lineages (Supplementary Note 5), so held-out values can be higher than external ones.`],
  ['95 % interval', 'Wilson interval for that probability; wide when the rule fires on few genomes.'],
  ['Base rate', 'The probability you would get by always guessing the same answer, on the same genomes as the probability next to it (external base rate next to an external probability, the training-set base rate next to a held-out one). A rule is only useful if it clearly beats this (“informative”).'],
  ['Coverage', 'Share of all training genomes with the trait (or without it, for “lacks” rules) that the rule catches — measured in-sample.'],
  ['Stability', 'In how many of the 5 cross-validation folds the identical rule was learned again from scratch.'],
  ['Trust level', 'Mechanistic: every informative condition is a pathway gene in the state that fits the rule (present for “has”, absent for “lacks”). Physiological: those plus aerobic/anaerobic lifestyle markers whose signal fits the trait. Contains lineage proxies: at least one gene has no known link, or is in the state that contradicts the pathway — likely to fail on new lineages. Post hoc (Addendum H; direction rule added in Addendum L).'],
], [1700, 7660], { left: [1] }));
C.push(spacer());


// ------------------------------------------------------------------ 2. biological interpretation
C.push(H('Supplementary Note 2. Extended biological interpretation', HeadingLevel.HEADING_1));
C.push(P('Interpretation, not a result: this section reads the rules in the light of established physiology. The numbers come from the analysis; the explanations are hypotheses from the literature on each pathway and were not tested here. Where a gene’s computed trust label (Supplementary Note 4) says “no known link”, the text says so too.', { size: 19, color: MUTED }));
const bio = (title, paras) => { C.push(H(title, HeadingLevel.HEADING_2)); paras.filter(x => x).forEach(x => C.push(P(x))); };
if (DA) {
  const s = bestRel(DA, 'acceptor_sulfate'), dk = fRule(DA, 'acceptor_sulfate', 'K27188', true);
  bio('Sulfate respiration', [
    `Sulfate reducers activate sulfate with ATP sulfurylase (sat), reduce it to sulfite with APS reductase (aprAB) and sulfite to sulfide with dissimilatory sulfite reductase (dsrAB); the membrane complex DsrMKJOP returns electrons from the quinone pool to the DsrC carrier. None of these genes is specific on its own: sat also serves sulfate assimilation, and dsrAB runs in reverse in sulfur-oxidising bacteria (reverse-Dsr).`,
    s ? `The reliable rule (${genesOf(s)}; ${pct(s.probability)} on ${s.probability_source}, trust level: ${s.interpretation_level})${hasK(s, 'K00176') || hasK(s, 'K00177') || hasK(s, 'K00174') ? ' combines the respiratory core with a subunit of the O₂-sensitive 2-oxoglutarate:ferredoxin oxidoreductase (kor), which is common in anaerobes (the four-subunit enzyme also occurs in some microaerophiles). A plausible reading is that requiring the pathway together with an anaerobic marker excludes aerobic sulfur oxidisers that carry dsrAB' : ''}. dsrB on its own reaches MCC ${f(markerMcc(DA, 'acceptor_sulfate'))} (the forest: ${f(modelMcc(DA, 'acceptor_sulfate'))}).` : '',
    dk ? `The forest’s favourite question, dsrK of the DsrMKJOP complex, raises the probability from ${pct(dk.external_base_rate ?? dk.base_rate)} to ${pct(pExt(dk))}.` : '']);
  const o = bestRel(DA, 'acceptor_oxygen');
  bio('Oxygen as terminal electron acceptor', [
    o ? `Aerobic respiration ends at a terminal oxidase that reduces O₂ to water; electrons usually reach it through the cytochrome bc1 complex and cytochrome c. The reliable rule is ${genesOf(o)} (${pct(o.probability)} on ${o.probability_source}; trust level: ${o.interpretation_level}).${hasK(o, 'K00413') && hasK(o, 'K05601') ? ' Read biologically: a respiratory chain able to feed c-type oxidases (cytochrome c1, petC), and no hybrid-cluster protein (hcp), an enzyme typical of anaerobic and facultative metabolism.' : ''}` : '',
    'The single oxidase genes perform worse (aa3 or cbb3 as a marker rule: MCC ' + f(markerMcc(DA, 'acceptor_oxygen')) + ', against ' + f(modelMcc(DA, 'acceptor_oxygen')) + ' for the forest). Two known reasons: aerobic respiration can use several oxidase families, and cytochrome bd — itself an O₂ reductase — is kept by many strict anaerobes to scavenge traces of O₂.']);
}
if (OX) {
  const ph = physRel(OX, 'uses_O2'), pa = physRel(OX, 'uses_O2', 'lacks trait'), broad = bestRel(OX, 'uses_O2');
  const hb = (OX.forest_rules || []).find(x => x.genes === 'K06886' && x.direction === 'raises');
  const oxLevels = [...new Set(OX.rules.filter(r => r.cls === 'has trait').map(r => r.interpretation_level))];
  bio('Oxygen use (aerobic versus anaerobic lifestyle)', [
    `Aerobes oxidise pyruvate with pyruvate dehydrogenase (aceE); many anaerobes use pyruvate formate-lyase, whose radical enzyme is destroyed by O₂. ${pa ? `The rule-out rule ${genesOf(pa)} identifies anaerobes with ${pct(pa.probability)} precision on held-out genera (trust level: ${pa.interpretation_level}).` : ''} ${ph ? `The rule ${genesOf(ph)} identifies O₂ users with ${pct(ph.probability)} (trust level: ${ph.interpretation_level}).` : `None of the rules that predict O₂ use is interpretable as physiology by the direction-aware check (trust levels: ${oxLevels.join(', ')}): each contains a gene with no curated link to O₂, so they may partly track lineage.`}`,
    broad ? `The rule that catches the most O₂ users is ${genesOf(broad)} (${pct(broad.probability)}, catching ${pct(broad.coverage)} of O₂ users in the training genomes; trust level: ${broad.interpretation_level}).` : '',
    hb ? `The forest’s most frequent first question is a truncated haemoglobin (K06886), an O₂-binding protein; carrying it raises the probability of O₂ use to ${pct(hb.p_heldout)}. O₂ use is written into many genes at once, so models can pick any of several correlated genes.` : '']);
}
if (DA) {
  const n = bestRel(DA, 'acceptor_nitrate');
  bio('Nitrate respiration', [
    n ? `The reliable rule is ${genesOf(n)} (${pct(n.probability)} on ${n.probability_source}; trust level: ${n.interpretation_level}). narI anchors the respiratory nitrate reductase NarGH in the membrane and connects it to the quinone pool.${hasK(n, 'K02276') || hasK(n, 'K02274') ? ' The rule also requires an aa3 oxidase subunit; one reading is that it picks organisms that can respire both O₂ and nitrate (facultative anaerobes), but the computed label treats this gene as unlinked to nitrate respiration, and the association may partly reflect lineage.' : ''}` : '',
    n ? `The rule catches only ${pct(n.coverage)} of nitrate respirers in the training genomes: others use the periplasmic Nap system or reduce nitrate only to nitrite.` : '']);
  const c = bestRel(DA, 'acceptor_CO2_methanogen_acetogen'), cdp = fRule(DA, 'acceptor_CO2_methanogen_acetogen', 'K00194', true);
  bio('CO₂ as acceptor: methanogens and acetogens', [
    c ? `The reliable rule is ${genesOf(c)} (${pct(c.probability)} on ${c.probability_source}; trust level: ${c.interpretation_level}). The CODH/ACS complex (cdhD) of the Wood–Ljungdahl pathway, which assembles acetyl-CoA, is shared by acetogens and methanogens.${hasK(c, 'K11261') ? ' K11261 is annotated as FwdE of formylmethanofuran dehydrogenase; in methanogens it belongs to the CO₂-to-methane pathway, but the bacteria this rule catches (acetogens) carry FwdE-like homologs with no methanofuran pathway, so for them it is a marker, not the mechanism.' : ''}${hasK(c, 'K00338') ? ' The rule also requires nuoI (a complex I subunit) to be absent; many of these organisms conserve energy with Ech, Rnf or Mtr complexes rather than complex I, but the computed label treats this condition as a lineage marker.' : ''}` : '',
    cdp ? `cdhD alone is the forest’s most frequent question: present → ${pct(pExt(cdp))}, absent → ${pct(pExt(fRule(DA, 'acceptor_CO2_methanogen_acetogen', 'K00194', false)))}.` : '']);
  const a = fRule(DA, 'donor_ammonia_oxidation', 'K10945', true);
  bio('Ammonia oxidation', [
    a ? `Ammonia monooxygenase (amoB) is the committed enzyme of ammonia oxidation; the forest rule “amoB present” is right ${pct(a.p_heldout)} of the time on held-out genera, and ${pct(a.p_external)} on the external genomes. Known cause of the drop: particulate methane monooxygenase (pmoB) has the same KEGG identifier, so methanotrophs look like ammonia oxidisers. The external set records this trait as “nitrification”, a broader label, and external precision is scored against it. Hydroxylamine oxidoreductase (hao) is a poor marker: HAO-like octaheme cytochromes also occur in anaerobic sulfate reducers.` : '']);
  const sx = fRule(DA, 'donor_sulfur_oxidation', 'K17222', true);
  bio('Sulfur oxidation', [
    sx ? `The Sox system (SoxAX, SoxYZ, SoxB) and sulfide:quinone oxidoreductase (sqr) are the rules’ genes. soxA present raises the probability of sulfur oxidation from ${pct(sx.base_rate)} to ${pct(sx.p_heldout)} on held-out genera, but only to ${pct(sx.p_external)} externally (base rate ${pct(sx.external_base_rate)}). Sulfur oxidisers also use other routes (reverse-Dsr, the tetrathionate pathway of acidophiles), and some external organisms known to oxidise sulfur may not be labelled so.` : '']);
  bio('Traits without a reliable rule — and why', [
    'H₂ oxidation: hydrogenases form a large, diverse family ([NiFe] groups 1–4, [FeFe]), and the KEGG identifiers cover only some of them; almost half of the H₂ oxidisers missed by the marker list are methanogens that use group-3 and group-4 hydrogenases. The models fall back on lineage markers (fumarate hydratase, a pyrophosphatase), which fail on new organisms.',
    'Iron oxidation and iron reduction: KEGG now has orthologs for MtrA/D, MtrB/E, MtrF (K28380, K28381, K28383) and a Cyc2-type iron:cytochrome-c reductase (K28382), but they are not in the 852-KO energy panel these rules were built from, and OmcS and many other multiheme cytochromes have no specific KO. The rules could not see the key iron proteins.',
    'Nitrite oxidation: nitrite oxidoreductase (nxrAB) shares KEGG identifiers with respiratory nitrate reductase (narGH), so no KO can tell the two apart.',
    'CO oxidation and the “other” categories: few positives, and the “other” categories group unrelated physiologies.']);
  bio('What the rules say overall', [
    'Rules whose genes are on the trait’s pathway are easiest to explain, and single pathway genes are most useful for ruling a trait out. But grounding in biology did not protect a rule from losing precision on new organisms (Supplementary Note 5): grounded rules dropped too. Biological grounding is therefore a reason to prefer a rule, not proof that it will transfer; the external test is the evidence.']);
}


// ------------------------------------------------------------------ 3. trusted rules, one by one
C.push(H('Supplementary Note 3. Trusted rules, one by one', HeadingLevel.HEADING_1));
C.push(P('Every rule that is reliable (at least 80 % correct, lower bound at least 60 %, and better than guessing from the base rate), read biologically: what each condition means, why the combination implies the trait, which organisms it actually picks out, and where it goes wrong. The genera are real examples from the data (numbers of genomes in brackets).'));
C.push(P('How broad is a rule? Its hits are compared with all organisms that have the trait. A rule is “as broad as the trait” when its phylum mix matches the trait’s; “narrower” when it draws mainly on one phylum — then it describes that lineage’s version of the physiology, and other rules for the same trait are needed to cover the rest.', { size: 19, color: MUTED }));
for (const [k, b] of doneTasks) {
  const tr = b.trusted || [];
  if (!tr.length) continue;
  C.push(H(b.name, HeadingLevel.HEADING_2));
  const pos = tr.filter(r => r.cls === 'has trait');
  const neg = tr.filter(r => r.cls === 'lacks trait');
  for (const r of pos) {
    C.push(P(nice(r.trait), { bold: true, size: 22, color: ACCENT, after: 60, keepNext: true }));
    C.push(ruleBox(r.rule, 'EEF6F0'));
    C.push(Rich([{ t: 'Probability ', b: true, s: 19 }, { t: `${pct(r.probability)}; covers ${pct(r.coverage)} of organisms with the trait; trust: ${r.interpretation_level}.`, s: 19 }], { after: 40 }));
    C.push(P('What each condition means:', { bold: true, size: 19, after: 20, keepNext: true }));
    r.conditions.split(' | ').forEach(c => C.push(bullet(c, { size: 19 })));
    C.push(Rich([{ t: 'Why the combination implies the trait: ', b: true, s: 19 }, { t: r.logic, s: 19 }], { after: 40 }));
    C.push(Rich([{ t: 'Organisms it picks out: ', b: true, s: 19 }, { t: `${r.captures} — ${r.n_genera_captured} genera in all. ${r.breadth ? r.breadth.charAt(0).toUpperCase() + r.breadth.slice(1) + '.' : ''}`, s: 19 }], { after: 40 }));
    const mis = [r.misfires_training ? `among training genomes: ${r.misfires_training}` : '', r.misfires_external ? `on the external genomes: ${r.misfires_external}` : ''].filter(x => x).join('; ');
    C.push(Rich([{ t: 'Where it goes wrong: ', b: true, s: 19 }, { t: mis || 'no misfires among the genomes checked.', s: 19 }], { after: 160 }));
  }
  if (neg.length) {
    C.push(P('Trusted rules that rule a trait out', { bold: true, after: 60 }));
    C.push(P('These fire when an organism lacks the genes the physiology needs. Because most traits are rare, they are only listed when they clearly beat the base rate.', { size: 19, color: MUTED }));
    C.push(table(['Trait', 'Rule (genes)', 'Probability', 'Why', 'Example organisms'],
      neg.map(r => [nice(r.trait), genesOf(r), pct(r.probability), r.logic, r.captures.split(', ').slice(0, 3).join(', ')]),
      [1500, 2000, 900, 3160, 1800], { left: [1, 3, 4] }));
    C.push(spacer());
  }
}


// ------------------------------------------------------------------ 3. rule cards
C.push(H('Supplementary Note 4. Rule cards: every rule that predicts a trait', HeadingLevel.HEADING_1));
for (const [k, b] of doneTasks) {
  C.push(H(b.name, HeadingLevel.HEADING_2));
  const pos = b.rules.filter(r => r.cls === 'has trait').sort((a, c) => a.trait.localeCompare(c.trait) || a.rank - c.rank);
  let last = null;
  for (const r of pos) {
    if (r.trait !== last) { C.push(P(nice(r.trait), { bold: true, size: 22, color: ACCENT, after: 60, keepNext: true })); last = r.trait; }
    const inf = r.informative_external ?? r.informative_heldout;
    const fill = r.high_probability_and_informative ? 'EEF6F0' : (inf ? 'FFF9EC' : 'FBEFEA');
    C.push(ruleBox(r.rule, fill));
    const ext = r.external_precision !== null && r.external_precision !== undefined;
    C.push(Rich([
      { t: 'Held-out genera: ', b: true, s: 19 }, { t: r.heldout_support ? `${pct(r.heldout_precision)} [${pct(r.heldout_ci_low)}–${pct(r.heldout_ci_high)}], n = ${r.heldout_support} (base rate ${pct(r.base_rate)}).  ` : 'never relearned in a held-out fold — no held-out estimate.  ', s: 19 },
      ...(ext ? [{ t: 'External genomes: ', b: true, s: 19 }, { t: `${pct(r.external_precision)} [${pct(r.external_ci_low)}–${pct(r.external_ci_high)}], n = ${r.external_support} (base rate ${pct(r.external_base_rate)}; ${r.external_n_pos ?? '–'} external genomes with the trait${(r.external_n_pos ?? 99) < 5 ? ' — too few for a reliable test' : ''}).  `, s: 19 }] : []),
      { t: 'Coverage ', b: true, s: 19 }, { t: `${pct(r.coverage)} of carriers in the training genomes; stability ${r.cv_stability}.`, s: 19 },
    ], { after: 40 }));
    C.push(Rich([{ t: 'Why it works: ', b: true, s: 19 }, { t: r.mechanism.split(' | ').join('; '), s: 19 }], { after: 40 }));
    C.push(Rich([{ t: `Trust — ${r.interpretation_level}: `, b: true, s: 19, c: r.interpretation_level === 'contains unexplained conditions' ? 'A15C00' : '2E7D32' },
      { t: r.interpretation, s: 19 }], { after: r.failure_note ? 40 : 140 }));
    if (r.failure_note) C.push(Rich([{ t: 'Why it fails on the external genomes: ', b: true, s: 19, c: 'B3261E' }, { t: r.failure_note, s: 19 }], { after: 140 }));
  }
}


// ------------------------------------------------------------------ 7. ruling out
C.push(H('Supplementary Table S1. Rules that rule a trait out', HeadingLevel.HEADING_1));
C.push(P('Because most traits are rare, “this organism lacks it” is right most of the time even by guessing. These rules are listed with their base rate; only those whose lower confidence bound beats the base rate are marked informative.'));
for (const [k, b] of doneTasks) {
  const neg = b.rules.filter(r => r.cls === 'lacks trait');
  if (!neg.length) continue;
  C.push(P(b.name, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Rule (genes)', 'Probability', 'Base rate', 'Informative', 'Trust'],
    neg.map(r => { const inf = r.informative_external ?? r.informative_heldout;
      const genes = r.rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, '').replace(/ AND /g, ' + ');
      return { cells: [nice(r.trait), genes, pct(r.probability), pct(baseFor(r)), inf ? 'yes' : 'no', r.interpretation_level],
        fills: [0, 0, 0, 0, inf ? GOOD : BAD, LEVEL_BG[r.interpretation_level]].map(x => x || undefined) }; }),
    [2000, 3500, 1000, 900, 900, 1060], { left: [1, 5] }));
  C.push(spacer());
}


// ------------------------------------------------------------------ 4. what the forest looks at (XAI)
C.push(H('Supplementary Table S2. Genes the random forest relies on', HeadingLevel.HEADING_1));
C.push(P('The genes the random forest relies on, measured on genera it was not trained on: first screened by SHAP, then confirmed by permutation importance (how much held-out accuracy drops when the gene is scrambled). Only genes in the forest’s top 10 in at least 3 of 5 folds are shown. Direction + means carrying the gene pushes the prediction towards the trait.'));
const TXD = JSON.parse(fs.readFileSync(`${ROOT}/results/trees_xai/report_data.json`, 'utf8'));
const ovl = ((TXD.tasks.da || {}).rf_nn_overlap || []).map(r => r.share_of_nn_in_rf).sort((a, b) => a - b);
const ovMed = ovl.length ? ovl[Math.floor((ovl.length - 1) / 2)] * 0.5 + ovl[Math.ceil((ovl.length - 1) / 2)] * 0.5 : null;
C.push(P(`Read with one caveat: a neural network that predicts as well as the forest relies on largely different genes (median ${pct(ovMed)} overlap of stable top-10 genes, donors/acceptors, all traits). When many genes carry the same information, equally good models explain themselves differently, so a gene named here is a lead, not proof — the rules in the rule cards (Supplementary Note 4) and the independent checks are the evidence.`, { size: 19, color: MUTED }));
for (const [k, b] of doneTasks) {
  const x = b.forest_xai || [];
  if (!x.length) continue;
  C.push(P(b.name, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Gene', 'Function (KEGG)', 'In top 10', 'Direction', 'Known marker'],
    x.map(r => [nice(r.trait), r.KO, r.name, `${Math.round(r.stability * 5)}/5 folds`, (r.shap_direction ?? 0) > 0 ? '+' : '−', r.known_marker ? 'yes' : 'no']),
    [2100, 800, 3700, 1000, 800, 960], { left: [1, 2] }));
  C.push(spacer());
}


// ------------------------------------------------------------------ 5. rules the forest itself uses
C.push(H('Supplementary Table S3. Rules extracted from the random forest', HeadingLevel.HEADING_1));
C.push(P('Every one of the forest’s 1 500 trees starts with a few yes/no gene questions. Counting which first one or two questions recur across the trees gives the rules the forest actually relies on (the SIRUS method; Bénard et al. 2021). Each rule gives the probability of the trait among genomes it covers, next to the base rate. Here both directions are useful: a rule that drops the probability far below the base rate reliably rules a trait out.'));
for (const [k, b] of doneTasks) {
  const fr = (b.forest_rules || []).filter(r => r.rank <= 4);
  if (!fr.length) continue;
  C.push(P(b.name, { bold: true, after: 60 }));
  const col = (r, ext) => { const hasE = ext && r.p_external !== undefined && r.p_external !== null;
    if (ext && !hasE) return { txt: '–', fill: undefined };
    const p = ext ? r.p_external : r.p_heldout, base = ext ? r.external_base_rate : r.base_rate;
    const inf = ext ? r.informative_external : r.informative_heldout;
    return { txt: `${pct(p)} (${pct(base)})`, fill: inf ? (p > base ? GOOD : 'E3EEF8') : undefined }; };
  C.push(table(['Trait', 'Forest rule (genes)', 'P held-out (base)', 'P external (base)', 'Trees', 'Stable', 'Trust'],
    fr.map(r => { const genes = r.rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, '').replace(/ AND /g, ' + ');
      const h = col(r, false), e = col(r, true);
      return { cells: [nice(r.trait), genes, h.txt, e.txt, pct(r.tree_frequency), r.stability, r.interpretation_level],
        fills: [undefined, undefined, h.fill, e.fill, undefined, undefined, LEVEL_BG[r.interpretation_level]] }; }),
    [1800, 2500, 1200, 1200, 700, 700, 1260], { left: [1, 6] }));
  C.push(P('Each probability is shown with the base rate of the same genomes in brackets. Green: significantly above the base rate; blue: significantly below it (the rule reliably rules the trait out). “Trees” = share of all trees that start with this question.', { size: 17, color: MUTED }));
  C.push(spacer());
}
if (DA && (DA.forest_top || []).length) { const ft = DA.forest_top, onP = ft.filter(x => x.on_pathway);
  C.push(P(`For ${onP.length} of ${ft.length} donor/acceptor traits the forest’s most frequent first question is a pathway gene (${onP.map(x => `${x.name} for ${nice(x.trait)}`).join(', ')}); for the rest it is a gene with no curated link, and such rules tend to fail on new organisms (for H₂ oxidation, fumarate hydratase), as the XAI-built rules did. No single question dominates: the most frequent one appears in ${pct(Math.min(...ft.map(x => x.tree_frequency)))}–${pct(Math.max(...ft.map(x => x.tree_frequency)))} of trees.`)); }


// ------------------------------------------------------------------ 6. single decision tree
C.push(H('Supplementary Table S4. Rules of the single decision tree', HeadingLevel.HEADING_1));
C.push(P('A single decision tree is a rule list by construction. Shown are its positive rules that were learned again in at least 3 of the 5 folds, with precision measured on the held-out genera of each fold.'));
for (const [k, b] of doneTasks) {
  const dr = b.dt_rules || [];
  if (!dr.length) continue;
  C.push(P(b.name, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Tree rule', 'Folds', 'Held-out support', 'Precision'],
    dr.map(r => [nice(r.trait), r.rule.replace(/ present/g, '+').replace(/ absent/g, '−').replace(/ AND /g, ', '), `${r.folds}/5`, r.support, pct(r.precision)]),
    [2200, 4300, 800, 1100, 960], { left: [1] }));
  C.push(P('K-number + = gene present, − = absent.', { size: 17, color: MUTED }));
  C.push(spacer());
}


// ------------------------------------------------------------------ 11. what transfers
C.push(H('Supplementary Note 5 and Table S5. Transfer to new lineages and the rule lists as classifiers', HeadingLevel.HEADING_1));
C.push(P(`The real test of a rule is organisms it has never seen. For donors and acceptors every rule was also applied to ${NEXT} external genomes whose genus is absent from training, labelled by a different curation. For ammonia oxidation the external label is “nitrification”, a broader trait.`));
const da = D.tasks.da;
if (da && da.transfer_by_level) {
  C.push(table(['Trust level', 'Positive rules', 'Held-out probability', 'External probability', 'Drop'],
    da.transfer_by_level.map(r => [r.interpretation_level, r.n, pct(r.heldout), pct(r.external), f(r.drop)]),
    [2800, 1400, 1800, 1800, 1560], { left: [0] }));
  C.push(spacer());
  const fails = da.rules.filter(r => r.cls === 'has trait' && r.external_precision !== null && r.external_precision !== undefined && r.external_precision < 0.5 && r.heldout_precision >= 0.8)
    .map(r => `${nice(r.trait)} from ${genesOf(r)} [${r.interpretation_level}] (held-out ${pct(r.heldout_precision)} → external ${pct(r.external_precision)}, n = ${r.external_support})${r.failure_note ? ` — ${r.failure_note}` : ''}`);
  if (fails.length) C.push(bullet(`Rules that were right at least 80 % of the time on held-out genera but below 50 % on the external genomes: ${fails.join('; ')}.`));
  C.push(P('This comparison was defined after the first rulebook run (Addendum H), the groups are small and the trust levels are themselves post hoc, so it is exploratory.', { size: 19, color: MUTED }));
}
for (const [k, b] of doneTasks) {
  if (!b.procedure) continue;
  const main = b.procedure.filter(r => (r.selection ?? 'main') === 'main'), strict = b.procedure.filter(r => r.selection === 'strict');
  C.push(P(`${b.name}: the whole positive rule list used as a classifier, on held-out genera. Recall counts every positive in every fold, including folds where no rule was learned. “Strict” re-runs the procedure with the candidate genes also chosen inside each training fold (a lineage-collapsed correlation screen), so the held-out genera played no part in choosing the genes (Addendum L). Strict and main estimates are close, so gene-selection leakage is not what makes held-out precision exceed external precision.`, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Flagged', 'Correct', 'Precision', 'Recall', 'Strict precision', 'Strict recall'],
    main.map(r => { const st = strict.find(x => x.trait === r.trait);
      return [nice(r.trait), r.fired, r.tp, r.fired ? pct(r.heldout_precision) : '–', pct(r.heldout_recall), st ? (st.fired ? `${pct(st.heldout_precision)} (${st.fired})` : '– (0)') : '–', st ? pct(st.heldout_recall) : '–']; }),
    [2600, 900, 900, 1100, 1000, 1500, 1360]));
  C.push(P('Strict precision is followed by the number of genomes flagged in brackets; – = no genome flagged.', { size: 17, color: MUTED }));
  C.push(spacer());
}


const doc = new Document({ creator: 'bio_pred', title: 'Supplementary Information',
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 360, hanging: 200 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [new TextRun({
      children: ['Page ', PageNumber.CURRENT, ' of ', PageNumber.TOTAL_PAGES], size: 18, color: MUTED, font: 'Calibri' })] })] }) },
    children: C }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(`${ROOT}/docs/rulebook_gene_rules_supplementary.docx`, b); console.log('written docs/rulebook_gene_rules_supplementary.docx'); });
