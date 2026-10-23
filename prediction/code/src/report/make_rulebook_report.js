// Manuscript: interpretable gene rules for microbial electron donors/acceptors and oxygen use.
// Every number is read from results/rulebook/rulebook_data.json (src/52_rulebook_report_data.py) and
// results/trees_xai/report_data.json (src/50_trees_xai_report_data.py); figures from src/58_rulebook_figures.py.
// Supplementary material: src/report/make_rulebook_supplement.js.
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell, ImageRun,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer, PageBreak,
  LineNumberRestartFormat,
} = require('docx');
const fs = require('fs');

const ROOT = '.';
const D = JSON.parse(fs.readFileSync(`${ROOT}/results/rulebook/rulebook_data.json`, 'utf8'));
const TX = JSON.parse(fs.readFileSync(`${ROOT}/results/trees_xai/report_data.json`, 'utf8'));
const TITLE = require('./rulebook_title.js');
const FIG = `${ROOT}/results/rulebook/figures`;
const HEAD_BG = 'EDF2FA', ALT_BG = 'F7F8FA';
const BODY_LINE = 480;          // double spacing for the main text (journal review format)

const f = (x, n = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–' : Number(x).toFixed(n);
const pct = x => (x === null || x === undefined || Number.isNaN(x)) ? '–' : `${Math.round(100 * x)}%`;
const PRETTY = {
  acceptor_CO2_methanogen_acetogen: 'CO₂ reduction (methanogens, acetogens)', acceptor_oxygen: 'O₂ respiration',
  acceptor_nitrate: 'nitrate reduction', acceptor_sulfate: 'sulfate reduction', acceptor_sulfur: 'sulfur reduction',
  acceptor_iron: 'iron reduction', acceptor_other: 'other acceptors', donor_sulfur_oxidation: 'sulfur oxidation', donor_ammonia_oxidation: 'ammonia oxidation',
  donor_nitrite_oxidation: 'nitrite oxidation', donor_hydrogen_oxidation: 'H₂ oxidation', donor_CO_oxidation: 'CO oxidation',
  donor_iron_oxidation: 'iron oxidation', donor_other: 'other donors', uses_O2: 'oxygen use',
};
const W = ['Zero', 'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven', 'Eight', 'Nine', 'Ten', 'Eleven', 'Twelve'];
const Num = k => (k >= 0 && k < W.length ? W[k] : String(k));   // a number that starts a sentence
const nice = t => PRETTY[t] || t.replace(/_/g, ' ');
const genesOf = r => r ? r.rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, '').replace(/ AND /g, ' + ') : '–';

// ------------------------------------------------------------------ paragraph helpers
const run = (t, o = {}) => new TextRun({ text: t, bold: o.b, italics: o.i, size: o.s ?? 22, font: 'Calibri', color: o.c });
const P = (text, o = {}) => new Paragraph({ spacing: { after: o.after ?? 120, line: o.line ?? BODY_LINE }, alignment: o.align,
  keepNext: o.keepNext, children: Array.isArray(text) ? text.map(x => typeof x === 'string' ? run(x, o) : run(x.t, { ...o, ...x })) : [run(text, o)] });
const H1 = t => new Paragraph({ heading: HeadingLevel.HEADING_1, keepNext: true, spacing: { before: 360, after: 160, line: 276 },
  children: [new TextRun({ text: t, bold: true, size: 28, font: 'Calibri', color: '000000' })] });
const H2 = t => new Paragraph({ heading: HeadingLevel.HEADING_2, keepNext: true, spacing: { before: 240, after: 120, line: 276 },
  children: [new TextRun({ text: t, bold: true, italics: true, size: 23, font: 'Calibri', color: '000000' })] });
const caption = (label, text) => new Paragraph({ spacing: { before: 80, after: 240, line: 276 },
  children: [run(`${label}${label ? ' ' : ''}`, { b: true, s: 19 }), run(text, { s: 19 })] });
const capAbove = (label, text) => new Paragraph({ keepNext: true, spacing: { before: 200, after: 80, line: 276 },
  children: [run(`${label} `, { b: true, s: 19 }), run(text, { s: 19 })] });
const bullet = text => new Paragraph({ numbering: { reference: 'bul', level: 0 }, spacing: { after: 60, line: BODY_LINE }, children: [run(text)] });
const ref = (k, text) => new Paragraph({ spacing: { after: 80, line: 276 }, indent: { left: 440, hanging: 440 },
  children: [run(`${k}.\t`, { s: 20 }), run(text, { s: 20 })] });

function table(headers, rows, widths, o = {}) {
  const mk = (cells, { bg, bold } = {}) => new TableRow({ cantSplit: true, tableHeader: bold, children: cells.map((t, i) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    shading: bg ? { type: ShadingType.CLEAR, fill: bg, color: 'auto' } : undefined,
    margins: { top: 40, bottom: 40, left: 70, right: 70 },
    children: [new Paragraph({ alignment: (o.left || [0]).includes(i) ? AlignmentType.LEFT : AlignmentType.CENTER,
      spacing: { after: 0, line: 240 }, children: [new TextRun({ text: String(t), bold, size: 17, font: 'Calibri' })] })],
  })) });
  return new Table({ columnWidths: widths, width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: { top: { style: BorderStyle.SINGLE, size: 6, color: '000000' }, bottom: { style: BorderStyle.SINGLE, size: 6, color: '000000' },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: 'D9D9D9' }, insideVertical: { style: BorderStyle.NONE } },
    rows: [mk(headers, { bg: HEAD_BG, bold: true }), ...rows.map((r, k) => mk(r, { bg: k % 2 ? ALT_BG : undefined }))] });
}
const figure = (file, wIn, hIn) => new Paragraph({ alignment: AlignmentType.CENTER, keepNext: true, spacing: { before: 200, after: 40 },
  children: [new ImageRun({ type: 'png', data: fs.readFileSync(`${FIG}/${file}`), transformation: { width: Math.round(wIn * 96), height: Math.round(hIn * 96) } })] });

// ------------------------------------------------------------------ numbers used in the text
const DA = D.tasks.da, OX = D.tasks.ox, NEXT = D.n_external_novel;
const macro = (task, arm, lvl = 'genus') => (TX.tasks[task].macro || []).find(r => r.arm === arm && r.split_level === lvl) || {};
const rfDA = macro('da', 'RF'), nullDA = macro('da', 'null_taxNN'), prevDA = macro('da', 'null_prevalence'), rfDAf = macro('da', 'RF', 'family');
const rfOX = macro('ox', 'RF'), nullOX = macro('ox', 'null_taxNN'), prevOX = macro('ox', 'null_prevalence');
const nPredDA = TX.tasks.da.n_predictable_best, nTraitsDA = new Set(DA.rules.map(r => r.trait)).size;
const reliable = b => b.rules.filter(r => r.cls === 'has trait' && r.high_probability_and_informative);
const relTraitsDA = [...new Set(reliable(DA).map(r => r.trait))];
const bestRel = (b, t) => reliable(b).filter(r => r.trait === t).sort((a, c) => c.coverage - a.coverage)[0];
const posDA = DA.rules.filter(r => r.cls === 'has trait'), negDA = DA.rules.filter(r => r.cls === 'lacks trait');
const infNeg = b => b.rules.filter(r => r.cls === 'lacks trait' && r.high_probability_and_informative);
const has = x => x !== null && x !== undefined;
const withExt = posDA.filter(r => has(r.external_precision) && has(r.heldout_precision) && r.external_support >= 5);
const median = a => { const s = [...a].sort((x, y) => x - y); const m = (s.length - 1) / 2; return s.length ? (s[Math.floor(m)] + s[Math.ceil(m)]) / 2 : null; };
const medH = median(withExt.map(r => r.heldout_precision)), medE = median(withExt.map(r => r.external_precision));
const tb = [...(DA.transfer_by_level || [])].sort((a, b) => a.drop - b.drop);
const proc = b => { const m = (b.procedure || []).filter(r => (r.selection ?? 'main') === 'main' && r.fired > 0);
  return m.map(r => ({ trait: r.trait, main: r.heldout_precision, strict: ((b.procedure || []).find(x => x.selection === 'strict' && x.trait === r.trait && x.fired > 0) || {}).heldout_precision }))
    .filter(x => has(x.strict)); };
const pr = [...proc(DA), ...proc(OX)], maxDiff = Math.max(...pr.map(x => Math.abs(x.main - x.strict)));
const medDiff = median(pr.map(x => x.strict - x.main));
const ft = DA.forest_top || [], ftOn = ft.filter(x => x.on_pathway);
const pred = id => D.predictions.find(p => p.id === id) || {};
const shifts = [...new Map(posDA.filter(r => has(r.external_base_rate) && Math.abs(r.base_rate - r.external_base_rate) >= 0.2)
  .map(r => [r.trait, { t: r.trait, a: r.base_rate, b: r.external_base_rate }])).values()];
const fRule = (b, trait, gene, present) => (b.forest_rules || []).find(r => r.trait === trait && r.genes === gene && r.rule.includes(`) ${present ? 'present' : 'absent'}`));
const pExt = r => r ? (r.p_external ?? r.p_heldout) : null;
const hasK = (r, ko) => !!r && r.rule.includes(`(${ko})`);
const s = bestRel(DA, 'acceptor_sulfate'), o = bestRel(DA, 'acceptor_oxygen'), n = bestRel(DA, 'acceptor_nitrate'), c = bestRel(DA, 'acceptor_CO2_methanogen_acetogen');
const amo = fRule(DA, 'donor_ammonia_oxidation', 'K10945', true), dsrKabs = fRule(DA, 'acceptor_sulfate', 'K27188', false), cdhDabs = fRule(DA, 'acceptor_CO2_methanogen_acetogen', 'K00194', false);
const oxRel = reliable(OX), oxNeg = infNeg(OX), oxPos = OX.rules.filter(r => r.cls === 'has trait');
const relLine = t => { const r = bestRel(DA, t); return `${nice(t)} (${genesOf(r)}; ${pct(r.probability)} on ${r.probability_source === 'external genomes' ? 'external genomes' : 'held-out genera'})`; };
const collapsed = withExt.filter(r => r.heldout_precision >= 0.99 && r.external_precision < 0.3);
const collapsedTraits = [...new Set(collapsed.map(r => r.trait))];
const knownFail = posDA.filter(r => r.failure_note && r.heldout_precision >= 0.8 && has(r.external_precision) && r.external_precision < 0.5);
const anaer = oxNeg.find(r => hasK(r, 'K00198') && hasK(r, 'K04070'));
const gp2 = pred('G-P2 (da)');
const amoPos = (posDA.find(r => r.trait === 'donor_ammonia_oxidation' && has(r.external_n_pos)) || {}).external_n_pos;

const C = [];
// ------------------------------------------------------------------ title page
C.push(new Paragraph({ spacing: { after: 240, line: 360 }, children: [new TextRun({ text: TITLE, bold: true, size: 32, font: 'Calibri' })] }));
C.push(P('Running title: Gene rules for microbial respiration across lineages', { line: 360 }));
C.push(P('[Author names]¹*', { line: 360, after: 40 }));
C.push(P('¹[Affiliation, city, country]', { line: 360, after: 40, s: 20 }));
C.push(P('*Corresponding author: [name, e-mail]', { line: 360, after: 240, s: 20 }));
C.push(P([{ t: 'Keywords: ', b: true }, 'electron acceptor; electron donor; oxygen; KEGG orthology; interpretable machine learning; rule learning; random forest; external validation; genome-based phenotype prediction'], { line: 360 }));
C.push(new Paragraph({ children: [new PageBreak()] }));

// ------------------------------------------------------------------ abstract
C.push(H1('Abstract'));
C.push(P([{ t: 'Background. ', b: true }, 'Short gene-presence rules can make genome-based predictions of microbial energy metabolism transparent, but whether they hold on lineages absent from training is rarely tested.']));
C.push(P([{ t: 'Methods. ', b: true }, `We learned IF–THEN rules of up to three KEGG orthologue conditions for ${nTraitsDA} electron-donor and electron-acceptor categories (358 genomes) and for oxygen use (2 927 genomes), with candidate genes from random-forest explanations. Rules were scored on held-out genera and, for donors and acceptors, on ${NEXT} independently labelled genomes from genera absent from training.`]));
C.push(P([{ t: 'Results. ', b: true }, `For donors and acceptors, a random forest reached a macro-averaged AUPRC of ${f(rfDA.macro_auprc, 3)} on held-out genera (taxonomic null ${f(nullDA.macro_auprc, 3)}). ${Num(relTraitsDA.length)} of ${nTraitsDA} categories had a reliable rule on the external genomes: ${relTraitsDA.map(t => `${nice(t)} ${pct(bestRel(DA, t).probability)}`).join(', ')}. Across the ${withExt.length} testable positive rules, median precision fell from ${pct(medH)} on held-out genera to ${pct(medE)} externally, and rules perfect on held-out genera fell as low as ${pct(Math.min(...withExt.filter(r => r.heldout_precision >= 0.99).map(r => r.external_precision)))}. Choosing candidate genes strictly inside training folds changed held-out precision by at most ${f(maxDiff * 100, 0)} percentage points, so the drop reflects lineage and curation shift, not gene-selection leakage. The pre-registered prediction that rules transfer within 0.10 was not supported (${gp2.evidence ? gp2.evidence.split(' (')[0] : '–'}).`]));
C.push(P([{ t: 'Conclusions. ', b: true }, 'A few gene rules predict electron acceptors precisely on new lineages. Held-out precision did not predict which rules transfer, and a biological reading of the rule was only a partial guide. Rules need an external test on lineage-novel genomes before they are trusted.']));

// ------------------------------------------------------------------ introduction
C.push(H1('Introduction'));
C.push(P('Which electron donors and acceptors an organism uses defines its role in biogeochemical cycles. Yet for most sequenced prokaryotes this has never been tested in culture. Genome-based prediction fills the gap: the enzymes of respiration and lithotrophy are encoded in the genome and can be recognised by orthology (1, 2). Curated trait compilations (3) and taxonomy-based functional annotation (4) make it possible to train and test such predictions at scale.'));
C.push(P('Two problems limit current approaches. First, accurate models are often opaque. A random forest (5) can predict a trait without revealing which genes it uses, and explanation methods such as SHAP (6) describe one fitted model rather than a rule a microbiologist can check (7). Second, genomes are phylogenetically clustered. A gene that marks one well-sampled lineage can look predictive without any functional link to the trait, so accuracy on held-out genomes from known lineages can overstate accuracy on new ones.'));
C.push(P('Short IF–THEN rules over gene presence address the first problem directly, because each condition can be read against known enzyme function (8). Whether they also address the second is an empirical question. Here we learn such rules for electron donors, electron acceptors and oxygen use. We measure their precision on held-out genera and on independently curated genomes from novel genera, and ask which rules transfer and why the others do not.'));

// ------------------------------------------------------------------ methods
C.push(H1('Methods'));
C.push(H2('Genomes, labels and features'));
C.push(P(`Donor/acceptor task. 358 genomes were labelled from a curated electron donor–acceptor database [reference to be added] in ${nTraitsDA} binary categories: 7 electron donors and 7 electron acceptors. Proteomes were retrieved from NCBI RefSeq and annotated with KOfam profiles (2). Features were the presence or absence of the 852 KOs of an energy-metabolism panel. A trait not recorded for an organism was treated as absent (closed-world labels).`));
C.push(P('Oxygen-use task. 2 927 genomes whose species has an oxygen-relationship record in the trait synthesis of Madin et al. (3). Species recorded as aerobic, obligate aerobic, facultative or microaerophilic were labelled O₂ users; anaerobic and obligate anaerobic species were labelled non-users. Features were the full KO profiles, restricted to KOs present in 1–99% of genomes (7 594 KOs).'));
C.push(P(`External set. ${NEXT} genomes of species whose genus is absent from the donor/acceptor training data, one RefSeq genome per species. They were labelled from the pathway annotations of Madin et al. (3), which combine 26 sources including FAPROTAX (4) and are curated independently of the training database. Twenty initially selected genomes whose genus occurs in training were removed. The set was enriched for positives: every species positive for a mapped trait, plus an equal random sample of negatives. Its base rates therefore differ from those of the training data. For ammonia oxidation the external label is nitrification, a broader category.`));
C.push(H2('Validation design'));
C.push(P('All evaluation used nested cross-validation with five frozen folds that hold out whole genera; a stricter family-level split was also used. Folds were fixed and hashed before any model was fitted, and hyper-parameters and decision thresholds were selected in inner folds only. Model performance was measured as the area under the precision–recall curve (AUPRC), with 95% confidence intervals from a cluster bootstrap over genera (1 000 resamples). Each model was compared with a prevalence null and with a taxonomic null. The taxonomic null predicts the mean label of the training genomes that share the lowest taxonomic rank with the test genome.'));
C.push(H2('Random forest and explanations'));
C.push(P('Random forests (5) of 500 trees were trained in three seeds with balanced class weights; the minimum leaf size was chosen in inner folds. On each held-out fold, KOs were screened by mean absolute TreeSHAP value (6). The 50 highest-ranked KOs, plus all known marker KOs, were confirmed by permutation importance (10 repeats, drop in AUPRC). A KO was called stable when it ranked among the ten most important in at least three of five folds. Associations were also tested without a model, by a Cochran–Mantel–Haenszel test (9) stratified by phylum and genome size, with Benjamini–Hochberg correction (10).'));
C.push(H2('Rule learning'));
C.push(P('For each trait, up to 16 candidate KOs were taken, in this order, from:'));
C.push(bullet('genes on which the forest and the model-free test agreed;'));
C.push(bullet('genes stable in the forest;'));
C.push(bullet('genes used by depth-3 surrogate trees of the forest;'));
C.push(bullet('known markers significant in the stratified test.'));
C.push(P('Rules were conjunctions of one to three conditions (a KO present or absent), with at most one condition per KO. They were learned by sequential covering (8). At each step, the rule with the highest Wilson 95% lower bound on precision (11) was added, provided it fired on at least 8 genomes from at least 5 genera and that bound was at least 0.60. Up to three rules predicting a trait and two ruling it out were learned per trait.'));
C.push(P('Validation was nested: rules were learned on four genus folds and scored on the fifth. A rule’s held-out precision pools the folds in which the identical rule was learned again, and its stability is the number of such folds. Final rules were learned on all genomes. A rule was called high-probability when its precision was at least 0.80 with a Wilson lower bound of at least 0.60 (pre-registered). Precision was measured on the external genomes where available and on held-out genera otherwise. A rule was called informative when that lower bound also exceeded the base rate of the same genomes (added after the first run). “Reliable” denotes both.'));
C.push(P('The candidate genes were chosen with information from all genera, which could make held-out precision optimistic. We therefore re-ran the whole procedure with the candidate genes chosen inside each training fold only (strictly nested). There, genes were ranked by the absolute correlation between genus-mean presence and genus-mean label across genera. We compared the held-out precision of each trait’s positive rule list under both procedures.'));
C.push(H2('Rules extracted from the forest'));
C.push(P('Following SIRUS (12), every root path of depth one or two in the 1 500 trees of the forest (three seeds × 500 trees) was counted. The ten most frequent paths per trait were kept as rules and validated like the learned rules.'));
C.push(H2('Biological reading of rules'));
C.push(P('Each condition was compared with a curated table of gene functions, compiled from KEGG orthology and module definitions (1) and standard physiology. A condition supports a rule in two cases:'));
C.push(bullet('it is a gene of the trait’s pathway in the state that fits the rule: present for a rule that predicts the trait, absent for a rule that rules it out. An absent gene counts only when at least 80% of organisms with the trait carry it;'));
C.push(bullet('it is a marker of oxygen lifestyle whose signal fits the trait’s oxygen relationship.'));
C.push(P('A rule was read as mechanistic when all its conditions are pathway genes. It was read as physiological when all are pathway genes or lifestyle markers. Otherwise it was read as containing unexplained conditions: at least one condition has no known link to the trait, is uninformative, or contradicts the pathway. This reading was defined after the first rule run and is used to interpret rules, not to select them.'));
C.push(H2('Pre-registration and disclosure'));
C.push(P('The procedure and its predictions were written down before rules were learned, as addenda to the project’s pre-registration. The following were added later and are labelled as such: the base-rate criterion, the biological reading, the strictly nested estimate, and the restriction of the external set to lineage-novel genomes. The time stamps of the pre-registration were recorded by hand. Because the project was not under version control, they cannot be verified independently, and we do not claim that the plan demonstrably preceded the analysis. All deviations are listed in the code repository.'));
C.push(H2('Software'));
C.push(P('Analyses used Python 3.12 with scikit-learn (13), shap (6), pandas and NumPy.'));

// ------------------------------------------------------------------ results
C.push(H1('Results'));
C.push(H2('Gene presence predicts donors, acceptors and oxygen use beyond lineage'));
C.push(P(`For the ${nTraitsDA} donor/acceptor categories, the random forest reached a macro-averaged AUPRC of ${f(rfDA.macro_auprc, 3)} on held-out genera (95% CI ${f(rfDA.macro_auprc_ci_low, 3)}–${f(rfDA.macro_auprc_ci_high, 3)}) and ${f(rfDAf.macro_auprc, 3)} on held-out families. The taxonomic null reached ${f(nullDA.macro_auprc, 3)} and the prevalence null ${f(prevDA.macro_auprc, 3)}. ${Num(nPredDA)} of ${nTraitsDA} categories passed the pre-registered gate, which requires an MCC lower bound above 0.2 on held-out families and a significant gain over the taxonomic null. For oxygen use the forest reached ${f(rfOX.macro_auprc, 3)}, against ${f(nullOX.macro_auprc, 3)} for the taxonomic null and ${f(prevOX.macro_auprc, 3)} for prevalence. Oxygen use is thus largely predictable from lineage alone.`));
C.push(H2('A few rules are reliable'));
C.push(P(`We learned ${DA.rules.length} rules for donors and acceptors (${posDA.length} predicting a trait, ${negDA.length} ruling one out) and ${OX.rules.length} for oxygen use. ${Num(relTraitsDA.length)} donor/acceptor categories have at least one reliable rule (Table 1).${o ? ` O₂ as terminal acceptor is predicted by ${genesOf(o)}, correct for ${pct(o.probability)} of the external genomes it fires on.` : ''}${s ? ` Sulfate respiration is predicted by ${genesOf(s)} (${pct(s.probability)}).` : ''}${c ? ` CO₂ reduction by methanogens and acetogens is predicted by ${genesOf(c)} (${pct(c.probability)}).` : ''}${n ? ` Nitrate respiration is predicted by ${genesOf(n)} (${pct(n.probability)}).` : ''} Each of these rules catches only part of the organisms with the trait (Table 1, coverage). No reliable rule was found for the other ${nTraitsDA - relTraitsDA.length} categories.`));
const better = (a, b) => !b || (a.high_probability_and_informative !== b.high_probability_and_informative ? a.high_probability_and_informative
  : (a.high_probability_and_informative ? a.coverage > b.coverage : (a.probability ?? 0) > (b.probability ?? 0)));
const bestRows = [];
for (const b of [DA, OX]) {
  const by = {};
  for (const r of b.rules.filter(x => x.cls === 'has trait' && (x.informative_external ?? x.informative_heldout))) if (better(r, by[r.trait])) by[r.trait] = r;
  for (const t of [...new Set(b.rules.map(r => r.trait))]) {
    const r = by[t];
    if (!r) { bestRows.push([nice(t), 'no informative rule', '–', '–', '–', '–', '–']); continue; }
    const ext = r.probability_source === 'external genomes';
    const lo = ext ? r.external_ci_low : r.heldout_ci_low, hi = ext ? r.external_ci_high : r.heldout_ci_high, nn = ext ? r.external_support : r.heldout_support;
    bestRows.push([nice(t), genesOf(r), `${pct(r.probability)} (${pct(lo)}–${pct(hi)}) ${ext ? 'E' : 'H'}${r.high_probability_and_informative ? '' : ' †'}`,
      pct(ext ? r.external_base_rate : r.base_rate), nn, pct(r.coverage), r.interpretation_level]);
  }
}
C.push(capAbove('Table 1.', 'The best informative rule for each trait. Among rules that predict a trait and exceed the base rate, the reliable rule with the highest coverage is shown; where no rule is reliable, the most precise informative rule is shown and marked †.'));
C.push(table(['Trait', 'Rule (genes)', 'Precision (95% CI)', 'Base rate', 'n', 'Coverage', 'Biological reading'], bestRows,
  [1800, 2700, 1500, 750, 500, 850, 1400], { left: [0, 1, 6] }));
C.push(caption('', 'E, measured on external genomes from novel genera; H, measured on held-out genera. Base rate: share of the same genomes with the trait. n: genomes the rule fired on. Coverage: share of training genomes with the trait that the rule catches. Biological reading: see Methods. The last row is oxygen use (2 927 genomes). Every rule with its evidence is given in Supplementary Note 4.'));

C.push(H2('Precision on held-out genera does not carry over to new lineages'));
C.push(P(`${Num(withExt.length)} positive donor/acceptor rules fired on at least five external genomes. Their median precision was ${pct(medH)} on held-out genera and ${pct(medE)} on the external genomes (Fig. 1). Every group of rules lost precision, by ${f(tb[0].drop * 100, 0)}–${f(tb[tb.length - 1].drop * 100, 0)} percentage points on average (Table 2); the groups are small, and the smallest loss was for ${tb[0].interpretation_level} rules (${tb[0].n} rules).${knownFail.some(r => r.trait === 'donor_ammonia_oxidation') ? ' The amoB rule failed because most external amoB carriers are methanotrophs, whose methane monooxygenase has the same KO as ammonia monooxygenase (20).' : ''}${knownFail.some(r => hasK(r, 'K27925')) ? ` The tetrathionate-hydrolase rule fired on only ${(knownFail.find(r => hasK(r, 'K27925')) || {}).external_support} external genomes, one of which is probably mislabelled.` : ''}${collapsedTraits.length ? ` Rules that were perfectly precise on held-out genera but contain unexplained conditions fell to ${pct(Math.min(...collapsed.map(r => r.external_precision)))}–${pct(Math.max(...collapsed.map(r => r.external_precision)))} externally (${collapsedTraits.map(nice).join(', ')}).` : ''}`));
C.push(figure('fig1_transfer.png', 3.5, 3.4));
C.push(caption('Figure 1.', `Precision of positive donor/acceptor rules on held-out genera (x) and on external genomes from novel genera (y), for the ${withExt.length} rules that fired on at least five external genomes. Marker size is proportional to the number of external genomes the rule fired on; colour and shape give the biological reading of the rule. Dashed line: equal precision.`));
C.push(capAbove('Table 2.', 'Transfer of positive donor/acceptor rules to external genomes, by biological reading (exploratory; groups defined after the first rule run).'));
C.push(table(['Biological reading', 'Rules', 'Mean precision, held-out genera', 'Mean precision, external genomes', 'Drop'],
  tb.map(r => [r.interpretation_level, r.n, pct(r.heldout), pct(r.external), f(r.drop)]), [2700, 800, 2100, 2100, 1100], { left: [0] }));
C.push(caption('', 'Positive rules with both a held-out and an external estimate are included.'));
C.push(P(`The drop is not explained by how candidate genes were chosen. When candidate genes were selected inside each training fold only, the held-out precision of each trait’s positive rule list changed by a median of ${(medDiff < 0 ? '−' : '+') + f(Math.abs(medDiff) * 100, 1)} percentage points and by at most ${f(maxDiff * 100, 0)} (Fig. 2). The external genomes differ from the training genomes in lineage, curation and the share of positives${shifts.length ? ` (${shifts.map(x => `${nice(x.t)}: ${pct(x.a)} of training genomes and ${pct(x.b)} of external genomes`).join('; ')})` : ''}.`));
C.push(figure('fig2_strict.png', 5.0, 3.7));
C.push(caption('Figure 2.', 'Held-out precision of the whole positive rule list per trait, for the main procedure (circles) and for the strictly nested procedure, in which candidate genes were chosen inside each training fold only (diamonds). Shown are the traits for which both procedures flagged at least one genome.'));

C.push(H2('The forest’s own rules are mostly lineage questions'));
C.push(P(`Counting the first one or two questions of all 1 500 trees gives the rules the forest relies on. For ${ftOn.length} of ${ft.length} donor/acceptor categories the most frequent first question is a gene of the trait’s pathway (${ftOn.map(x => `${x.name} for ${nice(x.trait)}`).join('; ')}). For the other ${ft.length - ftOn.length} it is a gene with no curated link to the trait. No single question dominates: the most frequent one appears in only ${pct(Math.min(...ft.map(x => x.tree_frequency)))}–${pct(Math.max(...ft.map(x => x.tree_frequency)))} of trees. Pathway questions were most useful for ruling a trait out. Without dsrK, the probability of sulfate respiration was ${pct(pExt(dsrKabs))}; without cdhD, the probability of a methanogenic or acetogenic lifestyle was ${pct(pExt(cdhDabs))} (Supplementary Table S3).`));
C.push(H2('Rules that rule a trait out'));
C.push(P(`Because most traits are rare, a rule that predicts absence is usually right even by chance. Of the ${negDA.length} donor/acceptor rules that rule a trait out, ${infNeg(DA).length} exceed the base rate of the same genomes; for oxygen use, ${oxNeg.length} of ${OX.rules.filter(r => r.cls === 'lacks trait').length} do.${anaer ? ` Anaerobes are identified by the anaerobic CO dehydrogenase (cooS) together with the pyruvate formate-lyase activating enzyme, with ${pct(anaer.probability)} precision on held-out genera.` : ''} All rule-out rules are listed in Supplementary Table S1.`));
C.push(H2('Pre-registered predictions'));
C.push(P('Table 3 scores the predictions made for the rules before they were learned. Verdicts are binary, as specified.'));
const PT = [['G-P1', 'da', 'At least one high-probability positive rule for at least 5 of the 7 predictable donor/acceptor traits'],
  ['G-P2', 'da', 'External precision of those rules within 0.10 of their held-out precision'],
  ['I-P1', 'da', 'The forest’s most frequent path shares a gene with the rules, for every predictable donor/acceptor trait'],
  ['I-P1', 'ox', 'As I-P1, for oxygen use']];
C.push(capAbove('Table 3.', 'Pre-registered predictions for the rules, scored.'));
const cleanEv = (id, e) => { if (!e) return '–';
  if (id === 'G-P1') { const m = e.match(/^(\d+) of (\d+).*alone (\d+) of (\d+)/); return m ? `${m[1]} of ${m[2]} traits on held-out genera; ${m[3]} of ${m[4]} on external genomes` : e; }
  if (id === 'G-P2') { const m = e.match(/^(\d+) of (\d+) rules within 0\.10 \(median \|difference\| ([0-9.]+)/); return m ? `${m[1]} of ${m[2]} rules within 0.10; median absolute difference ${m[3]}` : e; }
  return e.replace('predictable traits', 'predictable traits share a gene'); };
C.push(table(['ID', 'Prediction', 'Verdict', 'Evidence'], PT.map(([id, t, txt]) => { const p = pred(`${id} (${t})`); return [`${id} (${t === 'da' ? 'donors/acceptors' : 'oxygen use'})`, txt, p.verdict || '–', cleanEv(id, p.evidence)]; }),
  [1300, 3100, 1100, 4100], { left: [0, 1, 3] }));
C.push(caption('', 'G-P1 is scored on held-out genera, the reading implied by G-P2; on the external genomes alone it would not be supported. Predictions for model performance and explanations are reported in a companion analysis.'));

C.push(H2('Biological reading of the reliable rules'));
C.push(P('The following is interpretation based on established physiology; it was not tested here.', { i: true }));
C.push(P(`Sulfate reduction runs through ATP sulfurylase (sat), APS reductase (aprAB) and dissimilatory sulfite reductase (dsrAB) (14). None of these genes is specific on its own: sat also serves sulfate assimilation, and dsrAB runs in reverse in sulfur-oxidising bacteria.${s && (hasK(s, 'K00177') || hasK(s, 'K00176') || hasK(s, 'K00174')) ? ' The reliable rule adds a subunit of 2-oxoglutarate:ferredoxin oxidoreductase (kor), an O₂-sensitive enzyme common in anaerobes. Requiring this anaerobic marker together with the pathway may exclude aerobic sulfur oxidisers that carry dsrAB.' : ''}`));
C.push(P(`${o && hasK(o, 'K00413') ? 'For O₂ as terminal acceptor, the rule reads the respiratory chain through cytochrome c1 of the bc1 complex (petC)' + (hasK(o, 'K05601') ? ', together with the absence of the hybrid-cluster protein (hcp), which is made under anoxia. The absence of a single anaerobic enzyme is weak evidence on its own, so the biological reading counts this condition as unexplained.' : '.') : ''} Single oxidase genes are weaker markers, because aerobic respiration can use several oxidase families, and cytochrome bd is also kept by many anaerobes to scavenge traces of O₂ (15).`));
C.push(P(`${c ? 'For CO₂ reduction, the CODH/acetyl-CoA synthase complex (cdhD) of the Wood–Ljungdahl pathway is shared by acetogens and methanogens (16, 17).' + (hasK(c, 'K11261') ? ' The FwdE-annotated KO in the rule belongs to the methanofuran pathway in methanogens; acetogens carry a homolog outside that pathway.' : '') : ''}${n ? ` For nitrate, narI anchors the respiratory nitrate reductase in the membrane${hasK(n, 'K02276') ? '; the rule also requires an aa3 oxidase subunit, which fits organisms that respire both O₂ and nitrate but has no curated link to nitrate respiration itself' : ''}.` : ''}`));
C.push(P(`Where no reliable rule exists, annotation explains much of the failure. Hydrogenases form a large, diverse family that the KO set covers only in part (18). Nitrite oxidoreductase shares its KOs with respiratory nitrate reductase (19). amoB shares its KO with particulate methane monooxygenase (20): the amoB rule was right ${pct(amo && amo.p_heldout)} of the time on held-out genera and ${pct(amo && amo.p_external)} on the external genomes${has(amoPos) ? `, where only ${amoPos} genome${amoPos === 1 ? ' was an ammonia oxidiser' : 's were ammonia oxidisers'}` : ''}. The multiheme cytochromes of iron respiration (21) were not in the energy panel used here.`));

// ------------------------------------------------------------------ discussion
C.push(H1('Discussion'));
C.push(P(`Short gene rules predicted ${relTraitsDA.length} electron-acceptor categories with high precision on genomes from novel genera, and their conditions can be read against known physiology. The central result, however, is a warning: precision on held-out genera did not predict precision on new lineages. Some rules that were perfect on held-out genera lost most of their precision externally, and the strictly nested analysis shows that this is not an artefact of gene selection.`));
C.push(P('Three factors explain the gap. First, a genus-level hold-out keeps test genomes within families and phyla that are well represented in training, so a gene that marks such a lineage stays predictive within it. The external genomes come from other genera, often in other families. Second, the two curations differ in what they record, and the external set was enriched for positives, so base rates differ. Third, labels are closed-world in both sets: a trait never tested counts as absent, which lowers the measured precision for under-reported traits. A rule with few external positives, such as the ammonia-oxidation rule, cannot be judged externally at all.'));
C.push(P('Interpretability helped as a filter, not as a guarantee. Rules with unexplained conditions failed most often, but rules built only on pathway genes also failed when the annotation itself is ambiguous, as for amoB and pmoB. A biological reading is therefore a reason to prefer a rule; an external test is the evidence that it holds. The same caution applies to model explanations. The forest’s most frequent questions were pathway genes for only a minority of traits, and in a companion analysis a neural network of similar accuracy relied on largely different genes.'));
C.push(H2('Limitations'));
C.push(P('The donor/acceptor data contain 358 genomes, so several traits have few positives and wide confidence intervals. The external test depends on a single independent curation, maps its categories onto ours (nitrification for ammonia oxidation) and has no equivalent for three categories; oxygen use has no external test at all. KO annotation cannot separate some enzyme families. Several analysis choices were made after the first results: the base-rate criterion, the biological reading, the strictly nested estimate and the lineage-novel restriction of the external set. The time stamps of the pre-registration cannot be verified independently. Finally, the biological reading relies on a curated gene table that is necessarily incomplete.'));
C.push(H2('Conclusions'));
C.push(P('Gene-presence rules can make genome-based predictions of microbial energy metabolism transparent, and a few of them are precise on new lineages. Their precision must be measured on lineage-novel genomes from an independent curation. Neither genus-level cross-validation nor a plausible biological story substitutes for that test.'));

// ------------------------------------------------------------------ back matter
C.push(H1('Data and code availability'));
C.push(P('All code, frozen splits and rule tables, and the scripts that generate every number, table and figure in this article, are available at [repository URL or DOI to be added]. Genomes are public in NCBI RefSeq; trait labels are from the cited sources.'));
C.push(H1('Declarations'));
C.push(P([{ t: 'Funding: ', b: true }, '[to be added]. ', { t: 'Competing interests: ', b: true }, '[to be confirmed by the authors]. ', { t: 'Author contributions: ', b: true }, '[to be added]. ', { t: 'Acknowledgements: ', b: true }, '[to be added].']));
C.push(H1('References'));
const REFS = [
  'Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 2000;28:27–30.',
  'Aramaki T, Blanc-Mathieu R, Endo H, Ohkubo K, Kanehisa M, Goto S, Ogata H. KofamKOALA: KEGG Ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2020;36:2251–2252.',
  'Madin JS, Nielsen DA, Brbić M, Corkrey R, Danko D, Edwards K, et al. A synthesis of bacterial and archaeal phenotypic trait data. Sci Data. 2020;7:170.',
  'Louca S, Parfrey LW, Doebeli M. Decoupling function and taxonomy in the global ocean microbiome. Science. 2016;353:1272–1277.',
  'Breiman L. Random forests. Mach Learn. 2001;45:5–32.',
  'Lundberg SM, Erion G, Chen H, DeGrave A, Prutkin JM, Nair B, et al. From local explanations to global understanding with explainable AI for trees. Nat Mach Intell. 2020;2:56–67.',
  'Rudin C. Stop explaining black box machine learning models for high stakes decisions and use interpretable models instead. Nat Mach Intell. 2019;1:206–215.',
  'Fürnkranz J. Separate-and-conquer rule learning. Artif Intell Rev. 1999;13:3–54.',
  'Mantel N, Haenszel W. Statistical aspects of the analysis of data from retrospective studies of disease. J Natl Cancer Inst. 1959;22:719–748.',
  'Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B. 1995;57:289–300.',
  'Wilson EB. Probable inference, the law of succession, and statistical inference. J Am Stat Assoc. 1927;22:209–212.',
  'Bénard C, Biau G, da Veiga S, Scornet E. SIRUS: Stable and Interpretable RUle Set for classification. Electron J Stat. 2021;15:427–505.',
  'Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, et al. Scikit-learn: machine learning in Python. J Mach Learn Res. 2011;12:2825–2830.',
  'Pereira IAC, Ramos AR, Grein F, Marques MC, da Silva SM, Venceslau SS. A comparative genomic analysis of energy metabolism in sulfate reducing bacteria and archaea. Front Microbiol. 2011;2:69.',
  'Borisov VB, Gennis RB, Hemp J, Verkhovsky MI. The cytochrome bd respiratory oxygen reductases. Biochim Biophys Acta. 2011;1807:1398–1413.',
  'Ragsdale SW, Pierce E. Acetogenesis and the Wood–Ljungdahl pathway of CO₂ fixation. Biochim Biophys Acta. 2008;1784:1873–1898.',
  'Thauer RK, Kaster AK, Seedorf H, Buckel W, Hedderich R. Methanogenic archaea: ecologically relevant differences in energy conservation. Nat Rev Microbiol. 2008;6:579–591.',
  'Greening C, Biswas A, Carere CR, Jackson CJ, Taylor MC, Stott MB, et al. Genomic and metagenomic surveys of hydrogenase distribution indicate H₂ is a widely utilised energy source for microbial growth and survival. ISME J. 2016;10:761–777.',
  'Daims H, Lücker S, Wagner M. A new perspective on microbes formerly known as nitrite-oxidizing bacteria. Trends Microbiol. 2016;24:699–712.',
  'Holmes AJ, Costello A, Lidstrom ME, Murrell JC. Evidence that particulate methane monooxygenase and ammonia monooxygenase may be evolutionarily related. FEMS Microbiol Lett. 1995;132:203–208.',
  'Shi L, Dong H, Reguera G, Beyenal H, Lu A, Liu J, et al. Extracellular electron transfer mechanisms between microorganisms and minerals. Nat Rev Microbiol. 2016;14:651–662.',
];
REFS.forEach((r, i) => C.push(ref(i + 1, r)));
C.push(H1('Supplementary material'));
C.push(P('Supplementary Information (rulebook_gene_rules_supplementary.docx): Note 1, how to read a rule; Note 2, extended biological interpretation; Note 3, the reliable rules one by one; Note 4, rule cards for every rule; Table S1, rules that rule a trait out; Table S2, genes the random forest relies on; Table S3, rules extracted from the forest; Table S4, decision-tree rules; Note 5 and Table S5, transfer to new lineages and the rule lists as classifiers.', { line: 276 }));

const doc = new Document({ creator: '[Authors]', title: TITLE,
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 360, hanging: 240 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: 22 } } } },
  sections: [{ properties: { page: { size: { width: 11906, height: 16838 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } },
    lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({
      children: [PageNumber.CURRENT], size: 18, font: 'Calibri' })] })] }) },
    children: C }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(`${ROOT}/docs/rulebook_gene_rules.docx`, b); console.log('written docs/rulebook_gene_rules.docx'); });
