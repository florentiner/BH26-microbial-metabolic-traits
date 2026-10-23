// Interpretable models, XAI and gene–trait correlations — every number read from
// results/trees_xai/report_data.json (built by src/50_trees_xai_report_data.py). Pending tasks are marked, not guessed.
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer,
} = require('docx');
const fs = require('fs');

const ROOT = '.';
const R = JSON.parse(fs.readFileSync(`${ROOT}/results/trees_xai/report_data.json`, 'utf8'));
const T = R.tasks;
const ACCENT = '2A78D6', MUTED = '52514E', HEAD_BG = 'EDF2FA', ALT_BG = 'F7F8FA';
const OK_BG = 'DCEFE0', PART_BG = 'FFF3D6', NO_BG = 'FBE3D6', PEND_BG = 'EEEEEE';
const VBG = { supported: OK_BG, 'not supported': NO_BG, pending: PEND_BG, 'not run': PEND_BG };
const ADDENDA = [...fs.readFileSync(`${ROOT}/docs/preregistration_trees_xai.md`, 'utf8').matchAll(/^## Addendum ([A-Z])/gm)].map(m => m[1]);

const f = (x, n = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–' : Number(x).toFixed(n);
const pct = x => (x === null || x === undefined) ? '–' : `${Math.round(100 * x)} %`;
const q = x => (x === null || x === undefined) ? '–' : (x < 0.001 ? '< 0.001' : Number(x).toFixed(3));
const nice = t => t.replace(/^donor_/, 'donor: ').replace(/^acceptor_/, 'acceptor: ').replace(/_/g, ' ');

const P = (text, o = {}) => new Paragraph({ spacing: { after: o.after ?? 120, line: 276 },
  children: [new TextRun({ text, bold: o.bold, italics: o.italics, size: o.size ?? 21, color: o.color, font: 'Calibri' })] });
const Rich = (runs, o = {}) => new Paragraph({ spacing: { after: o.after ?? 120, line: 276 },
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i, size: 21, font: 'Calibri' })) });
const H = (text, level) => new Paragraph({ heading: level,
  spacing: { before: level === HeadingLevel.HEADING_1 ? 340 : 240, after: 130 },
  children: [new TextRun({ text, font: 'Calibri', bold: true, size: level === HeadingLevel.HEADING_1 ? 30 : 24,
    color: level === HeadingLevel.HEADING_1 ? ACCENT : '1A1A1A' })] });
const bullet = t => new Paragraph({ numbering: { reference: 'bul', level: 0 }, spacing: { after: 70, line: 276 },
  children: [new TextRun({ text: t, size: 21, font: 'Calibri' })] });
const rb = runs => new Paragraph({ numbering: { reference: 'bul', level: 0 }, spacing: { after: 70, line: 276 },
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, size: 21, font: 'Calibri' })) });
const spacer = () => new Paragraph({ spacing: { after: 100 }, children: [] });

function table(headers, rows, widths, o = {}) {
  const mk = (cells, { bg, bold, fills } = {}) => new TableRow({ children: cells.map((t, i) => new TableCell({
    width: { size: widths[i], type: WidthType.DXA },
    shading: (fills && fills[i]) || bg ? { type: ShadingType.CLEAR, fill: (fills && fills[i]) || bg, color: 'auto' } : undefined,
    margins: { top: 50, bottom: 50, left: 80, right: 80 },
    children: [new Paragraph({ alignment: i === 0 || (o.left || []).includes(i) ? AlignmentType.LEFT : AlignmentType.RIGHT,
      spacing: { after: 0, line: 240 }, children: [new TextRun({ text: String(t), bold, size: 17, font: 'Calibri' })] })],
  })) });
  return new Table({ columnWidths: widths, width: { size: widths.reduce((a, b) => a + b, 0), type: WidthType.DXA },
    borders: { top: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' }, bottom: { style: BorderStyle.SINGLE, size: 2, color: 'C9D4E4' },
      left: { style: BorderStyle.NONE }, right: { style: BorderStyle.NONE },
      insideHorizontal: { style: BorderStyle.SINGLE, size: 1, color: 'E4E7EC' }, insideVertical: { style: BorderStyle.NONE } },
    rows: [mk(headers, { bg: HEAD_BG, bold: true }), ...rows.map((r, k) => mk(r.cells || r, { bg: k % 2 ? ALT_BG : undefined, fills: r.fills }))] });
}

const ARM = { RF: 'random forest', ET: 'ExtraTrees', NN: 'neural network (GPU)', L1: 'L1 logistic', DT: 'single decision tree',
  FIGS: 'FIGS (sum of small trees)', B2: 'B2 — released model (XGBoost/logistic)', B2XL: 'B2XL — XGBoost + logistic',
  PB1: 'PB1 — 852-KO panel', PB2: 'PB2 — full KO profile', PB3: 'PB3 — KO + modules', PB4: 'PB4 — multi-label MLP',
  null_taxNN: 'taxonomic null (lineage only)', null_prevalence: 'prevalence null' };
const done = Object.entries(T).filter(([, b]) => b.status === 'done');
const pending = Object.entries(T).filter(([, b]) => b.status !== 'done').map(([k]) => k);
const TN = { da: 'Donors & acceptors', ox: 'Oxygen use', pw: 'Pathways' };

const C = [];
C.push(new Paragraph({ spacing: { after: 60 }, children: [new TextRun({ text: 'Interpretable models, XAI and gene–trait correlations',
  bold: true, size: 38, font: 'Calibri', color: '1A1A1A' })] }));
C.push(P('Decision tree, random forest and alternatives over KO profiles: what can be predicted, which genes drive it, and how rules compare with the best model', { size: 23, color: MUTED, after: 60 }));
C.push(P(`Generated from results/trees_xai/ by src/report/make_trees_xai_report.js. Pre-registration: docs/preregistration_trees_xai.md (with Addenda ${ADDENDA[0]}–${ADDENDA[ADDENDA.length - 1]} and the correction to Addendum F). Tasks complete: ${done.map(([k]) => TN[k]).join(', ') || 'none'}${pending.length ? `; still training: ${pending.map(k => TN[k]).join(', ')}` : ''}.`, { size: 18, color: MUTED, after: 200 }));

// ------------------------------------------------------------------ 1. summary
C.push(H('1. Summary', HeadingLevel.HEADING_1));
for (const [k, b] of done) {
  const g = b.macro.filter(r => r.split_level === 'genus');
  const get = a => g.find(r => r.arm === a);
  const best = get(b.best_model), dt = get('DT'), nul = get('null_taxNN');
  C.push(rb([{ t: `${TN[k]}: `, b: true }, { t: `best model ${ARM[b.best_model] || b.best_model}, macro AUPRC ${f(best.macro_auprc, 3)} on held-out genera (lineage-only null ${f(nul.macro_auprc, 3)}); ${b.n_predictable_best ?? '–'} of ${b.per_trait_best?.length ?? '–'} traits pass the predictability gate on held-out families.${dt ? ` A single decision tree reaches ${f(dt.macro_auprc, 3)}.` : ''}` }]));
}
const pr = R.predictions;
const cnt = v => pr.filter(p => p.preregistered !== false && p.verdict === v).length;
C.push(rb([{ t: 'Pre-registered predictions: ', b: true }, { t: `${cnt('supported')} supported, ${cnt('not supported')} not supported, ${cnt('pending')} pending (section 7; verdicts are yes/no as specified; post-hoc rows marked * are not counted).` }]));
C.push(spacer());

// ------------------------------------------------------------------ 2. methods
C.push(H('2. What was done', HeadingLevel.HEADING_1));
C.push(table(['Step', 'How'], [
  ['Models', 'Single decision tree (depth 2–6), random forest (500 trees), ExtraTrees, L1 logistic, FIGS; a neural network on the Apple GPU (MPS); compared with the project’s existing best models'],
  ['Validation', 'Nested cross-validation on frozen folds with whole genera (and separately whole families) held out; hyper-parameters chosen on inner folds only'],
  ['Cut-off', 'Per trait, the probability that maximises MCC on inner-fold predictions — never tuned on test genomes'],
  ['"Predictable"', 'On held-out families: MCC 95 % lower bound > 0.2 and better than the lineage-only null at BH q < 0.05'],
  ['XAI', 'Random forest: SHAP screen, then permutation importance on the top 50 + known markers; tree: rules scored on held-out genomes; network: Expected Gradients (post hoc, see section 8). Stability = share of folds in the top 10'],
  ['Independent checks', 'Model-free CMH test per gene, stratified by phylum × genome size; Boruta all-relevant selection; agreement with known marker genes'],
  ['Rules vs best model', 'Tree rules, textbook marker rules and a depth-3 surrogate of the best model, all on the same held-out genomes'],
], [2000, 7360], { left: [1] }));
C.push(spacer());

// ------------------------------------------------------------------ 3. prediction
C.push(H('3. Can it be predicted, and at what cut-off?', HeadingLevel.HEADING_1));
for (const [k, b] of done) {
  C.push(H(`3${'abc'[['da', 'ox', 'pw'].indexOf(k)]}. ${b.name}`, HeadingLevel.HEADING_2));
  const g = b.macro.filter(r => r.split_level === 'genus'), fm = b.macro.filter(r => r.split_level === 'family');
  C.push(table(['Model', 'Macro AUPRC, genus [95 % CI]', 'Family', 'Macro MCC', 'Predictable traits'],
    g.map(r => { const fr = fm.find(x => x.arm === r.arm) || {};
      return [ARM[r.arm] || r.arm, `${f(r.macro_auprc, 3)} [${f(r.macro_auprc_ci_low, 3)}, ${f(r.macro_auprc_ci_high, 3)}]`, f(fr.macro_auprc, 3), f(r.macro_mcc), r.arm.startsWith('null') ? '–' : `${fr.n_predictable ?? '–'}`]; }),
    [3300, 2400, 1000, 1100, 1560]));
  C.push(spacer());
  if (b.per_trait_best) {
    C.push(P(`Per trait — ${ARM[b.best_model] || b.best_model}, held-out families:`, { bold: true, after: 80 }));
    C.push(table(['Trait', 'Positives', 'AUPRC', 'Lineage null', 'MCC [low]', 'Cut-off', 'Predictable'],
      b.per_trait_best.map(r => ({ cells: [nice(r.trait), r.n_pos, f(r.auprc_family), f(r.null_auprc), `${f(r.mcc)} [${f(r.mcc_ci_low)}]`, f(r.threshold), r.predictable ? 'yes' : 'no'],
        fills: [0, 0, 0, 0, 0, 0, r.predictable ? OK_BG : NO_BG].map(x => x || undefined) })),
      [3000, 900, 900, 1100, 1300, 900, 1260]));
    C.push(P('The cut-off differs between traits, so one global threshold such as 0.5 would be wrong for most of them.', { size: 18, color: MUTED }));
  }
}
for (const k of pending) C.push(P(`${TN[k]}: models still training — this section will be filled when they finish.`, { italics: true, color: MUTED }));

// ------------------------------------------------------------------ 4. XAI
C.push(H('4. Which genes drive the predictions (XAI)', HeadingLevel.HEADING_1));
C.push(P('For each predictable trait, the genes that the random forest relies on in at least 60 % of held-out folds. Direction: + means carrying the gene pushes the prediction towards the trait. "Known" marks genes already tied to the trait by KEGG modules, curated marker groups or the reference’s functional groups — the check that the explanation recovers real biology.'));
for (const [k, b] of done) {
  const st = b.xai_rf_stable || [];
  if (!st.length) { C.push(P(`${TN[k]}: forest XAI still running.`, { italics: true, color: MUTED })); continue; }
  const pred = new Set((b.per_trait_best || []).filter(r => r.predictable).map(r => r.trait));
  const rows = st.filter(r => pred.has(r.trait)).slice(0, 60);
  C.push(P(`${TN[k]}:`, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Gene', 'Function (KEGG)', 'Stability', 'Direction', 'Known'],
    rows.map(r => [nice(r.trait), r.KO, (r.name || '').slice(0, 55), pct(r.stability), (r.shap_direction ?? 0) > 0 ? '+' : '−', r.known_marker ? 'yes' : 'new']),
    [2100, 800, 3800, 900, 800, 960], { left: [1, 2] }));
  C.push(spacer());
  const en = b.enrichment_rf || [];
  if (en.length) {
    C.push(P('Do the explanations recover known biology? Known markers among each trait’s stable genes vs chance:', { after: 60 }));
    C.push(table(['Trait', 'Stable genes', 'Known among them', 'Expected by chance', 'q'],
      en.filter(r => pred.has(r.trait)).map(r => ({ cells: [nice(r.trait), r.n_stable_top, r.n_stable_top_known, f(r.expected_by_chance), q(r.q_bh)],
        fills: [0, 0, 0, 0, r.q_bh < 0.05 ? OK_BG : undefined].map(x => x || undefined) })),
      [3400, 1300, 1600, 1700, 1360]));
    C.push(spacer());
  }
}

// Rashomon caveat: equally good models, different genes
{
  const median = a => { if (!a.length) return null; const s = [...a].sort((x, y) => x - y), m = Math.floor(s.length / 2);
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
  const ov = (T.da && T.da.rf_nn_overlap) || [];
  const med = median(ov.map(r => r.share_of_nn_in_rf));
  const medZ = median(((T.da && T.da.rf_nn_overlap_zero) || []).map(r => r.share_of_nn_in_rf));
  const nnRow = T.da && T.da.macro.find(r => r.arm === 'NN' && r.split_level === 'genus');
  const rfRow = T.da && T.da.macro.find(r => r.arm === 'RF' && r.split_level === 'genus');
  if (nnRow && rfRow && med !== null) {
    C.push(H('Equally good models, different explanations', HeadingLevel.HEADING_2));
    C.push(P(`The neural network predicts donors and acceptors as well as the forest (macro AUPRC ${f(nnRow.macro_auprc, 3)} vs ${f(rfRow.macro_auprc, 3)}), yet the genes each relies on overlap little: a median ${pct(med)} of the network’s stable genes are also stable for the forest (Expected Gradients; ${pct(medZ)} with the pre-registered Integrated Gradients, whose all-absent baseline favours housekeeping genes). When many genes carry the same information, models that are equally accurate can explain themselves differently. A gene named by one model is therefore a lead; a gene–trait link is treated as established only when independent methods agree (section 5) or when a rule built from it holds up on new organisms (the rulebook).`));
    C.push(spacer());
  }
}

// interpretation of the donor/acceptor explanations (written after reading them; labelled as interpretation)
if (T.da && T.da.xai_rf_stable && T.da.xai_rf_stable.length) {
  C.push(H('Reading the donor/acceptor explanations (interpretation)', HeadingLevel.HEADING_2));
  C.push(P('The enrichment test counts a gene as "known" only if it was on a list fixed before the analysis. Where a trait fails that test, the genes were read one by one to see why. This reading is interpretation, not part of the pre-registered scoring, which stays as reported in section 7.', { size: 19, color: MUTED }));
  C.push(rb([{ t: 'Textbook genes, recovered: ', b: true }, { t: 'sulfate reduction — dsrB and the DsrMK complex; sulfur oxidation — SoxAX, SoxB and sulfide:quinone oxidoreductase; ammonia oxidation — amoB and amoC; nitrate — NO reductase norBC and narJ, narI (nitrate is not predictable, so read these with care).' }]));
  C.push(rb([{ t: 'hao is not used. ', b: true }, { t: 'The marker the preliminary analysis showed to be misleading (it fires in anaerobic sulfate reducers) is absent from the forest’s stable genes for ammonia oxidation.' }]));
  C.push(rb([{ t: 'Real biology outside the pre-set list: ', b: true }, { t: 'for the O₂ acceptor, the forest reasons mostly from the absence of anaerobic machinery (heterodisulfide reductase hdrA2/B2, the methanogen hydrogenase mvhD, CODH/ACS), plus heme O synthase (ctaB), which is required to build aa3/bo3 oxidases. For the CO₂ acceptor it uses formylmethanofuran dehydrogenase (methanogenesis) and CODH/ACS δ (Wood–Ljungdahl pathway). These are core genes of those pathways that the list did not contain.' }]));
  C.push(rb([{ t: 'Likely lineage proxies — lower trust: ', b: true }, { t: 'for sulfur reduction the stable genes (ppsA, serC) are not sulfur metabolism, and for H₂ oxidation, apart from hyaB, they are NADH dehydrogenase genes (complex I subunits nuoI and nuoM, and the non-electrogenic type-II NADH dehydrogenase), a manganese-dependent pyrophosphatase and fumarate hydratase. These explanations should not be read as gene–trait biology.' }]));
  C.push(spacer());
}

// ------------------------------------------------------------------ 5. correlations
C.push(H('5. Gene–trait correlations confirmed by independent methods', HeadingLevel.HEADING_1));
C.push(P('A model explanation is one view. The same question was asked without any model (CMH test per gene, stratified by phylum and genome size) and with Boruta, which keeps every gene that beats its own shuffled copy. A correlation is called robust only when the forest, the model-free test and Boruta agree, with the same direction.'));
const removed = done.filter(([, b]) => b.cmh_size && b.cmh_phylum).map(([, b]) => 1 - b.cmh_size.n_significant / b.cmh_phylum.n_significant);
C.push(P(`Why genome size is controlled: organisms that use O₂ carry about 40 % more genes than anaerobes even within one phylum, so almost any gene looks "enriched" in them. Controlling for it removed ${removed.length ? `${Math.round(100 * Math.min(...removed))}–${Math.round(100 * Math.max(...removed))} %` : 'part'} of the model-free associations (Addendum B).`, { size: 19, color: MUTED }));
const cmhRows = done.filter(([, b]) => b.cmh_size).map(([k, b]) => [TN[k], b.cmh_phylum?.n_significant ?? '–', b.cmh_size.n_significant, `${b.cmh_size.n_positive} / ${b.cmh_size.n_negative}`, b.boruta ? b.boruta.n_confirmed : 'not run']);
if (cmhRows.length) C.push(table(['Task', 'CMH significant (phylum only)', 'CMH significant (phylum × size)', '+ / −', 'Boruta confirmed'], cmhRows, [2200, 2000, 2100, 1300, 1760]));
C.push(spacer());
for (const [k, b] of done) {
  const cs = b.consensus || [];
  if (!cs.length) { C.push(P(`${TN[k]}: consensus pending (needs the forest XAI).`, { italics: true, color: MUTED })); continue; }
  const robust = cs.filter(r => r.status === 'robust');
  C.push(P(`${TN[k]} — ${robust.length} robust gene–trait correlations (all three methods agree):`, { bold: true, after: 60 }));
  if (robust.length) C.push(table(['Trait', 'Gene', 'Function (KEGG)', 'Direction', 'MH odds ratio', 'Also in NN'],
    robust.slice(0, 40).map(r => [nice(r.trait), r.KO, (r.name || '').slice(0, 55), r.direction > 0 ? '+' : '−', f(r.mh_odds_ratio, 1), r.in_nn_stable === undefined ? '–' : (r.in_nn_stable ? 'yes' : 'no')]),
    [2100, 800, 3900, 800, 900, 860], { left: [1, 2] }));
  C.push(spacer());
}

// ------------------------------------------------------------------ 6. rules vs best model
C.push(H('6. Rules compared with the best model', HeadingLevel.HEADING_1));
C.push(P('Three kinds of rules, scored on the same held-out genomes as the best model (MCC at each one’s cut-off). Δ > 0 means the best model is better; a CI that excludes 0 means the difference is real.'));
for (const [k, b] of done) {
  const rv = b.rules_vs_best || [];
  if (!rv.length) { C.push(P(`${TN[k]}: comparison pending.`, { italics: true, color: MUTED })); continue; }
  C.push(H(`${TN[k]} — best model: ${ARM[b.best_model] || b.best_model}`, HeadingLevel.HEADING_2));
  C.push(table(['Trait', 'Rule set', 'Best model MCC', 'Rule MCC', 'Δ [95 % CI]', 'Agreement κ'],
    rv.map(r => { const sig = r.delta_ci_low > 0 ? NO_BG : (r.delta_ci_high < 0 ? OK_BG : PART_BG);
      return { cells: [nice(r.trait), r.rule_set.replace(/ \(K.*\)$/, ''), f(r.best_mcc), f(r.rule_mcc), `${f(r.delta_mcc_best_minus_rule)} [${f(r.delta_ci_low)}, ${f(r.delta_ci_high)}]`, f(r.kappa_rule_vs_model)],
        fills: [0, 0, 0, 0, sig, 0].map(x => x || undefined) }; }),
    [2400, 1700, 1100, 900, 1900, 1360], { left: [1] }));
  C.push(P('Green: rule significantly better; yellow: tied; red: best model significantly better.', { size: 18, color: MUTED }));
  const sf = b.surrogate_fidelity || [];
  if (sf.length) {
    C.push(P('Can the best model be written as a few rules? A depth-3 tree trained to imitate it, scored on held-out genera:', { bold: true, after: 60 }));
    C.push(table(['Trait', 'Fidelity (agreement)', 'Fidelity κ', 'Model positives reproduced', 'Surrogate MCC vs truth'],
      sf.map(r => ({ cells: [nice(r.trait), pct(r.fidelity_accuracy), f(r.fidelity_kappa), pct(r.fidelity_on_model_positives), f(r.surrogate_mcc_vs_labels)],
        fills: [0, r.fidelity_accuracy >= 0.9 ? OK_BG : PART_BG, 0, 0, 0].map(x => x || undefined) })),
      [3000, 1600, 1200, 1900, 1660]));
    C.push(spacer());
  }
  const sr = (b.surrogate_rules || []).filter(r => (b.per_trait_best || []).some(p => p.trait === r.trait && p.predictable));
  if (sr.length) {
    C.push(P('The surrogate’s main positive rules for predictable traits (support = genomes the rule covers):', { bold: true, after: 60 }));
    C.push(table(['Trait', 'Rule', 'Support', 'Agrees with model', 'Precision vs truth'],
      sr.map(r => [nice(r.trait), r.rule, r.support, pct(r.agrees_with_model), pct(r.precision_vs_labels)]),
      [2200, 4200, 900, 1000, 1060], { left: [1] }));
    C.push(spacer());
  }
}

// rules the forest itself uses, and recurring single-tree rules
C.push(H('Rules the random forest itself uses', HeadingLevel.HEADING_2));
C.push(P('The first one or two gene questions that recur most across the forest’s 1 500 trees (SIRUS method). Full tables, with the external genomes, are in the rulebook (docs/rulebook_gene_rules_supplementary.docx, Table S3).'));
for (const [k, b] of done) {
  const fr = (b.forest_rules || []).filter(r => (b.per_trait_best || []).some(p => p.trait === r.trait && p.predictable));
  if (!fr.length) continue;
  C.push(table(['Trait', 'Forest rule', 'P held-out', 'Base rate', 'Trees', 'Stable'],
    fr.map(r => [nice(r.trait), r.rule.replace(/^IF /, '').replace(/ THEN .*/, '').replace(/ \(K\d{5}\)/g, ''), pct(r.p_heldout), pct(r.base_rate), pct(r.tree_frequency), r.stability]),
    [2400, 3300, 1000, 1000, 800, 860], { left: [1] }));
  C.push(spacer());
}
C.push(H('Rules of the single decision tree (recurring in ≥ 3 of 5 folds)', HeadingLevel.HEADING_2));
for (const [k, b] of done) {
  const dr = b.dt_rules_recurring || [];
  if (!dr.length) continue;
  C.push(P(`${TN[k]}:`, { bold: true, after: 60 }));
  C.push(table(['Trait', 'Rule', 'Folds', 'Held-out support', 'Precision'],
    dr.map(r => [nice(r.trait), r.rule.replace(/ present/g, '+').replace(/ absent/g, '−').replace(/ AND /g, ', '), `${r.folds}/5`, r.support, pct(r.precision)]),
    [2200, 4300, 800, 1100, 960], { left: [1] }));
  C.push(spacer());
}

// ------------------------------------------------------------------ 7. predictions
C.push(H('7. Pre-registered predictions, scored', HeadingLevel.HEADING_1));
C.push(table(['ID', 'Prediction', 'Verdict', 'Evidence'],
  pr.map(p => ({ cells: [p.id, p.prediction, p.verdict, p.evidence || ''], fills: [0, 0, VBG[p.verdict], 0].map(x => x || undefined) })),
  [900, 3700, 1200, 3560], { left: [1, 3] }));
C.push(spacer());

// ------------------------------------------------------------------ 8. deviations and limits
C.push(H('8. Changes to the plan and limits', HeadingLevel.HEADING_1));
C.push(bullet('Addendum B (post hoc, disclosed): the model-free CMH test is stratified by phylum × genome gene-count, because gene count confounds presence/absence even within a phylum. It can only remove associations.'));
C.push(bullet('Addendum E (before any forest XAI output): permutation importance over every gene would have cost ~9–250 CPU-hours per task; it is run on the top 50 SHAP-screened genes plus all known markers instead.'));
C.push(bullet('Neural-network bug (found and fixed): inputs were scaled by (sd + 10⁻⁶), so a gene absent from — or present in — every training genome of a fold became a ±10⁶ input on test genomes. It corrupted the network’s predictions and explanations; all network results were re-run after the fix (macro AUPRC 0.573 → 0.605). Tree, forest and logistic models were never affected.'));
C.push(bullet('Addendum F (post hoc, disclosed): even after the fix, the pre-registered Integrated Gradients baseline — a genome with no genes — puts the network’s attribution on near-universal housekeeping genes for every trait. Expected Gradients, which uses real training genomes as the reference like TreeSHAP does for the forest, returns trait-specific genes. Both are reported; predictions marked * are post hoc.'));
C.push(bullet('The neural network runs on the Apple GPU (MPS). Its dropout draws random numbers on the GPU, so a GPU run behaves like a different random seed of a CPU run (checked for the pathway network, docs/deviations.md #22: AUPRC difference ≈ 0.005).'));
C.push(bullet('Addendum L (after an independent audit): verdicts are scored yes/no as pre-registered (no “partly”); every clause is scored for every task it covers (T2 both clauses and for oxygen use, A-P2 for oxygen use, T4 on SHAP, A-P1 on the unadjusted CMH test, C-P1 on predictable traits); post-hoc rows (*) are listed but not counted.'));
C.push(bullet('A single decision tree ranks genomes poorly (coarse probabilities), which is why it fails the AUPRC-based gate even where its yes/no calls match the forest.'));
C.push(bullet('Labels are closed-world: a trait not recorded for an organism counts as absent, which lowers precision for under-reported traits.'));
if (pending.length) C.push(bullet(`Not yet included: ${pending.map(k => TN[k]).join(', ')} — models still training; this report is regenerated when they finish.`));

// ------------------------------------------------------------------ 9. files
C.push(H('9. Where the results are', HeadingLevel.HEADING_1));
C.push(table(['File (results/trees_xai/<task>/)', 'Content'], [
  ['oof/<model>.parquet', 'held-out predictions of every model'],
  ['metrics_macro.csv, metrics_per_trait.csv, paired.csv', 'performance, the predictability gate, model-vs-model contrasts'],
  ['xai_importance_<model>.csv, xai_stable_<model>.csv', 'per-fold gene importance; genes stable across folds'],
  ['xai_enrichment_<model>.csv', 'known-marker enrichment of each model’s stable genes'],
  ['xai_rules_dt.csv', 'every decision-tree rule with held-out support and precision'],
  ['alt_cmh.csv, alt_cmh_size_adjusted.csv, alt_boruta.csv', 'model-free and all-relevant gene–trait associations'],
  ['consensus.csv, rf_nn_overlap.csv', 'robust correlations; agreement between forest and neural network'],
  ['rules_vs_best.csv, surrogate_fidelity.csv, surrogate_rules.csv', 'rules compared with the best model'],
  ['../predictions_scored.csv', 'pre-registered predictions with verdicts'],
], [4200, 5160], { left: [1] }));

const doc = new Document({ creator: 'bio_pred', title: 'Interpretable models, XAI and gene–trait correlations',
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 360, hanging: 200 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [new TextRun({
      children: ['Page ', PageNumber.CURRENT, ' of ', PageNumber.TOTAL_PAGES], size: 18, color: MUTED, font: 'Calibri' })] })] }) },
    children: C }] });
Packer.toBuffer(doc).then(b => { fs.writeFileSync(`${ROOT}/docs/trees_xai_correlation_report.docx`, b); console.log('written docs/trees_xai_correlation_report.docx'); });
