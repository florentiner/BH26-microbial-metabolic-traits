// Correlation / theory report for the marker-KO set in jaccard_similarity.html.
// Every number is read from results/jaccard/report_data.json (built by src/44_correlation_summary.py).
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, Table, TableRow, TableCell,
  WidthType, ShadingType, AlignmentType, BorderStyle, LevelFormat, PageNumber, Footer,
} = require('docx');
const fs = require('fs');

const ROOT = '.';
const D = JSON.parse(fs.readFileSync(`${ROOT}/results/jaccard/report_data.json`, 'utf8'));

const ACCENT = '2A78D6';
const MUTED = '52514E';
const HEAD_BG = 'EDF2FA';
const ALT_BG = 'F7F8FA';
const GRADE_BG = { 'High': 'DCEFE0', 'Moderate': 'FFF3D6', 'Low': 'FBE3D6', 'Not supported': 'EEEEEE' };

const f = (x, n = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–' : Number(x).toFixed(n);
const pct = x => (x === null || x === undefined) ? '–' : `${Math.round(100 * x)} %`;
const q = x => (x === null || x === undefined) ? '–' : (x < 0.001 ? '< 0.001' : Number(x).toFixed(3));

const P = (text, o = {}) => new Paragraph({
  spacing: { after: o.after ?? 120, line: 276 },
  children: [new TextRun({ text, bold: o.bold, italics: o.italics, size: o.size ?? 21, color: o.color, font: 'Calibri' })],
});
const Rich = (runs, o = {}) => new Paragraph({
  spacing: { after: o.after ?? 120, line: 276 },
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i, size: 21, color: r.c, font: 'Calibri' })),
});
const H = (text, level) => new Paragraph({
  heading: level,
  spacing: { before: level === HeadingLevel.HEADING_1 ? 340 : 240, after: 130 },
  children: [new TextRun({ text, font: 'Calibri', bold: true,
    size: level === HeadingLevel.HEADING_1 ? 30 : 24, color: level === HeadingLevel.HEADING_1 ? ACCENT : '1A1A1A' })],
});
const bullet = (text, o = {}) => new Paragraph({
  numbering: { reference: 'bul', level: 0 }, spacing: { after: 70, line: 276 },
  children: [new TextRun({ text, size: 21, font: 'Calibri', bold: o.bold })],
});
const richBullet = runs => new Paragraph({
  numbering: { reference: 'bul', level: 0 }, spacing: { after: 70, line: 276 },
  children: runs.map(r => new TextRun({ text: r.t, bold: r.b, italics: r.i, size: 21, font: 'Calibri' })),
});

function table(headers, rows, widths, o = {}) {
  const mkRow = (cells, { bg, bold, fills } = {}) => new TableRow({
    children: cells.map((t, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: (fills && fills[i]) || bg ? { type: ShadingType.CLEAR, fill: (fills && fills[i]) || bg, color: 'auto' }
        : undefined,
      margins: { top: 60, bottom: 60, left: 90, right: 90 },
      children: [new Paragraph({
        alignment: i === 0 || (o.leftCols || []).includes(i) ? AlignmentType.LEFT : AlignmentType.RIGHT,
        spacing: { after: 0, line: 240 },
        children: [new TextRun({ text: String(t), bold, size: 18, font: 'Calibri' })],
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
      ...rows.map((r, k) => mkRow(r.cells || r, { bg: k % 2 ? ALT_BG : undefined, fills: r.fills }))],
  });
}
const spacer = () => new Paragraph({ spacing: { after: 100 }, children: [] });

const T = Object.fromEntries(D.theories.map(t => [t.id, t]));
const S = D.summary;
const gv = name => D.group_vs_trait.find(r => r.functional_group === name);
const sk = name => D.single_ko.find(r => r.functional_group === name);
const oxd = (m, p) => D.oxygen_deep_dive.find(r => r.marker_set === m && r.phenotype === p);

// ------------------------------------------------------------------ content
const children = [];
children.push(new Paragraph({ spacing: { after: 60 }, children: [new TextRun({
  text: 'Marker-gene correlations in the 586-genome KO set', bold: true, size: 38, font: 'Calibri', color: '1A1A1A' })] }));
children.push(P('Which correlations hold, how sure we are, and the data behind each theory', { size: 24, color: MUTED, after: 60 }));
children.push(P(`Source: jaccard_similarity.html (${S.n_genomes} genomes, ${S.n_marker_kos} marker KOs, ${S.n_groups} functional groups) joined to curated donor/acceptor labels and to Madin et al. 2020 oxygen requirement and habitat. Generated from results/jaccard/ by src/report/make_correlation_report.js.`, { size: 18, color: MUTED, after: 200 }));

// ---- 1. verdicts
children.push(H('1. Verdicts at a glance', HeadingLevel.HEADING_1));
children.push(P('Each note from the preliminary analysis was turned into a testable claim. Every confidence grade is computed by an explicit rule, never assigned by judgement; the rules are in section 8.', { after: 100 }));
const vrows = D.theories.map(t => ({
  cells: [t.claim, t.verdict, t.grade],
  fills: [undefined, undefined, GRADE_BG[t.grade]],
}));
children.push(table(['Claim', 'Verdict', 'Confidence'], vrows, [3300, 4600, 1460], { leftCols: [1] }));
children.push(spacer());

// ---- 2. how to read the numbers
children.push(H('2. How to read the numbers', HeadingLevel.HEADING_1));
children.push(P('The genomes are phylogenetically clustered: several genomes of one genus that share a marker and a label are close to one observation, not several. Counting each genome separately makes any lineage-linked marker look far more certain than it is. Two levels of lineage control are therefore used, and they answer different questions:'));
children.push(richBullet([{ t: 'Genus-level test ', b: true }, { t: '(trait shuffled between whole genera). Removes repeated counting within a genus. It answers "does the marker detect the trait?" — the right test for claims 1 and 3, where a marker confined to one lineage is still a valid detector.' }]));
children.push(richBullet([{ t: 'Within-phylum test ', b: true }, { t: '(trait shuffled only among genera of the same phylum). Also removes confounding by phylum and domain. It answers "is there a functional link, or do both simply belong to the same lineage?" — the right test for enzyme ↔ environment correlations (claim 2, section 6e) and for co-occurrence (6c). It is conservative: a link that is real but confined to one phylum (methanogens and anoxia) cannot pass it, so failing means "not separable from lineage in this data", not "false".' }]));
children.push(table(['Metric', 'What it measures', 'Good value'], [
  ['Precision', 'Of the genomes the marker flags, the share that truly have the trait', 'high = few false alarms'],
  ['Recall', 'Of the genomes that truly have the trait, the share the marker catches', 'high = few misses'],
  ['MCC (per genome)', 'Balanced agreement of marker and trait, −1…+1; 0 = no better than chance', '> 0.5 strong'],
  ['MCC (per genus)', 'Same, one genome per genus (median of 200 draws); removes within-genus replication', '> 0.5 strong'],
  ['Genus-permutation p', 'Trait shuffled between whole genera: removes within-genus replication', '< 0.05'],
  ['Within-phylum p', 'Trait shuffled only within each phylum: removes confounding by phylum / domain', '< 0.05'],
  ['q (BH)', 'p corrected for testing many groups at once (Benjamini-Hochberg)', '< 0.05'],
], [2000, 5360, 2000], { leftCols: [1] }));
children.push(spacer());

// ---- 3. claim 1
const amo = gv('Ammonia oxidizing (AOB/AOA)'), nob = gv('Nitrite oxidizing (NOB)');
const kamo = sk('Ammonia oxidizing (AOB/AOA)'), knob = sk('Nitrite oxidizing (NOB)');
children.push(H('3. Claim 1 — ammonia-oxidation and nitrite-oxidation enzymes are predicted accurately', HeadingLevel.HEADING_1));
children.push(Rich([{ t: 'Verdict: ', b: true }, { t: `true for ammonia oxidation with the right marker; not true for nitrite oxidation. ` },
  { t: 'Both groups catch essentially every true oxidiser (high recall), which is what an inspection of known organisms sees. Precision — the share of flagged genomes that are real — is the part that was not visible.' }]));
const kn1 = T.C1a.key_numbers;
children.push(table(['Marker rule', 'Flagged', 'With trait', 'Precision', 'Recall', 'MCC genome', 'MCC genus', 'q'], [
  ['Ammonia: 4-KO group', amo.n_marker_present, amo.n_phenotype_pos, f(amo.precision), f(amo.recall), f(amo.mcc_per_genome), f(amo.mcc_genus_level), q(amo.q_bh_permutation)],
  [`Ammonia: amoB (${kamo.most_chosen_ko}) alone`, kn1.amoB_flagged, amo.n_phenotype_pos, f(kn1.amoB_precision), f(kn1.amoB_recall), f(kn1.amoB_mcc_per_genome), f(kn1.amoB_mcc_genus), q(kn1.amoB_q)],
  ['Nitrite: 5-KO group', nob.n_marker_present, nob.n_phenotype_pos, f(nob.precision), f(nob.recall), f(nob.mcc_per_genome), f(nob.mcc_genus_level), q(nob.q_bh_permutation)],
], [2700, 800, 650, 1000, 850, 1060, 1000, 1300]));
children.push(P(`amoB was also chosen on training genera only and scored on held-out genera: MCC ${f(kamo.mcc_single_ko_out_of_genus)} per genome, the same KO chosen in ${pct(kamo.selection_stability)} of folds. For nitrite oxidation no single KO does better out-of-genus (MCC ${f(knob.mcc_single_ko_out_of_genus)}).`, { size: 19, color: MUTED }));
children.push(H('Why it is true for ammonia oxidation', HeadingLevel.HEADING_2));
children.push(bullet(`amoB (${kamo.most_chosen_ko}) alone reaches MCC ${f(kamo.mcc_single_ko_out_of_genus)} on genera never used to choose it, and it was chosen in ${pct(kamo.selection_stability)} of cross-validation folds. The choice is a property of the biology, not of the sample.`));
children.push(bullet('The 4-KO group rule is weaker almost entirely because of hydroxylamine oxidoreductase, hao (K10535): among the labelled genomes it is detected in 60, 44 of which are not ammonia oxidisers — mostly sulfate reducers (Desulfosarcina, Thermodesulfobacterium, Desulfotignum, Desulfatitalea) and, across all 586 genomes, hydrothermal-vent anaerobes (Thermovibrio, Desulfurobacterium, Geoglobus, Nautilia). These carry HAO-like octaheme cytochromes that the KOfam hao model also matches. amoB and amoC make only 5 false calls each.'));
const amoOx = D.group_vs_oxygen_top.find(r => r.functional_group === 'Ammonia oxidizing (AOB/AOA)' && r.phenotype === 'oxygen:anaerobe');
if (amoOx) {
  children.push(bullet(`Independent confirmation from a second data source: the 4-KO ammonia group correlates with being an anaerobe (MCC ${f(amoOx.mcc_genus_level)} at genus level, odds ratio ${f(amoOx.odds_ratio, 1)}, q ${q(amoOx.q_bh_permutation)}), although ammonia oxidisers are obligate aerobes. The sign is biologically backwards, and it is exactly what the hao false positives in (anaerobic) sulfate reducers predict.`));
}
children.push(bullet('Two smaller definition effects: amoA/B/C are one KOfam family shared with particulate methane monooxygenase, so 3 verrucomicrobial methanotrophs (Methylacidimicrobium) are flagged; and amoA (K10944) is under-detected — found in 9 genomes against 32 for amoB — so a rule requiring amoA would miss most ammonia oxidisers.'));
children.push(H('Why it is not true for nitrite oxidation', HeadingLevel.HEADING_2));
children.push(bullet(`The group flags ${nob.n_marker_present} genomes for ${nob.n_phenotype_pos} true nitrite oxidisers (precision ${pct(nob.precision)}).`));
const nobDrop = D.shared_ko_test.find(r => r.functional_group === 'Nitrite oxidizing (NOB)' && r.variant !== 'all KOs');
children.push(bullet(`KOfam annotates nitrite oxidoreductase (nxrAB) with the same KOs as respiratory nitrate reductase (narGH, K00370/K00371). Removing those two KOs does not fix precision — it collapses recall to ${f(nobDrop.recall)} — because there is no NOB-specific KO. The ambiguity is in the annotation system, not the data.`));
children.push(bullet('What would work: placing the NxrA/NarG sequence on a reference tree (Nitrospira/Nitrospina/Nitrobacter NxrA clades are distinct from NarG), or requiring the periplasmic NxrA signal. Presence/absence of a KO cannot do it.'));

// ---- 4. claim 2
const all13 = oxd('Aerobic respiration (all 13 KOs)', 'Madin O2-capable');
const tru = oxd('aa3 or cbb3 (true respiratory oxidases)', 'Madin O2-capable');
const aa3 = oxd('aa3 cytochrome c oxidase only', 'Madin O2-capable');
const bd = oxd('cytochrome bd only', 'Madin O2-capable');
children.push(H('4. Claim 2 — aerobic-respiration enzymes correlate with oxygen level', HeadingLevel.HEADING_1));
children.push(Rich([{ t: 'Verdict: supported, strongly — once cytochrome bd is separated out. ', b: true },
  { t: `The hope in the note holds. Scored against the recorded oxygen requirement (Madin et al.: aerobic, facultative or microaerophilic vs anaerobic; n = ${all13.n} genomes), the true respiratory oxidases give MCC ${f(tru.mcc_genus_level)} at genus level (q ${q(tru.q_bh_permutation)}).` }]));
children.push(table(['Marker set', 'KOs', 'Precision', 'Recall', 'MCC genus', 'Odds ratio', 'q'], [
  ['All 13 "aerobic respiration" KOs', all13.n_kos, f(all13.precision), f(all13.recall), f(all13.mcc_genus_level), f(all13.odds_ratio, 0), q(all13.q_bh_permutation)],
  ['aa3 + cbb3 oxidases (coxABC, ccoNOP)', tru.n_kos, f(tru.precision), f(tru.recall), f(tru.mcc_genus_level), f(tru.odds_ratio, 0), q(tru.q_bh_permutation)],
  ['aa3 cytochrome c oxidase only', aa3.n_kos, f(aa3.precision), f(aa3.recall), f(aa3.mcc_genus_level), f(aa3.odds_ratio, 0), q(aa3.q_bh_permutation)],
  ['Cytochrome bd only (cydAB)', bd.n_kos, f(bd.precision), f(bd.recall), f(bd.mcc_genus_level), f(bd.odds_ratio, 1), q(bd.q_bh_permutation)],
], [3500, 600, 1000, 900, 1100, 1060, 1200]));
const envf = (grp, cond) => D.env_functional.find(r => r.functional_group === grp && r.condition === cond);
const e1 = envf('aa3 + cbb3 oxidases (respiratory O2 reduction)', 'oxygen: uses O2 (aerobe/facultative/microaerophile)');
if (e1) {
  children.push(Rich([{ t: 'It survives the strictest control. ', b: true }, { t: `Collapsing to one point per genus and shuffling oxygen use only among genera of the same phylum, aa3 + cbb3 oxidases still give φ ${f(e1.statistic)} (within-phylum q ${q(e1.q_bh_within_phylum)}): ${pct(e1.trait_rate_marker_pos)} of genera carrying them use O₂, against ${pct(e1.trait_rate_marker_neg)} of genera without them. The link is functional, not a matter of which phylum an organism belongs to. The 13-KO group as defined does not survive this test (section 6e) — cytochrome bd and the non-respiratory members blur it enough that its signal cannot be separated from lineage.` }]));
}
children.push(H('Why it is true — and why cytochrome bd dilutes it', HeadingLevel.HEADING_2));
const prev = D.oxidase_prevalence;
const pv = ko => prev.find(r => r.KO === ko);
children.push(bullet(`aa3-type cytochrome c oxidase (coxB) is in ${pct(pv('K02275').prevalence_in_O2_capable)} of oxygen-using organisms and ${pct(pv('K02275').prevalence_in_anaerobes)} of anaerobes — almost a clean switch. It is the canonical low-affinity terminal oxidase of aerobic respiration.`));
children.push(bullet(`Cytochrome bd (cydA) is in ${pct(pv('K00425').prevalence_in_O2_capable)} of oxygen users but also ${pct(pv('K00425').prevalence_in_anaerobes)} of anaerobes. It is a high-affinity oxidase that many strict anaerobes keep to scavenge and detoxify traces of O₂, not to respire. Including it in an "aerobic respiration" group adds false positives without adding information.`));
children.push(bullet('cbb3 oxidase sits in between: it is the microaerobic oxidase, so it tracks oxygen use but less sharply than aa3.'));
children.push(bullet(`The recorded oxygen requirement and the curated O₂-acceptor label agree at Cohen's κ ${f(T.C2.key_numbers.label_kappa)} (substantial). On the ${all13.n} genomes where both exist, the correlation holds against either source.`));

// ---- 5. claim 3
const h2 = gv('Hydrogen oxidizing'), kh2 = sk('Hydrogen oxidizing');
children.push(H('5. Claim 3 — the hydrogen-oxidation definition is insufficient', HeadingLevel.HEADING_1));
children.push(Rich([{ t: 'Verdict: supported. ', b: true },
  { t: `With ${h2.n_kos_in_group} KOs the group catches ${pct(h2.recall)} of known H₂ oxidisers (MCC ${f(h2.mcc_per_genome)} per genome, ${f(h2.mcc_genus_level)} per genus) — the weakest of the well-defined groups, and no single KO rescues it (best single KO chosen out-of-genus: MCC ${f(kh2.mcc_single_ko_out_of_genus)} per genome).` }]));
children.push(H('Why it is true', HeadingLevel.HEADING_2));
children.push(bullet('Hydrogen oxidation is carried by a very diverse enzyme family. The 4 KOs cover two [NiFe] group 1 large subunits (hyaB, hybC) and the [FeFe] hydA pair. Missing are, at least: the uptake hydrogenase small subunits (hyaA/hybO), group 1 hupSL of rhizobia and cyanobacteria, the NAD-reducing group 3d hoxHY of Knallgas bacteria (Cupriavidus), group 2 and group 3b/3c enzymes, and the membrane-bound respiratory complexes (group 4 / Ech).'));
children.push(bullet(`The misses are definitional, not biological: ${T.C3.key_numbers.n_misses} known H₂ oxidisers are missed across ${T.C3.key_numbers.n_genera} genera, ${T.C3.key_numbers.n_families} families and ${T.C3.key_numbers.n_phyla} phyla (no single genus above 9 %), and they are not explained by genome quality (section 7).`));
children.push(bullet('The largest single block tells exactly which enzymes are missing: 40 of the 94 misses (43 %) are hydrogenotrophic methanogens (Methanobacteriota). They oxidise H₂ with the F₄₂₀-reducing hydrogenase FrhABG ([NiFe] group 3a), MvhADG (group 3c) and the energy-converting Ech/Eha/Ehb complexes (group 4) — none of which is in the 4-KO list. The next block, 26 Bacillota, points to the electron-bifurcating [FeFe] HydABC complexes beyond hydA.'));
children.push(bullet('Recommended fix: classify hydrogenase large subunits with HydDB (Søndergaard et al. 2016) rather than a KO list; it assigns the group and therefore the direction (uptake vs evolving).'));

// ---- 6. new correlations
children.push(H('6. New correlations found', HeadingLevel.HEADING_1));
children.push(H('6a. The group definitions overlap — and two contain enzymes of another physiology', HeadingLevel.HEADING_2));
children.push(P(`${D.n_group_pairs_sharing} pairs of functional groups share KOs, so the groups are not independent tests. amoABC (K10944–K10946) alone sits in four groups. Two memberships are not about the named physiology at all:`));
children.push(richBullet([{ t: '"Nitrite oxidizing (NOB)" contains amoABC, ', b: true }, { t: 'ammonia monooxygenase. Classic nitrite oxidisers do not carry it (only comammox Nitrospira does both), so every ammonia oxidiser also lights up the NOB group.' }]));
children.push(richBullet([{ t: '"Aerobic respiration" contains five enzymes that are not respiratory: ', b: true }, { t: 'amoABC, soluble methane monooxygenase mmoX (K16157) and CO dehydrogenase coxL (K03520). They use O₂ as a substrate, not as the terminal electron acceptor; only 8 of the 13 KOs are respiratory oxidases.' }]));
children.push(P('What each problem costs, scored against the curated labels:', { after: 80 }));
children.push(table(['Group', 'Variant', 'Flagged', 'FP', 'Precision', 'Recall', 'MCC'],
  D.definition_variants.map(r => [r.group, r.variant, r.flagged, r.fp, f(r.precision), f(r.recall), f(r.mcc)]),
  [1900, 3400, 750, 650, 900, 800, 960], { leftCols: [1] }));
children.push(spacer());
const dv = n => D.definition_variants.find(r => r.variant.startsWith(n));
children.push(P(`Reading: in the NOB group amoABC adds ${dv('amoABC only').fp} false positives for ${dv('amoABC only').flagged - dv('amoABC only').fp} true ones, but most of the damage comes from narGH (${dv('narGH').fp} false positives); removing amoABC lifts precision only from ${f(dv('as defined (amoABC').precision)} to ${f(dv('narGH').precision)}. In aerobic respiration the five non-respiratory enzymes are conceptually misplaced but numerically harmless (MCC ${f(dv('as defined (13').mcc)} → ${f(dv('respiratory oxidases only').mcc)} without them); the real diluter is cytochrome bd (aa3 + cbb3 alone: MCC ${f(dv('aa3 + cbb3').mcc)}), the same conclusion section 4 reached from a different phenotype source. Full overlap list: results/jaccard/group_definition_overlaps.csv.`, { after: 160 }));
children.push(H('6b. Single-marker detectors that beat their own group', HeadingLevel.HEADING_2));
children.push(P('For these groups, one specific KO predicts the trait far better than "any KO of the group present". The KO is chosen on training genera only and scored on held-out genera, 20 repeats; stability is how often the same KO was chosen.'));
const robust = D.single_ko.filter(r => r.verdict === 'single KO better (robust)' && r.gain_out_of_genus >= 0.1)
  .sort((a, b) => b.gain_out_of_genus - a.gain_out_of_genus);
children.push(table(['Functional group', 'Group rule MCC', 'Single KO', 'Out-of-genus MCC', 'Gain', 'Stability'],
  robust.map(r => [r.functional_group, f(r.mcc_any_ko_rule), r.most_chosen_ko, f(r.mcc_single_ko_out_of_genus), `+${f(r.gain_out_of_genus)}`, pct(r.selection_stability)]),
  [3200, 1200, 1100, 1500, 1060, 1300]));
children.push(P('Reading: mcrA (K00399) for methanogens, dsrB (K11181) for sulfate reducers, amoB (K10945) for ammonia oxidisers and acsA/cooS (K00198) for acetogens are each the textbook diagnostic gene for that physiology, so the data recover the known biology. The reverse also holds: for hydrogen oxidation, denitrification, sulfur reduction and CO oxidation, the multi-KO rule is better, because those physiologies are carried by many alternative enzymes.', { after: 160 }));

children.push(H('6c. Physiologies that co-occur in the same genome beyond phylum identity', HeadingLevel.HEADING_2));
children.push(P(`Pairs of functional groups, one data point per genus, tested with the within-phylum null. ${D.cooc_beyond_phylum.length} pairs co-occur beyond phylum identity; ${D.cooc_lineage_bound} more are significant only across genera; ${D.cooc_artefact.length} involve groups that share KOs. The pairs with a clear biological reading:`));
const READ = new Set(['Sulfur oxidizing (SOB)|Aerobic respiration', 'Aerobic respiration|Nitrate-reducing (denitrifiers)', 'CO oxidizing (carboxydotroph)|Formate oxidizing', 'Hydrogen oxidizing|Sulfur-reducing', 'Hydrogen oxidizing|Aerobic respiration', 'Sulfate-reducing (SRB)|Sulfur-reducing', 'Sulfur oxidizing (SOB)|Nitrate-reducing (denitrifiers)']);
children.push(table(['Group A', 'Group B', 'Genera with both', 'φ', 'q within phylum', 'Biological reading'],
  D.cooc_beyond_phylum.filter(r => READ.has(`${r.group_a}|${r.group_b}`)).map(r => [r.group_a, r.group_b, r.n_genera_both, f(r.phi_genus_level), q(r.q_within_phylum),
    ({
      'Sulfur oxidizing (SOB)|Aerobic respiration': 'aerobic sulfur lithotrophs',
      'Aerobic respiration|Nitrate-reducing (denitrifiers)': 'facultative denitrifiers switch O₂ → NO₃⁻',
      'CO oxidizing (carboxydotroph)|Formate oxidizing': 'carboxydotrophs (CODH + formate dehydrogenase)',
      'Hydrogen oxidizing|Sulfur-reducing': 'H₂-driven sulfur respiration',
      'Sulfate-reducing (SRB)|Sulfur-reducing': 'dissimilatory sulfur-compound respiration',
      'Sulfur oxidizing (SOB)|Nitrate-reducing (denitrifiers)': 'sulfur-driven denitrification',
      'Hydrogen oxidizing|Aerobic respiration': 'Knallgas (aerobic H₂-oxidising) bacteria',
      'Hydrogen oxidizing|Sulfate-reducing (SRB)': 'H₂ as the classic SRB electron donor',
      'Sulfur-reducing|Methanogens (hydrogenotrophic)': 'shared H₂-driven anaerobic niches',
    })[`${r.group_a}|${r.group_b}`] || '–']),
  [2050, 2050, 1000, 700, 900, 2660], { leftCols: [1, 5] }));
children.push(P(`Artefacts (groups share KOs): ${D.cooc_artefact.slice(0, 5).map(r => `${r.group_a} + ${r.group_b} (φ ${f(r.phi_genus_level)})`).join('; ')}. NOB + denitrifiers is driven entirely by narGH; sulfur oxidation + sulfate reduction partly by dsrAB/aprAB, which reverse-Dsr sulfur oxidisers genuinely run in the oxidative direction.`, { size: 19, color: MUTED, after: 100 }));
const hm = D.hdr_check.find(r => r.genes.startsWith('mcr')), hh = D.hdr_check.find(r => r.genes.startsWith('hdr'));
children.push(Rich([{ t: 'Two surviving pairs are artefacts of non-specific genes, and one of them reverses sign. ', b: true },
  { t: `"Sulfate reducers + methanogens" co-occur only because the methanogen group includes heterodisulfide reductase hdrABC, which sulfate reducers also use (${hh.n_genera} genera carry it, ${hh.n_genera_also_dsrAB} of them with dsrAB). The real methyl-CoM reductase mcrABG ${hm.n_genera_also_dsrAB ? 'rarely' : 'never'} co-occurs with dsrAB in this data: ${hm.n_genera_also_dsrAB} of ${hm.n_genera} genera, φ ${f(hm.phi_with_dsrAB)} — methanogenesis and sulfate respiration are almost never combined in one organism. "Ammonia oxidisers + DNRA" is driven entirely by hao (φ ${f(D.hao_dnra_phi)}; amoB: φ ${f(D.amob_dnra_phi)}), the same hao artefact as in sections 3 and 6d.` }], { after: 160 }));

children.push(H('6d. Oxygen and habitat', HeadingLevel.HEADING_2));
children.push(P(`Across all functional groups, ${D.n_oxygen_sig} of ${D.n_oxygen_tests} group × oxygen-class tests and ${D.n_habitat_sig} of ${D.n_habitat_tests} group × habitat tests survive the genus-level test at q < 0.05. This section uses the genus-level test only; section 6e applies the within-phylum test to oxygen and the other conditions.`));
if (D.group_vs_oxygen_top.length) {
  children.push(table(['Functional group', 'Oxygen class', 'Precision', 'Recall', 'MCC genus', 'q'],
    D.group_vs_oxygen_top.slice(0, 8).map(r => [r.functional_group, r.phenotype.replace('oxygen:', ''), f(r.precision), f(r.recall), f(r.mcc_genus_level), q(r.q_bh_permutation)]),
    [3600, 1500, 1000, 900, 1160, 1200]));
  children.push(spacer());
}
if (D.group_vs_habitat_top.length) {
  children.push(table(['Functional group', 'Habitat', 'Precision', 'Recall', 'MCC genus', 'q'],
    D.group_vs_habitat_top.slice(0, 8).map(r => [r.functional_group, r.phenotype.replace('habitat:', ''), f(r.precision), f(r.recall), f(r.mcc_genus_level), q(r.q_bh_permutation)]),
    [3600, 1500, 1000, 900, 1160, 1200]));
  const hd = D.habitat_decomposition;
  const hg = Object.entries(D.habitat_top_genera).map(([g, n]) => `${g} ×${n}`).join(', ');
  children.push(P(`The single habitat hit is not a habitat effect. Split by marker, the whole marine signal comes from hao-only genomes (MCC ${f(hd[2].mcc_vs_marine)}), while the genuine ammonia monooxygenase genes show none (amoABC only: MCC ${f(hd[1].mcc_vs_marine)}). Apart from genuine marine ammonia oxidisers (Nitrosococcus), the flagged marine genomes are hydrothermal-vent anaerobes (top genera: ${hg}). This is the same hao defect seen in sections 3 and 4, now from a third, independent data source.`, { after: 120 }));
  children.push(table(['Ammonia marker', 'Flagged', 'Marine among flagged', 'MCC vs marine'],
    hd.map(r => [r.marker, r.flagged, `${r.marine_among_flagged} (${pct(r.share_marine)})`, f(r.mcc_vs_marine)]),
    [3600, 1600, 2400, 1760]));
  children.push(spacer());
  children.push(P('Beyond that artefact, no functional group × habitat association survives the genus-permutation null. Habitat as recorded (isolation source) is dominated by lineage: once whole genera are shuffled, nothing remains. This is expected — isolation source records where a strain was first isolated, not the conditions it grows under.', { after: 120 }));
} else {
  children.push(P('No functional group × habitat association survives the genus-permutation null. Habitat, as recorded (isolation source), is dominated by lineage: once whole genera are shuffled, nothing remains. This is expected — isolation source records where a strain was first cultured, not the conditions it grows under.', { after: 120 }));
}

children.push(H('6e. Enzyme ↔ environment correlations (the type of correlation in the preliminary notes)', HeadingLevel.HEADING_2));
const cnt = (fam, ms, v) => (D.env_counts.find(r => r.family === fam && r.marker_set === ms && r.verdict.startsWith(v)) || { n: 0 }).n;
children.push(P(`Every functional group — and ten refined markers from sections 3–6b — tested against every environmental condition Madin et al. record with usable coverage: oxygen use, strict anaerobiosis, growth and optimum temperature, thermophily, optimum pH and doubling time (salinity has only 35 genomes and was not tested). Cell traits (Gram stain, sporulation, motility) were a separate family. ${D.n_env_tests} tests in all, one data point per genus.`));
children.push(table(['Marker family', 'Tests', 'Functional link (survives within phylum)', 'Lineage-bound only'], [
  ['Groups as defined — environment', cnt('environment', 'group as defined', 'functional') + cnt('environment', 'group as defined', 'lineage') + cnt('environment', 'group as defined', 'not'), cnt('environment', 'group as defined', 'functional'), cnt('environment', 'group as defined', 'lineage')],
  ['Refined markers — environment', cnt('environment', 'refined marker', 'functional') + cnt('environment', 'refined marker', 'lineage') + cnt('environment', 'refined marker', 'not'), cnt('environment', 'refined marker', 'functional'), cnt('environment', 'refined marker', 'lineage')],
  ['Groups as defined — cell traits', cnt('cell trait', 'group as defined', 'functional') + cnt('cell trait', 'group as defined', 'lineage') + cnt('cell trait', 'group as defined', 'not'), cnt('cell trait', 'group as defined', 'functional'), cnt('cell trait', 'group as defined', 'lineage')],
  ['Refined markers — cell traits', cnt('cell trait', 'refined marker', 'functional') + cnt('cell trait', 'refined marker', 'lineage') + cnt('cell trait', 'refined marker', 'not'), cnt('cell trait', 'refined marker', 'functional'), cnt('cell trait', 'refined marker', 'lineage')],
], [3400, 1000, 2900, 2060]));
children.push(spacer());
children.push(P('Functional links — enzyme and condition are associated even among genera of the same phylum:', { bold: true, after: 80 }));
const envRead = {
  'aa3 + cbb3 oxidases (respiratory O2 reduction)': 'respiratory O₂ reduction',
  'aa3 cytochrome c oxidase only': 'canonical aerobic respiration',
  'nirK or nirS (NO-forming nitrite reductase)': 'denitrifiers are facultative aerobes',
  'narGH / nxrAB only': 'nitrate respiration in facultative organisms',
  'dsrB only (dissimilatory sulfite reductase)': 'sulfate reducers are anaerobes',
  'hao only (hydroxylamine oxidoreductase)': 'ARTEFACT: HAO-like cytochromes of anaerobes',
  'Nitrate-reducing (denitrifiers)': 'denitrifiers are facultative aerobes',
  'Nitrite oxidizing (NOB)': 'driven by narGH (see above)',
  'Ammonia oxidizing (AOB/AOA)': 'ARTEFACT: hao — ammonia oxidisers are aerobes',
  'Nitrate-reducing (DNRA)': 'DNRA favoured in strongly reducing settings',
};
const cond = c => c.replace('oxygen: uses O2 (aerobe/facultative/microaerophile)', 'uses O₂').replace('oxygen: strict anaerobe', 'strict anaerobe');
const seenEnv = new Set();
const envRows = D.env_functional.filter(r => { const k = r.functional_group; if (seenEnv.has(k)) return false; seenEnv.add(k); return true; });
children.push(table(['Marker', 'Condition', 'φ / ρ', 'q within phylum', 'Reading'],
  envRows.map(r => [r.functional_group, cond(r.condition), f(r.statistic), q(r.q_bh_within_phylum), envRead[r.functional_group] || '–']),
  [3100, 1300, 800, 1200, 2960], { leftCols: [1, 4] }));
children.push(P('One row per marker (its strongest condition). "Uses O₂" and "strict anaerobe" are near-complements, so each marker usually appears significant for both with opposite sign.', { size: 18, color: MUTED }));
children.push(P('Significant across genera but not separable from lineage — consistent with biology, not proven by this data:', { bold: true, after: 80 }));
const lb = D.env_lineage_bound.filter(r => ['aa3 + cbb3 oxidases (respiratory O2 reduction)', 'mcrA only (methyl-CoM reductase)', 'Methanogens (hydrogenotrophic)', 'Acetogens (homoacetogens)', 'dsrB only (dissimilatory sulfite reductase)', 'Aerobic respiration'].includes(r.functional_group)).slice(0, 7);
children.push(table(['Marker', 'Condition', 'φ / ρ', 'Marker+ vs marker−', 'q within phylum'],
  lb.map(r => [r.functional_group, cond(r.condition), f(r.statistic),
    r.kind === 'numeric' ? `${f(r.median_marker_pos, 1)} vs ${f(r.median_marker_neg, 1)}` : `${pct(r.trait_rate_marker_pos)} vs ${pct(r.trait_rate_marker_neg)}`, q(r.q_bh_within_phylum)]),
  [3300, 2000, 900, 1860, 1300], { leftCols: [1] }));
children.push(spacer());
children.push(bullet('Temperature: genera with aa3 + cbb3 oxidases grow at a median 30 °C against 52 °C without them — consistent with O₂ being far less soluble in hot water, but hyperthermophiles are concentrated in a few archaeal and bacterial phyla, so the data cannot separate the two.'));
children.push(bullet('Methanogenesis and anoxia: mcrA is found only in strict anaerobes, but every mcrA-carrying genus is in the same few archaeal lineages, so there is no within-phylum contrast to test. Real biology, unidentifiable here.'));
children.push(bullet('Acetogens and temperature (31 °C vs 69 °C) is an artefact of the group definition. It includes acetate kinase ack (K00925; in 67 % of bacteria here, 5 % of archaea) and the Rnf complex rnfACDEG (about 40 % of bacteria, 0 % of archaea) — energy-conservation genes used by many organisms that are not acetogens. With them the group is present in 90 % of bacteria and 40 % of archaea; without them, 29 % and 24 %. So "no acetogen marker" mostly means "archaeon", and hyperthermophiles here are largely archaea.'));
children.push(bullet('No correlation with optimum pH, doubling time or any cell trait survives the within-phylum test.'));

// ---- 7. causes
children.push(H('7. The three suggested causes of error, tested', HeadingLevel.HEADING_1));
children.push(table(['Suggested cause', 'Test', 'Result', 'Confidence'], [
  { cells: ['(a) KO definition incomplete or too broad', 'Group rule vs single KO, out-of-genus; removing KOs shared by two groups', T.Ca.verdict, T.Ca.grade], fills: [0, 0, 0, GRADE_BG[T.Ca.grade]] },
  { cells: ['(b) Curated values wrong (strain, "experimental" share)', 'Agreement between two independent oxygen sources', T.Cb.verdict, T.Cb.grade], fills: [0, 0, 0, GRADE_BG[T.Cb.grade]] },
  { cells: ['(c) Poor genome quality blocks detection', 'CheckM2 completeness / contamination of missed vs detected positives', T.Cc.verdict, T.Cc.grade], fills: [0, 0, 0, GRADE_BG[T.Cc.grade]] },
].map(r => ({ cells: r.cells, fills: r.fills.map(x => x || undefined) })), [2500, 2700, 2900, 1260], { leftCols: [1, 2] }));
children.push(spacer());
children.push(P(`Cause (a) explains most of what was seen. Cause (b) is real but modest: κ ≈ 0.7 between two curations means roughly one genome in seven is labelled differently by the two sources, which caps any correlation below 1. Cause (c) is unlikely to be the main explanation here because most genomes are near-complete reference assemblies (${T.Cc.verdict.replace(/^Not supported in this dataset \(/, '').replace(/\)$/, '')}); it will matter for metagenome-assembled genomes, which the pathway study tests separately by degrading proteomes to 90 / 70 / 50 %.`));

// ---- 8. confidence rule
children.push(H('8. How confidence was graded', HeadingLevel.HEADING_1));
children.push(P('Association claims ("marker X tracks trait Y": claims 1 and 2) are graded by the rule below, which was fixed before any grade was computed. The marker’s own genus-level MCC and genus-permutation q are used — for claim 1, amoB’s.', { after: 80 }));
children.push(table(['Grade', 'Criteria (all must hold)'], [
  { cells: ['High', 'Genus-level MCC ≥ 0.5; genus-permutation q < 0.01; where a KO was selected, the choice survives out-of-genus validation with ≥ 80 % stability; a known biological mechanism'], fills: [GRADE_BG.High] },
  { cells: ['Moderate', 'q < 0.05 and genus-level MCC 0.3–0.5; or robust out-of-genus but an unstable KO choice'], fills: [GRADE_BG.Moderate] },
  { cells: ['Low', 'q < 0.05 but genus-level MCC < 0.3; or the signal rests on KOs shared between two groups'], fills: [GRADE_BG.Low] },
  { cells: ['Not supported', 'q ≥ 0.05, or the effect disappears at genus level'], fills: [GRADE_BG['Not supported']] },
].map(r => ({ cells: r.cells, fills: [r.fills[0], undefined] })), [1700, 7660], { leftCols: [1] }));
children.push(spacer());
children.push(P('Four claims are not associations and need their own rules. These were written during this analysis — after the data had been seen, though before each grade was computed — and are applied mechanically:', { after: 80 }));
children.push(table(['Claim', 'High', 'Moderate', 'Not supported'], [
  ['Definition insufficient (claim 3)', 'recall < 0.7, misses in ≥ 10 genera, and one lineage with a nameable missing enzyme family holds ≥ 25 % of misses', 'recall < 0.7 only', 'recall ≥ 0.7'],
  ['(a) over-broad definitions', '≥ 3 groups where one KO beats the group by ≥ 0.3 MCC out-of-genus with ≥ 80 % stability', '1–2 such groups', 'none'],
  ['(b) curated values wrong', 'κ < 0.6 between two curations', '0.6 ≤ κ < 0.8', 'κ ≥ 0.8'],
  ['(c) poor genome quality', '—', 'missed positives significantly less complete (q < 0.05) in ≥ 1 group', 'otherwise; a difference in the other direction does not count'],
], [2300, 3100, 2000, 1960], { leftCols: [1, 2, 3] }));
children.push(spacer());

// ---- 9. recommendations
children.push(H('9. What to change in the marker set', HeadingLevel.HEADING_1));
children.push(richBullet([{ t: 'Use the diagnostic subunit, not "any KO": ', b: true }, { t: 'amoB for ammonia oxidation, mcrA for methanogenesis, dsrB for sulfate reduction, acsA/cooS for acetogenesis.' }]));
children.push(richBullet([{ t: 'Drop hao (K10535) as an ammonia-oxidation marker, ', b: true }, { t: 'or require it together with amoB: on its own it mostly detects HAO-like octaheme cytochromes of sulfate reducers and hydrothermal-vent anaerobes.' }]));
children.push(richBullet([{ t: 'Remove amoABC from the NOB group ', b: true }, { t: '(or keep it only as a separate "comammox" call), and remove amoABC, mmoX and coxL from "aerobic respiration" — they are O₂-using enzymes, not respiratory oxidases.' }]));
children.push(richBullet([{ t: 'Remove hdrABC from the methanogen group, ', b: true }, { t: 'and likewise ack and the Rnf complex from the acetogen group — or require them together with mcrA or acsB. On their own, hdrABC is carried by sulfate reducers and ack / Rnf by most bacteria, whatever their metabolism.' }]));
children.push(richBullet([{ t: 'Make groups disjoint or mark shared KOs explicitly, ', b: true }, { t: 'so that co-occurrence between groups is never built into the definitions.' }]));
children.push(richBullet([{ t: 'Split aerobic respiration: ', b: true }, { t: 'aa3 + cbb3 oxidases as the oxygen-respiration marker; report cytochrome bd separately as "oxygen tolerance", not respiration.' }]));
children.push(richBullet([{ t: 'Replace the hydrogenase KO list with HydDB classification ', b: true }, { t: 'so that group (and direction) is assigned from sequence.' }]));
children.push(richBullet([{ t: 'Resolve NxrA vs NarG by sequence placement, ', b: true }, { t: 'since KOfam cannot separate them; until then, report nitrite oxidation as "not determinable from KO presence".' }]));
children.push(richBullet([{ t: 'Keep multi-KO rules ', b: true }, { t: 'for denitrification, sulfur reduction and CO oxidation, where alternative enzymes do the same job.' }]));

// ---- 10. caveats
children.push(H('10. Limits of this analysis', HeadingLevel.HEADING_1));
children.push(bullet(`Phenotypes come from two curations; ${S.n_with_curated_labels} genomes have curated donor/acceptor labels and ${S.n_with_oxygen} have a recorded oxygen requirement. Species are matched by name, so strain-level mix-ups (cause b) cannot be excluded.`));
children.push(bullet('The within-phylum test is conservative. Links confined to one lineage (methanogenesis and anoxia) cannot pass it, and phyla with few genera add little power, so "lineage-bound" results are unresolved rather than refuted.'));
children.push(bullet('Closed-world labels: a trait absent from a record is treated as absent from the organism. This lowers precision for traits that are under-reported.'));
children.push(bullet('The habitat field is the isolation source, which records where a strain was first isolated rather than the conditions it lives in; habitat results should be read in that light.'));
children.push(bullet('Correction during the analysis: the first habitat/oxygen scan paired genomes by position, and the HTML stores genome metadata in a different order from the marker table (584 of 586 positions differ). All frames are now joined by genome ID, with an assertion that fails if they are not aligned; every number in this report comes from the corrected run.'));

// ------------------------------------------------------------------ document
const doc = new Document({
  creator: 'bio_pred',
  title: 'Marker-gene correlations — theory report',
  numbering: { config: [{ reference: 'bul', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•',
    alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 360, hanging: 200 } } } }] }] },
  styles: { default: { document: { run: { font: 'Calibri', size: 21 } } } },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 },
      margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.RIGHT,
      children: [new TextRun({ children: ['Page ', PageNumber.CURRENT, ' of ', PageNumber.TOTAL_PAGES],
        size: 18, color: MUTED, font: 'Calibri' })] })] }) },
    children,
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync(`${ROOT}/docs/marker_correlation_theory_report.docx`, b);
  console.log('written docs/marker_correlation_theory_report.docx');
});
