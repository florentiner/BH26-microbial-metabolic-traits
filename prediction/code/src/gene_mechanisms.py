"""What each gene does, for explaining rules.

Use: interpretation and grounding labels of rules. It does NOT choose rule genes in the rulebook (51) or forest rules
(53), but PATHWAY plus the lifestyle markers ARE the candidate set of the universal-rule search (56, Addendum J), and
some markers (korA/C/D, hcp, K04070, aceE, coxL, cooS) were added to this table after they appeared in XAI or rulebook
outputs on overlapping genomes - disclosed in Addendum L (audit M5).

MECH: KO -> (short gene name, one-line mechanism). PATHWAY: trait -> KOs that are part of that physiology, used to
flag rule genes that are *not* on the trait's pathway (possible lineage proxies). Sources: KEGG orthology and module
definitions, standard microbial physiology.
"""
MECH = {
    # dissimilatory sulfate / sulfite reduction
    "K00958": ("sat", "ATP sulfurylase: activates sulfate to APS, the entry step of sulfate respiration (also used for assimilation)"),
    "K00394": ("aprA", "APS reductase α: reduces APS to sulfite in sulfate respiration"),
    "K00395": ("aprB", "APS reductase β: electron-transfer subunit of APS reductase"),
    "K11180": ("dsrA", "dissimilatory sulfite reductase α: reduces sulfite to sulfide, the core of sulfate respiration"),
    "K11181": ("dsrB", "dissimilatory sulfite reductase β: partner of DsrA in the sulfite-reducing core"),
    "K27187": ("dsrM", "DsrMKJOP membrane complex: recycles DsrC and couples sulfite reduction to the quinone pool"),
    "K27188": ("dsrK", "DsrMKJOP membrane complex, catalytic subunit K: reduces the DsrC trisulfide"),
    "K16950": ("asrA", "anaerobic sulfite reductase A: reduces sulfite to sulfide (Firmicutes-type)"),
    "K16951": ("asrB", "anaerobic sulfite reductase B"),
    "K16952": ("asrC", "anaerobic sulfite reductase C"),
    "K08352": ("psrA", "polysulfide / thiosulfate reductase catalytic subunit: respires sulfur compounds"),
    "K17219": ("sreA", "sulfur reductase molybdopterin subunit: respires elemental sulfur"),
    # sulfur oxidation
    "K17218": ("sqr", "sulfide:quinone oxidoreductase: oxidises sulfide and feeds electrons to the quinone pool"),
    "K17222": ("soxA", "SoxAX cytochrome: attaches thiosulfate to the SoxYZ carrier - first step of the Sox system"),
    "K17223": ("soxX", "SoxAX cytochrome partner of SoxA"),
    "K17224": ("soxB", "SoxB hydrolase: releases sulfate from the SoxYZ carrier"),
    "K17225": ("soxC", "SoxCD sulfane dehydrogenase: oxidises the sulfane sulfur on SoxYZ"),
    "K17226": ("soxY", "SoxYZ carrier protein of the Sox thiosulfate-oxidising system"),
    "K17227": ("soxZ", "SoxYZ carrier protein of the Sox thiosulfate-oxidising system"),
    "K19713": ("tsdA", "thiosulfate dehydrogenase: oxidises thiosulfate to tetrathionate"),
    # nitrogen
    "K10944": ("amoA", "ammonia monooxygenase A (same KO family as methane monooxygenase pmoA)"),
    "K10945": ("amoB", "ammonia monooxygenase B: first enzyme of ammonia oxidation"),
    "K10946": ("amoC", "ammonia monooxygenase C"),
    "K10535": ("hao", "hydroxylamine oxidoreductase; HAO-like octaheme cytochromes also occur in anaerobes"),
    "K00370": ("narG", "respiratory nitrate reductase α (also nitrite oxidoreductase nxrA of nitrite oxidisers)"),
    "K00371": ("narH", "respiratory nitrate reductase β (also nxrB)"),
    "K00373": ("narJ", "chaperone that assembles the respiratory nitrate reductase"),
    "K00374": ("narI", "respiratory nitrate reductase γ: membrane anchor linking NarGH to quinols"),
    "K02567": ("napA", "periplasmic nitrate reductase: reduces nitrate outside the cytoplasm"),
    "K02568": ("napB", "periplasmic nitrate reductase electron-transfer subunit"),
    "K00368": ("nirK", "copper nitrite reductase: nitrite to NO, the committed step of denitrification"),
    "K15864": ("nirS", "cd1 nitrite reductase: nitrite to NO in denitrification"),
    "K04561": ("norB", "nitric oxide reductase B: NO to N2O in denitrification"),
    "K02305": ("norC", "nitric oxide reductase C"),
    "K00376": ("nosZ", "nitrous oxide reductase: N2O to N2, the last step of denitrification"),
    "K03385": ("nrfA", "cytochrome c nitrite reductase: nitrite to ammonium in respiratory DNRA"),
    "K15876": ("nrfH", "NrfH quinol dehydrogenase partner of NrfA"),
    "K00362": ("nirB", "NAD(P)H nitrite reductase (assimilatory / fermentative nitrite to ammonium)"),
    "K00363": ("nirD", "NAD(P)H nitrite reductase small subunit"),
    "K02588": ("nifH", "nitrogenase iron protein: nitrogen fixation"),
    # methanogenesis and Wood-Ljungdahl
    "K00399": ("mcrA", "methyl-CoM reductase α: final, methane-forming step of methanogenesis (also run in reverse by anaerobic methane oxidisers)"),
    "K00401": ("mcrB", "methyl-CoM reductase β"),
    "K00402": ("mcrG", "methyl-CoM reductase γ"),
    "K11261": ("fwdE", "formylmethanofuran dehydrogenase: first step of CO2 reduction in methanogenesis"),
    "K03388": ("hdrA2", "heterodisulfide reductase A: flavin-based electron bifurcation (methanogens; also sulfate reducers)"),
    "K03389": ("hdrB2", "heterodisulfide reductase B"),
    "K03390": ("hdrC2", "heterodisulfide reductase C"),
    "K14127": ("mvhD", "F420-non-reducing hydrogenase subunit D of methanogens"),
    "K14138": ("acsB", "acetyl-CoA synthase: joins CO, a methyl group and CoA in the Wood-Ljungdahl pathway"),
    "K00192": ("cdhA", "CODH/ACS α (archaeal Wood-Ljungdahl complex)"),
    "K00193": ("cdhC", "CODH/ACS β"),
    "K00194": ("cdhD", "CODH/ACS δ: corrinoid methyltransferase of the Wood-Ljungdahl complex"),
    "K00197": ("cdhE", "CODH/ACS γ: corrinoid iron-sulfur protein of the Wood-Ljungdahl complex"),
    "K00198": ("cooS", "anaerobic CO dehydrogenase: oxidises CO / reduces CO2 (Wood-Ljungdahl, carboxydotrophs)"),
    # oxygen
    "K02274": ("coxA", "aa3 cytochrome c oxidase subunit I: reduces O2 to water - low-affinity aerobic respiration"),
    "K02275": ("coxB", "aa3 cytochrome c oxidase subunit II"),
    "K02276": ("coxC", "aa3 cytochrome c oxidase subunit III"),
    "K00404": ("ccoN", "cbb3 cytochrome oxidase: high-affinity O2 respiration under microaerobic conditions"),
    "K00405": ("ccoO", "cbb3 cytochrome oxidase subunit II"),
    "K00406": ("ccoP", "cbb3 cytochrome oxidase subunit III"),
    "K00425": ("cydA", "cytochrome bd: high-affinity O2 scavenging; also kept by anaerobes for O2 protection"),
    "K00426": ("cydB", "cytochrome bd subunit II"),
    "K02257": ("ctaB", "heme O synthase: makes heme o, required to build heme-copper oxidases"),
    "K02259": ("ctaA", "heme A synthase: makes heme a for aa3 oxidase"),
    "K05919": ("dfx", "superoxide reductase: detoxifies O2 radicals in anaerobes"),
    "K00627": ("pdhC", "pyruvate dehydrogenase E2: aerobic pyruvate oxidation"),
    # hydrogen, CO, formate, photosynthesis, C1
    "K06281": ("hyaB", "[NiFe] uptake hydrogenase large subunit: oxidises H2"),
    "K05922": ("hybC", "[NiFe] hydrogenase-2 large subunit: oxidises H2"),
    "K03518": ("coxS", "aerobic CO dehydrogenase small subunit: carboxydotrophy"),
    "K03519": ("coxM", "aerobic CO dehydrogenase medium subunit"),
    "K03520": ("coxL", "aerobic CO dehydrogenase large subunit"),
    "K00123": ("fdoG", "formate dehydrogenase α: formate as electron donor"),
    "K00124": ("fdoH", "formate dehydrogenase iron-sulfur subunit"),
    "K00127": ("fdoI", "formate dehydrogenase cytochrome b subunit"),
    "K08928": ("pufL", "type-II reaction centre L: anoxygenic photosynthesis"),
    "K08929": ("pufM", "type-II reaction centre M: anoxygenic photosynthesis"),
    "K02689": ("psaA", "photosystem I: oxygenic photosynthesis"),
    "K02703": ("psbA", "photosystem II D1: oxygenic photosynthesis"),
    "K14028": ("mxaF", "methanol dehydrogenase: methylotrophy"),
    "K16157": ("mmoX", "soluble methane monooxygenase: methanotrophy"),
}

MECH.update({
    "K27190": ("dsrO", "DsrMKJOP membrane complex subunit O: part of the sulfate-respiration electron-transfer module"),
    "K17229": ("fccB", "flavocytochrome c sulfide dehydrogenase: oxidises sulfide in sulfur oxidisers"),
    "K27925": ("tth", "tetrathionate hydrolase: sulfur oxidation in acidophilic sulfur oxidisers (S4I pathway)"),
    "K00386": ("sorB", "sulfite dehydrogenase (cytochrome) subunit: oxidises sulfite in sulfur oxidisers"),
    "K02575": ("narK", "NNP nitrate/nitrite transporter: imports nitrate (for respiration or assimilation)"),
    "K00174": ("korA", "2-oxoglutarate:ferredoxin oxidoreductase α: O2-sensitive enzyme of anaerobic central metabolism"),
    "K00176": ("korD", "2-oxoglutarate:ferredoxin oxidoreductase δ: O2-sensitive; common in anaerobes (the four-subunit enzyme also occurs in some microaerophiles)"),
    "K00177": ("korC", "2-oxoglutarate:ferredoxin oxidoreductase γ: O2-sensitive; common in anaerobes (the four-subunit enzyme also occurs in some microaerophiles)"),
    "K00410": ("petB/C", "cytochrome bc1 complex (complex III): feeds electrons to cytochrome c and c-type oxidases"),
    "K00413": ("petC", "cytochrome c1 of the bc1 complex (complex III) of respiratory chains"),
    "K05601": ("hcp", "hybrid-cluster protein (hydroxylamine reductase): made under anoxia; found in anaerobes and facultative anaerobes"),
    "K07306": ("dmsA", "anaerobic DMSO reductase: anaerobic respiration"),
    "K00381": ("cysI", "assimilatory sulfite reductase: sulfur assimilation, not respiration"),
})
# added after the oxygen-use rulebook was first built (Addendum H: changes interpretation labels only)
MECH.update({
    "K02297": ("cyoA", "cytochrome bo3 ubiquinol oxidase subunit II: a terminal O2 reductase of aerobic respiration"),
    "K02298": ("cyoB", "cytochrome bo3 ubiquinol oxidase subunit I"),
    "K00163": ("aceE", "pyruvate dehydrogenase E1: aerobic pyruvate oxidation (anaerobes use O2-sensitive PFOR instead)"),
    "K04070": ("pflA-like", "putative pyruvate formate-lyase activating enzyme: PFL is destroyed by O2, so it suggests anaerobic fermentation"),
})
# genes that indicate an aerobic or anaerobic lifestyle rather than one specific pathway
LIFESTYLE = {"K00174", "K00176", "K00177", "K00410", "K00413", "K05601", "K07306", "K05919", "K00627", "K02257",
             "K02259", "K03388", "K03389", "K03390", "K14127",
             "K00163", "K04070", "K03520", "K00198"}   # aceE, PFL activase, aerobic / anaerobic CO dehydrogenase

_SULF_RED = {"K00958", "K00394", "K00395", "K11180", "K11181", "K27187", "K27188", "K27190", "K16950", "K16951",
             "K16952"}
_SULF_OX = {"K17218", "K17222", "K17223", "K17224", "K17225", "K17226", "K17227", "K19713", "K11180", "K11181",
            "K00394", "K00395", "K00958", "K17229", "K27925", "K00386"}
_DENIT = {"K00370", "K00371", "K00373", "K00374", "K02567", "K02568", "K00368", "K15864", "K04561", "K02305", "K00376",
          "K02575"}
_O2 = {"K02274", "K02275", "K02276", "K00404", "K00405", "K00406", "K00425", "K00426", "K02257", "K02259",
       "K02297", "K02298"}
_METH = {"K00399", "K00401", "K00402", "K11261", "K03388", "K03389", "K03390", "K14127"}
_WL = {"K14138", "K00192", "K00193", "K00194", "K00197", "K00198"}
PATHWAY = {
    "acceptor_sulfate": _SULF_RED, "sulfate_reduction": _SULF_RED,
    "acceptor_sulfur": {"K08352", "K17219", "K16950", "K16951", "K16952"}, "sulfur_reduction": {"K08352", "K17219"},
    "donor_sulfur_oxidation": _SULF_OX, "sulfur_oxidation": _SULF_OX,
    "donor_ammonia_oxidation": {"K10944", "K10945", "K10946", "K10535"}, "nitrification": {"K10944", "K10945", "K10946"},
    "donor_nitrite_oxidation": {"K00370", "K00371"},
    "acceptor_nitrate": _DENIT | {"K03385", "K15876"}, "denitrification": _DENIT, "nitrate_reduction": _DENIT,
    "dnra": {"K03385", "K15876", "K00362", "K00363"},
    "acceptor_oxygen": _O2, "uses_O2": _O2 | {"K00627"}, "aerobic_chemoheterotrophy": _O2,
    "acceptor_CO2_methanogen_acetogen": _METH | _WL, "methanogenesis": _METH, "acetogenesis": _WL,
    "donor_hydrogen_oxidation": {"K06281", "K05922"}, "hydrogen_oxidation": {"K06281", "K05922"},
    "donor_CO_oxidation": {"K03518", "K03519", "K03520", "K00198"},
    "nitrogen_fixation": {"K02588"}, "photoautotrophy": {"K08928", "K08929", "K02689", "K02703"},
    "photoheterotrophy": {"K08928", "K08929"}, "methylotrophy": {"K14028"}, "methanotrophy": {"K16157"},
}
# genes whose ABSENCE is informative, with what their absence means
ABSENCE_MEANING = {
    "K03388": "no methanogen/sulfate-reducer electron-bifurcation machinery", "K03389": "no heterodisulfide reductase",
    "K14127": "no methanogen hydrogenase", "K00197": "no Wood-Ljungdahl CODH/ACS", "K00194": "no Wood-Ljungdahl CODH/ACS",
    "K02274": "no aa3 oxidase", "K02275": "no aa3 oxidase", "K02276": "no aa3 oxidase", "K00404": "no cbb3 oxidase",
    "K02257": "cannot make heme o, so no heme-copper oxidase", "K05919": "no anaerobe-type superoxide defence",
    "K10535": "no hydroxylamine oxidoreductase", "K00399": "no methyl-CoM reductase, so not a methanogen",
}

# lifestyle markers split by what they indicate (used to explain rules; also the candidate set of src/56, Addendum J)
ANAEROBIC_MARKERS = {"K00174", "K00176", "K00177", "K05601", "K07306", "K05919", "K03388", "K03389", "K03390",
                     "K14127", "K04070", "K00198"}
AEROBIC_MARKERS = {"K00410", "K00413", "K00627", "K00163", "K02257", "K02259", "K03520"}
# the oxygen relationship each trait implies
TRAIT_LIFESTYLE = {
    "acceptor_sulfate": "anaerobic", "acceptor_CO2_methanogen_acetogen": "anaerobic", "acceptor_sulfur": "anaerobic",
    "acceptor_iron": "anaerobic", "acceptor_oxygen": "aerobic", "uses_O2": "aerobic",
    "acceptor_nitrate": "facultative", "donor_ammonia_oxidation": "aerobic", "donor_nitrite_oxidation": "aerobic",
    "donor_sulfur_oxidation": "mostly aerobic", "donor_hydrogen_oxidation": "either", "donor_CO_oxidation": "either",
    "donor_iron_oxidation": "mostly aerobic",
}

# markers for the pathway task, each KEGG name checked against ko_list before adding. Deliberately NOT added:
# K00244/K00245 (KEGG does not separate fumarate reductase from succinate dehydrogenase) and K17050/K17051
# (selenate/chlorate reductase, not arsenate reductase) - fumarate and arsenate reduction have no clean KO marker.
MECH.update({
    "K02586": ("nifD", "nitrogenase MoFe protein α: nitrogen fixation"),
    "K02591": ("nifK", "nitrogenase MoFe protein β: nitrogen fixation"),
    "K02690": ("psaB", "photosystem I P700 apoprotein A2: oxygenic photosynthesis"),
    "K02706": ("psbD", "photosystem II D2: oxygenic photosynthesis"),
    "K01601": ("rbcL", "RuBisCO large subunit: CO2 fixation in the Calvin cycle"),
    "K00855": ("prkB", "phosphoribulokinase: the Calvin-cycle enzyme that regenerates the RuBisCO substrate"),
    "K23995": ("xoxF", "lanthanide-dependent methanol dehydrogenase: methanol oxidation"),
    "K14029": ("mxaI", "methanol dehydrogenase small subunit: methanol oxidation"),
    "K06282": ("hyaA", "[NiFe] uptake hydrogenase small subunit: H2 oxidation"),
    "K00436": ("hoxH", "NAD-reducing [NiFe] hydrogenase large subunit: H2 oxidation (Knallgas bacteria)"),
    "K18007": ("hoxY", "NAD-reducing [NiFe] hydrogenase small subunit"),
    "K00656": ("pflB", "pyruvate formate-lyase: splits pyruvate to acetyl-CoA and formate in fermentation"),
    "K04072": ("adhE", "bifunctional acetaldehyde/alcohol dehydrogenase: ethanol fermentation"),
    "K00016": ("ldh", "L-lactate dehydrogenase: lactic fermentation"),
    "K03778": ("ldhA", "D-lactate dehydrogenase: lactic fermentation"),
    "K00925": ("ack", "acetate kinase: makes ATP from acetyl phosphate in fermentation"),
    "K00929": ("buk", "butyrate kinase: butyrate fermentation"),
})
PATHWAY["nitrogen_fixation"] = {"K02588", "K02586", "K02591"}
PATHWAY["photoautotrophy"] = {"K08928", "K08929", "K02689", "K02690", "K02703", "K02706", "K01601", "K00855"}
PATHWAY["methylotrophy"] = {"K14028", "K14029", "K23995"}
PATHWAY["hydrogen_oxidation"] = {"K06281", "K05922", "K06282", "K00436", "K18007"}
PATHWAY["donor_hydrogen_oxidation"] = PATHWAY["hydrogen_oxidation"]
PATHWAY["fermentation"] = {"K00656", "K04072", "K00016", "K03778", "K00925", "K00929", "K04070"}
TRAIT_LIFESTYLE.update({"sulfate_reduction": "anaerobic", "methanogenesis": "anaerobic", "acetogenesis": "anaerobic",
                        "fermentation": "either", "denitrification": "facultative", "nitrate_reduction": "facultative",
                        "dnra": "either", "nitrogen_fixation": "either", "hydrogen_oxidation": "either",
                        "photoautotrophy": "either", "photoheterotrophy": "either", "methylotrophy": "mostly aerobic",
                        "aerobic_chemoheterotrophy": "aerobic", "sulfur_oxidation": "mostly aerobic",
                        "sulfur_reduction": "anaerobic"})


# ---- direction-aware support of one rule condition (audit C5) ---------------------------------------------------
# what an informative ABSENCE argues against: it supports a rule only when the rule rules one of these traits out
_O2_TRAITS = {"acceptor_oxygen", "uses_O2", "aerobic_chemoheterotrophy"}
ABSENCE_AGAINST = {
    "K00399": {"methanogenesis"}, "K03388": {"methanogenesis"}, "K03389": {"methanogenesis"},
    "K14127": {"methanogenesis"},
    "K00194": {"acetogenesis", "methanogenesis", "acceptor_CO2_methanogen_acetogen"},
    "K00197": {"acetogenesis", "methanogenesis", "acceptor_CO2_methanogen_acetogen"},
    "K02274": _O2_TRAITS, "K02275": _O2_TRAITS, "K02276": _O2_TRAITS, "K00404": _O2_TRAITS, "K02257": _O2_TRAITS,
    "K10535": {"donor_ammonia_oxidation", "nitrification"},
}
# an absence that says something about the oxygen lifestyle: only heme O synthase, without which neither aa3 nor bo3
# oxidase can be built (audit round 2 n3: absences of single oxidases or single anaerobic enzymes are treated alike,
# as weak evidence -> neutral)
ABSENCE_LIFESTYLE = {"K02257": "anaerobic"}
# share of a trait's carriers that carry each gene, registered by the scripts that have the data (audit round 2 N3):
# an absent pathway gene argues against the trait only when it is near-necessary (>= NEAR_NECESSARY of carriers)
CARRIER_SHARE = {}
NEAR_NECESSARY = 0.80


def register_carriers(trait, share):
    """share: mapping KO -> fraction of the trait's carriers (training genomes) that carry the KO."""
    CARRIER_SHARE[trait] = dict(share)


def register_task(X, traits, Y):
    """Register carrier shares of every curated gene for every trait of a task (X: genomes x KOs, Y: labels)."""
    import numpy as np
    genes = [k for k in set().union(*PATHWAY.values(), ABSENCE_AGAINST) if k in X.columns]
    B = X[genes].to_numpy() > 0
    for j, t in enumerate(traits):
        y = np.asarray(Y)[:, j] == 1
        if y.any():
            register_carriers(t, dict(zip(genes, B[y].mean(axis=0))))


def _wanted(trait, cls):
    """The oxygen lifestyle a rule's markers should point to, or None when the trait occurs in both."""
    need = TRAIT_LIFESTYLE.get(trait, "")
    if need.startswith("aerobic") or need == "mostly aerobic":
        w = "aerobic"
    elif need == "anaerobic":
        w = "anaerobic"
    else:
        return None                      # "either", "facultative", unknown: oxygen markers carry no information
    if cls == "lacks trait":
        w = {"aerobic": "anaerobic", "anaerobic": "aerobic"}[w]
    return w


def condition_support(ko, present, trait, cls):
    """How one condition of a rule connects to the trait. Returns (kind, text), kind in
    pathway   - a gene of the trait's pathway, in the state that fits the rule (present -> has, absent -> lacks)
    absence   - an informative absence that argues against the trait, in a rule that rules the trait out
    lifestyle - an oxygen-lifestyle signal that fits the lifestyle the rule implies
    neutral   - carries no information here (an absent aerobic marker; any oxygen marker for 'either' traits)
    opposite  - a pathway gene in the state that does not fit (e.g. required absent in a 'has' rule), or a lifestyle
                signal at odds with the trait
    none      - no known connection to the trait
    """
    name = MECH.get(ko, (ko,))[0]
    has = cls == "has trait"
    def necessary():
        sh = CARRIER_SHARE.get(trait, {}).get(ko)
        return sh is None or sh >= NEAR_NECESSARY, sh
    if ko in PATHWAY.get(trait, set()):
        if present == has:
            if not present:            # absent pathway gene in a rule-out: evidence only if near-necessary
                ok, sh = necessary()
                if not ok:
                    return "neutral", (f"{name} absent: part of the pathway, but {1 - sh:.0%} of organisms with the "
                                       f"trait also lack it (alternative enzymes), so its absence says little")
                return "pathway", (f"{name} absent: a pathway gene carried by "
                                   f"{'most' if sh is None else f'{sh:.0%} of'} organisms with the trait")
            return "pathway", f"{name}: part of the pathway"
        return "opposite", (f"{name}: a pathway gene required {'absent' if has else 'present'} - the rule covers "
                            f"organisms without it (another enzyme or a lineage), not the pathway itself")
    if not present and not has and trait in ABSENCE_AGAINST.get(ko, set()):
        ok, sh = necessary()
        if ok:
            return "absence", f"{name} absent: {ABSENCE_MEANING.get(ko, 'key enzyme missing')}"
        return "neutral", (f"{name} absent: {ABSENCE_MEANING.get(ko, 'enzyme missing')}, but {1 - sh:.0%} of organisms "
                           f"with the trait also lack it, so its absence says little")
    want = _wanted(trait, cls)
    signal = None
    if not present and ko in ABSENCE_LIFESTYLE:
        signal = ABSENCE_LIFESTYLE[ko]
    elif ko in ANAEROBIC_MARKERS | AEROBIC_MARKERS:
        if not present:
            return "neutral", (f"{name} absent: the absence of one {'anaerobic' if ko in ANAEROBIC_MARKERS else 'aerobic'}"
                               f" enzyme is weak evidence about oxygen (alternatives exist)")
        signal = "anaerobic" if ko in ANAEROBIC_MARKERS else "aerobic"
    if signal is not None:
        if want is None:
            return "neutral", f"{name}: {signal}-lifestyle signal, uninformative for a trait found in both lifestyles"
        if signal == want:
            return "lifestyle", f"{name} {'present' if present else 'absent'}: {signal}-lifestyle signal consistent with the rule"
        return "opposite", f"{name} {'present' if present else 'absent'}: {signal}-lifestyle signal at odds with the rule"
    return "none", f"{name}: no known role for this trait"

# audit round 2 n12: anaerobic sulfur oxidisers (Chlorobi, nitrate-dependent), photoferrotrophs and anaerobic
# methylotrophs exist, so these traits carry no oxygen-lifestyle expectation
TRAIT_LIFESTYLE.update({"donor_sulfur_oxidation": "either", "sulfur_oxidation": "either",
                        "donor_iron_oxidation": "either", "methylotrophy": "either"})
