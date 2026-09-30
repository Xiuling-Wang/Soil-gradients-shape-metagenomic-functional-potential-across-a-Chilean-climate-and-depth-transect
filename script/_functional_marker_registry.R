# Canonical functional-marker registry for the Chile metagenomic project.
#
# The registry was frozen before running the redesigned module tests on
# 2026-07-13. KO identities were checked against KEGG Orthology and KEGG
# pathway modules. Marker presence indicates encoded potential only; it does
# not establish pathway completeness, direction, activity, or flux.

nitrogen_markers_curated <- tibble::tribble(
  ~module, ~tier, ~category, ~KO, ~symbol, ~strict_core, ~gene_filter, ~specificity_note, ~source_basis,
  "nitrogen_transformations", "primary", "N fixation", "K02588", "nifH", TRUE, "none", "Nitrogenase component", "KEGG M00175",
  "nitrogen_transformations", "primary", "N fixation", "K02586", "nifD", TRUE, "none", "Nitrogenase component", "KEGG M00175",
  "nitrogen_transformations", "primary", "N fixation", "K02591", "nifK", TRUE, "none", "Nitrogenase component", "KEGG M00175",
  "nitrogen_transformations", "primary", "Hydroxylamine oxidation", "K10535", "hao", TRUE, "none", "Hydroxylamine dehydrogenase marker; not a complete nitrification module", "KEGG M00528",
  "nitrogen_transformations", "primary", "Nitrate/nitrite redox", "K00370", "narG/nxrA", FALSE, "none", "KO does not distinguish nitrate reduction from nitrite oxidation", "KEGG M00529/M00530",
  "nitrogen_transformations", "primary", "Nitrate/nitrite redox", "K00371", "narH/nxrB", FALSE, "none", "KO does not distinguish nitrate reduction from nitrite oxidation", "KEGG M00529/M00530",
  "nitrogen_transformations", "primary", "Nitrate/nitrite redox", "K00374", "narI", TRUE, "none", "Respiratory nitrate-reductase component", "KEGG M00529/M00530",
  "nitrogen_transformations", "primary", "Nitrate/nitrite redox", "K02567", "napA", TRUE, "none", "Periplasmic nitrate-reductase component", "KEGG M00529/M00530",
  "nitrogen_transformations", "primary", "Nitrate/nitrite redox", "K02568", "napB", TRUE, "none", "Periplasmic nitrate-reductase component", "KEGG M00529/M00530",
  "nitrogen_transformations", "primary", "Denitrification", "K00368", "nirK", TRUE, "none", "NO-forming nitrite reductase", "KEGG M00529",
  "nitrogen_transformations", "primary", "Denitrification", "K15864", "nirS", TRUE, "none", "NO-forming nitrite reductase; KO also includes hydroxylamine-reductase homologs", "KEGG M00529",
  "nitrogen_transformations", "primary", "Denitrification", "K04561", "norB", TRUE, "none", "Nitric-oxide reductase component", "KEGG M00529",
  "nitrogen_transformations", "primary", "Denitrification", "K02305", "norC", TRUE, "none", "Nitric-oxide reductase component", "KEGG M00529",
  "nitrogen_transformations", "primary", "Denitrification", "K00376", "nosZ", TRUE, "none", "Nitrous-oxide reductase", "KEGG M00529",
  "nitrogen_transformations", "primary", "DNRA", "K00362", "nirB", TRUE, "none", "NADH nitrite-reductase component", "KEGG M00530",
  "nitrogen_transformations", "primary", "DNRA", "K00363", "nirD", TRUE, "none", "NADH nitrite-reductase component", "KEGG M00530",
  "nitrogen_transformations", "primary", "DNRA", "K03385", "nrfA", TRUE, "none", "Cytochrome-c nitrite reductase", "KEGG M00530",
  "nitrogen_transformations", "primary", "DNRA", "K15876", "nrfH", TRUE, "none", "Cytochrome-c nitrite-reductase component", "KEGG M00530",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K00367", "narB", TRUE, "none", "Assimilatory nitrate reductase", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K10534", "NR", TRUE, "none", "NAD(P)H nitrate reductase", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K00372", "nasA/nasC", TRUE, "none", "Assimilatory nitrate-reductase component", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K00360", "nasB", TRUE, "none", "Assimilatory nitrate-reductase component", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K00366", "nirA", TRUE, "none", "Ferredoxin nitrite reductase", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K17877", "NIT-6", TRUE, "none", "NAD(P)H nitrite reductase", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K26139", "nasD/nasB", TRUE, "none", "NAD(P)H nitrite-reductase component", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K26138", "nasE", TRUE, "none", "NAD(P)H nitrite-reductase component", "KEGG M00531",
  "nitrogen_transformations", "primary", "Assimilatory nitrate reduction", "K00361", "nasB", TRUE, "none", "NAD(P)H nitrite reductase", "KEGG M00531",
  "nitrogen_transformations", "primary", "Urea hydrolysis", "K01428", "ureC", TRUE, "none", "Urease component", "KEGG Orthology",
  "nitrogen_transformations", "primary", "Urea hydrolysis", "K01429", "ureB", TRUE, "none", "Urease component", "KEGG Orthology",
  "nitrogen_transformations", "primary", "Urea hydrolysis", "K01430", "ureA", TRUE, "none", "Urease component", "KEGG Orthology"
)

pedogenic_markers_curated <- tibble::tribble(
  ~module, ~tier, ~category, ~KO, ~symbol, ~strict_core, ~gene_filter, ~specificity_note, ~source_basis,
  "pedogenic_resource_acquisition", "primary", "Organic-acid weathering", "K00117", "gcd", TRUE, "none", "PQQ-dependent glucose dehydrogenase; encoded potential is not a weathering-rate measurement", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Organic-acid weathering", "K06135", "pqqA", TRUE, "none", "PQQ-biosynthesis component", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Organic-acid weathering", "K06136", "pqqB", TRUE, "none", "PQQ-biosynthesis component", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Organic-acid weathering", "K06137", "pqqC", TRUE, "none", "PQQ-biosynthesis component", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Organic-acid weathering", "K06138", "pqqD", TRUE, "none", "PQQ-biosynthesis component", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Organic-acid weathering", "K06139", "pqqE", TRUE, "none", "PQQ-biosynthesis component", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "P acquisition", "K01113", "phoD", TRUE, "none", "Alkaline phosphatase", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "P acquisition", "K01077", "phoA/ALP", TRUE, "none", "Alkaline phosphatase", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "P acquisition", "K02040", "pstS", FALSE, "none", "High-affinity phosphate-transport component; resource acquisition but not pedogenesis-specific", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "P acquisition", "K06163", "phnJ", TRUE, "none", "C-P lyase component", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "P acquisition", "K00937", "ppk1", FALSE, "none", "Polyphosphate kinase; resource storage but not pedogenesis-specific", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Catecholate siderophore synthesis", "K00216", "entA", TRUE, "none", "Representative catecholate-siderophore marker, not a universal siderophore inventory", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Catecholate siderophore synthesis", "K02361", "entC", TRUE, "none", "Representative catecholate-siderophore marker, not a universal siderophore inventory", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Catecholate siderophore synthesis", "K01252", "entB", TRUE, "none", "Representative catecholate-siderophore marker, not a universal siderophore inventory", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "Catecholate siderophore synthesis", "K02363", "entE", TRUE, "none", "Representative catecholate-siderophore marker, not a universal siderophore inventory", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "EPS and aggregation", "K01991", "wza/gfcE", FALSE, "none", "Representative polysaccharide-export marker; taxonomically non-universal", "Rodriguez 2025; KEGG Orthology",
  "pedogenic_resource_acquisition", "primary", "EPS and aggregation", "K11936", "pgaC/icaA", FALSE, "none", "Representative PNAG-synthesis marker; taxonomically non-universal", "Rodriguez 2025; KEGG Orthology"
)

alternative_energy_markers_curated <- tibble::tribble(
  ~module, ~tier, ~category, ~KO, ~symbol, ~strict_core, ~gene_filter, ~specificity_note, ~source_basis,
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17222", "soxA", TRUE, "none", "Sox-system component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17223", "soxX", TRUE, "none", "Sox-system component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17224", "soxB", TRUE, "none", "Sox-system component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17225", "soxC", TRUE, "none", "Sox-system component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17226", "soxY", TRUE, "none", "Sox-system component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17227", "soxZ", TRUE, "none", "Sox-system component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17218", "sqr", TRUE, "none", "Sulfide:quinone oxidoreductase", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur oxidation", "K17229", "fccB", TRUE, "none", "Sulfide-dehydrogenase component", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur redox, direction unresolved", "K00394", "aprA", FALSE, "none", "KO abundance alone cannot distinguish oxidative from reductive direction", "Horstmann 2025; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur redox, direction unresolved", "K00395", "aprB", FALSE, "none", "KO abundance alone cannot distinguish oxidative from reductive direction", "Horstmann 2025; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur redox, direction unresolved", "K11180", "dsrA", FALSE, "none", "KO abundance alone cannot distinguish oxidative from reductive direction", "Horstmann 2025; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur redox, direction unresolved", "K11181", "dsrB", FALSE, "none", "KO abundance alone cannot distinguish oxidative from reductive direction", "Horstmann 2025; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Sulfur redox, direction unresolved", "K08352", "phsA/psrA", FALSE, "none", "Thiosulfate/polysulfide-reductase homolog", "Horstmann 2025; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Hydrogenase homologs", "K06281", "hyaB/hybC", FALSE, "hydrogenase_name", "Hydrogenase large-subunit homolog; KO does not resolve affinity group or in-situ direction", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Hydrogenase homologs", "K06282", "hyaA/hybO", FALSE, "hydrogenase_name", "Hydrogenase small-subunit homolog; KO does not resolve affinity group or in-situ direction", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Hydrogenase homologs", "K05922", "hydB", FALSE, "hydrogenase_name", "Quinone-reactive NiFe-hydrogenase homolog; in-situ direction unresolved", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Hydrogenase homologs", "K00436", "hoxH", FALSE, "hydrogenase_name", "NAD-reducing hydrogenase catalytic-subunit homolog; in-situ direction unresolved", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Putative CO oxidation", "K03518", "coxS", TRUE, "coxS_specific", "Retain only annotations explicitly assigned to aerobic CO-dehydrogenase small subunits", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Putative CO oxidation", "K03519", "coxM/cutM", TRUE, "coxM_specific", "Retain only annotations explicitly assigned to aerobic CO-dehydrogenase medium subunits", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "Putative CO oxidation", "K03520", "coxL/cutL", TRUE, "coxL_specific", "Retain only explicit aerobic CO-dehydrogenase annotations; active-site motif could not be checked without catalog protein sequences", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "CBB carbon-fixation markers", "K01601", "rbcL/cbbL", TRUE, "rubisco_name", "RuBisCO large-chain homolog; pathway completeness requires accompanying CBB components", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "CBB carbon-fixation markers", "K01602", "rbcS/cbbS", TRUE, "rubisco_name", "RuBisCO small-chain homolog; pathway completeness requires accompanying CBB components", "Horstmann 2025; Medina Caro 2026; KEGG Orthology",
  "alternative_energy_chemolithotrophy", "primary", "CBB carbon-fixation markers", "K00855", "prk", TRUE, "prk_specific", "Retain only annotations explicitly named prk/cbbP/phosphoribulokinase, excluding uridine-kinase-like assignments", "Horstmann 2025; Medina Caro 2026; KEGG Orthology"
)

stress_markers_curated <- tibble::tribble(
  ~module, ~tier, ~category, ~KO, ~symbol, ~strict_core, ~gene_filter, ~specificity_note, ~source_basis,
  "stress_response", "supplementary", "Protein-folding stress", "K04077", "groEL", TRUE, "none", "General chaperonin; abundance is encoded capacity, not a stress-response measurement", "KEGG Orthology",
  "stress_response", "supplementary", "Protein-folding stress", "K04043", "dnaK", TRUE, "none", "General chaperone; abundance is encoded capacity, not a stress-response measurement", "KEGG Orthology",
  "stress_response", "supplementary", "Protein-folding stress", "K04079", "htpG", TRUE, "none", "General chaperone; abundance is encoded capacity, not a stress-response measurement", "KEGG Orthology",
  "stress_response", "supplementary", "Protein-folding stress", "K04078", "groES", TRUE, "none", "General chaperonin; abundance is encoded capacity, not a stress-response measurement", "KEGG Orthology",
  "stress_response", "supplementary", "Protein-folding stress", "K03695", "clpB", TRUE, "none", "General protein-disaggregation capacity", "KEGG Orthology",
  "stress_response", "supplementary", "Oxidative-stress protection", "K04564", "sodA/sodB", TRUE, "none", "Fe/Mn superoxide-dismutase family", "KEGG Orthology",
  "stress_response", "supplementary", "Oxidative-stress protection", "K04565", "sodC", TRUE, "none", "Cu/Zn superoxide-dismutase family", "KEGG Orthology",
  "stress_response", "supplementary", "Oxidative-stress protection", "K03781", "katE", TRUE, "none", "Catalase", "KEGG Orthology",
  "stress_response", "supplementary", "Oxidative-stress protection", "K00432", "gpx", TRUE, "none", "Glutathione peroxidase", "KEGG Orthology",
  "stress_response", "supplementary", "Osmoprotection", "K03762", "proP", TRUE, "none", "Proline/betaine transporter", "KEGG Orthology",
  "stress_response", "supplementary", "Osmoprotection", "K00697", "otsA", TRUE, "none", "Trehalose-6-phosphate synthase", "KEGG Orthology",
  "stress_response", "supplementary", "Osmoprotection", "K01087", "otsB", TRUE, "none", "Trehalose-6-phosphate phosphatase", "KEGG Orthology",
  "stress_response", "supplementary", "Cold-shock proteins", "K03704", "cspA/cspB", TRUE, "none", "Cold-shock protein family; KO does not distinguish csp paralogs", "KEGG Orthology"
)

functional_marker_registry <- dplyr::bind_rows(
  nitrogen_markers_curated,
  pedogenic_markers_curated,
  alternative_energy_markers_curated,
  stress_markers_curated
) %>%
  dplyr::distinct(module, category, KO, .keep_all = TRUE)

# Accurate replacements for the retired broad carbon panel. This object exists
# only so legacy scripts cannot silently regenerate known-wrong labels. It is
# not one of the three primary modules in the redesigned analysis.
legacy_carbon_markers_curated <- tibble::tribble(
  ~KO, ~symbol, ~category, ~gene_filter,
  "K01601", "rbcL/cbbL", "CBB markers", "rubisco_name",
  "K01602", "rbcS/cbbS", "CBB markers", "rubisco_name",
  "K00855", "prk", "CBB markers", "prk_specific",
  "K15230", "aclA", "rTCA marker homologs", "aclA_specific",
  "K15231", "aclB", "rTCA marker homologs", "aclB_specific",
  "K00421", "lccA", "Laccase homologs", "none",
  "K05909", "laccase", "Laccase homologs", "none",
  "K03518", "coxS", "Putative CO oxidation", "coxS_specific",
  "K03519", "coxM/cutM", "Putative CO oxidation", "coxM_specific",
  "K03520", "coxL/cutL", "Putative CO oxidation", "coxL_specific",
  "K16157", "mmoX", "Soluble methane monooxygenase", "none",
  "K16158", "mmoY", "Soluble methane monooxygenase", "none",
  "K10944", "pmoA/amoA homolog", "Methane/ammonia monooxygenase homologs", "none",
  "K10945", "pmoB/amoB homolog", "Methane/ammonia monooxygenase homologs", "none",
  "K10946", "pmoC/amoC homolog", "Methane/ammonia monooxygenase homologs", "none",
  "K14028", "mxaF", "Methanol oxidation", "none"
)

legacy_marker_corrections <- tibble::tribble(
  ~old_KO, ~old_label, ~correct_KO, ~correct_label, ~action, ~reason,
  "K14067", "cbbM", NA_character_, "mtkA", "exclude", "K14067 is malate-CoA ligase subunit beta, not RuBisCO",
  "K09709", "aclA", "K15230", "aclA", "replace", "K09709 is 3-methylfumaryl-CoA hydratase",
  "K09710", "aclB", "K15231", "aclB", "replace", "K09710 is a ribosome-associated protein",
  "K01961", "accA", NA_character_, "accC", "exclude_from_specific_C_fixation", "K01961 is non-specific acetyl-CoA carboxylase biotin-carboxylase subunit",
  "K01534", "ppc", "K01595", "ppc", "replace_but_exclude_from_specific_C_fixation", "K01534 is a Zn/Cd-exporting ATPase; PEP carboxylase is also not C-fixation-specific",
  "K00525", "lac", "K00421/K05909", "laccase", "replace", "K00525 is ribonucleotide reductase",
  "K00184", "coxS", "K03518", "coxS", "replace", "K00184 is dimethyl-sulfoxide reductase iron-sulfur subunit",
  "K16157", "pmoA", "K10944", "pmoA/amoA homolog", "relabel_and_replace", "K16157 is soluble methane-monooxygenase mmoX",
  "K16158", "pmoB", "K10945", "pmoB/amoB homolog", "relabel_and_replace", "K16158 is soluble methane-monooxygenase mmoY",
  "K14080", "mxaF", "K14028", "mxaF", "replace", "K14080 is methyltransferase mtaA",
  "K06135", "pqqB", "K06135", "pqqA", "relabel", "KEGG assigns K06135 to pqqA",
  "K06136", "pqqC", "K06136", "pqqB", "relabel", "KEGG assigns K06136 to pqqB; pqqC is K06137",
  "K07093", "phoX", NA_character_, "uncharacterized DUF839 protein", "exclude", "Neither KEGG nor local eggNOG supports a phoX assignment",
  "K03715", "proP", "K03762", "proP", "replace", "K03715 is a galactolipid glycosyltransferase",
  "K01193", "otsA", "K00697", "otsA", "replace", "K01193 is beta-fructofuranosidase",
  "K01239", "otsB", "K01087", "otsB", "replace", "K01239 is purine nucleosidase",
  "K03705", "cspB", "K03704", "cspA/cspB", "replace_and_deduplicate", "K03705 is hrcA; cspA/cspB homologs map to K03704",
  "K04565", "sodB", "K04565", "sodC", "relabel", "K04565 is the Cu/Zn superoxide-dismutase family",
  "K18005", "hybC", "K18005", "hoxF", "relabel", "K18005 is the diaphorase moiety of NAD-reducing NiFe hydrogenase",
  "K18016", "hupL", "K18016", "mbhL", "relabel", "K18016 is membrane-bound hydrogenase subunit alpha",
  "K18017", "hupS", "K18017", "mbhK", "relabel", "K18017 is membrane-bound hydrogenase subunit beta"
)

cazy_substrate_families <- list(
  Cellulose = paste0("GH", c(5, 6, 9, 12, 44, 45, 48)),
  Hemicellulose = paste0("GH", c(8, 10, 11, 26, 30, 43, 53, 54, 62, 67, 115)),
  Pectin = c("GH28", paste0("PL", c(1, 2, 3, 9, 10, 11)), "CE8"),
  Chitin = c("GH18", "GH19", "CE4"),
  `Starch/glycogen` = paste0("GH", c(13, 15, 31, 57, 4)),
  Peptidoglycan = paste0("GH", c(23, 24, 25, 73, 102, 103, 104, 108))
)

# Vectorized annotation-name filter shared by inferential and figure scripts.
# It is intentionally conservative for broad KOs whose orthology alone does
# not establish the requested functional identity.
functional_marker_filter_pass <- function(gene_filter, gene_name, description) {
  gene_name <- stringr::str_to_lower(tidyr::replace_na(gene_name, ""))
  description <- stringr::str_to_lower(tidyr::replace_na(description, ""))
  text <- stringr::str_c(gene_name, description, sep = " | ")
  token <- function(pattern) {
    stringr::str_detect(gene_name, stringr::regex(pattern, ignore_case = TRUE))
  }

  dplyr::case_when(
    gene_filter == "none" ~ TRUE,
    gene_filter == "hydrogenase_name" ~ stringr::str_detect(text, "hydrogenase"),
    gene_filter == "rubisco_name" ~
      token("(^|[^a-z0-9])(rbcl|rbcs|cbbl|cbbs)([^a-z0-9]|$)") |
      stringr::str_detect(
        description,
        "rubisco|ribulose.{0,20}bisphosphate.{0,20}carboxylase"
      ),
    gene_filter == "prk_specific" ~
      (
        token("(^|[^a-z0-9])(prk|cbbp)([^a-z0-9]|$)") |
          stringr::str_detect(description, "phosphoribulokinase")
      ) & !stringr::str_detect(text, "uridine kinase"),
    gene_filter == "coxS_specific" ~
      token("(^|[^a-z0-9])(coxs|cuts)([^a-z0-9]|$)") |
      stringr::str_detect(
        description,
        "carbon monoxide dehydrogenase.{0,40}small"
      ),
    gene_filter == "coxM_specific" ~
      token("(^|[^a-z0-9])(coxm|cutm)([^a-z0-9]|$)") |
      stringr::str_detect(
        description,
        "carbon monoxide dehydrogenase.{0,40}(medium|middle)"
      ),
    gene_filter == "coxL_specific" ~
      token("(^|[^a-z0-9])(coxl|cutl)([^a-z0-9]|$)") |
      stringr::str_detect(
        description,
        "carbon monoxide dehydrogenase.{0,40}large"
      ),
    gene_filter == "aclA_specific" ~
      token("(^|[^a-z0-9])acla([^a-z0-9]|$)") |
      stringr::str_detect(
        description,
        "atp citrate lyase.{0,40}(alpha|subunit a)"
      ),
    gene_filter == "aclB_specific" ~
      token("(^|[^a-z0-9])aclb([^a-z0-9]|$)") |
      stringr::str_detect(
        description,
        "atp citrate lyase.{0,40}(beta|subunit b)"
      ),
    TRUE ~ FALSE
  )
}
