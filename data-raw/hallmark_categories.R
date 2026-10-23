## Reproducible recipe for data/hallmark_categories.rda
## Source: hardcoded HALLMARK -> category mapping used in the lab's
## pw_hm() (rmd-ssgsea-functions.R), as replicated in the heatmap scripts
## used during the ssGSEA2.0 -> ssGSEA3.0 validation. Bundled here as real
## package data so
## it's no longer copy-pasted between scripts.

hallmark_categories <- c(
  HALLMARK_TNFA_SIGNALING_VIA_NFKB = "signaling", HALLMARK_HYPOXIA = "pathway",
  HALLMARK_CHOLESTEROL_HOMEOSTASIS = "metabolic", HALLMARK_MITOTIC_SPINDLE = "proliferation",
  HALLMARK_WNT_BETA_CATENIN_SIGNALING = "signaling", HALLMARK_TGF_BETA_SIGNALING = "signaling",
  HALLMARK_IL6_JAK_STAT3_SIGNALING = "immune", HALLMARK_DNA_REPAIR = "DNA damage",
  HALLMARK_G2M_CHECKPOINT = "proliferation", HALLMARK_APOPTOSIS = "pathway",
  HALLMARK_NOTCH_SIGNALING = "signaling", HALLMARK_ADIPOGENESIS = "development",
  HALLMARK_ESTROGEN_RESPONSE_EARLY = "signaling", HALLMARK_ESTROGEN_RESPONSE_LATE = "signaling",
  HALLMARK_ANDROGEN_RESPONSE = "signaling", HALLMARK_MYOGENESIS = "development",
  HALLMARK_PROTEIN_SECRETION = "pathway", HALLMARK_INTERFERON_ALPHA_RESPONSE = "immune",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "immune", HALLMARK_APICAL_JUNCTION = "cellular component",
  HALLMARK_APICAL_SURFACE = "cellular component", HALLMARK_HEDGEHOG_SIGNALING = "signaling",
  HALLMARK_COMPLEMENT = "immune", HALLMARK_UNFOLDED_PROTEIN_RESPONSE = "pathway",
  HALLMARK_PI3K_AKT_MTOR_SIGNALING = "signaling", HALLMARK_MTORC1_SIGNALING = "signaling",
  HALLMARK_E2F_TARGETS = "proliferation", HALLMARK_MYC_TARGETS_V1 = "proliferation",
  HALLMARK_MYC_TARGETS_V2 = "proliferation", HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "development",
  HALLMARK_INFLAMMATORY_RESPONSE = "immune", HALLMARK_XENOBIOTIC_METABOLISM = "metabolic",
  HALLMARK_FATTY_ACID_METABOLISM = "metabolic", HALLMARK_OXIDATIVE_PHOSPHORYLATION = "metabolic",
  HALLMARK_GLYCOLYSIS = "metabolic", HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "pathway",
  HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY = "pathway", HALLMARK_P53_PATHWAY = "proliferation",
  HALLMARK_UV_RESPONSE_UP = "DNA damage", HALLMARK_UV_RESPONSE_DN = "DNA damage",
  HALLMARK_ANGIOGENESIS = "development", HALLMARK_HEME_METABOLISM = "metabolic",
  HALLMARK_COAGULATION = "immune", HALLMARK_IL2_STAT5_SIGNALING = "signaling",
  HALLMARK_BILE_ACID_METABOLISM = "metabolic", HALLMARK_PEROXISOME = "cellular component",
  HALLMARK_ALLOGRAFT_REJECTION = "immune", HALLMARK_SPERMATOGENESIS = "development",
  HALLMARK_KRAS_SIGNALING_UP = "signaling", HALLMARK_KRAS_SIGNALING_DN = "signaling",
  HALLMARK_PANCREAS_BETA_CELLS = "development"
)

hallmark_category_colors <- c(
  "signaling" = "#88CCEE", "immune" = "#CC6677", "development" = "#999933",
  "proliferation" = "#44AA99", "cellular component" = "#DDDDDD", "metabolic" = "#DDCC77",
  "DNA damage" = "#FFAABB", "pathway" = "#AA4499"
)

usethis::use_data(hallmark_categories, hallmark_category_colors, overwrite = TRUE)
