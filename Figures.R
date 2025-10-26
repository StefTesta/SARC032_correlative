##### Code for images of SU2C-SARC032 Correlative TME Analysis Manuscript #####
# load packages ----
set.seed(123)

library(tidyverse)    
library(data.table)
library(readxl)
library(qs)
library(broom)
library(broom.mixed)
library(ggpubr)
library(ggsignif)
library(ggsci)
library(ggrepel)
library(ggbeeswarm)
library(scales)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(scCustomize)
library(lme4)
library(lmerTest)
library(emmeans)
library(survival)
library(survminer)
library(coxphf)
library(forestplot)
library(forestploter)
library(compositions)
library(tidyplots)
library(rstatix)
library(SingleCellExperiment)
library(scater)

# Load data from Supp. Tables ----
flow_results_filtered <- read_excel(".../data//Supplementary_Tables.xlsx", 
                                    sheet = 4)
colnames(flow_results_filtered) <- flow_results_filtered[2, ]
flow_results_filtered <- flow_results_filtered[-c(1,2) ,]
flow_results_filtered <- as.data.frame(flow_results_filtered)

SIC_and_SE_assignments <- read_excel(".../data//Supplementary_Tables.xlsx", 
                                     sheet = 1)
SIC_and_SE_assignments <- SIC_and_SE_assignments[-c(1), ]
colnames(SIC_and_SE_assignments) <- SIC_and_SE_assignments[1, ]
SIC_and_SE_assignments <- SIC_and_SE_assignments[-c(1), ]
SIC_and_SE_assignments <- as.data.frame(SIC_and_SE_assignments)
SIC_and_SE_assignments$`DFS Time (Days)` <- as.numeric(SIC_and_SE_assignments$`DFS Time (Days)`)
SIC_and_SE_assignments$`Percent Necrosis` <- as.numeric(SIC_and_SE_assignments$`Percent Necrosis`)

src_lookup <- c(
  "803-072"="SRC83",  "001-074"="SRC82",  "022-059"="SRC75",  "004-058"="SRC74",
  "004-045"="SRC68",  "022-046"="SRC67",  "022-044"="SRC66",  "001-041"="SRC65",
  "092-037"="SRC64",  "004-040"="SRC63",  "005-027"="SRC55",  "803-029"="SRC54",
  "802-103"="SRC535", "802-149"="SRC532", "802-148"="SRC531", "005-135"="SRC530",
  "091-018"="SRC53",  "076-102"="SRC529", "076-146"="SRC527", "801-151"="SRC523",
  "802-132"="SRC521", "802-120"="SRC520", "001-023"="SRC52",  "022-025"="SRC51",
  "022-026"="SRC50",  "029-021"="SRC48",  "005-015"="SRC47",  "004-107"="SRC464",
  "004-112"="SRC463", "001-134"="SRC462", "005-147"="SRC461", "029-020"="SRC46",
  "801-129"="SRC459", "802-117"="SRC458", "802-122"="SRC457", "092-144"="SRC456",
  "802-108"="SRC454", "802-138"="SRC452", "098-155"="SRC450", "048-016"="SRC45",
  "801-141"="SRC449", "092-140"="SRC447", "001-017"="SRC44",  "005-013"="SRC42",
  "071-009"="SRC40",  "091-010"="SRC39",  "022-012"="SRC38",  "001-131"="SRC364",
  "092-124"="SRC362", "091-006"="SRC36",  "015-002"="SRC35",  "001-004"="SRC34",
  "071-003"="SRC33",  "801-042"="SRC279","051-094"="SRC278","051-087"="SRC276",
  "802-086"="SRC275","802-019"="SRC274","802-051"="SRC273","802-050"="SRC271",
  "802-063"="SRC269","802-062"="SRC267","802-071"="SRC265","802-081"="SRC264",
  "802-075"="SRC263","802-080"="SRC261","802-085"="SRC260","801-097"="SRC258",
  "802-056"="SRC256","802-096"="SRC255","802-070"="SRC254","016-098"="SRC251",
  "051-092"="SRC250","005-109"="SRC248","016-114"="SRC247","004-110"="SRC246",
  "016-104"="SRC245","004-090"="SRC221","016-095"="SRC220","002-083"="SRC217",
  "005-111"=NA_character_,"016-043"=NA_character_,
  "022-152"="SRC526","802-100"="SRC455","092-142"="SRC451","001-014"="SRC43",
  "001-128"="SRC365","802-078"="SRC262","022-105"="SRC244","001-076"="SRC84",
  "803-069"="SRC81","001-066"="SRC80","001-038"="SRC62","015-028"="SRC56",
  "802-143"="SRC533","076-068"="SRC528","802-150"="SRC524","048-022"="SRC49",
  "034-137"="SRC374","092-126"="SRC363","022-127"="SRC361","048-119"="SRC360",
  "802-052"="SRC268","001-115"="SRC259","001-089"="SRC219","071-073"="SRC79",
  "001-067"="SRC78","092-049"="SRC73","091-055"="SRC72","071-034"="SRC60",
  "029-033"="SRC59","022-032"="SRC58","001-030"="SRC57","802-139"="SRC536",
  "802-153"="SRC534","802-121"="SRC522","802-130"="SRC453","005-007"="SRC41",
  "092-008"="SRC37","801-088"="SRC277","051-057"="SRC272","051-060"="SRC270",
  "051-077"="SRC266","016-116"="SRC257","802-079"="SRC253","001-101"="SRC223",
  "029-099"="SRC222"
)

SIC_and_SE_assignments$Patient <- as.character(SIC_and_SE_assignments$Patient)

if (!"SRC code" %in% names(SIC_and_SE_assignments)) {
  SIC_and_SE_assignments[["SRC code"]] <- unname(src_lookup[ SIC_and_SE_assignments$Patient ])
} else {
  na_idx <- is.na(SIC_and_SE_assignments[["SRC code"]])
  SIC_and_SE_assignments[["SRC code"]][na_idx] <-
    unname(src_lookup[ SIC_and_SE_assignments$Patient[na_idx] ])
}

# Figure 3A ----- ----
cd4_clusters <- c(
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count",
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD103+ | Count",    # ITGAE                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD152+ | Count",    # CTLA4                          
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD279 + | Count",   # PD-1                           
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD366+ | Count",    # CD366                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD39+ | Count",     # ENTPD1 
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q1: CD45RA BUV563- , CD197 AC7+ | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q2: CD38 PCPC55+ , HLADR BV480+ | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q2: CD45RA BUV563+ , CD197 AC7+ | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q3: CD45RA BUV563+ , CD197 AC7- | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q4: CD45RA BUV563- , CD197 AC7- | Count" 
)

flow_results_filtered <- flow_results_filtered %>%
  mutate(across(-1, as.numeric))

flow_results_filtered <- as.data.table(flow_results_filtered)

flow_results_macro_clusters <- flow_results_filtered[, colnames(flow_results_filtered) %in% 
                                                       unique(c("Patient", cd4_clusters)), with = F]
flow_results_macro_clusters <- flow_results_macro_clusters %>%
  mutate(
    CD103_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD103+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD152_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD152+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD279_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD279 + | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD366_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD366+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD39_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD39+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD38pos_DRpos_activated_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q2: CD38 PCPC55+ , HLADR BV480+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD45RAneg_CCR7pos_CM_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q1: CD45RA BUV563- , CD197 AC7+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD45RApos_CCR7pos_Naive_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q2: CD45RA BUV563+ , CD197 AC7+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD45RApos_CCR7neg_Eff_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q3: CD45RA BUV563+ , CD197 AC7- | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
    CD45RAneg_CCR7neg_EM_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q4: CD45RA BUV563- , CD197 AC7- | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count`,
  )

flow_results_macro_clusters <- flow_results_macro_clusters %>% 
  relocate(Patient, .before = CD103_frac)

flow_results_macro_clusters <- flow_results_macro_clusters[, 12:22]
flow_results_macro_clusters <- flow_results_macro_clusters %>% 
  left_join(SIC_and_SE_assignments)

macro_clusters_flow_heatmap_df <- as.matrix(flow_results_macro_clusters[, 1:11])
rownames(macro_clusters_flow_heatmap_df) <- macro_clusters_flow_heatmap_df[, 1]
macro_clusters_flow_heatmap_df <- macro_clusters_flow_heatmap_df[, -1]

macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df) %>% 
  mutate(across(1:10, as.numeric))

macro_clusters_flow_heatmap_df <- t(macro_clusters_flow_heatmap_df)

column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df))

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                               "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")

column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df), ]

column_annotation <- column_annotation %>% 
  relocate(Arm)

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
 col = list(`Tumor Grade` = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             Arm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))


col_custom <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),viridis(5, option = "C"))

macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df) %>%
  rownames_to_column(var = "Cluster") %>%
  mutate(Cluster = case_when(
    Cluster == "CD103_frac"    ~ "CD103+ CD4+ T Cells",
    Cluster == "CD152_frac"  ~ "CD152+ CD4+ T Cells",
    Cluster == "CD279_frac"  ~ "CD279+ CD4+ T Cells",
    Cluster == "CD366_frac"   ~ "CD366+ CD4+ T Cells",
    Cluster == "CD39_frac" ~ "CD39+ CD4+ T Cells",
    Cluster == "CD38pos_DRpos_activated_frac" ~ "CD38+ HLA-DR+ CD4+ T Cells",
    Cluster == "CD45RAneg_CCR7pos_CM_frac" ~ "CM CD4+ T Cells",
    Cluster == "CD45RApos_CCR7pos_Naive_frac" ~ "Naive CD4+ T Cells",
    Cluster == "CD45RApos_CCR7neg_Eff_frac" ~ "Effector CD4+ T Cells",
    Cluster == "CD45RAneg_CCR7neg_EM_frac" ~ "EM CD4+ T Cells",
    TRUE ~ Cluster
  )) %>%
  column_to_rownames(var = "Cluster") 

macro_clusters_flow_heatmap_df <- as.matrix(macro_clusters_flow_heatmap_df)

ht = Heatmap(macro_clusters_flow_heatmap_df, 
             col = col_custom, 
             column_split = as.factor(column_annotation$Arm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(10, "cm"), 
             heatmap_width = unit(20, "cm"),
             cluster_column_slices = F, 
             name = "% of total Lymphocytes",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../flow_CD4_clusters.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure 3B ----
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

data_for_analysis <- tumor_long %>%
  inner_join(sample_annotation, by = "Sample")

data_for_analysis <- data_for_analysis %>% 
  left_join(SIC_and_SE_assignments[, c("SRC code", "DFS Event", "DFS Time (Days)")], 
            by = c("Sample" = "SRC code"))

data_for_analysis_cd4 <- data_for_analysis

violin_boxplot <- ggplot(data_for_analysis[data_for_analysis$CellType %in% c("CD152+ CD4+ T Cells",
                                                                             "CD279+ CD4+ T Cells",
                                                                             "CD38+ HLA-DR+ CD4+ T Cells", 
                                                                             "Naive CD4+ T Cells"), ], aes(x = Arm, y = ZScore, fill = Arm)) + 
  geom_violin(scale = "width", 
              trim = F,
              alpha = 0.5) +
  geom_boxplot(width = 0.05, 
               #color = "darkblue", 
               alpha = .8) + 
  geom_quasirandom(method = "smiley") + 
  facet_wrap(~CellType, nrow = 2, 
             scales = "free") + 
  scale_fill_manual(values = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b")) + 
  labs(title = "", y = "% of Total CD3+CD4+", x = "Arm") + 
  stat_compare_means(comparisons = list(c("Control", "Experimental")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     hide.ns = F)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox") +
  theme_pubr(legend = "right", border = T) + 
  theme(strip.text = element_text(size = 13))

ggsave(violin_boxplot,
       width = 22.5, height = 18, 
       units = "cm", 
       filename = '.../flow_CD4_clusters_boxplot.pdf')

# Get FDRs 
test_results <- data_for_analysis_cd4 %>%
  group_by(CellType) %>%
  wilcox_test(ZScore ~ Arm) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  dplyr::select(CellType, group1, group2, p, p.adj, p.adj.signif)

# Figure 3C ----
SE1_patients <- SIC_and_SE_assignments %>% 
  filter(`Sarcoma Ecotype Assignment` == "SE1") %>% 
  pull(Patient)

SICE_patients <- SIC_and_SE_assignments %>% 
  filter(`Sarcoma Immune Class` == "E") %>% 
  pull(Patient)

macro_clusters_flow_heatmap_df_SE1_SICE <- macro_clusters_flow_heatmap_df[, colnames(macro_clusters_flow_heatmap_df) %in% unique(c(SE1_patients, SICE_patients))]
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df_SE1_SICE) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

# create  annotation 
column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df_SE1_SICE))
column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df_SE1_SICE), ]
column_annotation <- column_annotation %>% 
  relocate(Arm)
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

sample_annotation <- sample_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Sarcoma Ecotype Assignment", 
                                       "Sarcoma Immune Class")], 
            by = c("Sample" = "Patient"))

sample_annotation$SICE_vs_SE1_classifier <- NA
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Ecotype Assignment` == "SE1"] <- "SE1"
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Immune Class` == "E"] <- "SICE"
sample_annotation$split_factor <- paste(sample_annotation$Arm, sample_annotation$SICE_vs_SE1_classifier, sep = "_")
sample_annotation$split_factor <- as.factor(sample_annotation$split_factor)
sample_annotation$split_factor <- factor(sample_annotation$split_factor, 
                                         levels = c("Control_SE1", "Experimental_SE1",
                                                    "Control_SICE", "Experimental_SICE"))

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#D682FF", "Grade 3" = "#0095FF"), 
             Arm = c("Control" = "#9D1536", "Experimental" = "#0093AF"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))


col_custom <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),viridis(5, option = "C"))

# Heatmap of macro-clusters 
ht = Heatmap(macro_clusters_flow_heatmap_df_SE1_SICE, 
             col = col_custom, 
             column_split = as.factor(sample_annotation$split_factor),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(10, "cm"), 
             heatmap_width = unit(18, "cm"),
             cluster_column_slices = F, 
             name = "% of total CD4+",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../flow_CD4_non_t_regs_clusters_SE1_SICE.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure 3D ----
macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df)
macro_clusters_flow_heatmap_df_SE1_SICE <- macro_clusters_flow_heatmap_df[, colnames(macro_clusters_flow_heatmap_df) %in% unique(c(SE1_patients, SICE_patients))]
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df_SE1_SICE) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

# create  annotation 
column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df_SE1_SICE))
column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df_SE1_SICE), ]
column_annotation <- column_annotation %>% 
  relocate(Arm)
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

sample_annotation <- sample_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Sarcoma Ecotype Assignment", 
                                       "Sarcoma Immune Class")], 
            by = c("Sample" = "Patient"))

sample_annotation$SICE_vs_SE1_classifier <- NA
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$ecotype_assignment == "E1"] <- "SE1"
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Immune Class` == "E"] <- "SICE"

data_for_analysis <- tumor_long %>%
  inner_join(sample_annotation, by = "Sample")

data_for_analysis_macro <- data_for_analysis

macro_violin <- data_for_analysis %>% 
  filter(CellType %in% c("CM CD4+ T Cells", 
                         "EM CD4+ T Cells",
                         "CD279+ CD4+ T Cells",
                         "CD366+ CD4+ T Cells"
  )) %>% 
  ggplot(aes(x = Arm, y = ZScore, fill = Arm)) +
  geom_violin(scale = "width", 
              trim = F,
              alpha = 0.5) +
  geom_boxplot(width = 0.05, 
               alpha = .8) + 
  geom_quasirandom(method = "smiley") + 
  facet_grid(SICE_vs_SE1_classifier ~ CellType,
             scales="free_y",
             switch = "y") + 
  scale_fill_manual(values = c("Control" = "#9D1536", "Experimental" = "#0093AF")) + 
  labs(title = "", y = "Z-Scored Cell %", x = "Arm") + 
  stat_compare_means(comparisons = list(c("Control", "Experimental")), 
                     label = "p.signif", 
                     method = "wilcox.test", hide.ns = F) + 
  stat_compare_means(method = "wilcox") +
  theme_pubr(legend = "right", border = T) + 
  theme(strip.text = element_text(size = 13))

ggsave(macro_violin, units = "cm", 
       width = 27.5, height = 15, 
       filename = '.../violin_box_CD4_non_t_regs_SE1_vs_SICE.pdf'
)

# Get FDRs 
test_results <- data_for_analysis %>%
  group_by(CellType, SICE_vs_SE1_classifier) %>%
  wilcox_test(ZScore ~ Arm) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(CellType, SICE_vs_SE1_classifier, group1, group2, p, p.adj, p.adj.signif)


# Figure Supplemental 3A ---- 
cd8_clusters <- c(
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count",
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD103+ | Count",    # ITGAE                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD152+ | Count",    # CTLA4                          
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD279 + | Count",   # PD-1                           
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD366+ | Count",    # CD366                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+ | Count",     # ENTPD1 
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+/CD103+ | Count",                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+/CD279 + | Count",  
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q1: CD45RA BUV563- , CD197 AC7+ | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q2: CD38 PCPC55+ , HLADR BV480+ | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q2: CD45RA BUV563+ , CD197 AC7+ | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q3: CD45RA BUV563+ , CD197 AC7- | Count",     
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q4: CD45RA BUV563- , CD197 AC7- | Count" 
)

flow_results_filtered <- flow_results_filtered %>%
  mutate(across(-1, as.numeric))

flow_results_filtered <- as.data.table(flow_results_filtered)

flow_results_macro_clusters <- flow_results_filtered[, colnames(flow_results_filtered) %in% 
                                                       unique(c("Patient", cd8_clusters)), with = F]

flow_results_macro_clusters <- flow_results_macro_clusters %>%
  mutate(
    CD103_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD103+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD152_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD152+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD279_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD279 + | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD366_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD366+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD39_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    
    CD39_CD103_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+/CD103+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    
    CD39_CD279_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+/CD279 + | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    
    CD38pos_DRpos_activated_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q2: CD38 PCPC55+ , HLADR BV480+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD45RAneg_CCR7pos_CM_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q1: CD45RA BUV563- , CD197 AC7+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD45RApos_CCR7pos_Naive_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q2: CD45RA BUV563+ , CD197 AC7+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD45RApos_CCR7neg_Eff_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q3: CD45RA BUV563+ , CD197 AC7- | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
    CD45RAneg_CCR7neg_EM_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q4: CD45RA BUV563- , CD197 AC7- | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count`,
  )

flow_results_macro_clusters <- flow_results_macro_clusters %>% 
  relocate(Patient, .before = CD103_frac)

flow_results_macro_clusters <- flow_results_macro_clusters[, 14:26]
flow_results_macro_clusters <- flow_results_macro_clusters %>% 
  left_join(SIC_and_SE_assignments)

macro_clusters_flow_heatmap_df <- as.matrix(flow_results_macro_clusters[, 1:13])
rownames(macro_clusters_flow_heatmap_df) <- macro_clusters_flow_heatmap_df[, 1]
macro_clusters_flow_heatmap_df <- macro_clusters_flow_heatmap_df[, -1]

macro_clusters_flow_heatmap_df[is.na(macro_clusters_flow_heatmap_df)] <- 0

macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df) %>% 
  mutate(across(1:12, as.numeric))

macro_clusters_flow_heatmap_df <- t(macro_clusters_flow_heatmap_df)

column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df))

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")

column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df), ]

column_annotation <- column_annotation %>% 
  relocate(Arm)

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             Arm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))


col_custom <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),viridis(5, option = "C"))

macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df) %>%
  rownames_to_column(var = "Cluster") %>%
  mutate(Cluster = case_when(
    Cluster == "CD103_frac"    ~ "CD103+ CD8+ T Cells",
    Cluster == "CD152_frac"  ~ "CD152+ CD8+ T Cells",
    Cluster == "CD279_frac"  ~ "CD279+ CD8+ T Cells",
    Cluster == "CD366_frac"   ~ "CD366+ CD8+ T Cells",
    Cluster == "CD39_frac" ~ "CD39+ CD8+ T Cells",
    Cluster == "CD39_CD103_frac" ~ "CD39+ CD103+ CD8+ T Cells",
    Cluster == "CD39_CD279_frac" ~ "CD39+ CD279+ CD8+ T Cells",
    Cluster == "CD38pos_DRpos_activated_frac" ~ "CD38+ HLA-DR+ CD8+ T Cells",
    Cluster == "CD45RAneg_CCR7pos_CM_frac" ~ "CM CD8+ T Cells",
    Cluster == "CD45RApos_CCR7pos_Naive_frac" ~ "Naive CD8+ T Cells",
    Cluster == "CD45RApos_CCR7neg_Eff_frac" ~ "Effector CD8+ T Cells",
    Cluster == "CD45RAneg_CCR7neg_EM_frac" ~ "EM CD8+ T Cells",
    TRUE ~ Cluster
  )) %>%
  column_to_rownames(var = "Cluster") 

macro_clusters_flow_heatmap_df <- as.matrix(macro_clusters_flow_heatmap_df)

ht = Heatmap(macro_clusters_flow_heatmap_df, 
             col = col_custom, 
             column_split = as.factor(column_annotation$Arm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(10, "cm"), 
             heatmap_width = unit(20, "cm"),
             cluster_column_slices = F, 
             name = "% of total Lymphocytes",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../flow_CD8_clusters.pdf",
    width = w, height = h)
draw(ht)
dev.off()


# Figure Supplemental 3B ---- 
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

data_for_analysis <- tumor_long %>%
  inner_join(sample_annotation, by = "Sample")

data_for_analysis_cd8 <- data_for_analysis

violin_boxplot <- ggplot(data_for_analysis[data_for_analysis$CellType %in% c("CD152+ CD8+ T Cells",
                                                                             "CD279+ CD8+ T Cells",
                                                                             "EM CD8+ T Cells", 
                                                                             "Naive CD8+ T Cells"), ],
                         aes(x = Arm, y = ZScore, fill = Arm)) + 
  geom_violin(scale = "width", 
              trim = F,
              alpha = 0.5) +
  geom_boxplot(width = 0.05, 
               alpha = .8) + 
  geom_quasirandom(method = "smiley") + 
  facet_wrap(~CellType, nrow = 2, 
             scales = "free") + 
  scale_fill_manual(values = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b")) + 
  labs(title = "", y = "% of Total CD3+CD8+", x = "Arm") + 
  stat_compare_means(comparisons = list(c("Control", "Experimental")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     hide.ns = F)+ 
  stat_compare_means(method = "wilcox") +
  theme_pubr(legend = "right", border = T) + 
  theme(strip.text = element_text(size = 13))

ggsave(violin_boxplot,
       width = 22.5, height = 18, 
       units = "cm", 
       filename = '.../flow_CD8_clusters_boxplot.pdf')

# Get FDRs 
test_results <- data_for_analysis_cd8 %>%
  group_by(CellType) %>%
  wilcox_test(ZScore ~ Arm) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(CellType, group1, group2, p, p.adj, p.adj.signif)

# Figure 3E ----
macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df)
macro_clusters_flow_heatmap_df_SE1_SICE <- macro_clusters_flow_heatmap_df[, colnames(macro_clusters_flow_heatmap_df) %in% unique(c(SE1_patients, SICE_patients))]
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df_SE1_SICE) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

# create  annotation 
column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df_SE1_SICE))
column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 
column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df_SE1_SICE), ]
column_annotation <- column_annotation %>% 
  relocate(Arm)
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

sample_annotation <- sample_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Sarcoma Ecotype Assignment", 
                                       "Sarcoma Immune Class")], 
            by = c("Sample" = "Patient"))

sample_annotation$SICE_vs_SE1_classifier <- NA
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Ecotype Assignment` == "SE1"] <- "SE1"
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Immune Class` == "E"] <- "SICE"
sample_annotation$split_factor <- paste(sample_annotation$Arm, sample_annotation$SICE_vs_SE1_classifier, sep = "_")
sample_annotation$split_factor <- as.factor(sample_annotation$split_factor)
sample_annotation$split_factor <- factor(sample_annotation$split_factor, 
                                         levels = c("Control_SE1", "Experimental_SE1",
                                                    "Control_SICE", "Experimental_SICE"))

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#D682FF", "Grade 3" = "#0095FF"), 
             Arm = c("Control" = "#9D1536", "Experimental" = "#0093AF"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99")))


col_custom <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),viridis(5, option = "C"))

# Heatmap of macro-clusters 
ht = Heatmap(macro_clusters_flow_heatmap_df_SE1_SICE, 
             col = col_custom, 
             column_split = as.factor(sample_annotation$split_factor),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(10, "cm"), 
             heatmap_width = unit(18, "cm"),
             cluster_column_slices = F, 
             name = "% of total CD8+",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../flow_CD8_non_t_regs_clusters_SE1_SICE.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure 3F ----
macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df)
macro_clusters_flow_heatmap_df_SE1_SICE <- macro_clusters_flow_heatmap_df[, colnames(macro_clusters_flow_heatmap_df) %in% unique(c(SE1_patients, SICE_patients))]
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df_SE1_SICE) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

# create  annotation 
column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df_SE1_SICE))
column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 
column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df_SE1_SICE), ]
column_annotation <- column_annotation %>% 
  relocate(Arm)
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

sample_annotation <- sample_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Sarcoma Ecotype Assignment", 
                                       "Sarcoma Immune Class")], 
            by = c("Sample" = "Patient"))

sample_annotation$SICE_vs_SE1_classifier <- NA
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Ecotype Assignment` == "E1"] <- "SE1"
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Immune Class` == "E"] <- "SICE"

data_for_analysis <- tumor_long %>%
  inner_join(sample_annotation, by = "Sample")

macro_violin <- data_for_analysis %>% 
  filter(CellType %in% c("CM CD8+ T Cells", 
                         "EM CD8+ T Cells",
                         "CD279+ CD8+ T Cells",
                         "CD366+ CD8+ T Cells"
  )) %>% 
  ggplot(aes(x = Arm, y = ZScore, fill = Arm)) +
  geom_violin(scale = "width", 
              trim = F,
              alpha = 0.5) +
  geom_boxplot(width = 0.05, 
               #color = "darkblue", 
               alpha = .8) + 
  geom_quasirandom(method = "smiley") + 
  facet_grid(SICE_vs_SE1_classifier ~ CellType,
             #labeller=plot_labeller,
             scales="free_y",
             switch = "y") + 
  scale_fill_manual(values = c("Control" = "#9D1536", "Experimental" = "#0093AF")) + 
  labs(title = "", y = "Z-Scored Cell %", x = "Arm") + 
  stat_compare_means(comparisons = list(c("Control", "Experimental")), 
                     label = "p.signif", 
                     method = "wilcox.test", hide.ns = F)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox") +
  theme_pubr(legend = "right", border = T) + 
  theme(strip.text = element_text(size = 13))

ggsave(macro_violin, units = "cm", 
       width = 27.5, height = 15, 
       filename = '/.../violin_box_CD8_non_t_regs_SE1_vs_SICE.pdf'
)

# Get FDRs 
test_results <- data_for_analysis %>%
  group_by(CellType, SICE_vs_SE1_classifier) %>%
  wilcox_test(ZScore ~ Arm) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(CellType, SICE_vs_SE1_classifier, group1, group2, p, p.adj, p.adj.signif)

# Figure 3G ----
cox_flow_results <- fread(".../data/Cox_T2_Flow.csv")
cox_flow_results <- as.data.frame(cox_flow_results)
cox_flow_results$V1 <- NULL

cox_results_sig <- cox_flow_results[cox_flow_results$p_value < 0.05, ]

# Create data frame for plotting 
dt <- cox_flow_results
dt$est <- log(dt$HR)
dt$low <- log(dt$Lower_CI)
dt$hi <- log(dt$Upper_CI)
dt$est_formatted <- sprintf("%.2f", dt$est)
dt$ci_formatted <- sprintf("(%.2f - %.2f)", dt$low, dt$hi)
dt$Effect_CI <- paste0(dt$est_formatted, " ", dt$ci_formatted)
dt$Cell_Type_copy <- dt$Cell_Type
dt$Cell_Type <- paste(dt$Cell_Type, dt$Arm, sep = "_")
dt <- dt[dt$Cell_Type_copy %in% cox_results_sig$Cell_Type, ]
dt <- dt[, c("Cell_Type", "est", "low", "hi", "se", "p_value", "Effect_CI")]
dt <- dt %>% 
  arrange(Cell_Type)
plot_data <- dt[, c("Cell_Type", "Effect_CI", "p_value")]
colnames(plot_data) <- c("Cell Type", "Hazard Ratio (95% CI)", "P-value")
plot_data[, 4] <- ""
colnames(plot_data)[colnames(plot_data) == "V4"] <- "                                                     "
plot_data$`Cell Type` <- ifelse(grepl("Control|Experimental", plot_data$`Cell Type`), 
                                paste0("   ", plot_data$`Cell Type`), 
                                plot_data$`Cell Type`)

plot_data$`P-value` <- as.numeric(plot_data$`P-value`)

plot_data <- plot_data %>%
  mutate(
    cox_coef = as.numeric(str_extract(`Hazard Ratio (95% CI)`, "-?\\d+\\.?\\d*")),
    signed_log10p = sign(cox_coef) * (-log10(`P-value`))
  )

plot_data <- plot_data %>%
  mutate(
    Category = case_when(
      str_detect(`Cell Type`, "CD4\\+") ~ "CD4+ T cells",
      str_detect(`Cell Type`, "CD8\\+") ~ "CD8+ T cells",
      str_detect(`Cell Type`, "T-regs") ~ "CD4+ T cells",
      TRUE ~ "Other"
    )
  )


plot_data <- plot_data %>%
  mutate(
    Arm = case_when(
      grepl("_Control", `Cell Type`)      ~ "Control",
      grepl("_Experimental", `Cell Type`) ~ "Experimental",
      grepl("_All", `Cell Type`)          ~ "Both",
      TRUE                                ~ "Unknown"
    ),
    `Cell Type` = str_replace_all(`Cell Type`, "_Control|_Experimental|_All", "")
  )

plot_data$`Cell Type` <- gsub("   ", "", plot_data$`Cell Type`)

barplot_double_faceted_T2_flow <- ggplot(
  data = plot_data, 
  aes(y = `Cell Type`, x = signed_log10p, fill = signed_log10p)
) +
  geom_col(width = 0.7, color = "black") +
  
  scale_fill_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick2",
    midpoint = 0,
    name = "Signed -log10(p)\n(blue = better DFS,\n red = worse DFS)"
  ) +
  facet_grid(Category ~ Arm, scales = "free", space = "fixed", switch = "y") +
  
  geom_vline(xintercept = c(-1.3, 0, 1.3), 
             linetype = "dashed", 
             color = "black") +
  
  # ⬇️ move y-axis to the RIGHT
  scale_y_discrete(position = "right") +
  
  theme_pubr(border = TRUE, base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,  
    y = "–log10(p) Cox PH",
    title = "Cox DFS Flow Cytometry T2"
  )

ggsave(barplot_double_faceted_T2_flow, width = 15, height = 20, units = "cm", 
       filename = '.../T2_cox_barplot.pdf')

# Figure Supplemental 3C ----
tregs_clusters <- c(
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count", # denominator 
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD103+ | Count",    # ITGAE                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD152+ | Count",    # CTLA4                          
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD279 + | Count",   # PD-1                           
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD366+ | Count",    # CD366                         
  "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD39+ | Count"     # ENTPD1 
)

flow_results_filtered <- flow_results_filtered %>%
  mutate(across(-1, as.numeric))

flow_results_filtered <- as.data.table(flow_results_filtered)

flow_results_macro_clusters <- flow_results_filtered[, colnames(flow_results_filtered) %in% 
                                                       unique(c("Patient", tregs_clusters)), with = F]

flow_results_macro_clusters <- flow_results_macro_clusters %>%
  mutate(
    CD103_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD103+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count`,
    CD152_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD152+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count`,
    CD279_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD279 + | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count`,
    CD366_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD366+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count`,
    CD39_frac            = `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD39+ | Count` /
      `PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count`
  )

flow_results_macro_clusters <- flow_results_macro_clusters %>% 
  relocate(Patient, .before = CD103_frac)

flow_results_macro_clusters <- flow_results_macro_clusters[, 7:12]

flow_results_macro_clusters <- flow_results_macro_clusters %>% 
  left_join(SIC_and_SE_assignments)

macro_clusters_flow_heatmap_df <- as.matrix(flow_results_macro_clusters[, 1:6])
rownames(macro_clusters_flow_heatmap_df) <- macro_clusters_flow_heatmap_df[, 1]
macro_clusters_flow_heatmap_df <- macro_clusters_flow_heatmap_df[, -1]

macro_clusters_flow_heatmap_df[is.na(macro_clusters_flow_heatmap_df)] <- 0

macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df) %>% 
  mutate(across(1:5, as.numeric))

macro_clusters_flow_heatmap_df <- t(macro_clusters_flow_heatmap_df)

column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df))

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")

column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df), ]

column_annotation <- column_annotation %>% 
  relocate(Arm)

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
 col = list(`Tumor Grade` = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             Arm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))


col_custom <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),viridis(5, option = "C"))

macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df) %>%
  rownames_to_column(var = "Cluster") %>%
  mutate(Cluster = case_when(
    Cluster == "CD103_frac"    ~ "CD103+ T-regs",
    Cluster == "CD152_frac"  ~ "CD152+ T-regs",
    Cluster == "CD279_frac"  ~ "CD279+ T-regs",
    Cluster == "CD366_frac"   ~ "CD366+ T-regs",
    Cluster == "CD39_frac" ~ "CD39+ T-regs",
    TRUE ~ Cluster
  )) %>%
  column_to_rownames(var = "Cluster") 

macro_clusters_flow_heatmap_df <- as.matrix(macro_clusters_flow_heatmap_df)

# Heatmap of macro-clusters 
ht = Heatmap(macro_clusters_flow_heatmap_df, 
             col = col_custom, 
             column_split = as.factor(column_annotation$Arm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(10, "cm"), 
             heatmap_width = unit(20, "cm"),
             cluster_column_slices = F, 
             name = "% of total T-regs",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../flow_tregs_clusters.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure Supplemental 3D ----
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

data_for_analysis <- tumor_long %>%
  inner_join(sample_annotation, by = "Sample")

data_for_analysis_tregs <- data_for_analysis

violin_boxplot <- ggplot(data_for_analysis[data_for_analysis$CellType != "CD366+ T-regs", ],
                         aes(x = Arm, y = ZScore, fill = Arm)) +
  geom_violin(scale = "width", 
              trim = F,
              alpha = 0.5) +
  geom_boxplot(width = 0.05, 
               alpha = .8) + 
  geom_quasirandom(method = "smiley") + 
  facet_wrap(~CellType, nrow = 2, 
             scales = "free") + 
  scale_fill_manual(values = c("Experimental" = "#61a5c2", "Control" = "#9e2a2b")) + 
  labs(title = "", y = "% of Total T-regs", x = "Arm") + 
  stat_compare_means(comparisons = list(c("Control", "Experimental")), 
                     label = "p.signif", 
                     method = "wilcox.test", 
                     hide.ns = F)+ 
  stat_compare_means(method = "wilcox") +
  theme_pubr(legend = "right", border = T) +
  theme(strip.text = element_text(size = 13))

ggsave(violin_boxplot,
       width = 32.5, height =18, 
       units = "cm", 
       filename = '.../flow_T-regs_clusters_boxplot.pdf')

test_results <- data_for_analysis_tregs %>%
  group_by(CellType) %>%
  wilcox_test(ZScore ~ Arm) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(CellType, group1, group2, p, p.adj, p.adj.signif)

# Figure Supplemental 3E ----
macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df)
macro_clusters_flow_heatmap_df_SE1_SICE <- macro_clusters_flow_heatmap_df[, colnames(macro_clusters_flow_heatmap_df) %in% unique(c(SE1_patients, SICE_patients))]
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df_SE1_SICE) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

# create  annotation 
column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df_SE1_SICE))
column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 
column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df_SE1_SICE), ]
column_annotation <- column_annotation %>% 
  relocate(Arm)
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

sample_annotation <- sample_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Sarcoma Ecotype Assignment", 
                                       "Sarcoma Immune Class")], 
            by = c("Sample" = "Patient"))

sample_annotation$SICE_vs_SE1_classifier <- NA
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Ecotype Assignment` == "E1"] <- "SE1"
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Immune Class` == "E"] <- "SICE"
sample_annotation$split_factor <- paste(sample_annotation$Arm, sample_annotation$SICE_vs_SE1_classifier, sep = "_")
sample_annotation$split_factor <- as.factor(sample_annotation$split_factor)
sample_annotation$split_factor <- factor(sample_annotation$split_factor, 
                                         levels = c("Control_SE1", "Experimental_SE1",
                                                    "Control_SICE", "Experimental_SICE"))

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#D682FF", "Grade 3" = "#0095FF"), 
             Arm = c("Control" = "#9D1536", "Experimental" = "#0093AF"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))


col_custom <- colorRamp2(c(0, 0.25, 0.5, 0.75, 1),viridis(5, option = "C"))

# Heatmap of macro-clusters 
ht = Heatmap(macro_clusters_flow_heatmap_df_SE1_SICE, 
             col = col_custom, 
             column_split = as.factor(sample_annotation$split_factor),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(10, "cm"), 
             heatmap_width = unit(18, "cm"),
             cluster_column_slices = F, 
             name = "% of total T-regs",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../flow_tregs_clusters_SE1_SICE.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure Supplemental 3F ----
macro_clusters_flow_heatmap_df <- as.data.frame(macro_clusters_flow_heatmap_df)
macro_clusters_flow_heatmap_df_SE1_SICE <- macro_clusters_flow_heatmap_df[, colnames(macro_clusters_flow_heatmap_df) %in% unique(c(SE1_patients, SICE_patients))]
tumor_long <- as.data.frame(macro_clusters_flow_heatmap_df_SE1_SICE) %>%
  rownames_to_column(var = "CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "ZScore")

# create  annotation 
column_annotation <- data.frame(samples = colnames(macro_clusters_flow_heatmap_df_SE1_SICE))
column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% flow_results_macro_clusters[flow_results_macro_clusters$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 
column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Histology")],
            by = c("samples" = "Patient"))

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(macro_clusters_flow_heatmap_df_SE1_SICE), ]
column_annotation <- column_annotation %>% 
  relocate(Arm)
sample_annotation <- column_annotation
sample_annotation <- rownames_to_column(sample_annotation, var = "Sample")

sample_annotation <- sample_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", 
                                       "Sarcoma Ecotype Assignment", 
                                       "Sarcoma Immune Class")], 
            by = c("Sample" = "Patient"))

sample_annotation$SICE_vs_SE1_classifier <- NA
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Ecotype Assignment` == "E1"] <- "SE1"
sample_annotation$SICE_vs_SE1_classifier[sample_annotation$`Sarcoma Immune Class` == "E"] <- "SICE"

data_for_analysis <- tumor_long %>%
  inner_join(sample_annotation, by = "Sample")

macro_violin <- data_for_analysis %>% 
  filter(CellType != "CD103+ T-regs") %>%
  ggplot(aes(x = Arm, y = ZScore, fill = Arm)) +
  geom_violin(scale = "width", 
              trim = F,
              alpha = 0.5) +
  geom_boxplot(width = 0.05, 
               alpha = .8) + 
  geom_quasirandom(method = "smiley") + 
  facet_grid(SICE_vs_SE1_classifier ~ CellType,
             scales="free_y",
             switch = "y") + 
  scale_fill_manual(values = c("Control" = "#9D1536", "Experimental" = "#0093AF")) + 
  labs(title = "", y = "Z-Scored Cell %", x = "Arm") + 
  stat_compare_means(comparisons = list(c("Control", "Experimental")), 
                     label = "p.signif", 
                     method = "wilcox.test", hide.ns = F)+ 
  stat_compare_means(method = "wilcox") +
  theme_pubr(legend = "right", border = T) + 
  theme(strip.text = element_text(size = 13))

ggsave(macro_violin, units = "cm", 
       width = 27.5, height = 15, 
       filename = '...y/violin_box_t_regs_SE1_vs_SICE.pdf'
)

# Get FDRs 
test_results <- data_for_analysis %>%
  group_by(CellType, SICE_vs_SE1_classifier) %>%
  wilcox_test(ZScore ~ Arm) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(CellType, SICE_vs_SE1_classifier, 
         group1, group2, 
         p, p.adj, p.adj.signif)

# Figure 4A ----
cd45_umap <- readRDS(".../data/data.cd45.sce.umap.RDS")

umap_mat <- reducedDim(cd45_umap, "UMAP2_on_PCA")
colnames(umap_mat) <- c("UMAP1","UMAP2")

df <- as.data.frame(umap_mat)
df$cluster   <- colData(cd45_umap)$cellCluster
df$sample_id <- colData(cd45_umap)$sample_id

df_sub <- df %>%
  slice_sample(n = 20000, replace = FALSE)

clusters <- sort(unique(df_sub$cluster))
clusters <- clusters[!is.na(clusters)]

PBMCs_palette <- c("#4292C6", "#67000D", "#BCBDDC", 
                   "#255668", "#54278F", "#006D2C", 
                   "#EB9619B3", "#41AB5D", "#35978F")

PBMCs_palette_mapping <- setNames(PBMCs_palette, clusters)

PBMCs_Plot <- ggplot(
  df_sub[!is.na(df_sub$cluster), ],
  aes(x = UMAP1, y = UMAP2, colour = cluster)
) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(
    values = PBMCs_palette_mapping,
    na.value = "grey70"   
  ) +
  labs(title = "UMAP (random 20k cells)") +
  theme_void()

ggsave(PBMCs_Plot, filename = '.../PBMCs_umap.pdf',
       width = 15, height = 10, units = "cm")

# t-cell clusters 
tcell_umap <- readRDS(".../data/tcell.sce.umap.RDS")

umap_mat <- reducedDim(tcell_umap, "UMAP2_on_PCA")
colnames(umap_mat) <- c("UMAP1","UMAP2")

df <- as.data.frame(umap_mat)
df$cluster   <- colData(tcell_umap)$cellCluster
df$sample_id <- colData(tcell_umap)$sample_id
df$toplot <- colData(tcell_umap)$ToPlot

df_sub <- df %>%
  slice_sample(n = 20000, replace = FALSE)

clusters <- sort(unique(df_sub$cluster))
clusters <- clusters[!is.na(clusters)]

tcells_palette <- sample(gnomeR::gnomer_palettes$main, length(clusters))

tcells_palette <- c("#c35abc",
                    "#76b341",
                    "#7560cd",
                    "#d09c45",
                    "#7882c9",
                    "#7f7f37",
                    "#d0406d",
                    "#51a976",
                    "#c85d3f",
                    "#45b0cf",
                    "#be6c91")

tcells_palette_mapping <- setNames(tcells_palette, clusters)

tcells_Plot <- ggplot(
  df_sub[!is.na(df_sub$cluster), ],
  aes(x = UMAP1, y = UMAP2, colour = cluster)
) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_color_manual(
    values = tcells_palette_mapping,
    na.value = "grey70"   
  ) +
  labs(title = "UMAP (random 20k cells)") +
  theme_void()

ggsave(tcells_Plot, filename = '.../plots/tcells_umap.pdf',
       width = 15, height = 10, units = "cm")

# Figure 4B ----
percentmatlist <- readRDS('.../data/percent.mat.list.RDS')
cd45_expression <- cd45_umap@assays@data@listData$exprs

sce      <- cd45_umap
expr_mat <- assay(sce, "exprs")           
clusters <- colData(sce)$cellCluster       

cluster_ids <- sort(unique(clusters))
cluster_ids <- cluster_ids[!is.na(cluster_ids)]

# Average expression per marker × cluster
expr_avg_mat <- sapply(cluster_ids, function(cl) {
  sel <- clusters == cl
  rowMeans(expr_mat[, sel], na.rm = TRUE)
})

colnames(expr_avg_mat) <- paste0("Cluster", cluster_ids)

expr_avg_df <- expr_avg_mat %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  pivot_longer(
    cols      = -marker,
    names_to  = "cluster",
    values_to = "avg_expr"
  )

# Average % positive per marker × cluster
percent_avg_mat <- sapply(percentmatlist, function(mat) {
  rowMeans(mat, na.rm = TRUE)
})
rownames(percent_avg_mat) <- sub("^x", "", rownames(percent_avg_mat))

percent_avg_df <- percent_avg_mat %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  pivot_longer(
    cols      = -marker,
    names_to  = "cluster",
    values_to = "avg_percent"
  )

dot_df <- inner_join(percent_avg_df, expr_avg_df, by = c("marker","cluster"))

# z-score average expression per marker across clusters 
dot_df_scaled <- dot_df %>%
  group_by(marker) %>%
  mutate(
    expr_mean = mean(avg_expr, na.rm = TRUE),
    expr_sd   = sd(avg_expr, na.rm = TRUE),
    scaled_expr = if_else(expr_sd > 0,
                          (avg_expr - expr_mean) / expr_sd,
                          0)
  ) %>%
  ungroup()

dotplot_cytof <- ggplot(dot_df_scaled, aes(
  x    = marker,
  y    = cluster,
  size = avg_percent,
  fill = scaled_expr
)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +  
  scale_size_continuous(range = c(1,8), name = "Avg % positive") +
  scale_fill_distiller(
    palette   = "RdBu",
    direction = -1,
    name      = "Scaled expr",
    limits    = c(-max(abs(dot_df_scaled$scaled_expr), na.rm = TRUE),
                  max(abs(dot_df_scaled$scaled_expr), na.rm = TRUE)),
    oob       = scales::squish
  ) +
  cowplot::theme_minimal_grid(line_size = 0.1, color = "black") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank()
  ) +
  ggtitle("Marker positivity and (scaled) expression by cluster")

ggsave(dotplot_cytof, filename = '.../plots/dotplot_Figure4b.pdf', 
       width = 32.5, height = 15, units = "cm")

# Figure 4C ----
phenotype_CyTOF <- readRDS('.../data/CyTOF/phenoData.df.RDS')
cluster_all_cells <- readRDS('.../data/CyTOF/clusterAssignments.RDS')
percent_all_cells <- readRDS('.../data/CyTOF/percent.mat.RDS')

rownames(percent_all_cells) <- cluster_all_cells[ rownames(percent_all_cells) ]

phenotype_CyTOF$PatientID <- gsub("x", "", phenotype_CyTOF$PatientID)
phenotype_CyTOF$PatientID <- gsub("_", "-", phenotype_CyTOF$PatientID)

phenotype_CyTOF <- phenotype_CyTOF %>% 
  left_join(sarc32_outcome, by = c("PatientID" = "Subject"))

phenotype_CyTOF$Visit <- factor(
  phenotype_CyTOF$Visit,
  levels = c(
    "Prior to 1st Tx", # "Pre-treatment" 
    "1 week Post Initiating treatment",  # "1 week on-treatment" 
    "During RT",  # "Week 3 of RT" 
    "Before Surgery",  # "Post-RT/pre-surgery" 
    "3 mo Post-Surgery", # "3 mo post-surgery" 
    "12 mo Post-Surgery" # "12 mo post-surgery" 
  )
)

# remove non evaluable patients
phenotype_CyTOF <- phenotype_CyTOF %>% 
  filter(PatientID %in% SIC_and_SE_assignments$Patient)

# load blood draw timing information 
blood_draw_timing <- fread('.../data/SARC032_blood_draw_timing.csv')

# Adjust visit labels for joining
blood_draw_timing$collection_time_point[blood_draw_timing$collection_time_point == "Baseline"] <- "Prior to 1st Tx"
blood_draw_timing$collection_time_point[blood_draw_timing$collection_time_point == "30 days post-cycle 6 (3 months after surgery)"] <- "3 mo Post-Surgery"
blood_draw_timing$collection_time_point[blood_draw_timing$collection_time_point == "Prior to Cycle 2 (before week 3 of radiotherapy)"] <- "1 week Post Initiating treatment During RT"
blood_draw_timing$collection_time_point[blood_draw_timing$collection_time_point == "Prior to surgery"] <- "Before Surgery"
blood_draw_timing$collection_time_point[blood_draw_timing$collection_time_point == "Final Sample"] <- "12 mo Post-Surgery"

# 1) Blood draw timing (from enrollment -> blood draw)
blood_timing <- blood_draw_timing %>%
  transmute(
    PatientID = Subject,
    Visit_key = collection_time_point,
    landmark_days = as.numeric(time_from_enrollment_to_blood_draw_days)
  )

blood_resolved <- blood_timing %>%
  filter(!is.na(Visit_key)) %>%
  arrange(PatientID, Visit_key, landmark_days) %>%
  group_by(PatientID, Visit_key) %>%
  mutate(dup_rank = row_number(),
         chosen   = dup_rank == 1) %>%       
  ungroup()

blood_timing_one_per <- blood_resolved %>%
  filter(chosen) %>%
  select(PatientID, Visit_key, landmark_days)

df0 <- phenotype_CyTOF %>%
  mutate(
    Visit_key = Visit,
    event  = `DFS event` == "Event occurred",
    t_event = as.numeric(`DFS time (days)`)
  ) %>%
  left_join(blood_timing_one_per, by = c("PatientID","Visit_key"))

percent_all_cells <- percent_all_cells[, colnames(percent_all_cells) %in% phenotype_CyTOF$Sample]

# CLR normalization 
clr_mat <- compositions::clr(as.matrix(t(percent_all_cells)))     

percent_all_cells <- clr_mat %>% 
  as.data.frame() %>%
  rownames_to_column(var = "CellType")

# Pivot to long format
percent_all_cells_long <- percent_all_cells %>%
  pivot_longer(
    cols = -CellType,           
    names_to = "Sample",        
    values_to = "Abundance"       
  )

names(percent_all_cells_long) <- c("Sample", "CellType", "Abundance")

percent_all_cells_long <- percent_all_cells_long %>% 
  left_join(df0[, c("PatientID" , "Sample", 
                    "Visit_key", "event", "t_event", "landmark_days",
                    "Collection_Date", "Visit")], by = c("Sample"))

percent_all_cells_long <- percent_all_cells_long %>% 
  left_join(SIC_and_SE_assignments, by = c("PatientID" = "Patient"))

percent_all_cells_long <- as.data.frame(percent_all_cells_long)

percent_all_cells_long$Visit <- factor(
  percent_all_cells_long$Visit,
  levels = c(
    "Prior to 1st Tx", 
    "1 week Post Initiating treatment", 
    "During RT", 
    "Before Surgery", 
    "3 mo Post-Surgery", 
    "12 mo Post-Surgery"
  )
)

# PBMCs line plot 
percent_all_cells_long$Arm <- percent_all_cells_long$`Treatment Arm`
line_blood <- ggline(
  data = percent_all_cells_long,
  x = "Visit",
  y = "Abundance",
  color = "Arm",
  palette = c("Control" = "#9D1536", "Experimental" = "#0093AF"),
  facet.by = "CellType",
  add = "mean_se",      
  scales = "free_y", 
  x.text.angle = 45, 
  xlab = "Time Points", 
  ylab = "Fraction of total PBMCs"
) + 
  theme(strip.text = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14))

ggsave(line_blood, 
       filename = '/.../time_series_PBMCs_CLR.pdf', 
       units = "cm", width = 32, height = 25)

# Wilcoxon test - between Arms comparisons 
stats_table_between_arms <- compare_means(
  Abundance ~ Arm,
  data = percent_all_cells_long,
  group.by = c("CellType", "Visit"),      
  method = "wilcox.test",      
  paired = FALSE          
)

# Adjust p-values considering all the cell types tested at each time point 
stats_table_between_arms <- stats_table_between_arms %>% 
  dplyr::group_by(Visit) %>% 
  mutate(p.adj = p.adjust(p, method = "BH"))

# Wilcoxon between time points, within Arms
stats_table_between_visit_withinArms <- compare_means(
  Abundance ~ Visit,
  data = percent_all_cells_long,
  group.by = c("CellType", "Arm"),      
  method = "wilcox.test",      
  paired = FALSE               
)

# Adjust p-values considering all the cell types tested at each time point 
stats_table_between_visit_withinArms <- stats_table_between_visit_withinArms %>% 
  dplyr::group_by(Arm) %>% 
  mutate(p.adj = p.adjust(p, method = "BH"))

# Wilcoxon between time points - All patients 
stats_table_between_visit <- compare_means(
  Abundance ~ Visit,
  data = percent_all_cells_long,
  group.by = c("CellType"),      
  method = "wilcox.test",      
  paired = FALSE               
)

# Adjust p-values considering all the cell types tested at each time point 
stats_table_between_visit <- stats_table_between_visit %>% 
  dplyr::group_by(CellType) %>% 
  mutate(p.adj = p.adjust(p, method = "BH"))

# Figure 4D ----
## --- 0) Build landmark time and inclusion flag
percent_all_cells_long <- percent_all_cells_long %>%
  mutate(
    time_since_landmark = t_event - landmark_days,
    include = !is.na(landmark_days) & !is.na(t_event) & (t_event > landmark_days)
  )

# (Optional) QA: see how many are kept/excluded per visit/arm
qa_counts <- percent_all_cells_long %>%
  group_by(Visit, Arm) %>%
  summarize(
    n_total = n_distinct(Sample),
    n_kept  = n_distinct(Sample[include]),
    n_excl_preL = n_distinct(Sample[!include & event & !is.na(landmark_days) & !is.na(t_event) & t_event <= landmark_days]),
    .groups = "drop"
  )
print(qa_counts)

## --- 1) Cox helper (landmark-correct)
fit_cox <- function(dat) {
  dat <- dat %>% filter(include)
  if (nrow(dat) < 5 || dplyr::n_distinct(dat$Abundance) < 2 || sum(dat$event) < 1) return(NULL)
  fit <- coxph(Surv(time_since_landmark, event) ~ Abundance, data = dat)
  s   <- summary(fit)
  tibble(
    HR        = s$conf.int[1, "exp(coef)"],
    Lower_CI  = s$conf.int[1, "lower .95"],
    Upper_CI  = s$conf.int[1, "upper .95"],
    coef      = s$coefficients[1, "coef"],
    se        = s$coefficients[1, "se(coef)"],
    p_value   = s$coefficients[1, "Pr(>|z|)"],
    n_included = nrow(dat),
    n_events   = sum(dat$event)
  )
}

cell_types <- sort(unique(percent_all_cells_long$CellType))
timepoints <- levels(percent_all_cells_long$Visit) %||% unique(percent_all_cells_long$Visit)

## --- 2) Loop exactly like before, but using time_since_landmark
cox_results <- tibble(
  Cell_Type = character(), Arm = character(), Timepoint = character(),
  HR = numeric(), Lower_CI = numeric(), Upper_CI = numeric(),
  coef = numeric(), se = numeric(), p_value = numeric(),
  n_included = integer(), n_events = integer()
)

for (cell in cell_types) {
  for (time in timepoints) {
    subset_data <- percent_all_cells_long %>%
      filter(CellType == cell, Visit == time)
    
    # All patients (pooled)
    res_all <- fit_cox(subset_data)
    if (!is.null(res_all)) {
      cox_results <- bind_rows(
        cox_results,
        res_all %>% mutate(Cell_Type = cell, Arm = "All", Timepoint = time, .before = 1)
      )
    } else {
      message(sprintf("Skip ALL: %s @ %s (insufficient variation/events)", cell, time))
    }
    
    # By Arm
    for (arm in c("Control","Experimental")) {
      res_arm <- fit_cox(subset_data %>% filter(Arm == arm))
      if (!is.null(res_arm)) {
        cox_results <- bind_rows(
          cox_results,
          res_arm %>% mutate(Cell_Type = cell, Arm = arm, Timepoint = time, .before = 1)
        )
      } else {
        message(sprintf("Skip %s: %s @ %s (insufficient variation/events)", arm, cell, time))
      }
    }
  }
}

# BH-FDR within each (Visit, Arm) slice (same rule you used)
cox_results <- cox_results %>%
  group_by(Timepoint, Arm) %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

cox_results_sig <- cox_results %>% filter(p_value < 0.05)

## --- 3) Heatmap (unchanged, built from the new cox_results)
cox_results <- cox_results %>%
  mutate(
    signed_log10p = sign(coef) * (-log10(p_value)),
    time_arm = paste(Arm, Timepoint, sep = "_")
  )

wide_df <- cox_results %>%
  select(Cell_Type, time_arm, signed_log10p) %>%
  pivot_wider(names_from = time_arm, values_from = signed_log10p)

mat <- wide_df %>%
  column_to_rownames("Cell_Type") %>%
  as.matrix()

col_anno_df <- data.frame(colnames(mat)) %>%
  setNames("colname")
col_anno_df$Arm  <- stringr::word(col_anno_df$colname, 1, sep = "_")
col_anno_df$Time <- stringr::word(col_anno_df$colname, 2, sep = "_")
rownames(col_anno_df) <- col_anno_df$colname
col_anno_df <- col_anno_df[colnames(mat), , drop = FALSE]

col_anno_df$Time <- factor(
  col_anno_df$Time,
  levels = c("Prior to 1st Tx",
             "1 week Post Initiating treatment",
             "During RT",
             "Before Surgery",
             "3 mo Post-Surgery",
             "12 mo Post-Surgery")
)
col_anno_df$Arm <- factor(col_anno_df$Arm, levels = c("Control","Experimental","All"))

ord <- unlist(lapply(levels(col_anno_df$Arm), function(a){
  idx <- which(col_anno_df$Arm == a)
  idx[order(col_anno_df$Time[idx])]
}))
mat_ord         <- mat[, ord, drop = FALSE]
col_anno_df_ord <- col_anno_df[ord, , drop = FALSE]
col_anno_df_ord$colname <- NULL

annotation_column <- HeatmapAnnotation(
  df = col_anno_df_ord,
  col = list(
    Arm = c("Control" = "#9D1536", "Experimental" = "#0093AF", "All" = "grey"),
    Time = c("Prior to 1st Tx" = "#CC79A7",
             "1 week Post Initiating treatment" = "#0072B2",
             "Before Surgery" = "#56B4E9",
             "During RT" = "#009E73",
             "3 mo Post-Surgery"  = "#F5C710",
             "12 mo Post-Surgery" = "#E69F00")
  ),
  which = "column",
  show_legend = TRUE, show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_ord,
  top_annotation = annotation_column,
  column_split = as.factor(col_anno_df_ord$Arm),
  show_column_names = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  name = "signed_p",
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(15, "cm"),
  col = colorRamp2(c(min(mat, na.rm=TRUE), 0, max(mat, na.rm=TRUE)),
                   c("#4575b4", "white", "#d73027"))
)

ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../cox_heatmap_PBMCs_CLR.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure 4E ----
cluster_T_cells <- readRDS('.../data/CyTOF/tcell.clusterAssignments.RDS')
percent_T_cells <- readRDS('.../data/CyTOF/tcell.percent.mat.RDS')

rownames(percent_T_cells) <- cluster_T_cells[ rownames(percent_T_cells) ]

clr_mat <- compositions::clr(as.matrix(t(percent_T_cells)))        # rows=samples, cols=cell subsets

percent_T_cells <- t(clr_mat) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "CellType")

# Pivot to long format
percent_T_cells_long <- percent_T_cells %>%
  pivot_longer(
    cols = -CellType,            # All columns except 'CellType' 
    names_to = "Sample",         # Column to store old column names (sample IDs)
    values_to = "Abundance"        # Column to store the values
  )

percent_T_cells_long <- percent_T_cells_long %>% 
  left_join(df0[, c("PatientID" , "Sample", 
                    "Visit_key", "event", "t_event", "landmark_days",
                    "Collection_Date", "Visit")], by = c("Sample"))

percent_T_cells_long <- percent_T_cells_long %>% 
  left_join(SIC_and_SE_assignments, by = c("PatientID" = "Patient"))

percent_T_cells_long <- as.data.frame(percent_T_cells_long)

percent_T_cells_long$Visit <- factor(
  percent_T_cells_long$Visit,
  levels = c(
    "Prior to 1st Tx", 
    "1 week Post Initiating treatment", 
    "During RT", 
    "Before Surgery", 
    "3 mo Post-Surgery", 
    "12 mo Post-Surgery"
  )
)

percent_T_cells_long <- percent_T_cells_long %>% 
  filter(PatientID %in% SIC_and_SE_assignments$Patient)

percent_T_cells_long$Arm <- percent_T_cells_long$`Treatment Arm`

# Line plot T cells 
line_t_cells <- ggline(
  data = percent_T_cells_long,
  x = "Visit",
  y = "Abundance",
  color = "Arm",
  palette = c("Control" = "#9D1536", "Experimental" = "#0093AF"),
  facet.by = "CellType",
  add = "mean_se",      
  scales = "free_y", # Add error bars for mean ± SD
  x.text.angle = 45, 
  ncol = 3, 
  xlab = "Time Points", 
  ylab = "Fraction of total T-cells"
) + 
  theme(strip.text = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14))

ggsave(line_t_cells, 
       filename = '/.../time_series_T_cells_CLR.pdf', 
       units = "cm", width = 32.5, height = 35.5)

# Wilcoxon test between Arms 
stats_table_between_arms <- compare_means(
  Abundance ~ Arm,
  data = percent_T_cells_long,
  group.by = c("CellType", "Visit"),     
  method = "wilcox.test",      
  paired = FALSE               
)

stats_table_between_arms <- stats_table_between_arms %>% 
  dplyr::group_by(Visit) %>% 
  mutate(p.adj = p.adjust(p, method = "BH"))

# Wilcoxon test between time points, within Arms
stats_table_between_visit_withinArms <- compare_means(
  Abundance ~ Visit,
  data = percent_T_cells_long,
  group.by = c("CellType", "Arm"),      
  method = "wilcox.test",      
  paired = FALSE               
)

stats_table_between_visit_withinArms <- stats_table_between_visit_withinArms %>% 
  dplyr::group_by(Arm) %>% 
  mutate(p.adj = p.adjust(p, method = "BH"))

# Wilcoxon test between time points - All patients 
stats_table_between_visit <- compare_means(
  Abundance ~ Visit,
  data = percent_T_cells_long,
  group.by = c("CellType"),     
  method = "wilcox.test",      
  paired = FALSE              
)

# Adjust p-values considering all the cell types tested at each time point 
stats_table_between_visit <- stats_table_between_visit %>% 
  dplyr::group_by(CellType) %>% 
  mutate(p.adj = p.adjust(p, method = "BH"))


# Figure 4F ----
percent_T_cells_long <- percent_T_cells_long %>%
  mutate(
    time_since_landmark = t_event - landmark_days,
    include = !is.na(landmark_days) & !is.na(t_event) & (t_event > landmark_days)
  )

qa_counts <- percent_T_cells_long %>%
  group_by(Visit, Arm) %>%
  summarize(
    n_total = n_distinct(Sample),
    n_kept  = n_distinct(Sample[include]),
    n_excl_preL = n_distinct(Sample[!include & event & !is.na(landmark_days) & !is.na(t_event) & t_event <= landmark_days]),
    .groups = "drop"
  )
print(qa_counts)

# Run Cox analysis 
fit_cox <- function(dat) {
  dat <- dat %>% filter(include)
  if (nrow(dat) < 5 || dplyr::n_distinct(dat$Abundance) < 2 || sum(dat$event) < 1) return(NULL)
  fit <- coxph(Surv(time_since_landmark, event) ~ Abundance, data = dat)
  s   <- summary(fit)
  tibble(
    HR        = s$conf.int[1, "exp(coef)"],
    Lower_CI  = s$conf.int[1, "lower .95"],
    Upper_CI  = s$conf.int[1, "upper .95"],
    coef      = s$coefficients[1, "coef"],
    se        = s$coefficients[1, "se(coef)"],
    p_value   = s$coefficients[1, "Pr(>|z|)"],
    n_included = nrow(dat),
    n_events   = sum(dat$event)
  )
}

cell_types <- sort(unique(percent_T_cells_long$CellType))
timepoints <- levels(percent_T_cells_long$Visit) %||% unique(percent_T_cells_long$Visit)

cox_results <- tibble(
  Cell_Type = character(), Arm = character(), Timepoint = character(),
  HR = numeric(), Lower_CI = numeric(), Upper_CI = numeric(),
  coef = numeric(), se = numeric(), p_value = numeric(),
  n_included = integer(), n_events = integer()
)

for (cell in cell_types) {
  for (time in timepoints) {
    subset_data <- percent_T_cells_long %>%
      filter(CellType == cell, Visit == time)
    
    # All patients (pooled)
    res_all <- fit_cox(subset_data)
    if (!is.null(res_all)) {
      cox_results <- bind_rows(
        cox_results,
        res_all %>% mutate(Cell_Type = cell, Arm = "All", Timepoint = time, .before = 1)
      )
    } else {
      message(sprintf("Skip ALL: %s @ %s (insufficient variation/events)", cell, time))
    }
    
    # By Arm
    for (arm in c("Control","Experimental")) {
      res_arm <- fit_cox(subset_data %>% filter(Arm == arm))
      if (!is.null(res_arm)) {
        cox_results <- bind_rows(
          cox_results,
          res_arm %>% mutate(Cell_Type = cell, Arm = arm, Timepoint = time, .before = 1)
        )
      } else {
        message(sprintf("Skip %s: %s @ %s (insufficient variation/events)", arm, cell, time))
      }
    }
  }
}

# BH-FDR within each (Visit, Arm)
cox_results <- cox_results %>%
  group_by(Timepoint, Arm) %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "BH")) %>%
  ungroup()

cox_results_sig <- cox_results %>% filter(p_value < 0.05)

cox_results <- cox_results %>%
  mutate(
    signed_log10p = sign(coef) * (-log10(p_value)),
    time_arm = paste(Arm, Timepoint, sep = "_")
  )

wide_df <- cox_results %>%
  select(Cell_Type, time_arm, signed_log10p) %>%
  pivot_wider(names_from = time_arm, values_from = signed_log10p)

mat <- wide_df %>%
  column_to_rownames("Cell_Type") %>%
  as.matrix()

col_anno_df <- data.frame(colnames(mat)) %>%
  setNames("colname")
col_anno_df$Arm  <- stringr::word(col_anno_df$colname, 1, sep = "_")
col_anno_df$Time <- stringr::word(col_anno_df$colname, 2, sep = "_")
rownames(col_anno_df) <- col_anno_df$colname
col_anno_df <- col_anno_df[colnames(mat), , drop = FALSE]

col_anno_df$Time <- factor(
  col_anno_df$Time,
  levels = c("Prior to 1st Tx",
             "1 week Post Initiating treatment",
             "During RT",
             "Before Surgery",
             "3 mo Post-Surgery",
             "12 mo Post-Surgery")
)
col_anno_df$Arm <- factor(col_anno_df$Arm, levels = c("Control","Experimental","All"))

ord <- unlist(lapply(levels(col_anno_df$Arm), function(a){
  idx <- which(col_anno_df$Arm == a)
  idx[order(col_anno_df$Time[idx])]
}))
mat_ord         <- mat[, ord, drop = FALSE]
col_anno_df_ord <- col_anno_df[ord, , drop = FALSE]
col_anno_df_ord$colname <- NULL

annotation_column <- HeatmapAnnotation(
  df = col_anno_df_ord,
  col = list(
    Arm = c("Control" = "#9D1536", "Experimental" = "#0093AF", "All" = "grey"),
    Time = c("Prior to 1st Tx" = "#CC79A7",
             "1 week Post Initiating treatment" = "#0072B2",
             "Before Surgery" = "#56B4E9",
             "During RT" = "#009E73",
             "3 mo Post-Surgery"  = "#F5C710",
             "12 mo Post-Surgery" = "#E69F00")
  ),
  which = "column",
  show_legend = TRUE, show_annotation_name = FALSE
)

ht <- Heatmap(
  mat_ord,
  top_annotation = annotation_column,
  column_split = as.factor(col_anno_df_ord$Arm),
  show_column_names = FALSE,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  name = "signed_p",
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(15, "cm"),
  col = colorRamp2(c(min(mat, na.rm=TRUE), 0, max(mat, na.rm=TRUE)),
                   c("#4575b4", "white", "#d73027"))
)

ht <- draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../cox_heatmap_Tcells_CLR.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# CIBERSORTx TR4-correction of scRNA-seq deconvolution -----
evaluable_patients <- SIC_and_SE_assignments$Patient

TR4 <- fread('.../data/TR4_CPM_final.csv')
TR4$V1 <- NULL
TR4$V1 <- NULL
cibersortx_homegrown_UMI_new_cd45neg <- fread('...data/scRNA_seq_CIBERSORTx_non_immune_cells_final.csv')
cibersortx_homegrown_UMI_new_cd45neg$V1 <- NULL
cibersortx_homegrown_UMI_new_cd45pos <- fread('...data/scRNA_seq_CIBERSORTx_immune_cells_final.csv')
cibersortx_homegrown_UMI_new_cd45pos$V1 <- NULL

TR4 <- TR4 %>% select(-`P-value`, -Correlation, -RMSE)
TR4$CD45neg <- (TR4$CD10 + TR4$CD31 + TR4$EPCAM)
cibersortx_homegrown_UMI_new_cd45neg <- cibersortx_homegrown_UMI_new_cd45neg %>% select(-`P-value`, -Correlation, -RMSE)
cibersortx_homegrown_UMI_new_cd45pos <- cibersortx_homegrown_UMI_new_cd45pos %>% select(-`P-value`, -Correlation, -RMSE)

combo <- TR4 %>% inner_join(cibersortx_homegrown_UMI_new_cd45neg, by = "Mixture")
combo <- combo %>% inner_join(cibersortx_homegrown_UMI_new_cd45pos, by = "Mixture")

combo_adjusted <- combo %>%
  mutate(across(all_of(setdiff(names(cibersortx_homegrown_UMI_new_cd45pos), c("Mixture", "timepoint", "patient_id"))), ~ .x * CD45)) %>% 
  mutate(across(all_of(setdiff(names(cibersortx_homegrown_UMI_new_cd45neg), c("Mixture", "timepoint", "patient_id"))), ~ .x * CD45neg))

combo_adjusted <- combo_adjusted %>% select(-c(CD45 ,CD45neg, 
                                               CD10, CD31, 
                                               EPCAM))

# retain only evaluable patients 
combo_adjusted <- combo_adjusted %>% 
  filter(patient_id %in% SIC_and_SE_assignments$Patient)

combo_adjusted$timepoint.x <- NULL
combo_adjusted$patient_id.x <- NULL
combo_adjusted$timepoint.y <- NULL
combo_adjusted$patient_id.y <- NULL

# Separate the data by timepoints
T1 <- combo_adjusted %>% filter(timepoint == "pre_treatment")
T2 <- combo_adjusted %>% filter(timepoint == "surgical_resection")

combined_rejoined_cpm_tr4_init <- rbind(T1, T2)

combined_rejoined_cpm_tr4_init <- as.data.frame(combined_rejoined_cpm_tr4_init)

#rownames(combined_rejoined_cpm_tr4_init) <- paste(combined_rejoined_cpm_tr4_init$patient, combined_rejoined_cpm_tr4_init$timepoint, sep = "_")
rownames(combined_rejoined_cpm_tr4_init) <- combined_rejoined_cpm_tr4_init$Mixture

# Figure 5B ----
combined_rejoined_cpm_tr4_init$timepoint[combined_rejoined_cpm_tr4_init$timepoint == "pre_treatment"] <- "T1"
combined_rejoined_cpm_tr4_init$timepoint[combined_rejoined_cpm_tr4_init$timepoint == "surgical_resection"] <- "T2"

T1_CLR <- combined_rejoined_cpm_tr4_init[combined_rejoined_cpm_tr4_init$timepoint == "T1", ]

#rownames(T1_CLR) <- paste(T1_CLR$patient, T1_CLR$timepoint, sep = "_")

T1_CLR$Mixture <- NULL

T1_CLR <- T1_CLR[, -c(42, 43)]
T1_CLR <- compositions::clr(as.matrix(T1_CLR)) 


df <- as.data.frame(T1_CLR) %>%
  rownames_to_column(var = "Sample_TP")

df <- df %>%
  separate(
    col    = Sample_TP,
    into   = c("Sample", "Time_Point"),
    sep    = "_",        
    remove = TRUE      
  )

cibersortx_long_T1 <- df %>%
  pivot_longer(
    cols       = -c(Sample, Time_Point),  
    names_to   = "Cell_Type",
    values_to  = "Value"
)

src2id <- setNames(names(src_lookup), unname(src_lookup))
cibersortx_long_T1$Sample <- src2id[cibersortx_long_T1$Sample]

cibersortx_long_T1 <- cibersortx_long_T1 %>% 
  left_join(SIC_and_SE_assignments,
            by = c("Sample" = "Patient"))

cell_types <- unique(cibersortx_long_T1$Cell_Type)

# Initialize a data frame to store Cox PH results
cox_results <- data.frame(Cell_Type = character(),
                          Arm = character(),
                          HR = numeric(),
                          Lower_CI = numeric(),
                          Upper_CI = numeric(),
                          coef = numeric(),
                          se = numeric(),
                          p_value = numeric(),
                          adjusted_p_value = numeric(),
                          stringsAsFactors = FALSE)

# Loop through each cell type
for (cell in cell_types) {
  subset_data <- subset(cibersortx_long_T1, Cell_Type == cell)
  
  if (length(unique(subset_data$Value)) > 1) {
    cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Value, data = subset_data)
    
    summary_cox <- summary(cox_model)
    HR <- summary_cox$coefficients[1, "exp(coef)"]
    Lower_CI <- summary_cox$conf.int[1, "lower .95"]
    Upper_CI <- summary_cox$conf.int[1, "upper .95"]
    coef <- summary_cox$coefficients[1, "coef"]
    se <- summary_cox$coefficients[1, "se(coef)"]
    p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
    
    # Store results for combined analysis
    cox_results <- rbind(cox_results, data.frame(Cell_Type = cell,
                                                 Arm = "All",  
                                                 HR = HR,
                                                 Lower_CI = Lower_CI,
                                                 Upper_CI = Upper_CI,
                                                 coef = coef,
                                                 se = se,
                                                 p_value = p_value,
                                                 stringsAsFactors = FALSE))
  } else {
    message(paste("Skipping combined analysis for cell type", cell, "- no variability in 'Value'"))
  }
  
  # Separate analysis for Control and Experimental arms
  for (arm in c("Control", "Experimental")) {
    arm_data <- subset(subset_data, `Treatment Arm` == arm)
    
    if (length(unique(arm_data$Value)) > 1) {
      cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Value, data = arm_data)
      
      summary_cox <- summary(cox_model)
      HR <- summary_cox$coefficients[1, "exp(coef)"]
      Lower_CI <- summary_cox$conf.int[1, "lower .95"]
      Upper_CI <- summary_cox$conf.int[1, "upper .95"]
      coef <- summary_cox$coefficients[1, "coef"]
      se <- summary_cox$coefficients[1, "se(coef)"]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      cox_results <- rbind(cox_results, data.frame(Cell_Type = cell,
                                                   Arm = arm,
                                                   HR = HR,
                                                   Lower_CI = Lower_CI,
                                                   Upper_CI = Upper_CI,
                                                   coef = coef,
                                                   se = se,
                                                   p_value = p_value,
                                                   stringsAsFactors = FALSE))
    } else {
      message(paste("Skipping cell type", cell, "in arm", arm, "- no variability in 'Value'"))
    }
  }
}

cells_sel <- c("NK cells CD56 low", "CD36+ CAFs", "M2 TREM2+ SPP1+",
               "Monocytes EREG+",
               "CD8+ Activated",
               "CD8+ GZMK+ Tm",
               "STS Mesenchymal-like ECM Remodeling",
               "Immunoregulatroy PMNs"
              
)

cox_results_sig <- cox_results[cox_results$p_value < 0.05 | cox_results$Cell_Type %in% cells_sel, ]
dt <- cox_results
dt$est <- log(dt$HR)
dt$low <- log(dt$Lower_CI)
dt$hi <- log(dt$Upper_CI)
dt$est_formatted <- sprintf("%.2f", dt$est)
dt$ci_formatted <- sprintf("(%.2f - %.2f)", dt$low, dt$hi)
dt$Effect_CI <- paste0(dt$est_formatted, " ", dt$ci_formatted)
dt$Cell_Type_copy <- dt$Cell_Type
dt$Cell_Type <- paste(dt$Cell_Type, dt$Arm, sep = "_")
dt <- dt[dt$Cell_Type_copy %in% cox_results_sig$Cell_Type, ]
dt <- dt[, c("Cell_Type", "est", "low", "hi", "se", "p_value", "Effect_CI")]
dt <- dt %>% 
  arrange(Cell_Type)
plot_data <- dt[, c("Cell_Type", "Effect_CI", "p_value")]
colnames(plot_data) <- c("Cell Type", "Hazard Ratio (95% CI)", "P-value")
plot_data[, 4] <- ""
colnames(plot_data)[colnames(plot_data) == "V4"] <- "                                                     "
plot_data$`Cell Type` <- ifelse(grepl("Control|Experimental", plot_data$`Cell Type`), 
                                paste0("   ", plot_data$`Cell Type`), 
                                plot_data$`Cell Type`)

plot_data$`P-value` <- as.numeric(plot_data$`P-value`)

plot_data <- plot_data %>%
  mutate(
    cox_coef = as.numeric(str_extract(`Hazard Ratio (95% CI)`, "-?\\d+\\.?\\d*")),
    signed_log10p = sign(cox_coef) * (-log10(`P-value`))
  )

# define broad categories of cell types 
plot_data <- plot_data %>%
  mutate(
    Category = case_when(
      str_detect(`Cell Type`, "CD4\\+") ~ "Lymphoid Cells",
      str_detect(`Cell Type`, "NK cells CD56 low") ~ "Lymphoid Cells",
      str_detect(`Cell Type`, "CD8\\+") ~ "Lymphoid Cells",
      str_detect(`Cell Type`, "B Cells & Plasma Cells") ~ "Lymphoid Cells",
      str_detect(`Cell Type`, "M2|Monocytes|Macro") ~ "Myeloid Cells",
      str_detect(`Cell Type`, "Immunoregulatroy PMNs") ~ "Myeloid Cells",
      str_detect(`Cell Type`, "CAF") ~ "CAFs",
      str_detect(`Cell Type`, "STS") ~ "Sarcoma Cells (STS)",
      TRUE ~ "Other"
    )
  )

plot_data <- plot_data %>%
  mutate(
    Arm = case_when(
      grepl("_Control", `Cell Type`)      ~ "Control",
      grepl("_Experimental", `Cell Type`) ~ "Experimental",
      grepl("_All", `Cell Type`)          ~ "Both",
      TRUE                                ~ "Unknown"
    ),
    `Cell Type` = str_replace_all(`Cell Type`, "_Control|_Experimental|_All", "")
  )

plot_data$`Cell Type` <- gsub("   ", "", plot_data$`Cell Type`)

barplot_double_faceted_T1 <- ggplot(
  data = plot_data, 
  aes(x = `Cell Type`, y = signed_log10p, fill = signed_log10p)
) +
  geom_col(width = 0.7, color = "black") +
  
  scale_fill_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick2",
    midpoint = 0,
    name = "Signed -log10(p)\n(blue = better DFS,\n red = worse DFS)"
  ) +
  
  facet_grid(Arm ~ Category, scales = "free_x", space = "free_x") +
  
  geom_hline(yintercept = c(-1.3, 0, 1.3), 
             linetype = "dashed", 
             color = "black") +
  theme_pubr(border = TRUE, base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL, 
    y = "–log10(p) Cox PH",
    title = "Cox DFS T1"
  )

ggsave(barplot_double_faceted_T1, width = 22, height = 20, units = "cm", 
       filename = '.../T1_cibersottx_DFS_by_macrogroup_triple_faceted_CLR.pdf')

# Figure 5C ----
T1_CLR <- combined_rejoined_cpm_tr4_init[combined_rejoined_cpm_tr4_init$timepoint == "T1", ]

#rownames(T1_CLR) <- paste(T1_CLR$patient, T1_CLR$timepoint, sep = "_")

T1_CLR$Mixture <- NULL

T1_CLR <- T1_CLR[, -c(42, 43)]
T1_CLR <- compositions::clr(as.matrix(T1_CLR)) 

tumor_long <- as.data.frame(T1_CLR) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "CLR")

tumor_long$Sample <- gsub("_T1", "", tumor_long$Sample)

src2id <- setNames(names(src_lookup), unname(src_lookup))

tumor_long$Sample <- src2id[tumor_long$Sample]

data_for_analysis <- tumor_long %>%
  inner_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Immune Class")],
             by = c("Sample" = "Patient"))

colnames(data_for_analysis)[colnames(data_for_analysis) == "Sarcoma Immune Class"] <- "Sarcoma Immune Class"

# Kruskal-Wallis test 
kw_results <- data_for_analysis %>%
  group_by(CellType) %>%
  kruskal_test(CLR ~ `Sarcoma Immune Class`) %>%
  adjust_pvalue(method = "BH") %>%
  rename(p_value = p, adj_p_value = p.adj)

significant_cells_sics <- kw_results %>%
  filter(adj_p_value < 0.05) %>%
  pull(CellType)

# Post-hoc pairwise Dunn's test 
names(data_for_analysis)[names(data_for_analysis) == "Sarcoma Immune Class"] <- "Sarcoma_Immune_Class"
pairwise_np <- data_for_analysis %>%
  group_by(CellType) %>%
  dunn_test(CLR ~ Sarcoma_Immune_Class, 
            p.adjust.method = "BH") %>%
  mutate(signif = ifelse(p.adj < 0.05, "yes", "no"))

sig_np <- pairwise_np %>% filter(signif == "yes")

selected_cells <- union(significant_cells_sics, unique(sig_np$CellType))

plot_cells <- setdiff(selected_cells, c("Vascular CAFs", "CD4+ Stressed", 
                                        "CD4+ Th17", 
                                        "Capillary Venous Endothelial Cells", 
                                        "STS Immunomodulatory"))

my_comparisons <- list( c("E", "D"), c("E", "C"), c("E", "B"), c("E", "A"))

Sarcoma_Immune_Class_palette = c("A" = "#1f78b4", "B" = "#a6cee3", "C" = "#33a02c", "D" = "#fdbf6f", "E" = "#e31a1c")

swarm_violin <- ggplot(data_for_analysis %>% filter(CellType %in% plot_cells), 
                       aes(x = Sarcoma_Immune_Class, y = CLR, fill = Sarcoma_Immune_Class)) +
  geom_violin(alpha = 0.2, trim = T, scale = "width") +
  geom_quasirandom(shape = 21, size = 1.5, stroke = 0.5, alpha = 1, method = "smiley") +
  scale_fill_manual(values = Sarcoma_Immune_Class_palette) + 
  labs(y = "CLR Normalized Cell %", x = "Sarcoma Immune Class") + 
  stat_compare_means(comparisons = my_comparisons, 
                     label = "p.signif", 
                     method = "t.test", 
                     hide.ns = F)+ 
  facet_wrap(~ CellType, ncol = 4, scale = "free_y") + 
  stat_compare_means(method = "anova") +
  theme_pubr(legend = "top", border = T) +
  theme(
    axis.title = element_text(size = 18),    
    axis.text = element_text(size = 16),     
    strip.text = element_text(size = 14),      
    legend.title = element_text(size = 12),    
    legend.text = element_text(size = 10)      
  )

ggsave(swarm_violin, 
       units = "cm", width = 35, height = 25, 
       filename = '.../swarm_violin_selected_CLR_SIC.pdf')

# Figure 5D -----
tumor_long <- as.data.frame(T1_CLR) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "CLR")

tumor_long$Sample <- gsub("_T1", "", tumor_long$Sample)

src2id <- setNames(names(src_lookup), unname(src_lookup))
tumor_long$Sample <- src2id[tumor_long$Sample]

data_for_analysis <- tumor_long %>%
  inner_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Ecotype Assignment")],
             by = c("Sample" = "Patient"))

data_for_analysis <- data_for_analysis[!(is.na(data_for_analysis$`Sarcoma Ecotype Assignment`)), ]

# Kruskal-Wallis test 
kw_results <- data_for_analysis %>%
  group_by(CellType) %>%
  kruskal_test(CLR ~ `Sarcoma Ecotype Assignment`) %>%
  adjust_pvalue(method = "BH") %>%
  rename(p_value = p, adj_p_value = p.adj)

# keep cell types with significant global differences
significant_cells_ecotypes <- kw_results %>%
  filter(adj_p_value < 0.05) %>%
  pull(CellType)

# Post-hoc pairwise Dunn's test 
data_for_analysis$ecotype_assignment <- data_for_analysis$`Sarcoma Ecotype Assignment`
pairwise_np <- data_for_analysis %>%
  group_by(CellType) %>%
  dunn_test(CLR ~ ecotype_assignment, p.adjust.method = "BH") %>%
  mutate(signif = ifelse(p.adj < 0.05, "yes", "no"))

sig_np <- pairwise_np %>% filter(signif == "yes")
selected_cells <- union(significant_cells_ecotypes, unique(sig_np$CellType))
plot_cells <- setdiff(selected_cells, c("CD4+ Th17", "STS Immunomodulatory", "STS Immunomodulatory"))
plot_cells <- union(plot_cells, c("M2 TREM2+ SPP1+", "CD4+ Tn"))

# Specify the comparisons you want
my_comparisons <- list( c("SE2", "SE1"), c("SE3", "SE2"), c("SE3", "SE1"))

ecotyper_palette = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")

swarm_violin_eco <- ggplot(
  data_for_analysis %>% 
    filter(CellType %in% plot_cells) %>% 
    filter(`Sarcoma Ecotype Assignment` != "NA"),
  aes(x = ecotype_assignment, y = CLR, fill = ecotype_assignment)
) +
  geom_violin(alpha = 0.2, trim = TRUE, scale = "width") +
  ggbeeswarm::geom_quasirandom(shape = 21, size = 1.5, stroke = 0.5, alpha = 1, method = "smiley") +
  scale_fill_manual(values = ecotyper_palette) +
  labs(y = "CLR Normalized Cell %", x = "Sarcoma Ecotype") +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "wilcox.test",
    hide.ns = FALSE
  ) +
  stat_compare_means(method = "kruskal.test") +
  facet_wrap(~ CellType, ncol = 3, scale = "free_y") +
  theme_pubr(legend = "top", border = TRUE) +
  theme(
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave(swarm_violin_eco, 
       units = "cm", width = 25, height = 25, 
       filename = ".../swarm_violin_selected_CLR_SEs.pdf")

# Figure 6B ----
T2_CLR <- combined_rejoined_cpm_tr4_init[combined_rejoined_cpm_tr4_init$timepoint == "T2", ]

#rownames(T2_CLR) <- paste(T2_CLR$patient, T2_CLR$timepoint, sep = "_")

T2_CLR$Mixture <- NULL

T2_CLR <- T2_CLR[, -c(42, 43)]
T2_CLR <- compositions::clr(as.matrix(T2_CLR)) 

df <- as.data.frame(T2_CLR) %>%
  rownames_to_column(var = "Sample_TP")

df <- df %>%
  separate(
    col    = Sample_TP,
    into   = c("Sample", "Time_Point"),
    sep    = "_",        
    remove = TRUE        
  )

cibersortx_long_T2 <- df %>%
  pivot_longer(
    cols       = -c(Sample, Time_Point), 
    names_to   = "Cell_Type",
    values_to  = "Value"
  )

src2id <- setNames(names(src_lookup), unname(src_lookup))
cibersortx_long_T2$Sample <- src2id[cibersortx_long_T2$Sample]

cibersortx_long_T2 <- cibersortx_long_T2 %>% 
  left_join(SIC_and_SE_assignments,
            by = c("Sample" = "Patient"))

surgery_timing <- fread('/Users/stefanotesta/Desktop/Moding Lab/SARC_32/Figures/last_draft/updated_time_cytof/SARC032_surgery_timing.csv')

cibersortx_long_T2 <- cibersortx_long_T2 %>% 
  left_join(surgery_timing, by = c("Sample" = "Subject"))

# Adjust DFS time to subtract time from enrollment to surgery 
cibersortx_long_T2$new_dfs_time <- cibersortx_long_T2$`DFS Time (Days)` - cibersortx_long_T2$time_from_enrollment_to_surgery_days

cell_types <- unique(cibersortx_long_T2$Cell_Type)

# Initialize a data frame to store Cox PH results
cox_results <- data.frame(Cell_Type = character(),
                          Arm = character(),
                          HR = numeric(),
                          Lower_CI = numeric(),
                          Upper_CI = numeric(),
                          coef = numeric(),
                          se = numeric(),
                          p_value = numeric(),
                          adjusted_p_value = numeric(),
                          stringsAsFactors = FALSE)

# Loop through each cell type
for (cell in cell_types) {
  # Subset data by cell type
  subset_data <- subset(cibersortx_long_T2, Cell_Type == cell)
  
  # Combined analysis (all samples together)
  if (length(unique(subset_data$Value)) > 1) {
    cox_model <- coxph(Surv(new_dfs_time, `DFS Event` == "Event occurred") ~ Value, data = subset_data)
    
    summary_cox <- summary(cox_model)
    HR <- summary_cox$coefficients[1, "exp(coef)"]
    Lower_CI <- summary_cox$conf.int[1, "lower .95"]
    Upper_CI <- summary_cox$conf.int[1, "upper .95"]
    coef <- summary_cox$coefficients[1, "coef"]
    se <- summary_cox$coefficients[1, "se(coef)"]
    p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
    
    # Store results for combined analysis
    cox_results <- rbind(cox_results, data.frame(Cell_Type = cell,
                                                 Arm = "All",  # Label for combined analysis
                                                 HR = HR,
                                                 Lower_CI = Lower_CI,
                                                 Upper_CI = Upper_CI,
                                                 coef = coef,
                                                 se = se,
                                                 p_value = p_value,
                                                 stringsAsFactors = FALSE))
  } else {
    message(paste("Skipping combined analysis for cell type", cell, "- no variability in 'Value'"))
  }
  
  # Separate analysis for Control and Experimental arms
  for (arm in c("Control", "Experimental")) {
    arm_data <- subset(subset_data, `Treatment Arm` == arm)
    
    if (length(unique(arm_data$Value)) > 1) {
      cox_model <- coxph(Surv(new_dfs_time, `DFS Event` == "Event occurred") ~ Value, data = arm_data)
      
      summary_cox <- summary(cox_model)
      HR <- summary_cox$coefficients[1, "exp(coef)"]
      Lower_CI <- summary_cox$conf.int[1, "lower .95"]
      Upper_CI <- summary_cox$conf.int[1, "upper .95"]
      coef <- summary_cox$coefficients[1, "coef"]
      se <- summary_cox$coefficients[1, "se(coef)"]
      p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
      
      # Store results for each arm
      cox_results <- rbind(cox_results, data.frame(Cell_Type = cell,
                                                   Arm = arm,
                                                   HR = HR,
                                                   Lower_CI = Lower_CI,
                                                   Upper_CI = Upper_CI,
                                                   coef = coef,
                                                   se = se,
                                                   p_value = p_value,
                                                   stringsAsFactors = FALSE))
    } else {
      message(paste("Skipping cell type", cell, "in arm", arm, "- no variability in 'Value'"))
    }
  }
}

cells_sel <- c("STS Epithelial-like",
               "Mastocytes",
               "STS Myofibroblastic",
               "Matrix CAFs",
               "CD8+ Stressed",
               "STS CTA high",
               "Proinflammatory PMNs",
               "Vascular CAFs")

cox_results_sig <- cox_results[cox_results$p_value < 0.05 | cox_results$Cell_Type %in% cells_sel, ]
# Create the data frame 'dt' for the forest plot
#dt <- cox_results_sig
dt <- cox_results
# Rename columns to match 'forest()' function requirements
dt$est <- log(dt$HR)
dt$low <- log(dt$Lower_CI)
dt$hi <- log(dt$Upper_CI)
# Optional: Round the estimates and confidence intervals for display
dt$est_formatted <- sprintf("%.2f", dt$est)
dt$ci_formatted <- sprintf("(%.2f - %.2f)", dt$low, dt$hi)
# Create a column combining the estimate and CI for display
dt$Effect_CI <- paste0(dt$est_formatted, " ", dt$ci_formatted)
dt$Cell_Type_copy <- dt$Cell_Type
dt$Cell_Type <- paste(dt$Cell_Type, dt$Arm, sep = "_")
dt <- dt[dt$Cell_Type_copy %in% cox_results_sig$Cell_Type, ]
#Rearrange columns as needed
dt <- dt[, c("Cell_Type", "est", "low", "hi", "se", "p_value", "Effect_CI")]
dt <- dt %>% 
  arrange(Cell_Type)
plot_data <- dt[, c("Cell_Type", "Effect_CI", "p_value")]
# Rename columns for better display
colnames(plot_data) <- c("Cell Type", "Hazard Ratio (95% CI)", "P-value")
# Convert p-values to formatted strings
#plot_data$`P-value` <- sprintf("%.3f", plot_data$`P-value`)
plot_data[, 4] <- ""
colnames(plot_data)[colnames(plot_data) == "V4"] <- "                                                     "
plot_data$`Cell Type` <- ifelse(grepl("Control|Experimental", plot_data$`Cell Type`), 
                                paste0("   ", plot_data$`Cell Type`), 
                                plot_data$`Cell Type`)

plot_data$`P-value` <- as.numeric(plot_data$`P-value`)

plot_data <- plot_data %>%
  mutate(
    cox_coef = as.numeric(str_extract(`Hazard Ratio (95% CI)`, "-?\\d+\\.?\\d*")),
    signed_log10p = sign(cox_coef) * (-log10(`P-value`))
  )

# define broad categories of cell types 
plot_data <- plot_data %>%
  mutate(
    Category = case_when(
      str_detect(`Cell Type`, "CD4\\+") ~ "Lymphoid Cells",
      str_detect(`Cell Type`, "CD8\\+") ~ "Lymphoid Cells",
      #str_detect(`Cell Type`, "DCs") ~ "Myeloid Cells",
      str_detect(`Cell Type`, "Mastocytes") ~ "Myeloid Cells",
      str_detect(`Cell Type`, "Proinflammatory PMNs") ~ "Myeloid Cells",
      str_detect(`Cell Type`, "CAF") ~ "CAFs",
      str_detect(`Cell Type`, "Vascular CAFs") ~ "CAFs",
      str_detect(`Cell Type`, "STS") ~ "Sarcoma Cells (STS)",
      TRUE ~ "Other"
    )
  )

plot_data <- plot_data %>%
  mutate(
    Arm = case_when(
      grepl("_Control", `Cell Type`)      ~ "Control",
      grepl("_Experimental", `Cell Type`) ~ "Experimental",
      grepl("_All", `Cell Type`)          ~ "Both",
      TRUE                                ~ "Unknown"
    ),
    # Then remove the suffix from Cell Type so that it doesn’t clutter the facet
    `Cell Type` = str_replace_all(`Cell Type`, "_Control|_Experimental|_All", "")
  )

plot_data$`Cell Type` <- gsub("   ", "", plot_data$`Cell Type`)

barplot_double_faceted_T2 <- ggplot(
  data = plot_data, 
  aes(x = `Cell Type`, y = signed_log10p, fill = signed_log10p)
) +
  geom_col(width = 0.7, color = "black") +
  
  scale_fill_gradient2(
    low = "dodgerblue3", mid = "white", high = "firebrick2",
    midpoint = 0,
    name = "Signed -log10(p)\n(blue = better DFS,\n red = worse DFS)"
  ) +
  
  # Suppose you want rows = Arm, columns = Category
  # so that you see, for each Category, subpanels for
  # Control/Experimental/Both in separate rows (or vice versa).
  facet_grid(Arm ~ Category, scales = "free_x", space = "free_x") +
  
  geom_hline(yintercept = c(-1.3, 0, 1.3), 
             linetype = "dashed", 
             color = "black") +
  theme_pubr(border = TRUE, base_size = 14) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = NULL,  # we often omit x label if it's just "Cell Type"
    y = "–log10(p) Cox PH",
    title = "Cox DFS T2 - All in One Plot"
  ) + 
  ylim(c(-2,2))


ggsave(barplot_double_faceted_T2, width = 22, height = 16, units = "cm", 
       filename = '.../T2_cibersottx_DFS_by_macrogroup_triple_faceted_CLR.pdf')

# Figure 6E ----
combined_rejoined_cpm_tr4_CLR <- compositions::clr(as.matrix(combined_rejoined_cpm_tr4_init[, -c(1, 43, 44)])) 

combined_rejoined_cpm_tr4_CLR <-  combined_rejoined_cpm_tr4_CLR %>%
  as.data.frame(combined_rejoined_cpm_tr4_CLR) %>% 
  rownames_to_column(var = "Mixture") %>% 
  mutate(patient = str_extract(Mixture, "^[^_]+"),
         timepoint = str_extract(Mixture, "(?<=_).+")) %>% 
  select(-Mixture)

cibersortx_long_cpm_tr4 <- melt(combined_rejoined_cpm_tr4_CLR, 
                                id.vars = c("patient", "timepoint"),
                                variable.name = "Cell_Type", 
                                value.name = "Value")
colnames(cibersortx_long_cpm_tr4)[colnames(cibersortx_long_cpm_tr4) == "patient"] <- "Sample"
colnames(cibersortx_long_cpm_tr4)[colnames(cibersortx_long_cpm_tr4) == "timepoint"] <- "Time_Point"
cibersortx_long_cpm_tr4$Sample[cibersortx_long_cpm_tr4$Sample == "SRC77_b"] <- "SRC77"

src2id <- setNames(names(src_lookup), unname(src_lookup))
cibersortx_long_cpm_tr4$Sample <- src2id[cibersortx_long_cpm_tr4$Sample]

cibersortx_long_cpm_tr4 <- cibersortx_long_cpm_tr4 %>% 
  left_join(SIC_and_SE_assignments, by = c("Sample" = "Patient"))

cibersortx_long_cpm_tr4 <- cibersortx_long_cpm_tr4[!(cibersortx_long_cpm_tr4$Sample %in% "SRC102b"), ]

paired_df <- cibersortx_long_cpm_tr4 %>%
  filter(Time_Point %in% c("T1","T2")) %>%
  pivot_wider(names_from = Time_Point, values_from = Value) %>%
  drop_na(T1, T2)

# Compute delta and fold‐change
paired_df <- paired_df %>%
  mutate(
    Delta       = T2 - T1,
    mean_delta  = Delta,            
    fold_change = exp(mean_delta)
  )

# Paired test *within* each arm
within_arm <- paired_df %>%
  group_by(Cell_Type, `Treatment Arm`) %>%
  summarize(
    p_within   = wilcox.test(T2, T1, paired=TRUE)$p.value, 
    delta      = mean(Delta),
    fc         = exp(delta),
    .groups    = "drop"
  ) %>%
  mutate(adj_p_within = p.adjust(p_within, "fdr"))

# Unpaired test *between* arms on Δ
between_arm <- paired_df %>%
  group_by(Cell_Type) %>%
  summarize(
    p_between  = wilcox.test(Delta ~ `Treatment Arm`)$p.value, 
    mean_delta_control = mean(Delta[`Treatment Arm`=="Control"]),
    mean_delta_exp     = mean(Delta[`Treatment Arm`=="Experimental"]),
    .groups    = "drop"
  ) %>%
  mutate(adj_p_between = p.adjust(p_between, "fdr"))

# Classify
results <- left_join(within_arm, between_arm, by="Cell_Type")

# which cell states change in each group 
cell_increase_exp <- within_arm %>% 
  filter(`Treatment Arm` == "Experimental") %>% 
  filter(adj_p_within < 0.05 & delta > 0) %>% 
  pull(Cell_Type)

cell_increase_control <- within_arm %>% 
  filter(`Treatment Arm` == "Control") %>% 
  filter(adj_p_within < 0.05 & delta > 0) %>% 
  pull(Cell_Type)

cell_decrease_exp <- within_arm %>% 
  filter(`Treatment Arm` == "Experimental") %>% 
  filter(adj_p_within < 0.05 & delta < 0) %>% 
  pull(Cell_Type)

cell_decrease_control <- within_arm %>% 
  filter(`Treatment Arm` == "Control") %>% 
  filter(adj_p_within < 0.05 & delta < 0) %>% 
  pull(Cell_Type)

nonimmunecells_selected_toplot <- c("Matrix CAFs", "STS Epithelial-like", 
                                    "STS Mesenchymal-like ECM Remodeling")
non_immune_selected <- ggplot(
  cibersortx_long_cpm_tr4[cibersortx_long_cpm_tr4$Cell_Type %in% nonimmunecells_selected_toplot, ],
  aes(x = Time_Point, y = Value, fill = Time_Point)
) +
  geom_violin(alpha = 0.3, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) + 
  geom_line(aes(group = Sample), color = "black", alpha = 0.2) + 
  geom_point(aes(group = Sample), size = 0.2) + 
  labs(
    title = "",
    x = "",
    y = "CIBERSORTx Fraction"
  ) +
  theme_pubr(border = T) +
  geom_signif(comparisons = list(c("T1", "T2")), test = "wilcox.test") +
  labs(fill = "Time Point") +
  facet_grid(`Treatment Arm` ~ Cell_Type,
             scales="free_y",
             switch = "y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("T1" = "#396CAE", "T2" = "#FBBF85"))

ggsave(non_immune_selected, 
       filename = ".../T2_T1_Non_immune_selected_CLR.pdf",
       units = "cm", width = 15, height = 15)

immunecells_selected_toplot <- c("CD8+ Early Activation",
                                 "CD4+ T memory", "CD4+ Tfh", "Macro SELENOP+ SLC40A1+")

immune_selected <- ggplot(
  cibersortx_long_cpm_tr4[cibersortx_long_cpm_tr4$Cell_Type %in% immunecells_selected_toplot, ],
  aes(x = Time_Point, y = Value, fill = Time_Point)
) +
  geom_violin(alpha = 0.3, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) + # Avoid outliers overlapping violins
  geom_line(aes(group = Sample), color = "black", alpha = 0.2) + # Add connecting lines
  geom_point(aes(group = Sample), size = 0.2) + # Add points for T1 and T2
  labs(
    title = "",
    x = "",
    y = "CIBERSORTx Fraction"
  ) +
  theme_pubr(border = T) +
  geom_signif(comparisons = list(c("T1", "T2")), test = "wilcox.test") +
  labs(fill = "Time Point") +
  facet_grid(`Treatment Arm` ~ Cell_Type,
             scales="free_y",
             switch = "y") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("T1" = "#396CAE", "T2" = "#FBBF85"))

ggsave(immune_selected, 
       filename = "/.../T2_T1_immune_selected_CLR.pdf",
       units = "cm", width = 24, height = 15)

# Figure 6A ----
necrosis_evaluable_wilcox <- SIC_and_SE_assignments %>% 
  filter(!(is.na(SIC_and_SE_assignments$`Percent Necrosis`))) %>% 
  mutate(necrosis = `Percent Necrosis`) %>% 
  ggplot(
    aes(x = `Treatment Arm`, y = necrosis, fill = `Treatment Arm`)
  ) + 
  geom_boxplot(width = .50, outlier.shape = NA, alpha = .5) +
  geom_point(aes(group = `Patient`), shape = 21, size = 1.5, stroke = .5, alpha = 1) +
  scale_fill_manual(values = c("Experimental" = "#0093af", "Control" = "#9d1536")) + 
  stat_compare_means(
    comparisons = list( c("Control", "Experimental")),
    label = "p.format",
    method = "wilcox.test",
    hide.ns = FALSE
  ) +
  theme_pubr(border = TRUE)

ggsave(necrosis_evaluable_wilcox, 
       units = "cm", width = 8, height = 10, 
       filename = ".../necrosis_box_evaluable.pdf")

# Figure S6A ----
path_data_filtered <- SIC_and_SE_assignments %>%
  mutate(necrosis = `Percent Necrosis`) %>% 
  mutate(
    DFS_event = ifelse(`DFS Event` == "Event occurred", 1, 0),
    necrosis10    = necrosis / 10,                 
    necrosis_high = ifelse(necrosis >= 95, 1, 0)   
  )

df_rt  <- path_data_filtered %>% filter(`Treatment Arm` == "Control")        
df_exp <- path_data_filtered %>% filter(`Treatment Arm` == "Experimental")   

cox_one <- function(dat, time_col, event_col, predictor, covariate_label, arm_label, outcome_label){
  f <- as.formula(paste0("Surv(`", time_col, "`, ", event_col, ") ~ ", predictor))
  fit <- coxph(f, data = dat)
  
  broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    transmute(
      outcome   = outcome_label,
      arm       = arm_label,
      covariate = covariate_label,
      hr        = estimate,
      lo        = conf.low,
      hi        = conf.high,
      p_wald    = p.value,
      n         = nrow(dat),
      events    = sum(dat[[event_col]], na.rm = TRUE)
    )
}

# DFS
dfs_rt_results <- bind_rows(
  cox_one(df_rt,  "DFS Time (Days)", "DFS_event", "necrosis10",
          "Necrosis after resection (per 10% increase)", "RT Only", "DFS"),
  cox_one(df_rt,  "DFS Time (Days)", "DFS_event", "necrosis_high",
          "% Necrosis \u226595% (vs <95%)",               "RT Only", "DFS")
)

dfs_exp_results <- bind_rows(
  cox_one(df_exp, "DFS Time (Days)", "DFS_event", "necrosis10",
          "Necrosis after resection (per 10% increase)", "RT + Pembro", "DFS"),
  cox_one(df_exp, "DFS Time (Days)", "DFS_event", "necrosis_high",
          "% Necrosis \u226595% (vs <95%)",               "RT + Pembro", "DFS")
)

all_results <- bind_rows(dfs_rt_results, dfs_exp_results)

all_results %>% arrange(outcome, covariate, arm)

arm_cols <- c("RT Only" = "#9d1536", "RT + Pembro" = "#0093af")
square_forest <- function(dat, xlab = "Multivariable Cox hazard ratio") {
  dat <- dat %>%
    mutate(
      arm = factor(arm, levels = c("RT Only", "RT + Pembro")),
      row = paste0(ifelse(arm == "RT Only","RT","RT+P"), " — ", covariate)
    )

  dat <- dat %>%
    mutate(row = factor(row, levels = rev(unique(row))))

  x_min <- max(0.05, min(dat$lo, na.rm=TRUE) * 0.8)
  x_max <- min(100, max(dat$hi, na.rm=TRUE) * 1.25)
  
  ggplot(dat, aes(x = hr, y = row, fill = arm)) +
    geom_vline(xintercept = 1, linetype = "22", color = "grey50", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5) +
    geom_point(shape = 22, size = 5, stroke = 0.5, color = "black") +
    scale_fill_manual(values = arm_cols, guide = "none") +
    scale_x_log10(
      limits = c(x_min, x_max),
      breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50),
      minor_breaks = NULL
    ) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(face = "bold"),
      axis.text.y   = element_text(size = 12)
    )
}

dfs <- all_results %>% filter(outcome == "DFS")
p_dfs_sq <- square_forest(dfs, xlab = "Univariate Cox DFS hazard ratio") +
  ggpubr::theme_pubr()
p_dfs_sq <- p_dfs_sq + theme_pubr(border = T)



# Figure 6F ----
ecotype_abundance <- read_excel('.../data/Supplementary_Tables.xlsx', 
                                sheet = 12)

colnames(ecotype_abundance) <- ecotype_abundance[2, ]
ecotype_abundance <- ecotype_abundance[-c(1,2), ]
ecotype_abundance <- as.data.frame(ecotype_abundance)
ecotype_abundance$`Sarcoma Ecotype Abundance` <- as.numeric(ecotype_abundance$`Sarcoma Ecotype Abundance`)

df_long <- ecotype_abundance %>%
  mutate(
    patient    = Patient,  # or: patient = unname(src_lookup[Patient])
    Time_Point = `Time Point`,
    Ecotype    = `Sarcoma Ecotype`,
    abundance  = `Sarcoma Ecotype Abundance`
  ) %>%
  select(patient, Time_Point, Ecotype, abundance, `Treatment Arm`)

ecotyper_palette = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")

# plot T1 vs T2 ecotype abundance for patients with matched time points 
violin_change <- df_long %>% 
  group_by(patient) %>%
  filter(n_distinct(Time_Point) == 2) %>%
  ungroup() %>% 
  ggplot(
    aes(x = Time_Point, y = abundance, fill = Ecotype)
  ) +
  geom_violin(alpha = 0.3, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA) + 
  geom_line(aes(group = patient), color = "black", alpha = 0.2) + 
  geom_point(aes(group = patient), size = 0.2) + 
  labs(
    title = "",
    x = "",
    y = "Sarcoma Ecotype Abundance"
  ) +
  theme_pubr(border = T, legend = "right") +
  labs(fill = "SEs") +
  facet_grid(`Treatment Arm` ~ Ecotype,
             scales="free_y",
             switch = "y") + 
  geom_signif(comparisons = list(c("Pre-treatment", "Post-treatment")), test = "wilcox.test") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 20), 
        axis.text =  element_text(size = 18), 
        axis.title.y = element_text(size = 20)) +
  scale_fill_manual(values = ecotyper_palette)

ggsave(violin_change, units = "cm", 
       width = 22.5, height = 17.5, 
       filename = '.../violin_abundance_bySE_change.pdf')

#run wilcoxon signed-rank
df_wide_test <- df_long %>% 
  group_by(patient) %>%
  filter(n_distinct(Time_Point) == 2) %>%
  ungroup() %>%  
  select(Ecotype, patient, Time_Point, abundance, `Treatment Arm`) %>%
  pivot_wider(
    names_from  = Time_Point,
    values_from = abundance
  ) %>%
  filter(!is.na(`Pre-treatment`) & !is.na(`Post-treatment`))  

# Run wilcoxon signed-rank per ecotype
df_wide_test_results <- df_wide_test %>%
  group_by(Ecotype, `Treatment Arm`) %>%
  summarise(
    t_statistic = wilcox.test(`Post-treatment`, `Pre-treatment`, paired = TRUE)$statistic,
    p_value     = wilcox.test(`Post-treatment`, `Pre-treatment`, paired = TRUE)$p.value,
    .groups     = "drop"
  )

df_wide_test_results$fdr <- p.adjust(df_wide_test_results$p_value, method = "BH")


# Figure Supplemental 2D ----
tcell_inflamed_data <- read_excel('.../data//Supplementary_Tables.xlsx', 
                                  sheet = 2)
colnames(tcell_inflamed_data) <- tcell_inflamed_data[2,]

tcell_inflamed_data <- tcell_inflamed_data[-c(1,2), ]
tcell_inflamed_data <- as.data.frame(tcell_inflamed_data)
tcell_inflamed_data$`T cell-inflamed GEP` <- as.numeric(tcell_inflamed_data$`T cell-inflamed GEP`)

tcell_inflamed_data$Patient <- as.character(tcell_inflamed_data$Patient)

tcell_inflamed_data <- tcell_inflamed_data %>% 
  left_join(SIC_and_SE_assignments, by = c("Patient", 
                                           "Treatment Arm"))

cox_results <- data.frame(Signature = character(),
                          Arm = character(),
                          HR = numeric(),
                          Lower_CI = numeric(),
                          Upper_CI = numeric(),
                          coef = numeric(),
                          se = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

sigs <- "T cell-inflamed GEP"

# Loop through each cell type
for (signature in sigs) {
  # Subset data by cell type
  tcell_inflamed_data$Time_Point <- tcell_inflamed_data$`Time Point`
  subset_data <- tcell_inflamed_data[tcell_inflamed_data$Time_Point == "Pre-treatment", c(signature, 
                                                                               "Treatment Arm", 
                                                                               "DFS Event", 
                                                                               "DFS Time (Days)")]
  str(subset_data)
  print(paste("Signature:", signature))
  print(length(subset_data[[signature]]))
  print(dim(subset_data))
  
  # Combined analysis (all samples together)
  cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ subset_data[[signature]], 
                     data = subset_data)
  
  summary_cox <- summary(cox_model)
  HR <- summary_cox$coefficients[1, "exp(coef)"]
  Lower_CI <- summary_cox$conf.int[1, "lower .95"]
  Upper_CI <- summary_cox$conf.int[1, "upper .95"]
  coef <- summary_cox$coefficients[1, "coef"]
  se <- summary_cox$coefficients[1, "se(coef)"]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
  # Store results for combined analysis
  cox_results <- rbind(cox_results, data.frame(Signature = signature,
                                               Arm = "All",  # Label for combined analysis
                                               HR = HR,
                                               Lower_CI = Lower_CI,
                                               Upper_CI = Upper_CI,
                                               coef = coef,
                                               se = se,
                                               p_value = p_value,
                                               stringsAsFactors = FALSE))
  
  # Separate analysis for Control and Experimental arms
  for (arm in c("Control", "Experimental")) {
    arm_data <- subset(subset_data, `Treatment Arm` == arm)
    str(arm_data)
    
    cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ arm_data[[signature]],
                       data = arm_data)
    
    summary_cox <- summary(cox_model)
    HR <- summary_cox$coefficients[1, "exp(coef)"]
    Lower_CI <- summary_cox$conf.int[1, "lower .95"]
    Upper_CI <- summary_cox$conf.int[1, "upper .95"]
    coef <- summary_cox$coefficients[1, "coef"]
    se <- summary_cox$coefficients[1, "se(coef)"]
    p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
    
    # Store results for each arm
    cox_results <- rbind(cox_results, data.frame(Signature = signature,
                                                 Arm = arm,
                                                 HR = HR,
                                                 Lower_CI = Lower_CI,
                                                 Upper_CI = Upper_CI,
                                                 coef = coef,
                                                 se = se,
                                                 p_value = p_value,
                                                 stringsAsFactors = FALSE))
  }
}

dt <- cox_results
dt$est <- log(dt$HR)
dt$low <- log(dt$Lower_CI)
dt$hi <- log(dt$Upper_CI)
dt$est_formatted <- sprintf("%.2f", dt$est)
dt$ci_formatted <- sprintf("(%.2f - %.2f)", dt$low, dt$hi)
dt$Effect_CI <- paste0(dt$est_formatted, " ", dt$ci_formatted)
dt$signature_copy <- dt$Signature
dt$Signature <- paste(dt$Signature, dt$Arm, sep = "_")
dt <- dt[, c("Signature", "est", "low", "hi", "se", "p_value", "Effect_CI")]
dt <- dt %>% 
  arrange(Signature)
plot_data <- dt[, c("Signature", "Effect_CI", "p_value")]
colnames(plot_data) <- c("Signature", "Hazard Ratio (95% CI)", "P-value")
plot_data$`P-value` <- sprintf("%.3f", plot_data$`P-value`)
plot_data[, 4] <- ""
colnames(plot_data)[colnames(plot_data) == "V4"] <- "                                                     "
plot_data$Signature <- ifelse(grepl("Control|Experimental", plot_data$Signature), 
                              paste0("   ", plot_data$Signature), 
                              plot_data$Signature)

# Create the forest plot
p <- forest(
  plot_data,
  est = dt$est,
  lower = dt$low,
  upper = dt$hi,
  ci_column = 4,  # Column number where the CI will be plotted
  ref_line = 0, 
  arrow_lab = c("Lower risk of progression", "Higher risk of progression")
)
p

# Figure Supplemental 2E ----
Sarcoma_Immune_Class_palette = c("A" = "#1f78b4", "B" = "#a6cee3", "C" = "#33a02c", "D" = "#fdbf6f", "E" = "#e31a1c")
my_comparisons <- list( c("E", "D"), c("E", "C"), c("E", "B"), c("E", "A"))

ifng_expanded_genes <- ggplot(tcell_inflamed_data %>% filter(!(is.na(`Sarcoma Immune Class`)) & Time_Point == "Pre-treatment"), 
                              aes(x = `Sarcoma Immune Class`, y = `T cell-inflamed GEP`,
                                  fill = `Sarcoma Immune Class`)) +
  geom_violin(alpha = 0.5, trim = T, scale = "width") +
  geom_quasirandom(shape = 21, size = 1.5, stroke = 0.5, alpha = 1, method = "smiley") +
  scale_fill_manual(values = Sarcoma_Immune_Class_palette) + 
  labs(title = "T cell-inflamed GEP by SIC", y = "T cell-inflamed GEP", x = "SIC") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = F)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test" ) +
  theme_pubr(legend = "right")

tcell_inflamed_data$Sarcoma_Immune_Class <- tcell_inflamed_data$`Sarcoma Immune Class` 
tcell_inflamed_data$IFNy_expanded_avg <- tcell_inflamed_data$`T cell-inflamed GEP`

pairwise_dunn <- tcell_inflamed_data %>% 
  filter(!(is.na(`Sarcoma Immune Class` )) & Time_Point == "Pre-treatment") %>% 
  dunn_test(IFNy_expanded_avg ~ Sarcoma_Immune_Class, p.adjust.method = "BH") %>%
  mutate(signif = ifelse(p.adj < 0.05, "yes", "no"))

# Figure Supplemental 2F ----
ecotyper_palette = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")
my_comparisons <- list( c("SE2", "SE1"), c("SE3", "SE2"), c("SE3", "SE1"))

tcell_inflamed_data$`Sarcoma Ecotype Assignment`[tcell_inflamed_data$`Sarcoma Ecotype Assignment` == "NA"] <- NA

ifng_expanded_genes <- ggplot(tcell_inflamed_data %>% filter(!(is.na(`Sarcoma Ecotype Assignment`)) & Time_Point == "Pre-treatment"), 
                              aes(x = `Sarcoma Ecotype Assignment`, 
                                  y = `T cell-inflamed GEP`, fill = `Sarcoma Ecotype Assignment`)) +
  geom_violin(alpha = 0.5, trim = T, scale = "width") +
  geom_quasirandom(shape = 21, size = 1.5, stroke = 0.5, alpha = 1, method = "smiley") +
  scale_fill_manual(values = ecotyper_palette) + 
  labs(title = "T cell-inflamed GEP by SIC", y = "T cell-inflamed GEP", x = "Sarcoma Ecotypes") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", hide.ns = F) + 
  stat_compare_means(method = "kruskal.test" ) +
  theme_pubr(legend = "right")

tcell_inflamed_data$ecotype_assignment <- tcell_inflamed_data$`Sarcoma Ecotype Assignment`

pairwise_dunn <- tcell_inflamed_data %>% 
  filter(!(is.na(`Sarcoma Ecotype Assignment`)) & Time_Point == "Pre-treatment") %>% 
  dunn_test(IFNy_expanded_avg ~ ecotype_assignment, p.adjust.method = "BH") %>%
  mutate(signif = ifelse(p.adj < 0.05, "yes", "no"))


# Figure 6G and Figure Supplemental 6I ----
# ensure T1 before T2
tcell_inflamed_data <- tcell_inflamed_data %>%
  mutate(Time_Point = factor(Time_Point, levels = c("Pre-treatment",
                                                    "Post-treatment")))
# paired Wilcoxon by Arm
paired_wsr <- function(data, value, strata = NULL) {
  val_quo   <- rlang::enquo(value)
  strat_quo <- rlang::enquo(strata)
  
  # columns to keep for grouping the final test
  grp_cols <- c("Treatment Arm", if (!rlang::quo_is_null(strat_quo)) rlang::as_name(rlang::ensym(strata)) else NULL)
  
  data %>%
    filter(Time_Point %in% c("Pre-treatment",
                             "Post-treatment"), !is.na(!!val_quo)) %>%
    group_by(across(all_of(c(grp_cols, "Patient")))) %>%
    filter(n_distinct(Time_Point) == 2) %>%
    ungroup() %>%
    transmute(!!!syms(grp_cols), Patient, Time_Point, value = !!val_quo) %>%
    pivot_wider(names_from = Time_Point, values_from = value) %>%
    group_by(across(all_of(grp_cols))) %>%
    summarise(
      n_pairs      = n(),
      median_T1    = median(`Pre-treatment`, na.rm = TRUE),
      median_T2    = median(`Post-treatment`, na.rm = TRUE),
      median_delta = median(`Post-treatment` - `Pre-treatment`, na.rm = TRUE),
      test         = list(wilcox.test(`Post-treatment`, `Pre-treatment`, 
                                      paired = TRUE, exact = FALSE)),
      .groups = "drop"
    ) %>%
    mutate(
      W        = map_dbl(test, ~ unname(.x$statistic)),
      p_value  = map_dbl(test, ~ .x$p.value),
      p_adj_BH = p.adjust(p_value, method = "BH")
    ) %>%
    select(-test)
}

by_arm <- paired_wsr(tcell_inflamed_data, `T cell-inflamed GEP`)

by_SIC_arm <- paired_wsr(
  filter(tcell_inflamed_data, !is.na(`Sarcoma Immune Class`)),
  IFNy_expanded_avg,
  strata = `Sarcoma Immune Class`
)
by_SIC_arm

by_Eco_arm <- paired_wsr(
  filter(tcell_inflamed_data, !is.na(ecotype_assignment)),
  IFNy_expanded_avg,
  strata = ecotype_assignment
)
by_Eco_arm

# Figure Supplemental 5I ---
SIC_data <- tcell_inflamed_data %>% 
  filter(!is.na(`Sarcoma Immune Class`)) %>% 
  filter(`Sarcoma Immune Class` != "NA")

# y-positions per panel for clean annotation
SIC_ypos <- SIC_data %>%
  group_by(`Treatment Arm`, `Sarcoma Immune Class`) %>%
  summarise(y.position = max(IFNy_expanded_avg, na.rm = TRUE) * 1.05, .groups = "drop")

SIC_pvals <- by_SIC_arm %>%
  left_join(SIC_ypos, by = c("Treatment Arm","Sarcoma Immune Class")) %>%
  mutate(group1 = "Pre-treatment", group2 = "Post-treatment",
         label  = ifelse(is.na(p_value), "n.s.", scales::pvalue(p_value, accuracy = 0.001)))

SIC_change2 <-
  SIC_data %>% 
  mutate(Sample = paste(Patient, Time_Point, sep = "_")) %>% 
  ggplot(
    aes(x = Time_Point, y = IFNy_expanded_avg, fill = `Sarcoma Immune Class`)
  ) +
  geom_violin(alpha = .5, trim = FALSE, scale = "width") +
  geom_boxplot(width = .10, outlier.shape = NA, fill = "grey", alpha = .5) +
  geom_line(aes(group = Patient), colour = "black", alpha = .20) +
  geom_point(aes(group = Sample), shape = 21, size = 1.5, stroke = .5, alpha = 1) +
  scale_fill_manual(values = Sarcoma_Immune_Class_palette) +
  facet_grid(`Treatment Arm` ~ `Sarcoma Immune Class`) +
  ggpubr::stat_pvalue_manual(
    SIC_pvals,
    label = "label", y.position = "y.position",
    xmin = "group1", xmax = "group2", tip.length = 0.01
  ) +
  theme_pubr(border = TRUE)

# Figure 6G ---
Eco_data <- tcell_inflamed_data %>% 
  filter(!is.na(ecotype_assignment)) %>% 
  filter(ecotype_assignment != "NA")

Eco_ypos <- Eco_data %>%
  group_by(`Treatment Arm`, ecotype_assignment) %>%
  summarise(y.position = max(IFNy_expanded_avg, na.rm = TRUE) * 1.05, .groups = "drop")

Eco_pvals <- by_Eco_arm %>%
  left_join(Eco_ypos, by = c("Treatment Arm","ecotype_assignment")) %>%
  mutate(group1 = "Pre-treatment", group2 = "Post-treatment",
         label  = ifelse(is.na(p_value), "n.s.", scales::pvalue(p_value, accuracy = 0.001)))

SE_change2 <-
  Eco_data %>% 
  mutate(Sample = paste(Patient, Time_Point, sep = "_")) %>% 
  ggplot(
    aes(x = Time_Point, y = IFNy_expanded_avg, fill = ecotype_assignment)
  ) +
  geom_violin(alpha = .2, trim = FALSE, scale = "width") +
  geom_boxplot(width = .10, outlier.shape = NA, fill = "grey", alpha = .5) +
  geom_line(aes(group = Patient), alpha = 0.5, linewidth = 0.5) +
  geom_point(aes(group = Sample), shape = 21, size = 1.5, stroke = .5, alpha = 1) +
  scale_fill_manual(values = ecotyper_palette) +
  #scale_color_manual(values = c("Experimental" = "#0093af", "Control" = "#9d1536")) + 
  facet_grid(`Treatment Arm` ~ ecotype_assignment) +
  ggpubr::stat_pvalue_manual(
    Eco_pvals,
    label = "label", y.position = "y.position",
    xmin = "group1", xmax = "group2", tip.length = 0.01
  ) +
  theme_pubr(border = TRUE) +
  theme(strip.text = element_text(size = 14))

# Figure Supplemental 2G ----
# Load RNA seq data 
T1_mixture <- fread('...data/SARC32_tpm_T1only.txt')
T1_mixture <- as.data.frame(T1_mixture)
T2_mixture <- fread('...data/SARC32_tpm_T2only.txt')
T2_mixture <- as.data.frame(T2_mixture)
T1_T2_mixture <- T1_mixture %>% 
  left_join(T2_mixture)

# Log transform to to Log2(TPM +1)
T1_T2_mixture <- T1_T2_mixture %>%
  dplyr::mutate(across(where(is.numeric), ~log2(. + 1)))

IC_genes <- T1_T2_mixture %>% 
  filter(Gene %in% c("PDCD1", "CD274", "P4HA1",
                     "PDCD1LG2", "CTLA4", 
                     "LAG3", "HAVCR2", "TIGIT", 
                     "BTLA", "VSIR", "CD276", 
                     "VTCN1", "CD28", "ICOS", 
                     "TNFRSF4", "TNFRSF9", 
                     "CD40", "CD27", 
                     "TNFSF14", "CD200", "CD200R", 
                     "LGALS9", "SIGLEC7", "SIGLEC9", 
                     "ADORA2A", "IDO1", "IDO2", 
                     "NCR3LG1", "HLA-G", "B2M"))
IC_genes <- t(IC_genes)
IC_genes <- as.data.frame(IC_genes)
colnames(IC_genes) <- IC_genes[1,]
IC_genes <- IC_genes[-1, ]
IC_genes <- IC_genes %>%
  dplyr::mutate(across(where(is.character), as.numeric))
IC_genes <- rownames_to_column(IC_genes, var = "Mixture")

IC_genes <- IC_genes %>%
  mutate(patient = str_extract(Mixture, "^[^_]+"))
IC_genes <- IC_genes %>%
  mutate(timepoint = str_extract(Mixture, "(?<=_).+"))

src2id <- setNames(names(src_lookup), unname(src_lookup))
IC_genes$patient <- src2id[IC_genes$patient]

IC_genes <- IC_genes %>% 
  left_join(SIC_and_SE_assignments, by = c("patient" = "Patient"))

IC_genes <- IC_genes %>% 
  left_join(surgery_timing, by = c("patient" = "Subject"))

# Adjust DFS time for post-treatment DFS correlations
IC_genes$new_dfs_time <- IC_genes$`DFS Time (Days)` - IC_genes$time_from_enrollment_to_surgery_days

IC_genes <- IC_genes[IC_genes$patient %in% SIC_and_SE_assignments$Patient, ]

# Initialize a data frame to store Cox PH results
cox_results <- data.frame(Gene = character(),
                          Arm = character(),
                          HR = numeric(),
                          Lower_CI = numeric(),
                          Upper_CI = numeric(),
                          coef = numeric(),
                          se = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

IC_targets <- intersect(c("PDCD1", "CD274", 
                          "PDCD1LG2", "CTLA4", 
                          "LAG3", "HAVCR2", "TIGIT", 
                          "BTLA", "VSIR", "CD276", 
                          "VTCN1", "CD28", "ICOS", 
                          "TNFRSF4", "TNFRSF9", 
                          "CD40", "CD27", 
                          #"P4HA1",
                          "TNFSF14", "CD200", "CD200R", 
                          "LGALS9", "SIGLEC7", "SIGLEC9", 
                          "ADORA2A", "IDO1", "IDO2", 
                          "NCR3LG1", "HLA-G", "B2M"), colnames(IC_genes))

# Loop through each cell type
for (gene in IC_targets) {
  # Subset data by cell type
  subset_data <- IC_genes[IC_genes$timepoint == "T1", c("Mixture", gene, 
                                                        "Treatment Arm", 
                                                        "DFS Event", 
                                                        "DFS Time (Days)", 
                                                        "Sarcoma Immune Class", 
                                                        "timepoint", "patient")]
  str(subset_data)
  print(paste("Gene:", gene))
  print(length(subset_data[[gene]]))
  print(dim(subset_data))
  
  # Combined analysis (all samples together)
  cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ subset_data[[gene]], 
                     data = subset_data)
  
  summary_cox <- summary(cox_model)
  HR <- summary_cox$coefficients[1, "exp(coef)"]
  Lower_CI <- summary_cox$conf.int[1, "lower .95"]
  Upper_CI <- summary_cox$conf.int[1, "upper .95"]
  coef <- summary_cox$coefficients[1, "coef"]
  se <- summary_cox$coefficients[1, "se(coef)"]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
  # Store results for combined analysis
  cox_results <- rbind(cox_results, data.frame(Gene = gene,
                                               Arm = "All",  # Label for combined analysis
                                               HR = HR,
                                               Lower_CI = Lower_CI,
                                               Upper_CI = Upper_CI,
                                               coef = coef,
                                               se = se,
                                               p_value = p_value,
                                               stringsAsFactors = FALSE))
  
  # Separate analysis for Control and Experimental arms
  for (arm in c("Control", "Experimental")) {
    arm_data <- subset(subset_data, `Treatment Arm` == arm)
    str(arm_data)
    
    cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ arm_data[[gene]], data = arm_data)
    
    summary_cox <- summary(cox_model)
    HR <- summary_cox$coefficients[1, "exp(coef)"]
    Lower_CI <- summary_cox$conf.int[1, "lower .95"]
    Upper_CI <- summary_cox$conf.int[1, "upper .95"]
    coef <- summary_cox$coefficients[1, "coef"]
    se <- summary_cox$coefficients[1, "se(coef)"]
    p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
    
    # Store results for each arm
    cox_results <- rbind(cox_results, data.frame(Gene = gene,
                                                 Arm = arm,
                                                 HR = HR,
                                                 Lower_CI = Lower_CI,
                                                 Upper_CI = Upper_CI,
                                                 coef = coef,
                                                 se = se,
                                                 p_value = p_value,
                                                 stringsAsFactors = FALSE))
  }
}

cox_results_sig <- cox_results[cox_results$p_value < 0.05, ]

# Create the data frame 'dt' for the forest plot
dt <- cox_results
# Rename columns to match 'forest()' function requirements
dt$est <- log(dt$HR)
dt$low <- log(dt$Lower_CI)
dt$hi <- log(dt$Upper_CI)
# Optional: Round the estimates and confidence intervals for display
dt$est_formatted <- sprintf("%.2f", dt$est)
dt$ci_formatted <- sprintf("(%.2f - %.2f)", dt$low, dt$hi)
# Create a column combining the estimate and CI for display
dt$Effect_CI <- paste0(dt$est_formatted, " ", dt$ci_formatted)
dt$Gene_copy <- dt$Gene
dt$Gene <- paste(dt$Gene, dt$Arm, sep = "_")
dt <- dt[dt$Gene_copy %in% cox_results_sig$Gene, ]
#Rearrange columns as needed
dt <- dt[, c("Gene", "est", "low", "hi", "se", "p_value", "Effect_CI")]
dt <- dt %>% 
  arrange(Gene)
plot_data <- dt[, c("Gene", "Effect_CI", "p_value")]
# Rename columns for better display
colnames(plot_data) <- c("Gene", "Hazard Ratio (95% CI)", "P-value")
# Convert p-values to formatted strings
plot_data$`P-value` <- sprintf("%.3f", plot_data$`P-value`)
plot_data[, 4] <- ""
colnames(plot_data)[colnames(plot_data) == "V4"] <- "                                                     "
plot_data$Gene <- ifelse(grepl("Control|Experimental", plot_data$Gene), 
                         paste0("   ", plot_data$Gene), 
                         plot_data$Gene)

# Create the forest plot
p <- forest(
  plot_data,
  est = dt$est,
  lower = dt$low,
  upper = dt$hi,
  ci_column = 4,  # Column number where the CI will be plotted
  ref_line = 0, 
  arrow_lab = c("Lower risk of progression", "Higher risk of progression"),
  #xlim = c(0, 4),
  #ticks_at = c(0.5, 1, 2, 3)
  #footnote = "Forest plot of hazard ratios for different cell types."
)
p

# Figure Supplemental 2H ----
# T1 -- CTLA4
subset_data <- IC_genes[IC_genes$timepoint == "T1", ]
threshold <- quantile(IC_genes$CTLA4, 0.5, na.rm = TRUE)
subset_data$Fraction <- ifelse(subset_data$CTLA4 > threshold, "High", "Low")
subset_data$Arm <- subset_data$`Treatment Arm`
surv_fit <- survfit(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction + Arm, 
                    data = subset_data)
ctla4_survplot <- ggsurvplot(surv_fit, data = subset_data,
                             pval = T, 
                             pval.method = T,
                             risk.table = TRUE, 
                             palette = c("#ee9b00", "#9b2226", "#90e0ef", "#0077b6"),
                             xlab = "Days", ylab = "Survival Probability", 
                             title = "Disease Free Survival",
                             legend.title = "CTLA4 Status")

coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction * Arm, data = subset_data)

ggsave(ctla4_survplot$plot,
       width = 8, height = 10, units = "cm",
       filename = '.../CTLA4_T1.pdf')
ggsave(ctla4_survplot$table,
       width = 20, height = 4, units = "cm",
       filename = '.../CTLA4_T1_table.pdf')

# Figure Supplemental 2I ----
# T1 -- TIGIT
subset_data <- IC_genes[IC_genes$timepoint == "T1", ]
threshold <- quantile(IC_genes$TIGIT, 0.5, na.rm = TRUE)
subset_data$Fraction <- ifelse(subset_data$TIGIT > threshold, "High", "Low")
subset_data$Arm <- subset_data$`Treatment Arm`
surv_fit <- survfit(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction + Arm, data = subset_data)
tigit_survplot <- ggsurvplot(surv_fit, data = subset_data,
                             pval = T, 
                             pval.method = T,
                             risk.table = TRUE, 
                             palette = c("#ee9b00", "#9b2226", "#90e0ef", "#0077b6"),
                             xlab = "Days", ylab = "Survival Probability", 
                             title = "Disease Free Survival",
                             legend.title = "TIGIT Status")

cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction * Arm, data = subset_data)
summary_cox <- summary(cox_model)
summary_cox$coefficients

ggsave(tigit_survplot$plot,
       width = 8, height = 10, units = "cm",
       filename = '.../Immune_checkpoints/TIGIT_T1.pdf')
ggsave(tigit_survplot$table,
       width = 20, height = 4, units = "cm",
       filename = '.../TIGIT_T1_table.pdf')

# Figure Supplemental 6B ----
cox_results <- data.frame(Gene = character(),
                          Arm = character(),
                          HR = numeric(),
                          Lower_CI = numeric(),
                          Upper_CI = numeric(),
                          coef = numeric(),
                          se = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)


for (gene in IC_targets) {
  subset_data <- IC_genes[IC_genes$timepoint == "T2", c("Mixture", gene, 
                                                        "Treatment Arm", 
                                                        "DFS Event", 
                                                        "new_dfs_time", 
                                                        "timepoint", "patient")]
  str(subset_data)
  print(paste("Gene:", gene))
  print(length(subset_data[[gene]]))
  print(dim(subset_data))
  
  cox_model <- coxph(Surv(new_dfs_time, `DFS Event` == "Event occurred") ~ subset_data[[gene]], data = subset_data)
  
  summary_cox <- summary(cox_model)
  HR <- summary_cox$coefficients[1, "exp(coef)"]
  Lower_CI <- summary_cox$conf.int[1, "lower .95"]
  Upper_CI <- summary_cox$conf.int[1, "upper .95"]
  coef <- summary_cox$coefficients[1, "coef"]
  se <- summary_cox$coefficients[1, "se(coef)"]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
  # Store results for combined analysis
  cox_results <- rbind(cox_results, data.frame(Gene = gene,
                                               Arm = "All",  # Label for combined analysis
                                               HR = HR,
                                               Lower_CI = Lower_CI,
                                               Upper_CI = Upper_CI,
                                               coef = coef,
                                               se = se,
                                               p_value = p_value,
                                               stringsAsFactors = FALSE))
  
  for (arm in c("Control", "Experimental")) {
    arm_data <- subset(subset_data, `Treatment Arm` == arm)
    str(arm_data)
    
    cox_model <- coxph(Surv(new_dfs_time, `DFS Event` == "Event occurred") ~ arm_data[[gene]], data = arm_data)
    
    summary_cox <- summary(cox_model)
    HR <- summary_cox$coefficients[1, "exp(coef)"]
    Lower_CI <- summary_cox$conf.int[1, "lower .95"]
    Upper_CI <- summary_cox$conf.int[1, "upper .95"]
    coef <- summary_cox$coefficients[1, "coef"]
    se <- summary_cox$coefficients[1, "se(coef)"]
    p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
    
    cox_results <- rbind(cox_results, data.frame(Gene = gene,
                                                 Arm = arm,
                                                 HR = HR,
                                                 Lower_CI = Lower_CI,
                                                 Upper_CI = Upper_CI,
                                                 coef = coef,
                                                 se = se,
                                                 p_value = p_value,
                                                 stringsAsFactors = FALSE))
  }
}

cox_results_sig <- cox_results[cox_results$p_value < 0.05, ]

cox_results_sig <- cox_results_sig[cox_results_sig$Gene != "TNFSF14", ]
dt <- cox_results
dt$est <- log(dt$HR)
dt$low <- log(dt$Lower_CI)
dt$hi <- log(dt$Upper_CI)
dt$est_formatted <- sprintf("%.2f", dt$est)
dt$ci_formatted <- sprintf("(%.2f - %.2f)", dt$low, dt$hi)
dt$Effect_CI <- paste0(dt$est_formatted, " ", dt$ci_formatted)
dt$Gene_copy <- dt$Gene
dt$Gene <- paste(dt$Gene, dt$Arm, sep = "_")
dt <- dt[dt$Gene_copy %in% cox_results_sig$Gene, ]
dt <- dt[, c("Gene", "est", "low", "hi", "se", "p_value", "Effect_CI")]
dt <- dt %>% 
  arrange(Gene)
plot_data <- dt[, c("Gene", "Effect_CI", "p_value")]
colnames(plot_data) <- c("Gene", "Hazard Ratio (95% CI)", "P-value")
plot_data$`P-value` <- sprintf("%.3f", plot_data$`P-value`)
plot_data[, 4] <- ""
colnames(plot_data)[colnames(plot_data) == "V4"] <- "                                                     "
plot_data$Gene <- ifelse(grepl("Control|Experimental", plot_data$Gene), 
                         paste0("   ", plot_data$Gene), 
                         plot_data$Gene)

p <- forest(
  plot_data,
  est = dt$est,
  lower = dt$low,
  upper = dt$hi,
  ci_column = 4,  
  ref_line = 0, 
  arrow_lab = c("Lower risk of progression", "Higher risk of progression"),
)
p

# Figure Supplemental 6C ----
subset_data <- IC_genes[IC_genes$timepoint == "T2", ]
threshold <- quantile(IC_genes$PDCD1, 0.5, na.rm = TRUE)
subset_data$Fraction <- ifelse(subset_data$PDCD1 > threshold, "High", "Low")
subset_data$Arm <- subset_data$`Treatment Arm`
surv_fit <- survfit(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction + Arm, data = subset_data)
PDCD1_survplot <- ggsurvplot(surv_fit, data = subset_data,
                             pval = T, 
                             pval.method = T,
                             risk.table = TRUE, 
                             palette = c("#ee9b00", "#9b2226", "#90e0ef", "#0077b6"),
                             xlab = "Days", ylab = "Survival Probability", 
                             title = "Disease Free Survival",
                             legend.title = "PDCD1 Status")

cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction * Arm, data = subset_data)
summary_cox <- summary(cox_model)
summary_cox$coefficients

ggsave(PDCD1_survplot$plot,
       width = 8, height = 10, units = "cm",
       filename = '.../PD1_T2.pdf')
ggsave(PDCD1_survplot$table,
       width = 20, height = 4, units = "cm",
       filename = '.../PD1_T2_table.pdf')

subset_data$Fraction <- factor(subset_data$Fraction, levels = c("Low","High"))
fit_int <- coxph(Surv(new_dfs_time, `DFS Event` == "Event occurred") ~ Fraction * Arm, data = subset_data)
sm <- summary(fit_int)

#p_int <- sm$coefficients["FractionHigh:ArmExperimental","Pr(>|z|)"]

b <- coef(fit_int); V <- vcov(fit_int)

lincomb <- function(L) {
  logHR <- sum(L * b)
  se    <- sqrt(as.numeric(t(L) %*% V %*% L))
  z     <- logHR / se
  p     <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tibble(HR = exp(logHR),
         CI_low = exp(logHR - 1.96*se),
         CI_high= exp(logHR + 1.96*se),
         p_value = p)
}

L_ctrl <- c("FractionHigh"=1, "ArmExperimental"=0, "FractionHigh:ArmExperimental"=0)
res_ctrl <- lincomb(L_ctrl) %>%
  mutate(Arm = "Control",
         label = sprintf("High vs Low: HR=%.2f (%.2f–%.2f), p=%.3f",
                         HR, CI_low, CI_high, p_value))
L_exp  <- c("FractionHigh"=1, "ArmExperimental"=0, "FractionHigh:ArmExperimental"=1)
res_exp <- lincomb(L_exp) %>%
  mutate(Arm = "Experimental",
         label = sprintf("High vs Low: HR=%.2f (%.2f–%.2f), p=%.3f",
                         HR, CI_low, CI_high, p_value))

# HRs High vs low in control vs experimental 
res_arm <- bind_rows(res_ctrl, res_exp)


# Figure Supplemental 6D ----
subset_data <- IC_genes[IC_genes$timepoint == "T2", ]
threshold <- quantile(IC_genes$NCR3LG1, 0.5, na.rm = TRUE)
subset_data$Fraction <- ifelse(subset_data$NCR3LG1 > threshold, "High", "Low")
subset_data$Arm <- subset_data$`Treatment Arm`
surv_fit <- survfit(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction + Arm, data = subset_data)
b7h6_survplot <- ggsurvplot(surv_fit, data = subset_data,
                            pval = T, 
                            pval.method = T,
                            risk.table = TRUE, 
                            palette = c("#ee9b00", "#9b2226", "#90e0ef", "#0077b6"),
                            xlab = "Days", ylab = "Survival Probability", 
                            title = "Disease Free Survival",
                            legend.title = "NCR3LG1 Status")

cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ Fraction * Arm, data = subset_data)
summary_cox <- summary(cox_model)
summary_cox$coefficients

ggsave(b7h6_survplot$plot,
       width = 8, height = 10, units = "cm",
       filename = '.../B7H6_T2.pdf')
ggsave(b7h6_survplot$table,
       width = 20, height = 4, units = "cm",
       filename = '.../B7H6_T2_table.pdf')

subset_data$Fraction <- factor(subset_data$Fraction, levels = c("Low","High"))
fit_int <- coxph(Surv(new_dfs_time, `DFS Event` == "Event occurred") ~ Fraction * Arm, data = subset_data)
sm <- summary(fit_int)

#p_int <- sm$coefficients["FractionHigh:ArmExperimental","Pr(>|z|)"]

b <- coef(fit_int); V <- vcov(fit_int)

lincomb <- function(L) {
  logHR <- sum(L * b)
  se    <- sqrt(as.numeric(t(L) %*% V %*% L))
  z     <- logHR / se
  p     <- 2 * pnorm(abs(z), lower.tail = FALSE)
  tibble(HR = exp(logHR),
         CI_low = exp(logHR - 1.96*se),
         CI_high= exp(logHR + 1.96*se),
         p_value = p)
}

L_ctrl <- c("FractionHigh"=1, "ArmExperimental"=0, "FractionHigh:ArmExperimental"=0)
res_ctrl <- lincomb(L_ctrl) %>%
  mutate(Arm = "Control",
         label = sprintf("High vs Low: HR=%.2f (%.2f–%.2f), p=%.3f",
                         HR, CI_low, CI_high, p_value))
L_exp  <- c("FractionHigh"=1, "ArmExperimental"=0, "FractionHigh:ArmExperimental"=1)
res_exp <- lincomb(L_exp) %>%
  mutate(Arm = "Experimental",
         label = sprintf("High vs Low: HR=%.2f (%.2f–%.2f), p=%.3f",
                         HR, CI_low, CI_high, p_value))

# HRs High vs low in control vs experimental 
res_arm <- bind_rows(res_ctrl, res_exp)


# Figure Supplemental 6E/6F ----
subset_B7H6 <- IC_genes[IC_genes$timepoint == "T2", c("Mixture", "NCR3LG1", 
                                                      "Treatment Arm", "DFS Event", 
                                                      "DFS Time (Days)",
                                                      "timepoint", "patient")]

subset_B7H6_exp <- subset_B7H6[subset_B7H6$`Treatment Arm` == "Experimental", ]
subset_B7H6_exp <- subset_B7H6_exp[!(is.na(subset_B7H6_exp$Mixture)), ]
subset_B7H6_exp$B7H6_cat <- ifelse(subset_B7H6_exp$NCR3LG1 >= median(subset_B7H6_exp$NCR3LG1), 
                                   "High", "Low")


NK_cell_data <- T1_T2_mixture %>% 
  filter(Gene %in% c("B2M", "HLA-E"))
NK_cell_data <- t(NK_cell_data)
NK_cell_data <- as.data.frame(NK_cell_data)
colnames(NK_cell_data) <- NK_cell_data[1,]
NK_cell_data <- NK_cell_data[-1, ]
NK_cell_data <- NK_cell_data %>%
  dplyr::mutate(across(where(is.character), as.numeric))
NK_cell_data <- rownames_to_column(NK_cell_data, var = "Mixture")
NK_cell_data <- NK_cell_data %>%
  mutate(patient = str_extract(Mixture, "^[^_]+"), 
         timepoint = str_extract(Mixture, "(?<=_).+"))

subset_B7H6_exp <- subset_B7H6_exp %>% 
  select(-c(timepoint, patient)) %>%
  left_join(NK_cell_data[NK_cell_data$Mixture %in% subset_B7H6_exp$Mixture, ], 
            join_by("Mixture"))

plot_B2M <- ggplot(subset_B7H6_exp, 
               aes(x = B7H6_cat, y = B2M, fill = B7H6_cat)) +
  geom_violin(alpha=0.5, trim = F, scale = "width") + 
  geom_boxplot(width = 0.1, outliers = F) + 
  labs(title = paste(resist, "Expression in B7-H6 Low vs High"),
       x = "B7-H6 Expression",
       y = paste(resist, "Expression - Log2(TPM+1)")) +
  theme_pubr(legend = "right") +
  geom_signif(
    comparisons = list(c("Low", "High"))
  ) +
  labs(fill = "B7-H6 Expression") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Low" = "#66bd63", "High" = "#f46d43")) 

ggsave(plot_B2M, filename = '.../data//b2m.pdf',
       width = 12.5, height = 10, units = "cm")

plot_HLAE <- ggplot(subset_B7H6_exp, 
                   aes(x = B7H6_cat, y = `HLA-E`, fill = B7H6_cat)) +
  geom_violin(alpha=0.5, trim = F, scale = "width") + 
  geom_boxplot(width = 0.1, outliers = F) + 
  labs(title = paste(resist, "Expression in B7-H6 Low vs High"),
       x = "B7-H6 Expression",
       y = paste(resist, "Expression - Log2(TPM+1)")) +
  theme_pubr(legend = "right") +
  geom_signif(
    comparisons = list(c("Low", "High"))
  ) +
  labs(fill = "B7-H6 Expression") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Low" = "#66bd63", "High" = "#f46d43")) 

ggsave(plot_HLAE, filename = '.../data//HLAE.pdf',
       width = 12.5, height = 10, units = "cm")


# Figure Supplemental 6G ----
T_cell_GEP_genes <- c( "CD3D",  "IL2RG",
                          "IDO1", "NKG7",
                          "CIITA", "HLA-E",
                          "CD3E", "CXCR6",
                          "CCL5", "LAG3",
                          "GZMK", "TAGAP",
                          "CD2", "CXCL10",
                          "HLA-DRA", "STAT1",
                          "CXCL13", "GZMB")

sample_annotation <- IC_genes[, c("Mixture", "timepoint", "Treatment Arm", "Sarcoma Immune Class")]
#sample_annotation <- sample_annotation[!(is.na(sample_annotation$Arm)), ]
sample_annotation <- sample_annotation[sample_annotation$`Sarcoma Immune Class` != "NA", ]

T1_T2_mixture_filtered <- T1_T2_mixture[T1_T2_mixture$Gene %in% unique(c(T_cell_GEP_genes)), ]
rownames(T1_T2_mixture_filtered) <- T1_T2_mixture_filtered$Gene
T1_T2_mixture_filtered$Gene <- NULL

T1_T2_mixture_filtered <- T1_T2_mixture_filtered[, colnames(T1_T2_mixture_filtered) %in% sample_annotation$Mixture]

z_score_rows <- function(x) {
  t(apply(x, 1, function(row) {
    (row - mean(row)) / sd(row)
  }))
}
T1_T2_mixture_filtered_zscored <- z_score_rows(as.matrix(T1_T2_mixture_filtered))

sample_annotation <- sample_annotation[match(colnames(T1_T2_mixture_filtered),
                                             sample_annotation$Mixture), ]

rownames(sample_annotation) <- sample_annotation$Mixture
sample_annotation$Mixture <- NULL

# Define a column annotation
col_annotation <- HeatmapAnnotation(
  df = sample_annotation,
  col = list(timepoint = c("T1" = "#386cb0", "T2" = "#fdc086"), 
             `Treatment Arm` = c("Experimental" = "#0093AF", "Control" = "#9E1C39"), 
             `Sarcoma Immune Class` = c("A" = "#1f78b4", "B" = "#a6cee3", "C" = "#33a02c", "D" = "#fdbf6f", "E" = "#e31a1c"))
)

#annotation_vector <- sample_annotation$Arm
annotation_vector <- sample_annotation$`Sarcoma Immune Class`
annotation_vector <- as.factor(annotation_vector)

# Plot the heatmap
IFNG_score_heatmap <- Heatmap(
  T1_T2_mixture_filtered_zscored,
  name = "Expression",
  cluster_column_slices = F,
  top_annotation = col_annotation,
  cluster_columns = TRUE,  # Disable clustering for columns
  cluster_rows = TRUE,      # Enable clustering for rows
  column_split = annotation_vector,  # Group columns by Sarcoma Immune Class
  show_column_names = F,
  show_row_names = T,
  #cluster_column_slices = FALSE,
  heatmap_width = unit(200, "mm"), 
  heatmap_height = unit(120, "mm"),
  #row_title = "Cell Types",
  #column_title = "T2 Samples SIC vs sc Cell States",
  heatmap_legend_param = list(title = "z-score")
)

ht <- draw(IFNG_score_heatmap, heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../heatmap_arm_vs_timepoint_IFNG_sig_genes_SIC.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure Supplemental 6H ----
sample_annotation <- IC_genes[, c("Sample", "Time_Point", "Arm", "Sarcoma Ecotype Assignment")]
#sample_annotation <- sample_annotation[!(is.na(sample_annotation$Arm)), ]
sample_annotation <- sample_annotation[sample_annotation$`Sarcoma Ecotype Assignment` != "NA", ]

T1_T2_mixture_filtered <- T1_T2_mixture[T1_T2_mixture$Gene %in% unique(c(T_cell_GEP_genes)), ]
rownames(T1_T2_mixture_filtered) <- T1_T2_mixture_filtered$Gene
T1_T2_mixture_filtered$Gene <- NULL

T1_T2_mixture_filtered <- T1_T2_mixture_filtered[, colnames(T1_T2_mixture_filtered) %in% 
                                                   sample_annotation$Mixture]

z_score_rows <- function(x) {
  t(apply(x, 1, function(row) {
    (row - mean(row)) / sd(row)
  }))
}
T1_T2_mixture_filtered_zscored <- z_score_rows(as.matrix(T1_T2_mixture_filtered))

sample_annotation <- sample_annotation[match(colnames(T1_T2_mixture_filtered),
                                             sample_annotation$Mixture), ]

rownames(sample_annotation) <- sample_annotation$Mixture
sample_annotation$Mixture <- NULL

# Define a column annotation
col_annotation <- HeatmapAnnotation(
  df = sample_annotation,
  col = list(timepoint = c("T1" = "#386cb0", "T2" = "#fdc086"), 
             `Treatment Arm` = c("Experimental" = "#0093AF", "Control" = "#9E1C39"), 
             `Sarcoma Ecotype Assignment` = c("SE1" = "#F16280", "SE2" = "#179ABE", "SE3" = "#63C2A6"))
)

#annotation_vector <- sample_annotation$Arm
annotation_vector <- sample_annotation$`Sarcoma Ecotype Assignment`
annotation_vector <- as.factor(annotation_vector)

# Plot the heatmap
IFNG_score_heatmap <- Heatmap(
  T1_T2_mixture_filtered_zscored,
  name = "Expression",
  cluster_column_slices = F,
  top_annotation = col_annotation,
  cluster_columns = TRUE,  # Disable clustering for columns
  cluster_rows = TRUE,      # Enable clustering for rows
  column_split = annotation_vector,  # Group columns by Sarcoma Immune Class
  show_column_names = F,
  show_row_names = T,
  #cluster_column_slices = FALSE,
  heatmap_width = unit(200, "mm"), 
  heatmap_height = unit(120, "mm"),
  #row_title = "Cell Types",
  #column_title = "T2 Samples SIC vs sc Cell States",
  heatmap_legend_param = list(title = "z-score")
)

ht <- draw(IFNG_score_heatmap, heatmap_legend_side = "left", annotation_legend_side = "left", merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../heatmap_arm_vs_timepoint_IFNG_sig_genes_SEs.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure 6H ----
simpson_diversity <- fread('...data/gini_simpson_diversity.csv')

simpson_diversity$V1 <- NULL

simpson_diversity <- simpson_diversity %>%
  left_join(SIC_and_SE_assignments, by = c("Sample" = "Patient"))

ecotyper_palette = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")

# Ecotype plot considering only Pre-treatment time point  
simpson_diversity |> 
  filter(!is.na(`Sarcoma Ecotype Assignment`))|> 
  filter(`Sarcoma Ecotype Assignment` != "NA") |> 
  filter(Time_Point == "T1") |>
  filter(!(is.na(`Treatment Arm`))) |> 
  mutate(ecotype_assignment = `Sarcoma Ecotype Assignment`) |>
  tidyplot(x = ecotype_assignment, 
           y = gini_simpson, 
           color = ecotype_assignment) |> 
  add_mean_bar(alpha = 0.4) |> 
  add_sem_errorbar() |> 
  add_data_points_beeswarm() |> 
  add_test_pvalue(ref.group = 1, 
                  method = "wilcox_test", 
                  p.adjust.method = "BH", 
                  label.size = 5) |>
  adjust_colors(ecotyper_palette) |>
  adjust_size(width = 10, height = 10, unit = "cm") |> 
  adjust_font(fontsize = 14) |> 
  adjust_x_axis_title(title = "Sarcoma Ecotypes", fontsize = 16) |>
  adjust_y_axis_title(title = "TCR Diversity (Gini-Simpson)", fontsize = 16) 

# Figure Supplemental 6J ----
simpson_diversity |> 
  filter(!is.na(`Sarcoma Immune Class`))|> 
  filter(`Sarcoma Immune Class` != "NA")|> 
  filter(Time_Point == "T1")|> 
  filter(!(is.na(`Treatment Arm`))) |> 
  tidyplot(x = `Sarcoma Immune Class`, 
           y = gini_simpson, 
           color = `Sarcoma Immune Class`) |> 
  add_mean_bar(alpha = 0.4) |> 
  add_sem_errorbar() |> 
  add_data_points_beeswarm() |> 
  add_test_pvalue(ref.group = 1, 
                  method = "wilcox_test", 
                  p.adjust.method = "BH", 
                  label.size = 5) |>
  adjust_colors(Sarcoma_Immune_Class_palette) |>
  adjust_size(width = 10, height = 10, unit = "cm") |> 
  adjust_font(fontsize = 14) |> 
  adjust_x_axis_title(title = "Sarcoma Immune Classes", fontsize = 16) |>
  adjust_y_axis_title(title = "TCR Diversity (Gini-Simpson)", fontsize = 16) 

# Figure 6I ----
time_arm <- simpson_diversity |>
  filter(!(is.na(`Treatment Arm`))) |> 
  tidyplot(
    x     = Time_Point,
    y     = gini_simpson,
    color = Time_Point
  ) |>
  add_mean_bar(alpha = 0.4) |>
  add_sem_errorbar() |>
  add_data_points_beeswarm() |>
  add_test_pvalue(
    method = "wilcox_test", 
    ref.group       = "T1",
    p.adjust.method = "BH",
    label.size      = 5
  ) |>
  adjust_colors(c("T1" = "#396CAE", "T2" = "#FBBF85")) |>
  adjust_size(width = 12.5, height = 10, unit = "cm") |>
  adjust_font(fontsize = 14) |>
  adjust_x_axis_title(title = "Time Point", fontsize = 16) |>
  adjust_y_axis_title(title = "TCR Diversity (Gini-Simpson)", fontsize = 16) +
  facet_wrap(~ `Treatment Arm`, scales = "free_x")


# Flow Cytometry vs CIBERSORTx cell fractions correlation Figure S4 ----
# Figure S4B ----
#Total CD4 correlation 
cd4_corr <- read.csv(file = '.../data/cd4_flow_corr.csv')
cd4_cibersortx_homegrown <- ggscatter(cd4_corr, x = "Flow_CD3_CD4", y = "total_cibersort_cd4", 
                                      xlab = "% Singlets CD3+ CD4+ Cells (Flow Cytometry)", 
                                      ylab = "Total CD4+ Fractions (CIBERSORTx)",
                                      add = "reg.line", 
                                      add.params = list(color = "#cb181d", fill = "lightgray"), # Customize reg. line
                                      conf.int = TRUE, # Add confidence interval
                                      color = "#ec7014") + 
  stat_cor(method = "spearman", label.x = 0.05, label.y = 0.25)

#Total CD8 correlation 
cd8_corr <- read.csv(file = '.../data/cd8_flow_corr.csv')
cd8_cibersortx_homegrown <- ggscatter(cd8_corr, x = "Flow_CD3_CD8", y = "total_cibersort_cd8", 
                                      xlab = "% Singlets CD3+ CD8+ Cells (Flow Cytometry)", 
                                      ylab = "Total CD8+ Fractions (CIBERSORTx)",
                                      add = "reg.line", 
                                      add.params = list(color = "#cb181d", fill = "lightgray"), # Customize reg. line
                                      conf.int = TRUE, # Add confidence interval
                                      color = "#0570b0") + 
  stat_cor(method = "spearman", label.x = 0.05, label.y = 0.32)

# B cells 
bcells_corr <- read.csv(file = '.../data/bcell_flow_corr.csv')
cd19_cibersortx_homegrown <- ggscatter(bcells_corr, x = "Flow_CD19", y = "total_cibersort_bcells", 
                                       xlab = "% Singlets CD19+ Cells (Flow Cytometry)", 
                                       ylab = "Total B-cells Fractions (CIBERSORTx)",
                                       add = "reg.line", 
                                       add.params = list(color = "#cb181d", fill = "lightgray"), # Customize reg. line
                                       conf.int = TRUE, # Add confidence interval
                                       color = "#238b45") + 
  stat_cor(method = "spearman", label.x = 0.01, label.y = 0.04)

# NKcells total 
nkcells_corr <- read.csv(file = '.../data/nkcells_corr_corr_flow_corr.csv')
nkcells_cibersortx_homegrown <- ggscatter(nkcells_corr, x = "Flow_nkcells", y = "total_cibersort_nkcells",  
                                          xlab = "% Singlets NK cells (Flow Cytometry)", 
                                          ylab = "Total NK cells Fractions (CIBERSORTx)",
                                          add = "reg.line", 
                                          add.params = list(color = "#cb181d", fill = "lightgray"), # Customize reg. line
                                          conf.int = TRUE, # Add confidence interval
                                          color = "#4a1486") + 
  stat_cor(method = "spearman", label.x = 0.005, label.y = 0.07)

# Total myeloid cells 
myeloid_corr <- read.csv(file = '.../data/myeloid_corr_flow_corr.csv')
total_myeloid_cibersortx_homegrown <- ggscatter(myeloid_corr, x = "Flow_myeloid", y = "total_cibersort_myeloid", 
                                                xlab = "% Singlets Total Myeloid Cells (Flow Cytometry)", 
                                                ylab = "Total Myeloid Cells (CIBERSORTx)",
                                                add = "reg.line",  
                                                add.params = list(color = "#cb181d", fill = "lightgray"), # Customize reg. line
                                                conf.int = TRUE, # Add confidence interval
                                                color = "#c51b7d") + 
  stat_cor(method = "spearman", label.x = 0.05, label.y = 0.6)



# Figure Supplemental 5A ----
combined_rejoined_cpm_tr4_init$Mixture <- NULL
cibersortx_long_cpm_tr4 <- melt(combined_rejoined_cpm_tr4_init, 
                                id.vars = c("patient_id", "timepoint"),
                                variable.name = "Cell_Type", 
                                value.name = "Value")
colnames(cibersortx_long_cpm_tr4)[colnames(cibersortx_long_cpm_tr4) == "patient_id"] <- "Sample"
colnames(cibersortx_long_cpm_tr4)[colnames(cibersortx_long_cpm_tr4) == "timepoint"] <- "Time_Point"
cibersortx_long_cpm_tr4 <- cibersortx_long_cpm_tr4 %>% 
  left_join(SIC_and_SE_assignments, by = c("Sample" = "Patient"))

cibersortx_long_cpm_tr4$Cell_Type <- factor(cibersortx_long_cpm_tr4$Cell_Type)
cibersortx_long_cpm_tr4$Sample <- factor(cibersortx_long_cpm_tr4$Sample)

heatmap_matrix <- cibersortx_long_cpm_tr4 %>%
  mutate(Sample_TP = paste0(Sample, "_", Time_Point)) %>%
  select(Sample_TP, Cell_Type, Value) %>%
  pivot_wider(
    names_from = Sample_TP, 
    values_from = Value
  ) %>%
  column_to_rownames("Cell_Type") %>%
  as.matrix()

heatmap_matrix <- heatmap_matrix[!(grepl("total_", rownames(heatmap_matrix))),]

# arcsin + z-scoring
# 1. perform arcsin normalization on the fractions 
data_asin <- asin(sqrt(heatmap_matrix))
# 2. Perform row-based z-score scaling
scaled_data <- t(scale(t(data_asin)))
scaled_data_T1 <- scaled_data[, grepl("_T1", colnames(scaled_data))]

row_annotation <- data.frame(cells = rownames(scaled_data_T1))
row_annotation$cell_types <- NA
row_annotation <- row_annotation %>%
  mutate(
    cell_types = case_when(
      grepl("CD4",         cells) ~ "CD4+ T-cells",
      grepl("CD8",         cells) ~ "CD8+ T-cells",
      grepl("NK cells",    cells) ~ "NK cells",
      grepl("DCs",         cells) ~ "Dendritic cells",
      grepl("PMNs",        cells) ~ "PMNs",
      grepl("Monocytes|M2|Macro", cells) ~ "Monocytes & Macrophages",
      grepl("Mastocytes",  cells) ~ "Mastocytes",
      grepl("Pericytes",   cells) ~ "Pericytes",
      grepl("CAF",         cells) ~ "CAFs",
      grepl("Endothelial", cells) ~ "Endothelial cells",
      grepl("B Cells & Plasma Cells", cells) ~ "B Cells & Plasma Cells",
      grepl("STS",         cells) ~ "Sarcoma cells",
      TRUE ~ cells  # fallback: leave as-is
    )
  )

# T1 only 
row_annotation <- column_to_rownames(row_annotation, var = "cells")
row_annotation$cell_types <- factor(
  row_annotation$cell_types,
  levels = c("CD8+ T-cells", "CD4+ T-cells", "NK cells", 
             "B Cells & Plasma Cells", "Monocytes & Macrophages", 
             "Dendritic cells", "PMNs", "Mastocytes", 
             "Sarcoma cells", "CAFs", "Endothelial cells", 
             "Pericytes"))

column_annotation <- data.frame(samples = colnames(scaled_data_T1))
#rownames(column_annotation) <- column_annotation$samples

column_annotation <- column_annotation %>%
  mutate(timepoint = str_extract(samples, "(?<=_).+"), 
         patient = str_extract(samples, "^[^_]+"))

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$patient %in% cibersortx_long_cpm_tr4[cibersortx_long_cpm_tr4$`Treatment Arm` == "Experimental", ]$Sample] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Immune Class", 
                               "Tumor Grade")], by = c("patient" = "Patient")) 

column_annotation$`Sarcoma Immune Class`[is.na(column_annotation$`Sarcoma Immune Class`)] <- "Unassigned"
column_annotation <- column_to_rownames(column_annotation, var = "samples")

column_annotation$patient <- NULL
column_annotation$timepoint <- NULL

column_annotation <- column_annotation[colnames(scaled_data_T1), ]

column_annotation <- column_annotation %>% 
  relocate(`Sarcoma Immune Class`)

column_annotation$`Sarcoma Immune Class` <- factor(
  column_annotation$`Sarcoma Immune Class`,
  levels = c("A", "B", "C", "D", "E")
)

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#a6dba0", "Grade 3" = "#984ea3"), 
             Arm = c("Experimental" = "#0093AF", "Control" = "#9E1C39"), 
             `Sarcoma Immune Class` = c("A" = "#1f78b4", "B" = "#a6cee3", 
                                        "C" = "#33a02c", "D" = "#fdbf6f", 
                                        "E" = "#e31a1c")))

col_custom <- colorRamp2(c(-4, 0, 4), c("#2166ac", "#ffffff", "#b2182b"))

ht = Heatmap(scaled_data_T1, 
             split = as.factor(row_annotation$cell_types), 
             column_split = column_annotation$`Sarcoma Immune Class`,
             cluster_column_slices = F,
             cluster_row_slices = F,
             row_title_rot = 0,
             border = T,
             name = "Row Z-score", 
             heatmap_height = unit(250, "mm"), 
             heatmap_width = unit(330, "mm"), 
             top_annotation = col_annotation,
             show_column_names = F, 
             col = col_custom)

ht = draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", 
          merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf('/.../scaled_T1_absolute_fractions.pdf', width = w, height = h)
draw(ht)
dev.off()

# Figure Supplemental 5B -----
cibersortx_long_cpm_tr4 <- melt(combined_rejoined_cpm_tr4_init, 
                                id.vars = c("patient_id", "timepoint"),
                                variable.name = "Cell_Type", 
                                value.name = "Value")
colnames(cibersortx_long_cpm_tr4)[colnames(cibersortx_long_cpm_tr4) == "patient_id"] <- "Sample"
colnames(cibersortx_long_cpm_tr4)[colnames(cibersortx_long_cpm_tr4) == "timepoint"] <- "Time_Point"
cibersortx_long_cpm_tr4 <- cibersortx_long_cpm_tr4 %>% 
  left_join(SIC_and_SE_assignments, by = c("Sample" = "Patient"))

cibersortx_long_cpm_tr4$Cell_Type <- factor(cibersortx_long_cpm_tr4$Cell_Type)
cibersortx_long_cpm_tr4$Sample <- factor(cibersortx_long_cpm_tr4$Sample)

heatmap_matrix <- cibersortx_long_cpm_tr4 %>%
  mutate(Sample_TP = paste0(Sample, "_", Time_Point)) %>%
  select(Sample_TP, Cell_Type, Value) %>%
  pivot_wider(
    names_from = Sample_TP, 
    values_from = Value
  ) %>%
  column_to_rownames("Cell_Type") %>%
  as.matrix()

heatmap_matrix <- heatmap_matrix[!(grepl("total_", rownames(heatmap_matrix))),]

# arcsin + z-scoring
# 1. perform arcsin normalization on the fractions 
data_asin <- asin(sqrt(heatmap_matrix))
# 2. Perform row-based z-score scaling
scaled_data <- t(scale(t(data_asin)))
scaled_data_T1 <- scaled_data[, grepl("_T1", colnames(scaled_data))]

row_annotation <- data.frame(cells = rownames(scaled_data_T1))
row_annotation$cell_types <- NA
row_annotation <- row_annotation %>%
  mutate(
    cell_types = case_when(
      grepl("CD4",         cells) ~ "CD4+ T-cells",
      grepl("CD8",         cells) ~ "CD8+ T-cells",
      grepl("NK cells",    cells) ~ "NK cells",
      grepl("DCs",         cells) ~ "Dendritic cells",
      grepl("PMNs",        cells) ~ "PMNs",
      grepl("Monocytes|M2|Macro", cells) ~ "Monocytes & Macrophages",
      grepl("Mastocytes",  cells) ~ "Mastocytes",
      grepl("Pericytes",   cells) ~ "Pericytes",
      grepl("CAF",         cells) ~ "CAFs",
      grepl("Endothelial", cells) ~ "Endothelial cells",
      grepl("B Cells & Plasma Cells", cells) ~ "B Cells & Plasma Cells",
      grepl("STS",         cells) ~ "Sarcoma cells",
      TRUE ~ cells  # fallback: leave as-is
    )
  )

# T1 only 
row_annotation <- column_to_rownames(row_annotation, var = "cells")
row_annotation$cell_types <- factor(
  row_annotation$cell_types,
  levels = c("CD8+ T-cells", "CD4+ T-cells", "NK cells", 
             "B Cells & Plasma Cells", "Monocytes & Macrophages", 
             "Dendritic cells", "PMNs", "Mastocytes", 
             "Sarcoma cells", "CAFs", "Endothelial cells", 
             "Pericytes"))

column_annotation <- data.frame(samples = colnames(scaled_data_T1))
column_annotation <- column_annotation %>%
  mutate(timepoint = str_extract(samples, "(?<=_).+"), 
         patient = str_extract(samples, "^[^_]+"))

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$patient %in% cibersortx_long_cpm_tr4[cibersortx_long_cpm_tr4$`Treatment Arm` == "Experimental", ]$Sample] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Ecotype Assignment", 
                               "Tumor Grade")], by = c("patient" = "Patient")) 

column_annotation$ecotype_assignment <- column_annotation$`Sarcoma Ecotype Assignment`

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation$patient <- NULL
column_annotation$timepoint <- NULL
column_annotation <- column_annotation[colnames(scaled_data_T1), ]
column_annotation <- column_annotation %>% 
  relocate(`ecotype_assignment`)

column_annotation$ecotype_assignment <- factor(
  column_annotation$ecotype_assignment,
  levels = c("SE1", "SE2", "SE3")
)

column_annotation <- column_annotation[!(is.na(column_annotation$ecotype_assignment)), ]
scaled_data_T1 <- scaled_data_T1[, colnames(scaled_data_T1) %in% rownames(column_annotation)]
column_annotation$`Sarcoma Ecotype Assignment` <- NULL

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#a6dba0", "Grade 3" = "#984ea3"), 
             Arm = c("Experimental" = "#0093AF", "Control" = "#9E1C39"), 
             `ecotype_assignment` = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")))

col_custom <- colorRamp2(c(-4, 0, 4), c("#2166ac", "#ffffff", "#b2182b"))

ht = Heatmap(scaled_data_T1, 
             split = as.factor(row_annotation$cell_types), 
             column_split = column_annotation$ecotype_assignment,
             cluster_column_slices = F,
             cluster_row_slices = F,
             row_title_rot = 0,
             border = T,
             name = "Row Z-score", 
             heatmap_height = unit(250, "mm"), 
             heatmap_width = unit(330, "mm"), 
             top_annotation = col_annotation,
             show_column_names = F, 
             col = col_custom)

ht = draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", 
          merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf('.../scaled_T1_absolute_fractions_EcoTyper.pdf', width = w, height = h)
draw(ht)
dev.off()


