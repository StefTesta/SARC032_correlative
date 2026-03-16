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
library(ggalluvial)
library(factoextra)
library(cluster)
library(NMF)
library(tidytext)
library(nnls)
library(scPCA)
library(ggdendro)
library(patchwork)
library(zCompositions)
library(limma)
library(Matrix)
library(ggridges)
library(LymphoSeq)
library(immunarch)
library(cobalt)

readr::local_edition(1)

# Load data from Supp. Tables ----
flow_results_filtered <- read_excel(".../data//Supplementary_Tables.xlsx", 
                                    sheet = 12)
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

# Load EcoType Abundance data ---- 
ecotype_abundance <- read_excel('.../Supplementary_Tables_flow_updated.xlsx', 
                                sheet = 8)

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
  dplyr::select(patient, Time_Point, Ecotype, abundance, `Treatment Arm`)

ecotyper_palette = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")


# Max Ecotype assignment for pre-treatment samples ----
pre_ecotype_calls_wide <- df_long %>%
  filter(Time_Point == "Pre-treatment",
         Ecotype %in% c("SE1","SE2","SE3")) %>%
  dplyr::select(patient, `Treatment Arm`, Ecotype, abundance) %>%
  pivot_wider(names_from = Ecotype, values_from = abundance) %>%
  mutate(
    Assigned_Max_Ecotype = c("SE1","SE2","SE3")[max.col(cbind(SE1, SE2, SE3), ties.method = "first")],
    Max_SE = pmax(SE1, SE2, SE3, na.rm = TRUE)
  )

###### Alluvial plot SIC vs SE  Figure S2E ###### -----

# 1) Keep only pre-treatment rows
df_pre <- df_long %>%
  dplyr::filter(Time_Point == "Pre-treatment")

# 2) Clean up column names 
SIC_and_SE_clean <- SIC_and_SE_assignments %>%
  dplyr::rename(
    patient      = Patient,
    TreatmentArm = `Treatment Arm`,
    SarcomaImmuneClass = `Sarcoma Immune Class`,
    SarcomaEcotypeAssignment = `Sarcoma Ecotype Assignment`,
    TumorGrade  = `Tumor Grade`,
    DFS_event   = `DFS Event`,
    DFS_time    = `DFS Time (Days)`,
    PercentNecrosis = `Percent Necrosis`,
    SRC_code    = `SRC code`
  )

df_pre_clean <- df_pre %>%
  dplyr::rename(
    TreatmentArm = `Treatment Arm`
  )

df_merged <- df_pre_clean %>%
  dplyr::inner_join(
    SIC_and_SE_clean,
    by = c("patient", "TreatmentArm")
  )

df_merged_max = df_merged %>% 
  left_join(pre_ecotype_calls_wide[, c("patient", "Assigned_Max_Ecotype")])

river_df <- df_merged_max %>%
  distinct(patient, SarcomaImmuneClass, Assigned_Max_Ecotype) %>%
  mutate(
    SIC = ifelse(is.na(SarcomaImmuneClass) | SarcomaImmuneClass == "",
                 "SIC NA", SarcomaImmuneClass),
    SE  = ifelse(is.na(Assigned_Max_Ecotype) | Assigned_Max_Ecotype == "",
                 "SE mixed", Assigned_Max_Ecotype)
  ) %>%
  group_by(SIC, SE) %>%
  summarise(n = n(), .groups = "drop")

river_df$SE[river_df$SE == "NA"] <- "SE mixed"
river_df$SE[is.na(river_df$SE)] <- "SE mixed"

river_df <- river_df %>%
  mutate(
    SIC = factor(SIC, levels = c("A", "B", "C", "D", "E", "SIC NA")),
    SE  = factor(SE,  levels = c("SE1", "SE2", "SE3"))
  )

Sarcoma_Immune_Class_palette = c("A" = "#1f78b4", 
                                 "B" = "#a6cee3", 
                                 "C" = "#33a02c", 
                                 "D" = "#fdbf6f", 
                                 "E" = "#e31a1c")

ecotyper_palette = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")

extra_palette <- c(
  "SIC NA"  = "grey50",
  "SE mixed" = "grey60"
)

label_palette <- c(Sarcoma_Immune_Class_palette,
                   ecotyper_palette,
                   extra_palette)

alluvional_detailed = ggplot(river_df,
                             aes(axis1 = SIC, axis2 = SE, y = n)) +
  geom_alluvium(aes(fill = SIC), width = 0, alpha = 0.6) +
  geom_stratum(width = 0.15, fill = "grey92", color = "grey40") +
  geom_text(
    stat  = "stratum",
    aes(label = after_stat(stratum),
        colour = after_stat(stratum)),
    size = 3.2,
    vjust = -0.1
  ) +
  scale_x_discrete(
    limits = c("Sarcoma immune class", "Sarcoma ecotype"),
    expand = c(.1, .1)
  ) +
  scale_fill_manual(values = Sarcoma_Immune_Class_palette, drop = FALSE) +
  scale_colour_manual(values = label_palette, guide = "none") +
  labs(
    x = NULL,
    y = "Number of patients",
    fill = "SIC",
    title = "Correspondence between sarcoma immune classes (SIC) and sarcoma ecotypes (SE)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    panel.grid    = element_blank(),
    legend.position = "right",
    plot.title    = element_text(face = "bold")
  )

ggsave(alluvional_detailed, width = 12.5, height = 12.5, units = 'cm',
       filename = '.../Alluvional_plot_Fig_S2E.pdf')

###### Figure S2d Heat map SE abundance vs SIC ###### ------
se_by_sic <- df_merged %>%
  mutate(
    SIC = ifelse(is.na(SarcomaImmuneClass) | SarcomaImmuneClass == "",
                 "SIC NA", SarcomaImmuneClass)
  ) %>%
  group_by(SIC, Ecotype) %>%
  summarise(
    mean_abundance   = mean(abundance, na.rm = TRUE),
    median_abundance = median(abundance, na.rm = TRUE),
    n_patients       = n_distinct(patient),
    .groups = "drop"
  )

se_by_sic_mtx <- se_by_sic %>%
  select(SIC, Ecotype, median_abundance) %>%
  pivot_wider(names_from = Ecotype,
              values_from = median_abundance)

se_by_sic_mtx <- column_to_rownames(se_by_sic_mtx, var = 'SIC')

row_annotation = HeatmapAnnotation(df = data.frame(SIC = c('A', 'B', 'C', 'D', 'E')), 
                                   col = list(SIC = Sarcoma_Immune_Class_palette),  
                                   which = "row", gp = gpar(col = NA, lwd = 0), 
                                   show_legend = T, show_annotation_name = F)

column_annotation = HeatmapAnnotation(df = data.frame(SEs = c('SE1', 'SE2', 'SE3')), 
                                      col = list(SEs = ecotyper_palette),  
                                      which = "column", gp = gpar(col = NA, lwd = 0), 
                                      show_legend = T, show_annotation_name = F)

ht = Heatmap(as.matrix(se_by_sic_mtx), 
             cluster_rows = F, 
             border = T,
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 0.5),
             col = colorRamp2(quantile(as.matrix(se_by_sic_mtx)), 
                              viridis(5, option = "C")), 
             top_annotation = column_annotation, 
             name = 'Median SE Abundance',
             heatmap_height =  unit(8, "cm"), 
             heatmap_width =  unit(5, "cm"), 
             left_annotation = row_annotation)

ht <- draw(ht, heatmap_legend_side = "left", annotation_legend_side = "left", 
           merge_legend = TRUE)
w = ComplexHeatmap:::width(ht)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(ht)
h = convertY(h, "inch", valueOnly = TRUE)
pdf(".../FigureS2d_Median_ecotype_abundance_per_SIC.pdf", width = w, height = h)
draw(ht)
dev.off()

###### Figure S2f Barplot SIC vs SE assignment ###### ------
sic_se_patient <- SIC_and_SE_clean %>%
  filter(SarcomaImmuneClass != "NA") %>% # these are patients with no pre-treatment sample
  mutate(
    SIC = ifelse(is.na(SarcomaImmuneClass) | SarcomaImmuneClass == "",
                 "SIC NA", SarcomaImmuneClass),
    SE_assign = ifelse(is.na(SarcomaEcotypeAssignment) | SarcomaEcotypeAssignment == "NA",
                       "SE mixed", SarcomaEcotypeAssignment)
  ) %>%
  select(patient, TreatmentArm, SIC, SE_assign, TumorGrade, Histology, DFS_event, DFS_time)

table_SIC <- sic_se_patient %>%
  dplyr::count(SIC) %>%
  mutate(prop = n / sum(n))

table_SE <- sic_se_patient %>%
  dplyr::count(SE_assign) %>%
  mutate(prop = n / sum(n))

# cross-tab SIC x SE_assign
xtab_SIC_SE <- sic_se_patient %>%
  dplyr::count(SIC, SE_assign) %>%
  group_by(SIC) %>%
  mutate(row_prop = n / sum(n)) %>%
  ungroup()

se_levels <- c("SE1", "SE2", "SE3", "SE mixed")
xtab_SIC_SE_plot <- xtab_SIC_SE %>%
  mutate(
    SE_assign = factor(SE_assign, levels = se_levels),
    SIC       = factor(SIC, levels = c("A", "B", "C", "D", "E"))
  )

ecotyper_palette_SEassign <- c(
  "SE1"      = "#F06180",
  "SE2"      = "#0F9ABE",
  "SE3"      = "#60C1A5",
  "SE mixed" = "grey70"
)

p_SIC_SE_stacked <- ggplot(
  xtab_SIC_SE_plot,
  aes(x = SIC, y = row_prop, fill = SE_assign)
) +
  geom_col(color = "black", width = 0.7, linewidth = 0.2) +
  scale_fill_manual(values = ecotyper_palette_SEassign) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Sarcoma immune class (SIC)",
    y = "Proportion of tumors",
    fill = "Dominant ecotype",
    title = "Distribution of EcoTyper ecotype assignments within each SIC"
  ) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p_SIC_SE_stacked
ggsave(p_SIC_SE_stacked, filename = '.../Figure_S2E_barplot_SEs_assignment_vs_SIC.pdf', 
       width = 12.5, height = 10, units = "cm")

###### Figure Supplemental 2g - DFS Forrest Plot T cell Inflamed GEP ###### ------
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
  tcell_inflamed_data$Time_Point <- tcell_inflamed_data$`Time Point`
  subset_data <- tcell_inflamed_data[tcell_inflamed_data$Time_Point == "Pre-treatment", c(signature, 
                                                                                          "Treatment Arm", 
                                                                                          "DFS Event", 
                                                                                          "DFS Time (Days)")]
  str(subset_data)
  print(paste("Signature:", signature))
  print(length(subset_data[[signature]]))
  print(dim(subset_data))
  
  cox_model <- coxph(Surv(`DFS Time (Days)`, `DFS Event` == "Event occurred") ~ subset_data[[signature]], 
                     data = subset_data)
  
  summary_cox <- summary(cox_model)
  HR <- summary_cox$coefficients[1, "exp(coef)"]
  Lower_CI <- summary_cox$conf.int[1, "lower .95"]
  Upper_CI <- summary_cox$conf.int[1, "upper .95"]
  coef <- summary_cox$coefficients[1, "coef"]
  se <- summary_cox$coefficients[1, "se(coef)"]
  p_value <- summary_cox$coefficients[1, "Pr(>|z|)"]
  
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
  ci_column = 4,
  ref_line = 0, 
  arrow_lab = c("Lower risk of progression", "Higher risk of progression")
)
p

###### Figure Supplemental 2h - T cell Inflamed GEP by SIC ###### -----
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

###### Figure Supplemental 2i - T cell Inflamed GEP by SE ###### -----
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


###### Figure 3B - Dot plot correlation SE abundance and CIBERSORTx cell states ###### -----
# Load data from supp tables
cibersortx_processed <- read_excel('.../Supplementary_Tables_flow_updated.xlsx', 
                                   sheet = 7)

colnames(cibersortx_processed) <- cibersortx_processed[2, ]
cibersortx_processed <- cibersortx_processed[-c(1,2), ]
cibersortx_processed <- as.data.frame(cibersortx_processed)

cibersortx_processed <- cibersortx_processed %>% 
  mutate(across(-c(1,2), as.numeric))

# Filter for Pre-treatment samples only
cibersortx_processed_pre <- cibersortx_processed %>% 
  filter(`Time Point` == "pre-treatment") %>% 
  select(-`Time Point`)

ciber_long <- cibersortx_processed_pre %>%
  pivot_longer(
    cols      = -Patient,
    names_to  = "cell_type",
    values_to = "fraction"
  )

# Add metadata 
SIC_and_SE_clean <- SIC_and_SE_assignments %>%
  dplyr::rename(
    Patient      = Patient,
    TreatmentArm = `Treatment Arm`,
    SarcomaImmuneClass       = `Sarcoma Immune Class`,
    SarcomaEcotypeAssignment = `Sarcoma Ecotype Assignment`,
    TumorGrade   = `Tumor Grade`,
    DFS_event    = `DFS Event`,
    DFS_time     = `DFS Time (Days)`,
    PercentNecrosis = `Percent Necrosis`,
    SRC_code     = `SRC code`
  )

ciber_pre_annot <- ciber_long %>%
  left_join(
    SIC_and_SE_clean %>%
      select(Patient, TreatmentArm, SarcomaImmuneClass,
             SarcomaEcotypeAssignment, TumorGrade, Histology),
    by = "Patient"
  ) %>%
  mutate(
    SIC = ifelse(
      is.na(SarcomaImmuneClass) | SarcomaImmuneClass == "",
      "SIC NA", SarcomaImmuneClass
    ),
    SE_assign = ifelse(
      is.na(SarcomaEcotypeAssignment) | SarcomaEcotypeAssignment == "NA",
      "SE mixed", SarcomaEcotypeAssignment
    )
  )

ciber_mat_wide <- ciber_pre_annot %>%
  dplyr::select(Patient, cell_type, fraction) %>%
  group_by(Patient, cell_type) %>%        
  summarize(fraction = mean(fraction), .groups = "drop") %>%
  pivot_wider(names_from = cell_type, values_from = fraction)

ciber_meta <- ciber_pre_annot %>%
  distinct(Patient, SIC, SE_assign, Histology, TumorGrade, TreatmentArm)

ciber_mat_wide <- ciber_mat_wide %>%
  dplyr::slice(match(ciber_meta$Patient, Patient))

ciber_features <- ciber_mat_wide %>%
  dplyr::select(-Patient) %>%
  as.data.frame()

rownames(ciber_features) <- ciber_mat_wide$Patient

# Handle zero fractions before doing CLR using Bayesian-multiplicative zero replacement
ciber_repl <- zCompositions::cmultRepl(
  ciber_features,
  label = 0,
  output = "prop",
  z.delete  = FALSE,
  method = "CZM"  
)

# Apply CLR normalization 
ciber_features_CLR <- compositions::clr(as.matrix(ciber_repl)) 

# Correlation between pre-treatment SE abundance and CIBERSROTx cell state abundance 
ecotype_pre_wide <- ecotype_abundance %>%
  filter(`Time Point` == "Pre-treatment") %>%
  dplyr::select(Patient, `Sarcoma Ecotype`, `Sarcoma Ecotype Abundance`) %>%
  pivot_wider(
    names_from  = `Sarcoma Ecotype`,
    values_from = `Sarcoma Ecotype Abundance`
  )

ciber_features_CLR_df <- as.data.frame(ciber_features_CLR) %>%
  mutate(Patient = rownames(.))

merged_corr <- ciber_features_CLR_df %>%
  inner_join(ecotype_pre_wide, by = "Patient")

se_cols   <- c("SE1", "SE2", "SE3")
cell_cols <- setdiff(colnames(merged_corr), c("Patient", se_cols))

# All pairwise (cell type, SE) correlations
cor_long <- expand_grid(
  CellType = cell_cols,
  Ecotype  = se_cols
) %>%
  rowwise() %>%
  mutate(
    test = list(
      suppressWarnings(
        cor.test(
          x = merged_corr[[CellType]],
          y = merged_corr[[Ecotype]],
          method = "spearman",
          exact  = FALSE
        )
      )
    ),
    rho = as.numeric(test$estimate),
    p   = test$p.value
  ) %>%
  ungroup() %>%
  dplyr::select(-test) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH")
  ) %>%
  arrange(p_adj)

#Build Dot plot 
cor_long_annot <- cor_long %>%
  mutate(
    cell_group = case_when(
      grepl("CD4", CellType) ~ "CD4+ T-cells",
      grepl("CD8", CellType) ~ "CD8+ T-cells",
      grepl("NK cells", CellType) ~ "NK cells",
      grepl("DCs", CellType) ~ "Dendritic cells",
      grepl("PMNs", CellType) ~ "PMNs",
      grepl("Monocytes|M2 |Macro ", CellType) ~ "Monocytes & Macrophages",
      grepl("Mastocytes", CellType) ~ "Mastocytes",
      grepl("Pericytes", CellType) ~ "Pericytes",
      grepl("CAF", CellType) ~ "CAFs",
      grepl("Endothelial", CellType) ~ "Endothelial cells",
      grepl("B Cells & Plasma Cells", CellType) ~ "B Cells & Plasma Cells",
      grepl("STS", CellType) ~ "Sarcoma cells",
      TRUE ~ "Other"
    )
  )

sig_pairs <- cor_long_annot %>%
  filter(!is.na(rho)) %>%
  filter(p_adj < 0.05, abs(rho) >= 0.3)

sig_pairs <- sig_pairs %>%
  mutate(
    Ecotype = factor(Ecotype, levels = c("SE1", "SE2", "SE3"))
  )

cell_order <- sig_pairs %>%
  group_by(CellType) %>%
  summarise(min_padj = min(p_adj, na.rm = TRUE)) %>%
  arrange(min_padj) %>%
  pull(CellType)

sig_pairs <- sig_pairs %>%
  mutate(
    CellType  = factor(CellType, levels = rev(cell_order)),  
    cell_group = factor(
      cell_group,
      levels = c("CD8+ T-cells", "CD4+ T-cells", "NK cells",
                 "B Cells & Plasma Cells", "Monocytes & Macrophages",
                 "Dendritic cells", "PMNs", "Mastocytes",
                 "Sarcoma cells", "CAFs", "Endothelial cells",
                 "Pericytes", "Other")
    )
  )

# 4) Add -log10(FDR) with a cap 
sig_pairs <- sig_pairs %>%
  mutate(
    neg_log10_FDR     = -log10(p_adj),
    neg_log10_FDR_cap = pmin(neg_log10_FDR, 10)  # cap at 10
  )

#vertical
p_corr_dot_grouped_vertical <- ggplot(sig_pairs, aes(x = Ecotype, y = CellType)) +
  geom_point( aes(size = neg_log10_FDR_cap, color = rho), alpha = 1 ) + 
  scale_size_continuous( name = expression(-log[10]("FDR")), range = c(2, 8) ) + 
  scale_color_gradient2( name = "Spearman \u03c1", low = "#2166ac",  mid = "white", high = "#b2182b", midpoint = 0 ) +
  facet_grid( rows = vars(cell_group), scales = "free_y", space = "free_y" ) + 
  labs( x = "Sarcoma ecotype", y = "CIBERSORTx cell state", 
        title = "Associations between pre-treatment SE abundance and deconvoluted cell states", 
        subtitle = "All significant correlations (FDR < 0.05, |ρ| ≥ 0.3)\nCLR-transformed CIBERSORTx fractions" ) + 
  theme_minimal() + 
  theme(
    axis.text.x  = element_text(angle = 45, 
                                size = 15,
                                color = "black",
                                hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15,
                               color = "black",
                               hjust = 1, vjust = 1),
    plot.title   = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold")
  )

ggsave(p_corr_dot_grouped_vertical, filename = '.../Figure_3b.pdf', 
       width = 22, height = 32, units = "cm")




###### Figure 3C - Dot plot correlation SE abundance and CIBERSORTx cell states ###### -----
# Plotting only cell states with FDR < 0.025 for KW test and pairwise comparisons with FDR < 0.025 for Dunn Test
tumor_long <- as.data.frame(ciber_features_CLR) %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "CLR")

data_for_analysis <- tumor_long %>%
  inner_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Ecotype Assignment")],
             by = c("Sample" = "Patient"))

# Add Max SE assignment
data_for_analysis <- data_for_analysis %>% 
  left_join(pre_ecotype_calls_wide[, c("patient", "Assigned_Max_Ecotype")], 
            by = c("Sample" = "patient"))

#Remove unassigned SE patients 
data_for_analysis_SE_assignment <- data_for_analysis %>% 
  filter(`Sarcoma Ecotype Assignment` != "NA")

# Kruskal-Wallis test for pre-treatment SE assignments 
kw_results_pre_assign <- data_for_analysis_SE_assignment %>%
  group_by(CellType) %>%
  kruskal_test(CLR ~ `Sarcoma Ecotype Assignment`) %>%
  adjust_pvalue(method = "BH") %>%
  dplyr::rename(p_value = p, adj_p_value = p.adj)

# keep cell types with significant global differences - pre-assignment 
# FDR at 0.025
sig_celltypes <- kw_results_pre_assign %>%
  filter(adj_p_value < 0.025) %>%
  pull(CellType) %>%
  unique()

# Post-hoc pairwise Dunn's test / FDR at 0.025  
data_for_analysis_SE_assignment$ecotype_assignment <- data_for_analysis_SE_assignment$`Sarcoma Ecotype Assignment`
pairwise_np_assign <- data_for_analysis_SE_assignment %>%
  group_by(CellType) %>%
  dunn_test(CLR ~ ecotype_assignment, p.adjust.method = "BH") %>%
  mutate(signif = ifelse(p.adj < 0.025, "yes", "no"))

se_levels <- c("SE1", "SE2", "SE3")
contrast_levels <- combn(se_levels, 2, FUN = \(x) paste0(x[1], "-", x[2]))

# 2) Median CLR per CellType x SE (for effect sizes)
med_by_SEs <- data_for_analysis_SE_assignment %>%
  filter(CellType %in% sig_celltypes) %>%
  mutate(ecotype_assignment = factor(ecotype_assignment, levels = se_levels)) %>%
  group_by(CellType, ecotype_assignment) %>%
  summarise(med_CLR = median(CLR, na.rm = TRUE), .groups = "drop")

pair_sig <- pairwise_np_assign %>%
  ungroup() %>%
  filter(CellType %in% sig_celltypes, p.adj < 0.025) %>%      
  mutate(
    g1 = as.character(group1),
    g2 = as.character(group2),
    ref  = if_else(match(g1, se_levels) < match(g2, se_levels), g1, g2),
    comp = if_else(ref == g1, g2, g1),
    contrast = factor(paste0(ref, "-", comp), levels = contrast_levels),
    FDR = p.adj
  ) %>%
  select(CellType, ref, comp, contrast, FDR)

# 4) Effect size: delta_med_CLR = median(comp) - median(ref)
plot_df <- pair_sig %>%
  left_join(med_by_SEs %>% dplyr::rename(ref  = ecotype_assignment, med_ref  = med_CLR),
            by = c("CellType","ref")) %>%
  left_join(med_by_SEs %>% dplyr::rename(comp = ecotype_assignment, med_comp = med_CLR),
            by = c("CellType","comp")) %>%
  mutate(
    delta_med_CLR = med_comp - med_ref,
    neglog10_FDR  = -log10(pmax(FDR, 1e-300))
  )

cap <- quantile(plot_df$neglog10_FDR, 0.99, na.rm = TRUE)
plot_df <- plot_df %>% mutate(neglog10_FDR_cap = pmin(neglog10_FDR, cap))

max_abs <- max(abs(plot_df$delta_med_CLR), na.rm = TRUE)

med_wide <- med_by_SEs %>%
  mutate(ecotype_assignment = factor(ecotype_assignment, levels = se_levels)) %>%
  pivot_wider(names_from = ecotype_assignment, values_from = med_CLR) %>%
  arrange(CellType)

Delta_mat <- sapply(strsplit(contrast_levels, "-"), function(x) {
  ref  <- x[1]
  comp <- x[2]
  med_wide[[comp]] - med_wide[[ref]]
})

Delta_mat <- as.matrix(Delta_mat)
rownames(Delta_mat) <- med_wide$CellType
colnames(Delta_mat) <- contrast_levels

# Clustering
d_row <- as.dist(1 - cor(t(Delta_mat), use = "pairwise.complete.obs"))
hc_row <- hclust(d_row, method = "average")

d_col <- as.dist(1 - cor(Delta_mat, use = "pairwise.complete.obs"))
hc_col <- hclust(d_col, method = "average")

row_order <- rownames(Delta_mat)[hc_row$order]
col_order <- colnames(Delta_mat)[hc_col$order]

plot_df2 <- plot_df %>%
  mutate(
    CellType  = factor(as.character(CellType), levels = rev(row_order)),
    contrast  = factor(as.character(contrast), levels = col_order)
  )

max_abs <- max(abs(plot_df2$delta_med_CLR), na.rm = TRUE)

# Plot 
bubble <- ggplot(plot_df2, aes(x = contrast, y = CellType)) +
  geom_point(aes(size = neglog10_FDR_cap, fill = delta_med_CLR),
             shape = 21, colour = "black", stroke = 0.25, alpha = 1) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-max_abs, max_abs),
                       name = "Δ median CLR\n(comp - ref)") +
  scale_size_area(max_size = 10,
                  name = expression(-log[10]("Dunn FDR"))) +
  scale_y_discrete(position = "right") +      
  scale_x_discrete(position = "bottom") +     
  labs(x = "SE pair (ref-comp)", y = NULL
  ) +
  theme_pubr(border = TRUE) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.title.y      = element_blank(),
    axis.text.x.top   = element_blank(),
    axis.ticks.x.top  = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey92"), 
    plot.margin       = margin(5, 5, 5, 5)
  )

n_row <- length(row_order)
n_col <- length(col_order)

dd_row <- ggdendro::dendro_data(hc_row, type = "rectangle") 
row_dendro <- ggplot() +
  geom_segment(data = dd_row$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +                                             
  scale_x_continuous(limits = c(0.5, n_row + 0.5), expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

dd_col <- ggdendro::dendro_data(hc_col, type = "rectangle") 
col_dendro <- ggplot() +
  geom_segment(data = dd_col$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_x_continuous(limits = c(0.5, n_col + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 0, 5))

final_panel <- (plot_spacer() | col_dendro) / (row_dendro | bubble) +
  plot_layout(widths = c(0.15, 1), 
              heights = c(0.15, 1), 
              guides = "collect") 
final_panel

ggsave(final_panel, 
       filename = '.../Figure_3C_dotplot_SE_pre_treatment_Dunn_KW_sig_0025.pdf',       
       units = "cm", 
       width = 30, height = 27.5)

###### Figure 3D - Dot plot correlation SE abundance and CIBERSORTx cell states ###### -----
sic_levels <- c("A","B","C","D","E")
contrast_levels <- combn(sic_levels, 2, FUN = \(x) paste0(x[1], "-", x[2]))

sig_celltypes <- kw_results %>%
  filter(adj_p_value < 0.025) %>%
  pull(CellType) %>%
  unique()

med_by_sic <- data_for_analysis %>%
  filter(CellType %in% sig_celltypes) %>%
  mutate(Sarcoma_Immune_Class = factor(Sarcoma_Immune_Class, levels = sic_levels)) %>%
  group_by(CellType, Sarcoma_Immune_Class) %>%
  summarise(med_CLR = median(CLR, na.rm = TRUE), .groups = "drop")

pair_sig <- pairwise_np %>%
  ungroup() %>%
  filter(CellType %in% sig_celltypes, p.adj < 0.025) %>%        
  mutate(
    g1 = as.character(group1),
    g2 = as.character(group2),
    ref  = if_else(match(g1, sic_levels) < match(g2, sic_levels), g1, g2),
    comp = if_else(ref == g1, g2, g1),
    contrast = factor(paste0(ref, "-", comp), levels = contrast_levels),
    FDR = p.adj
  ) %>%
  select(CellType, ref, comp, contrast, FDR)

plot_df <- pair_sig %>%
  left_join(med_by_sic %>% dplyr::rename(ref  = Sarcoma_Immune_Class, med_ref  = med_CLR),
            by = c("CellType","ref")) %>%
  left_join(med_by_sic %>% dplyr::rename(comp = Sarcoma_Immune_Class, med_comp = med_CLR),
            by = c("CellType","comp")) %>%
  mutate(
    delta_med_CLR = med_comp - med_ref,
    neglog10_FDR  = -log10(pmax(FDR, 1e-300))
  )

cap <- quantile(plot_df$neglog10_FDR, 0.99, na.rm = TRUE)
plot_df <- plot_df %>% mutate(neglog10_FDR_cap = pmin(neglog10_FDR, cap))

max_abs <- max(abs(plot_df$delta_med_CLR), na.rm = TRUE)

med_wide <- med_by_sic %>%
  mutate(Sarcoma_Immune_Class = factor(Sarcoma_Immune_Class, levels = sic_levels)) %>%
  pivot_wider(names_from = Sarcoma_Immune_Class, values_from = med_CLR) %>%
  arrange(CellType)

Delta_mat <- sapply(strsplit(contrast_levels, "-"), function(x) {
  ref  <- x[1]
  comp <- x[2]
  med_wide[[comp]] - med_wide[[ref]]
})

Delta_mat <- as.matrix(Delta_mat)
rownames(Delta_mat) <- med_wide$CellType
colnames(Delta_mat) <- contrast_levels

d_row <- as.dist(1 - cor(t(Delta_mat), use = "pairwise.complete.obs"))
hc_row <- hclust(d_row, method = "average")

d_col <- as.dist(1 - cor(Delta_mat, use = "pairwise.complete.obs"))
hc_col <- hclust(d_col, method = "average")

row_order <- rownames(Delta_mat)[hc_row$order]
col_order <- colnames(Delta_mat)[hc_col$order]

plot_df2 <- plot_df %>%
  mutate(
    CellType  = factor(as.character(CellType), levels = rev(row_order)),
    contrast  = factor(as.character(contrast), levels = col_order)
  )

max_abs <- max(abs(plot_df2$delta_med_CLR), na.rm = TRUE)

bubble <- ggplot(plot_df2, aes(x = contrast, y = CellType)) +
  geom_point(aes(size = neglog10_FDR_cap, fill = delta_med_CLR),
             shape = 21, colour = "black", stroke = 0.25, alpha = 1) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-max_abs, max_abs),
                       name = "Δ median CLR\n(comp - ref)") +
  scale_size_area(max_size = 10,
                  name = expression(-log[10]("Dunn FDR"))) +
  scale_y_discrete(position = "right") +    
  scale_x_discrete(position = "bottom") +   
  labs(x = "SIC pair (ref-comp)", y = NULL
      
  ) +
  theme_pubr(border = TRUE) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.text.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.title.y      = element_blank(),
    axis.text.x.top   = element_blank(),
    axis.ticks.x.top  = element_blank(), 
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey92"), 
    plot.margin       = margin(5, 5, 5, 5)
  )

n_row <- length(row_order)
n_col <- length(col_order)

dd_row <- ggdendro::dendro_data(hc_row, type = "rectangle") 
row_dendro <- ggplot() +
  geom_segment(data = dd_row$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +                                             
  scale_x_continuous(limits = c(0.5, n_row + 0.5), expand = c(0, 0)) +
  scale_y_reverse(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

dd_col <- ggdendro::dendro_data(hc_col, type = "rectangle") 
col_dendro <- ggplot() +
  geom_segment(data = dd_col$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_x_continuous(limits = c(0.5, n_col + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 0, 5))

final_panel <- (plot_spacer() | col_dendro) / (row_dendro | bubble) +
  plot_layout(widths = c(0.15, 1), 
              heights = c(0.15, 1), 
              guides = "collect") 
final_panel

ggsave(final_panel, filename = '.../Figure_3D_dotplot_SIC_pre_treatment_Dunn_KW_sig_0025.pdf', 
       units = "cm", 
       width = 40, height = 25)

###### Figure S3B - CIBERSORTx cell state abundance vs Flow Cytometry Abundance ###### ------
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






## Post-Treatment Samples Analysis - Helpers ## ----
norm_timepoint <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[ -]", "_") %>%
    str_replace_all("pretreatment|pre_treatment|pre", "pre_treatment") %>%
    str_replace_all("surgical_resection|post_treatment|post", "post_treatment")
}

make_sample_id <- function(patient, timepoint) {
  paste0(patient, "__", timepoint)
}

paired_only <- function(meta_df, patient_col = "Patient", time_col = "TimePoint",
                        pre_label = "pre_treatment", post_label = "post_treatment") {
  keep <- meta_df %>%
    distinct(.data[[patient_col]], .data[[time_col]]) %>%
    dplyr::count(.data[[patient_col]]) %>%
    dplyr::filter(n >= 2) %>%
    pull(.data[[patient_col]])
  
  meta_df %>% dplyr::filter(.data[[patient_col]] %in% keep) %>%
    dplyr::filter(.data[[time_col]] %in% c(pre_label, post_label))
}

# Build pre + post CIBERSORTx CLR matrices
SIC_and_SE_clean <- SIC_and_SE_assignments %>%
  dplyr::rename(
    Patient      = Patient,
    TreatmentArm = `Treatment Arm`,
    SIC_pre      = `Sarcoma Immune Class`,
    SE_pre       = `Sarcoma Ecotype Assignment`,
    TumorGrade   = `Tumor Grade`,
    Histology    = Histology
  ) %>%
  mutate(
    SIC_pre = ifelse(is.na(SIC_pre) | SIC_pre == "", "SIC NA", SIC_pre),
    SE_pre  = ifelse(is.na(SE_pre)  | SE_pre  == "NA", "SE mixed", SE_pre)
  )

ciber_long_all <- cibersortx_processed %>%
  mutate(TimePoint = norm_timepoint(`Time Point`)) %>%
  dplyr::filter(TimePoint %in% c("pre_treatment", "post_treatment")) %>%
  pivot_longer(
    cols      = -c(Patient, `Time Point`, TimePoint),
    names_to  = "cell_type",
    values_to = "fraction"
  ) %>%
  left_join(SIC_and_SE_clean %>% dplyr::select(Patient, TreatmentArm, SIC_pre, SE_pre, TumorGrade, Histology),
            by = "Patient") %>%
  mutate(
    SampleID = make_sample_id(Patient, TimePoint)
  )

ciber_long_all = ciber_long_all %>% 
  left_join(pre_ecotype_calls_wide[, c("patient", "Assigned_Max_Ecotype")], 
            by = c("Patient" = "patient"))

# Build wide fractions per SampleID
ciber_wide_all <- ciber_long_all %>%
  group_by(SampleID, Patient, TimePoint, TreatmentArm, 
           SIC_pre, SE_pre, Assigned_Max_Ecotype,
           TumorGrade, Histology, cell_type) %>%
  summarize(fraction = mean(fraction), .groups = "drop") %>%
  pivot_wider(names_from = cell_type, values_from = fraction, values_fill = 0)

ciber_frac_wide_all <- ciber_wide_all %>%
  column_to_rownames("SampleID") %>%
  dplyr::select(-Patient, -TimePoint, -TreatmentArm, -SIC_pre, -SE_pre, 
                -Assigned_Max_Ecotype, 
                -TumorGrade, -Histology)

ciber_meta_all <- ciber_wide_all %>%
  distinct(SampleID, Patient, TimePoint, TreatmentArm, Assigned_Max_Ecotype, 
           SIC_pre, SE_pre, TumorGrade, Histology) %>%
  arrange(Patient, TimePoint)

# CLR matrix aligned to meta
ciber_wide_all <- ciber_wide_all %>% dplyr::slice(match(ciber_meta_all$SampleID, SampleID))

ciber_wide_all <- ciber_wide_all %>% 
  column_to_rownames(var = "SampleID")

ciber_wide_all = ciber_wide_all[, -c(1:8)]

repl <- zCompositions::cmultRepl(
  ciber_wide_all,
  label = 0,
  output = "prop",
  z.delete = FALSE,
  method = "CZM"
)
clr_mat <- compositions::clr(as.matrix(repl))

ciber_clr_all = clr_mat

# Keep paired patients for paired analyses
ciber_meta_paired <- paired_only(ciber_meta_all, time_col = "TimePoint") %>%
  arrange(Patient, TimePoint)

ciber_clr_paired <- ciber_clr_all[ciber_meta_paired$SampleID, , drop = FALSE]


###### Figure 4D - Dot Plot Post- vs Pre-Treatment abundance change by treatment arm ##### ----- 
stopifnot(all(rownames(ciber_clr_paired) %in% ciber_meta_paired$SampleID))

meta <- ciber_meta_paired %>%
  dplyr::mutate(
    SampleID = as.character(SampleID),
    Patient  = factor(Patient),
    TimePoint = factor(TimePoint, levels = c("pre_treatment", "post_treatment")),
    TreatmentArm = factor(TreatmentArm)
  ) %>%
  dplyr::arrange(match(SampleID, rownames(ciber_clr_paired)))

Y <- ciber_clr_paired[meta$SampleID, , drop = FALSE]

Y <- as.matrix(unclass(Y))

Y_df <- as.data.frame(Y) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "CellState", values_to = "CLR") %>%
  dplyr::left_join(meta %>% dplyr::select(SampleID, Patient, TimePoint, 
                                          TreatmentArm, SIC_pre, SE_pre, 
                                          Assigned_Max_Ecotype),
                   by = "SampleID") %>%
  dplyr::select(-SampleID) %>%
  tidyr::pivot_wider(names_from = TimePoint, values_from = CLR) %>%
  dplyr::mutate(delta = post_treatment - pre_treatment)

Wilcoxon_paired_pre_post_global <- Y_df %>%
  group_by(CellState, TreatmentArm) %>%
  summarize(
    p_within   = wilcox.test(post_treatment, pre_treatment, paired=TRUE)$p.value, 
    delta      = mean(delta),
    .groups    = "drop"
  ) %>%
  mutate(adj_p_within = p.adjust(p_within, "fdr")) %>% #adjust globally
  dplyr::arrange((adj_p_within)) 

#Cells to plot 
wilc_sig <- Wilcoxon_paired_pre_post_global %>%
  dplyr::mutate(
    neglog10FDR = -log10(pmax(adj_p_within, 1e-300)),
    TreatmentArm = factor(TreatmentArm)
  ) %>%
  dplyr::filter(adj_p_within < 0.05) %>%
  dplyr::arrange(TreatmentArm, dplyr::desc(neglog10FDR))

# Order cell states by strongest significance anywhere
state_order <- Wilcoxon_paired_pre_post_global %>%
  dplyr::group_by(CellState) %>%
  dplyr::summarise(bestFDR = min(adj_p_within, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(bestFDR) %>%
  dplyr::pull(CellState)

wilc_sig$CellState <- factor(wilc_sig$CellState, levels = rev(state_order))

wilc_sig_dotplot <- ggplot(wilc_sig, aes(y = TreatmentArm, x = CellState)) +
  geom_point(aes(size = neglog10FDR, color = delta), alpha = 1) +
  scale_color_distiller(palette = "RdBu", direction = -1, name = "Mean ΔCLR\n(post-pre)") +
  scale_size_continuous(name = "-log10(FDR)") +
  theme_minimal() +
  labs(x = NULL, y = NULL, title = "Pre→post change within each arm (paired Wilcoxon)") + 
  theme(axis.text.x = element_text(size = 15, angle = 70, hjust = 1), 
        axis.text.y = element_text(size = 16))

ggsave(wilc_sig_dotplot, 
       width = 27.5, height = 14, 
       filename = '.../Figure_4d_dot_plot_wilcoxon_pre_vs_post_global.pdf',
       units = "cm")

###### Figure 4E Pre- vs Post- cell state abundance by treatment arm and SE ecotype Max ###### ------
Wilcoxon_paired_pre_post_SE_Max <- Y_df %>%
  group_by(CellState, TreatmentArm, Assigned_Max_Ecotype) %>%
  summarize(
    p_within   = wilcox.test(post_treatment, pre_treatment, paired=TRUE)$p.value, 
    delta      = mean(delta),
    .groups    = "drop"
  ) %>%
  mutate(adj_p_within = p.adjust(p_within, "fdr")) %>% 
  dplyr::arrange((adj_p_within)) 

wilcox_tbl  <- Wilcoxon_paired_pre_post_SE_Max
stratum_col <- "Assigned_Max_Ecotype"   
fdr_cutoff  <- 0.05

n_tbl <- Y_df %>%
  distinct(Patient, TreatmentArm, .data[[stratum_col]]) %>%
  dplyr::count(TreatmentArm, .data[[stratum_col]], name = "n_patients")

df <- wilcox_tbl %>%
  left_join(n_tbl, by = c("TreatmentArm", stratum_col))

df <- df %>%
  mutate(
    Stratum = .data[[stratum_col]],
    Stratum = factor(Stratum, levels = c("SE1","SE2","SE3")),
    neglog10FDR = -log10(pmax(adj_p_within, 1e-300)),
    neglog10FDR = ifelse(is.finite(neglog10FDR), neglog10FDR, NA_real_),
    sig_fdr = adj_p_within < fdr_cutoff
  )

cap <- quantile(df$neglog10FDR, 0.99, na.rm = TRUE)
df <- df %>% mutate(neglog10FDR_cap = pmin(neglog10FDR, cap))

keep_states <- df %>%
  group_by(CellState) %>%
  summarise(any_hit = any(sig_fdr, na.rm = TRUE), .groups = "drop") %>%
  filter(any_hit) %>%
  pull(CellState)

df <- df %>% filter(CellState %in% keep_states)

state_order <- df %>%
  group_by(CellState) %>%
  summarise(order_val = mean(delta, na.rm = TRUE), .groups = "drop") %>%
  arrange(order_val) %>%
  pull(CellState)

df <- df %>%
  mutate(CellState = factor(CellState, levels = state_order))

p <- ggplot(df, aes(x = Stratum, y = CellState)) +
  geom_point(
    data = df %>% filter(!sig_fdr),
    aes(size = neglog10FDR_cap, fill = delta),
    alpha = 0.5, stroke = 0,  shape = 21
  ) +
  geom_point(
    data = df %>% filter(sig_fdr),
    aes(size = neglog10FDR_cap, fill = delta),
    alpha = 1, stroke = 1, shape = 21, colour = "black"
  ) +
  facet_wrap(~ TreatmentArm, nrow = 1) +
  scale_size_continuous(
    name = expression(-log[10]("BH-FDR")),
    range = c(1.2, 8)
  ) +
  scale_fill_distiller(
    name = expression(Delta*CLR),
    palette = "RdBu",
    direction = -1
  ) +
  labs(x = stratum_col, y = NULL) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(
    strip.text = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey93")
  )

p

ggsave(p, unit = "cm",  width = 25, height = 10, 
       filename = '.../Figure_4E_ SE_max_signed_rank_control_exp_bubblePlot.pdf')



###### Figure 4F SICE vs SE1 cell state abundance change in Control vs Experimental ###### ------
between_arm_SE1_vs_SICE <- Y_df %>%
  mutate(SE_pre_simple = ifelse(SE_pre == "SE1", "SE1", "other"), 
         SIC_pre_simple = ifelse(SIC_pre == "E", "E", "other"), 
         SIC_SE_comb = case_when(
           SE_pre == "SE1" & SIC_pre != "E" ~ "SE1", 
           SE_pre != "SE1" & SIC_pre == "E" ~ "E", 
           SE_pre == "SE1" & SIC_pre == "E" ~ "Remove"
         )) %>%
  dplyr::filter(SIC_SE_comb != "Remove")  %>%
  group_by(CellState, TreatmentArm) %>%
  summarize(
    p_between  = wilcox.test(delta ~ SIC_SE_comb)$p.value, 
    mean_delta_SICE = mean(delta[SIC_SE_comb == "E"], na.rm = T),
    mean_delta_SE1     = mean(delta[SIC_SE_comb == "SE1"], na.rm = T),
    .groups    = "drop"
  ) %>% 
  mutate(adj_p_between = p.adjust(p_between, "fdr")) %>% #adjust globally
  dplyr::arrange((p_between)) 

sig_states <- between_arm_SE1_vs_SICE %>%
  filter(p_between < 0.05) %>%                
  pull(CellState) %>%
  unique()

df_ridge <- Y_df %>%
  mutate(
    SIC_SE_comb = case_when(
      SE_pre == "SE1" & SIC_pre != "E" ~ "SE1",
      SE_pre != "SE1" & SIC_pre == "E" ~ "E",
      SE_pre == "SE1" & SIC_pre == "E" ~ "Remove",
      TRUE ~ NA_character_
    ),
    SIC_SE_comb = factor(SIC_SE_comb, levels = c("SE1", "E"))
  ) %>%
  filter(SIC_SE_comb %in% c("SE1","E")) %>%
  filter(CellState %in% sig_states)

state_order <- df_ridge %>%
  group_by(CellState) %>%
  summarise(
    ord = median(delta[SIC_SE_comb == "E"],   na.rm = TRUE) -
      median(delta[SIC_SE_comb == "SE1"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ord) %>%
  pull(CellState)

df_ridge <- df_ridge %>%
  mutate(CellState = factor(CellState, levels = state_order))

p_ridge <- ggplot(df_ridge, aes(x = delta, y = CellState, fill = SIC_SE_comb)) +
  ggridges::geom_density_ridges(
    aes(group = interaction(CellState, SIC_SE_comb)),
    position = "identity",
    alpha = 0.55,
    scale = 1.1,
    rel_min_height = 0.01,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  facet_wrap(~ TreatmentArm, nrow = 1) +
  scale_fill_manual(
    values = c("SE1" = "grey80", "E" = "#e21f26"),
    name = "Baseline group"
  ) +
  labs(
    x = "Δ CLR (post - pre)",
    y = NULL,
    title = "Δ CLR distributions: baseline SE1 vs baseline SIC E",
    subtitle = "Ridgeline densities per cell state; faceted by Treatment Arm"
  ) +
  theme_pubr(border = TRUE, legend = "top") +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor = element_blank(), 
    strip.text = element_text(size = 18), 
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20), 
    axis.title.x = element_text(size = 20), 
    axis.text.y = element_text(size = 18)
  )

p_ridge

ggsave(p_ridge, unit = "cm",  width = 27.5, height = 15, 
       filename = '.../Figure_4F_ridgeplot_sice_se1_NOT_MAX.pdf')

###### Figure 4G - Box-Violin plot for Pre- vs Post-treatment SE abundance change by treatment arm ###### -----
ecotype_abundance <- read_excel('.../data/Supplementary_Tables.xlsx', 
                                sheet = 8)

colnames(ecotype_abundance) <- ecotype_abundance[2, ]
ecotype_abundance <- ecotype_abundance[-c(1,2), ]
ecotype_abundance <- as.data.frame(ecotype_abundance)
ecotype_abundance$`Sarcoma Ecotype Abundance` <- as.numeric(ecotype_abundance$`Sarcoma Ecotype Abundance`)

df_long <- ecotype_abundance %>%
  mutate(
    patient    = Patient,  
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
       filename = '.../Figure_4G_violin_abundance_bySE_change.pdf')

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

###### Figure Supplemental 5A and Supplemental 5B ###### ------
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

IC_genes <- IC_genes[IC_genes$patient %in% SIC_and_SE_assignments$Patient, ]

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
pdf(".../Figure_S5A_heatmap_arm_vs_timepoint_IFNG_sig_genes_SIC.pdf",
    width = w, height = h)
draw(ht)
dev.off()

# Figure Supp 5B 
sample_annotation <- IC_genes[, c("Sample", "Time_Point", "Arm", "Sarcoma Ecotype Assignment")]
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
pdf(".../Figure_S5B_heatmap_arm_vs_timepoint_IFNG_sig_genes_SEs.pdf",
    width = w, height = h)
draw(ht)
dev.off()

###### Figure Supplemental 5C and Figure Supplemental 5D ###### ------
tcell_inflamed_data <- tcell_inflamed_data %>%
  mutate(Time_Point = factor(Time_Point, levels = c("Pre-treatment",
                                                    "Post-treatment")))
paired_wsr <- function(data, value, strata = NULL) {
  val_quo   <- rlang::enquo(value)
  strat_quo <- rlang::enquo(strata)
  
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

# Figure S5C 
SIC_data <- tcell_inflamed_data %>% 
  filter(!is.na(`Sarcoma Immune Class`)) %>% 
  filter(`Sarcoma Immune Class` != "NA")

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

# Figure S5D 
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
  facet_grid(`Treatment Arm` ~ ecotype_assignment) +
  ggpubr::stat_pvalue_manual(
    Eco_pvals,
    label = "label", y.position = "y.position",
    xmin = "group1", xmax = "group2", tip.length = 0.01
  ) +
  theme_pubr(border = TRUE) +
  theme(strip.text = element_text(size = 14))

###### Figure Supplemental 4A ###### -----
# Filter for Pre-treatment samples only
cibersortx_processed_pre <- cibersortx_processed %>% 
  filter(`Time Point` == "pre-treatment") %>% 
  select(-`Time Point`)

ciber_long <- cibersortx_processed_pre %>%
  pivot_longer(
    cols      = -Patient,
    names_to  = "cell_type",
    values_to = "fraction"
  )

ciber_long <- ciber_long %>% 
  left_join(SIC_and_SE_assignments)

ciber_long$cell_type <- factor(ciber_long$cell_type)
ciber_long$Patient <- factor(ciber_long$Patient)

heatmap_matrix <- ciber_long %>%
  #mutate(Sample_TP = paste0(Sample, "_", Time_Point)) %>%
  select(Patient, cell_type, fraction) %>%
  pivot_wider(
    names_from = Patient, 
    values_from = fraction
  ) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

heatmap_matrix <- heatmap_matrix[!(grepl("total_", rownames(heatmap_matrix))),]

# arcsin + z-scoring
# 1. perform arcsin normalization on the fractions 
data_asin <- asin(sqrt(heatmap_matrix))
# 2. Perform row-based z-score scaling
scaled_data <- t(scale(t(data_asin)))

row_annotation <- data.frame(cells = rownames(scaled_data))
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

row_annotation <- column_to_rownames(row_annotation, var = "cells")
row_annotation$cell_types <- factor(
  row_annotation$cell_types,
  levels = c("CD8+ T-cells", "CD4+ T-cells", "NK cells", 
             "B Cells & Plasma Cells", "Monocytes & Macrophages", 
             "Dendritic cells", "PMNs", "Mastocytes", 
             "Sarcoma cells", "CAFs", "Endothelial cells", 
             "Pericytes"))

column_annotation <- data.frame(samples = colnames(scaled_data))

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% ciber_long[ciber_long$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Ecotype Assignment", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation$ecotype_assignment <- column_annotation$`Sarcoma Ecotype Assignment`

column_annotation <- column_to_rownames(column_annotation, var = "samples")
column_annotation <- column_annotation[colnames(scaled_data), ]
column_annotation <- column_annotation %>% 
  relocate(`ecotype_assignment`)

column_annotation$ecotype_assignment <- factor(
  column_annotation$ecotype_assignment,
  levels = c("SE1", "SE2", "SE3")
)

column_annotation <- column_annotation[!(is.na(column_annotation$ecotype_assignment)), ]
scaled_data <- scaled_data[, colnames(scaled_data) %in% rownames(column_annotation)]
column_annotation$`Sarcoma Ecotype Assignment` <- NULL

col_annotation <- HeatmapAnnotation(
  df = column_annotation,
  col = list(`Tumor Grade` = c("Grade 2" = "#a6dba0", "Grade 3" = "#984ea3"), 
             Arm = c("Experimental" = "#0093AF", "Control" = "#9E1C39"), 
             `ecotype_assignment` = c("SE1" = "#F06180", "SE2" = "#0F9ABE", "SE3" = "#60C1A5")))

col_custom <- colorRamp2(c(-4, 0, 4), c("#2166ac", "#ffffff", "#b2182b"))

ht = Heatmap(scaled_data, 
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
pdf('.../Figure_S4B_absolute_fractions_EcoTyper.pdf', width = w, height = h)
draw(ht)
dev.off()




###### Figure Supplemental 4B ###### -----
# Filter for Pre-treatment samples only
cibersortx_processed_pre <- cibersortx_processed %>% 
  filter(`Time Point` == "pre-treatment") %>% 
  select(-`Time Point`)

ciber_long <- cibersortx_processed_pre %>%
  pivot_longer(
    cols      = -Patient,
    names_to  = "cell_type",
    values_to = "fraction"
  )

ciber_long <- ciber_long %>% 
  left_join(SIC_and_SE_assignments)

ciber_long$cell_type <- factor(ciber_long$cell_type)
ciber_long$Patient <- factor(ciber_long$Patient)

heatmap_matrix <- ciber_long %>%
  #mutate(Sample_TP = paste0(Sample, "_", Time_Point)) %>%
  select(Patient, cell_type, fraction) %>%
  pivot_wider(
    names_from = Patient, 
    values_from = fraction
  ) %>%
  column_to_rownames("cell_type") %>%
  as.matrix()

heatmap_matrix <- heatmap_matrix[!(grepl("total_", rownames(heatmap_matrix))),]

# arcsin + z-scoring
# 1. perform arcsin normalization on the fractions 
data_asin <- asin(sqrt(heatmap_matrix))
# 2. Perform row-based z-score scaling
scaled_data <- t(scale(t(data_asin)))

row_annotation <- data.frame(cells = rownames(scaled_data))
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

column_annotation <- data.frame(samples = colnames(scaled_data))
#rownames(column_annotation) <- column_annotation$samples

column_annotation$Arm <- NA
column_annotation$Arm[column_annotation$samples %in% ciber_long[ciber_long$`Treatment Arm` == "Experimental", ]$Patient] <- "Experimental"
column_annotation$Arm[is.na(column_annotation$Arm)] <- "Control"

column_annotation <- column_annotation %>% 
  left_join(SIC_and_SE_assignments[, c("Patient", "Sarcoma Immune Class", 
                                       "Tumor Grade")], by = c("samples" = "Patient")) 

column_annotation$`Sarcoma Immune Class`[is.na(column_annotation$`Sarcoma Immune Class`)] <- "Unassigned"
column_annotation <- column_to_rownames(column_annotation, var = "samples")

column_annotation <- column_annotation[colnames(scaled_data), ]

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

ht = Heatmap(scaled_data, 
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
pdf('/.../Figure_S4A_scaled_absolute_fractions.pdf', width = w, height = h)
draw(ht)
dev.off()

###### Figure 4C - TCR Repetoire Analysis // TRUST-4 ###### -----
trust4_out_report <- repLoad(.path = ".../TRUST4/report/non_zero")
names(trust4_out_report$data) <- gsub(".trust4.out_report", "", names(trust4_out_report$data))
trust4_out_report$meta$Sample <- gsub(".trust4.out_report", "", trust4_out_report$meta$Sample)

# filter results - remove non coding sequences
trust4_out_report_coding <- repFilter(trust4_out_report, .method = "by.clonotype", 
                                      .query = list(CDR3.aa = exclude("partial", "out_of_frame")))
# filter results - remove IG genes
trust4_out_report_coding_Tcells <- repFilter(trust4_out_report_coding, .method = "by.clonotype", 
                                             .query = list(V.name = exclude("IG"),
                                                           D.name = exclude("IG"),
                                                           J.name = exclude("IG")), .match="substring")
# filter results - remove TCR genes
trust4_out_report_coding_Bcells <- repFilter(trust4_out_report_coding, .method = "by.clonotype", 
                                             .query = list(V.name = exclude("TR"),
                                                           D.name = exclude("TR"),
                                                           J.name = exclude("TR")), .match="substring")

# Add time point in the meta
trust4_out_report_coding_Tcells$meta$timepoint <- NA
trust4_out_report_coding_Tcells$meta$timepoint[grepl("T1", trust4_out_report_coding_Tcells$meta$Sample)] <- "T1"
trust4_out_report_coding_Tcells$meta$timepoint[grepl("T2", trust4_out_report_coding_Tcells$meta$Sample)] <- "T2"
trust4_out_report_coding_Tcells$meta$sample_simp <- str_extract(trust4_out_report_coding_Tcells$meta$Sample, "SRC\\d+")

#Load processed data 
trust4_out_report_coding_Tcells = qread('.../trust4_out_processed.qs')

patient_lookup <- setNames(names(src_lookup), src_lookup)

trust4_out_report_coding_Tcells$meta$Patient <- unname(
  patient_lookup[trust4_out_report_coding_Tcells$meta$sample_simp]
)

tcr_long <- imap_dfr(trust4_out_report_coding_Tcells$data, function(tbl, sample_id) {
  
  count_col <- intersect(c("Read.count","read_count","readCount","count","Clones"), colnames(tbl))[1]
  aa_col    <- intersect(c("CDR3.aa","CDR3_amino_acids","cdr3aa","CDR3_aa"), colnames(tbl))[1]
  c_col     <- intersect(c("C.name","C","c","chain"), colnames(tbl))[1]
  v_col     <- intersect(c("V.name","V","v"), colnames(tbl))[1]
  j_col     <- intersect(c("J.name","J","j"), colnames(tbl))[1]
  
  if (is.na(count_col) || is.na(aa_col)) {
    stop("Missing expected columns in sample: ", sample_id,
         "\nFound: ", paste(colnames(tbl), collapse = ", "))
  }
  
  out <- tbl %>%
    transmute(
      Sample = sample_id,
      count  = as.numeric(.data[[count_col]]),
      cdr3aa = as.character(.data[[aa_col]]),
      Cgene  = if (!is.na(c_col)) as.character(.data[[c_col]]) else NA_character_,
      Vgene  = if (!is.na(v_col)) as.character(.data[[v_col]]) else NA_character_,
      Jgene  = if (!is.na(j_col)) as.character(.data[[j_col]]) else NA_character_
    ) %>%
    filter(!is.na(cdr3aa), cdr3aa != "", !is.na(count), count > 0) %>%
    mutate(
      chain = case_when(
        !is.na(Cgene) & str_detect(Cgene, "^TRB") ~ "TRB",
        !is.na(Cgene) & str_detect(Cgene, "^TRA") ~ "TRA",
        !is.na(Cgene) & str_detect(Cgene, "^TRD") ~ "TRD",
        !is.na(Cgene) & str_detect(Cgene, "^TRG") ~ "TRG",
        TRUE ~ NA_character_
      )
    )
  
  out
})

meta_trust4 <- trust4_out_report_coding_Tcells$meta %>%
  mutate(
    Sample = as.character(Sample),
    TimePoint = case_when(
      str_detect(Sample, "T1") ~ "pre_treatment",
      str_detect(Sample, "T2") ~ "post_treatment",
      TRUE ~ NA_character_
    ),
    Patient = str_extract(Sample, "SRC\\d+")
  )

inv_lookup <- setNames(names(src_lookup), unname(src_lookup))
meta_trust4$Patient <- inv_lookup[meta_trust4$Patient]

meta_clin <- SIC_and_SE_clean %>%
  distinct(patient, TreatmentArm, SarcomaImmuneClass, SarcomaEcotypeAssignment, TumorGrade, Histology, 
           PercentNecrosis)

colnames(meta_clin)[colnames(meta_clin) == "patient"] = "Patient"

meta_all <- meta_trust4 %>%
  left_join(meta_clin)

chain_keep <- "TRB" 

tcr_counts <- tcr_long %>%
  filter(is.na(chain) | chain %in% chain_keep) %>%
  group_by(Sample, cdr3aa) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(
    productive_reads = sum(count),
    freq = count / productive_reads
  ) %>%
  ungroup()

tcr_metrics_sample <- tcr_counts %>%
  group_by(Sample) %>%
  summarise(
    productive_reads = unique(productive_reads),
    richness = n(),  # Unique clonotypes
    
    # Shannon Entropy 
    shannon = -sum(freq * log2(freq)),
    
    # Simpson 
    simpson_lambda = sum(freq^2),            
    gini_simpson   = 1 - simpson_lambda,    
    inv_simpson    = 1 / simpson_lambda,     
    
    # Evenness + clonality 
    pielou    = ifelse(richness > 1, shannon / log2(richness), NA_real_),
    clonality = 1 - pielou,
    
    top10_frac = sum(sort(freq, decreasing = TRUE)[seq_len(min(10, length(freq)))]),
    top1_frac  = max(freq),
    .groups = "drop"
  ) %>%
  left_join(meta_all, by = "Sample")

tcr_metrics_sample <- tcr_metrics_sample %>% 
  filter(Patient %in% SIC_and_SE_clean$patient)

SE_levels = c("SE1","SE2","SE3")

Max_ecotype_calls_wide <- df_long %>%
  dplyr::select(patient, Time_Point, `Treatment Arm`, Ecotype, abundance) %>%
  tidyr::pivot_wider(
    id_cols    = c(patient, Time_Point, `Treatment Arm`),
    names_from = Ecotype,
    values_from = abundance,
    values_fn  = mean,         
    values_fill = 0         
  ) %>%
  mutate(
    .mat = as.matrix(dplyr::select(., all_of(SE_levels))),
    .mat2 = .mat,
    .mat2[is.na(.mat2)] <- -Inf,
    .idx = max.col(.mat2, ties.method = "first"),    
    .idx = ifelse(rowSums(is.finite(.mat2)) == 0, NA_integer_, .idx),
    Assigned_Max_Ecotype = SE_levels[.idx],
    Max_SE = ifelse(rowSums(!is.na(.mat)) == 0, NA_real_, apply(.mat, 1, max, na.rm = TRUE))
  ) %>%
  dplyr::select(-.mat, -.mat2, -.idx)

tcr_metrics_sample = tcr_metrics_sample[tcr_metrics_sample$Patient %in% SIC_and_SE_assignments$Patient, ]

colnames(Max_ecotype_calls_wide)[colnames(Max_ecotype_calls_wide) == "Time_Point"] = "TimePoint"
Max_ecotype_calls_wide$TimePoint[Max_ecotype_calls_wide$TimePoint == "Pre-treatment"] = "pre_treatment"
Max_ecotype_calls_wide$TimePoint[Max_ecotype_calls_wide$TimePoint == "Post-treatment"] = "post_treatment"

colnames(Max_ecotype_calls_wide)[colnames(Max_ecotype_calls_wide) == "patient"] = "Patient"
colnames(Max_ecotype_calls_wide)[colnames(Max_ecotype_calls_wide) == "Treatment Arm"] = "TreatmentArm"

tcr_metrics_sample = tcr_metrics_sample %>% 
  left_join(Max_ecotype_calls_wide[, c("Patient", "TimePoint", "TreatmentArm",
                                       "Assigned_Max_Ecotype")])

Max_ecotype_calls_wide_pre = Max_ecotype_calls_wide[Max_ecotype_calls_wide$TimePoint == "pre_treatment", ]
colnames(Max_ecotype_calls_wide_pre)[colnames(Max_ecotype_calls_wide_pre) == "Assigned_Max_Ecotype"] = "Assigned_Max_Ecotype_pre"

tcr_metrics_sample = tcr_metrics_sample %>% 
  left_join(Max_ecotype_calls_wide_pre[, c("Patient", "TreatmentArm",
                                           "Assigned_Max_Ecotype_pre")])

keep <- tcr_metrics_sample %>%
  distinct(Patient, TimePoint) %>%
  dplyr::count(Patient) %>%
  dplyr::filter(n >= 2) %>%
  pull(Patient)

tcr_metrics_sample_paired_pre_post <- tcr_metrics_sample %>% 
  dplyr::filter(Patient %in% keep) %>%
  dplyr::filter(TimePoint %in% c("pre_treatment", "post_treatment"))

metrics <- c(
  "productive_reads",
  "richness",
  "shannon",
  "gini_simpson",
  "inv_simpson",
  "clonality",
  "top10_frac",
  "top1_frac"
)

colnames(tcr_metrics_sample_paired_pre_post)[colnames(tcr_metrics_sample_paired_pre_post) == "SarcomaImmuneClass"] = "SIC_pre"
colnames(tcr_metrics_sample_paired_pre_post)[colnames(tcr_metrics_sample_paired_pre_post) == "SarcomaEcotypeAssignment"] = "SE_pre"

paired <- tcr_metrics_sample_paired_pre_post %>%
  pivot_wider(
    id_cols = c(Patient, TreatmentArm, SIC_pre, SE_pre, Assigned_Max_Ecotype_pre),
    names_from = TimePoint,
    values_from = all_of(metrics),
    names_sep = "__"
  )

delta_long <- paired %>%
  pivot_longer(
    cols = matches("__pre_treatment$|__post_treatment$"),
    names_to = c("metric", "time"),
    names_pattern = "^(.*)__([a-z_]+)$",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = time, values_from = value) %>%
  mutate(delta = post_treatment - pre_treatment) %>%
  filter(metric %in% metrics)

within_arm_tests <- delta_long %>%
  group_by(TreatmentArm, metric) %>%
  summarise(
    n_pairs = sum(!is.na(pre_treatment) & !is.na(post_treatment)),
    p = wilcox.test(post_treatment, pre_treatment, paired = TRUE)$p.value,
    med_pre  = median(pre_treatment, na.rm = TRUE),
    med_post = median(post_treatment, na.rm = TRUE),
    med_delta = median(delta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(TreatmentArm) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(TreatmentArm, p_adj)

selected_metrics <- c("inv_simpson","shannon", "richness")

sig_tbl <- within_arm_tests %>% 
  filter(metric %in% selected_metrics)

plot_dat <- delta_long %>%
  semi_join(sig_tbl, by = c("TreatmentArm","metric")) %>%
  filter(complete.cases(pre_treatment, post_treatment)) %>%
  pivot_longer(c(pre_treatment, post_treatment),
               names_to = "TimePoint", values_to = "value") %>%
  mutate(
    TimePoint = factor(TimePoint, levels = c("pre_treatment", "post_treatment")),
    TreatmentArm = factor(TreatmentArm, levels = c("Control", "Experimental")),
    metric = factor(metric, levels = selected_metrics) 
  )

p <- ggplot(plot_dat, aes(x = TimePoint, y = value)) +
  geom_boxplot(
    aes(fill = TimePoint),
    width = 0.18, outlier.shape = NA,
    alpha = 0.45, 
    linewidth = 0.35
  ) +
  geom_line(
    aes(group = Patient),
    alpha = 0.2, linewidth = 0.35, color = "black"
  ) +
  
  geom_point(
    aes(group = Patient),
    size = 0.2, alpha = 0.5,
  ) +
  stat_summary(aes(group = 1), fun = median, geom = "line",
               linewidth = 0.5, color = "black") +
  
  facet_grid(metric ~ TreatmentArm, scales = "free_y") +
  
  scale_fill_manual(
    values = c("pre_treatment" = "#396CAE", "post_treatment" = "#FBBF85"),
    name = NULL
  ) +
  labs(
    x = NULL, y = "Metric value",
    title = "Within-arm paired changes (metrics shown)"
  ) +
  theme_pubr(border = TRUE, legend = "top") +
  scale_y_continuous(position = "right") +
  geom_signif(comparisons = list(c("pre_treatment", "post_treatment")), 
              test = "wilcox.test") +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    strip.text  = element_text(size = 12),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(p, width = 8, height = 16, 
       units = "cm", filename = '.../Figure_4C_boxplot_pre_post_control_exp.pdf')

###### Figure 5A - PBMCs Clusters CyTOF ###### ------
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

###### Figure 5C - T Cells Clusters CyTOF ###### ------
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



###### Figure 5B - Marker Expression PBMC clusters ###### -----
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

###### Figure 5D - Marker Expression PBMC clusters  ###### -----
percentmatlist_tcells <- readRDS('.../tcell.percent.mat.list.RDS')
tcells_expression <- tcell_umap@assays@data@listData$exprs

sce      <- tcell_umap
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
percent_avg_mat <- sapply(percentmatlist_tcells, function(mat) {
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

# 1) Build marker x cluster matrix of the values you're coloring by
mat_wide <- dot_df_scaled %>%
  select(marker, cluster, scaled_expr) %>%
  pivot_wider(names_from = cluster, values_from = scaled_expr)

mat <- mat_wide %>%
  column_to_rownames("marker") %>%
  as.matrix()

# If any NAs exist, replace with 0 (or marker mean) so cor/dist won't break
mat[is.na(mat)] <- 0

# 2) Correlation distance + hierarchical clustering
# Markers (rows): cluster by similarity across clusters
d_rows <- as.dist(1 - cor(t(mat), use = "pairwise.complete.obs", method = "pearson"))
hc_rows <- hclust(d_rows, method = "average")

# Clusters (cols): cluster by similarity across markers
d_cols <- as.dist(1 - cor(mat, use = "pairwise.complete.obs", method = "pearson"))
hc_cols <- hclust(d_cols, method = "average")

marker_order  <- hc_rows$labels[hc_rows$order]
cluster_order <- hc_cols$labels[hc_cols$order]

# 3) Apply clustered ordering to your plotting data
dot_df_scaled <- dot_df_scaled %>%
  mutate(
    marker  = factor(marker,  levels = marker_order),
    cluster = factor(cluster, levels = cluster_order)
  )

dotplot_cytof_tcells <- ggplot(dot_df_scaled, aes(
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

ggsave(dotplot_cytof_tcells, filename = '.../dotplot_t_cells.pdf', 
       width = 32.5, height = 15, units = "cm")

###### Figure 5E - Trend of PBMC clusters abundance over time ###### ----- 
cluster_T_cells <- readRDS('.../tcell.clusterAssignments.RDS')
percent_T_cells <- readRDS('.../tcell.percent.mat.RDS')

new_cluster_tcells_names = c(
  Cluster4 = "CD4+ Tcm", 
  Cluster2 = "CD4+ Tem", 
  Cluster5 = "Naive CD8+", 
  Cluster1 = "Naive CD4+", 
  Cluster10 = "Activated CD4+", 
  Cluster3 = "Tregs", 
  Cluster9 = "CD8+ CD57+ Teff", 
  Cluster8 = "CD8+ TEMRA", 
  Cluster7 = "Activated CD8+", 
  Cluster6 = "CD8+ T-BET+ Teff", 
  Cluster11 = "CD8+ Trm/Tex"
)

phenotype_CyTOF <- readRDS('.../phenoData.df.RDS')
cluster_all_cells <- readRDS('.../clusterAssignments.RDS')
percent_all_cells <- readRDS('.../percent.mat.RDS')

rownames(percent_all_cells) <- cluster_all_cells[ rownames(percent_all_cells) ]

phenotype_CyTOF$PatientID <- gsub("x", "", phenotype_CyTOF$PatientID)
phenotype_CyTOF$PatientID <- gsub("_", "-", phenotype_CyTOF$PatientID)

phenotype_CyTOF <- phenotype_CyTOF %>% 
  left_join(SIC_and_SE_assignments, by = c("PatientID" = "Patient"))

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
blood_draw_timing <- fread('.../SARC032_blood_draw_timing.csv')

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
    event  = `DFS Event` == "Event occurred",
    t_event = as.numeric(`DFS Time (Days)`)
  ) %>%
  left_join(blood_timing_one_per, by = c("PatientID","Visit_key"))

percent_all_cells <- percent_all_cells[, colnames(percent_all_cells) %in% phenotype_CyTOF$Sample]

# CLR normalization with c mult repl 
repl <- zCompositions::cmultRepl(
  percent_all_cells,
  label = 0,
  output = "prop",
  z.delete = FALSE,
  method = "CZM"
)
clr_mat_cytof <- compositions::clr(as.matrix(repl))

percent_all_cells <- clr_mat_cytof %>% 
  as.data.frame() %>%
  rownames_to_column(var = "CellType")

# Pivot to long format
percent_all_cells_long <- percent_all_cells %>%
  pivot_longer(
    cols = -CellType,           
    names_to = "Sample",        
    values_to = "Abundance"       
  )

percent_all_cells_long <- percent_all_cells_long %>% 
  left_join(df0[, c("PatientID" , "Sample", 
                    "Visit_key", "Visit")], by = c("Sample"))

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

percent_all_cells_long = percent_all_cells_long %>% 
  left_join(pre_ecotype_calls_wide[, c("patient", "Assigned_Max_Ecotype")], 
            by = c("PatientID" = "patient"))


baseline_visit <- "Prior to 1st Tx"

visit_levels <- c(
  "Prior to 1st Tx",
  "1 week Post Initiating treatment",
  "During RT",
  "Before Surgery",
  "3 mo Post-Surgery",
  "12 mo Post-Surgery"
)

df <- percent_all_cells_long %>%
  mutate(
    Arm   = factor(`Treatment Arm`, levels = c("Control","Experimental")),
    Visit = factor(Visit, levels = visit_levels)
  )

# baseline classes per patient
base_class <- df %>%
  filter(Visit == baseline_visit) %>%
  distinct(
    PatientID,
    SE_base  = Assigned_Max_Ecotype,
    SIC_base = `Sarcoma Immune Class`
  )

df <- df %>%
  left_join(base_class, by = "PatientID")

plot_pbmc_trends <- function(df_in, title = "PBMC CLR abundance over time") {
  
  # ensure 1 value per PatientID x Visit x CellType (in case duplicates exist)
  df_pt <- df_in %>%
    group_by(PatientID, Arm, Visit, CellType) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  # mean +/- SE per Arm x Visit x CellType
  df_sum <- df_pt %>%
    group_by(Arm, Visit, CellType) %>%
    summarise(
      n = sum(!is.na(Abundance)),
      mean = mean(Abundance, na.rm = TRUE),
      se   = sd(Abundance, na.rm = TRUE) / sqrt(n),
      ymin = mean - se,
      ymax = mean + se,
      .groups = "drop"
    )
  
  ggplot() +
    # SE band
    geom_ribbon(
      data = df_sum,
      aes(x = Visit, ymin = ymin, ymax = ymax, fill = Arm, group = Arm),
      alpha = 0.18, color = NA
    ) +
    geom_line(
      data = df_sum,
      aes(x = Visit, y = mean, group = Arm, color = Arm),
      linewidth = 1.05
    ) +
    geom_point(
      data = df_sum,
      aes(x = Visit, y = mean, color = Arm),
      size = 1.2
    ) +
    facet_wrap(~CellType, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("Control" = "#9D1536", "Experimental" = "#0093AF")) +
    scale_fill_manual(values  = c("Control" = "#9D1536", "Experimental" = "#0093AF")) +
    theme_pubr(border = T) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 14),
      legend.position = "top"
    ) +
    labs(
      title = title,
      x = "Time point",
      y = "CLR abundance (centered log-ratio)"
    )
}

p_all <- plot_pbmc_trends(df, title = "PBMC CLR abundance over time (all patients)")
p_all

ggsave(
  filename = ".../Figure_5E_time_series_PBMCs_CLR_spaghetti_band.pdf",
  plot = p_all,
  units = "cm", width = 22.5, height = 22.5
)

#Between arm comparisons
stats_between_arms <- compare_means(
  Abundance ~ Arm,
  data = df,
  group.by = c("CellType", "Visit"),
  method = "wilcox.test",
  paired = FALSE
) %>%
  group_by(Visit) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup()

#Vs baseline comparisons 
stats_vs_baseline_within_arm <- compare_means(
  Abundance ~ Visit,
  data = df,
  group.by = c("CellType", "Arm"),
  ref.group = baseline_visit,
  method = "wilcox.test",
  paired = FALSE
) %>%
  group_by(Arm) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup()

run_wilcox_tables <- function(dsub) {
  list(
    between_arms = compare_means(
      Abundance ~ Arm, data = dsub,
      group.by = c("CellType","Visit"),
      method = "wilcox.test", paired = FALSE
    ) %>% group_by(Visit) %>% mutate(p.adj = p.adjust(p, "BH")) %>% ungroup(),
    
    vs_baseline_within_arm = compare_means(
      Abundance ~ Visit, data = dsub,
      group.by = c("CellType","Arm"),
      ref.group = baseline_visit,
      method = "wilcox.test", paired = FALSE
    ) %>% group_by(Arm) %>% mutate(p.adj = p.adjust(p, "BH")) %>% ungroup()
  )
}
tab_all  <- run_wilcox_tables(df)
tab_se1  <- run_wilcox_tables(df %>% filter(SE_base == "SE1"))
tab_sice <- run_wilcox_tables(df %>% filter(SIC_base == "E"))



###### Figure 5F - Trend of T Cell clusters abundance over time ###### ----- 
rownames(percent_T_cells) <- new_cluster_tcells_names[ rownames(percent_T_cells) ]

repl <- zCompositions::cmultRepl(
  percent_T_cells,
  label = 0,
  output = "prop",
  z.delete = FALSE,
  method = "CZM"
)
clr_mat_tcells <- compositions::clr(as.matrix(repl))

percent_T_cells <- t(clr_mat_tcells) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "CellType")

# Pivot to long format
percent_T_cells_long <- percent_T_cells %>%
  pivot_longer(
    cols = -CellType,            # All columns except 'CellType' 
    names_to = "Sample",         # Column to store old column names (sample IDs)
    values_to = "Abundance"        # Column to store the values
  )

names(percent_T_cells_long) = c("Sample", "CellType", "Abundance")

percent_T_cells_long <- percent_T_cells_long %>% 
  left_join(df0[, c("PatientID" , "Sample", 
                    "Visit_key", "Visit")], by = c("Sample"))

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

percent_T_cells_long = percent_T_cells_long %>% 
  left_join(pre_ecotype_calls_wide[, c("patient", "Assigned_Max_Ecotype")], 
            by = c("PatientID" = "patient"))

percent_T_cells_long$Arm <- percent_T_cells_long$`Treatment Arm`

#remove non eligibile patients 
percent_T_cells_long = percent_T_cells_long %>% 
  filter(!(is.na(percent_T_cells_long$Arm)))

baseline_visit <- "Prior to 1st Tx"

visit_levels <- c(
  "Prior to 1st Tx",
  "1 week Post Initiating treatment",
  "During RT",
  "Before Surgery",
  "3 mo Post-Surgery",
  "12 mo Post-Surgery"
)

df <- percent_T_cells_long %>%
  mutate(
    Arm   = factor(`Treatment Arm`, levels = c("Control","Experimental")),
    Visit = factor(Visit, levels = visit_levels)
  )

# baseline classes per patient
base_class <- df %>%
  filter(Visit == baseline_visit) %>%
  distinct(
    PatientID,
    SE_base  = Assigned_Max_Ecotype,
    SIC_base = `Sarcoma Immune Class`
  )

df <- df %>%
  left_join(base_class, by = "PatientID")

plot_tcells_trends <- function(df_in, title = "T Cell Clusters CLR abundance over time") {
  
  # ensure 1 value per PatientID x Visit x CellType (in case duplicates exist)
  df_pt <- df_in %>%
    group_by(PatientID, Arm, Visit, CellType) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")
  
  # mean +/- SE per Arm x Visit x CellType
  df_sum <- df_pt %>%
    group_by(Arm, Visit, CellType) %>%
    summarise(
      n = sum(!is.na(Abundance)),
      mean = mean(Abundance, na.rm = TRUE),
      se   = sd(Abundance, na.rm = TRUE) / sqrt(n),
      ymin = mean - se,
      ymax = mean + se,
      .groups = "drop"
    )
  
  ggplot() +
    # SE band
    geom_ribbon(
      data = df_sum,
      aes(x = Visit, ymin = ymin, ymax = ymax, fill = Arm, group = Arm),
      alpha = 0.18, color = NA
    ) +
    geom_line(
      data = df_sum,
      aes(x = Visit, y = mean, group = Arm, color = Arm),
      linewidth = 1.05
    ) +
    geom_point(
      data = df_sum,
      aes(x = Visit, y = mean, color = Arm),
      size = 1.2
    ) +
    facet_wrap(~CellType, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("Control" = "#9D1536", "Experimental" = "#0093AF")) +
    scale_fill_manual(values  = c("Control" = "#9D1536", "Experimental" = "#0093AF")) +
    theme_pubr(border = T) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 14),
      legend.position = "top"
    ) +
    labs(
      title = title,
      x = "Time point",
      y = "CLR abundance (centered log-ratio)"
    )
}

p_all <- plot_tcells_trends(df, title = "T cell CLR abundance over time (all patients)")
p_all

#Between arm comparisons
stats_between_arms_tcells <- compare_means(
  Abundance ~ Arm,
  data = df,
  group.by = c("CellType", "Visit"),
  method = "wilcox.test",
  paired = FALSE
) %>%
  group_by(Visit) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup()

#Vs baseline comparisons 
stats_vs_baseline_within_arm_tcells <- compare_means(
  Abundance ~ Visit,
  data = df,
  group.by = c("CellType", "Arm"),
  ref.group = baseline_visit,
  method = "wilcox.test",
  paired = FALSE
) %>%
  group_by(Arm) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup()

ggsave(
  filename = ".../Figure_5F_time_series_Tcells_CLR_spaghetti_band.pdf",
  plot = p_all,
  units = "cm", width = 19.5, height = 23.5
)




## Flow Cytometry Anaysis module -----
pick_col <- function(nms, pattern, label = pattern) {
  hit <- grep(pattern, nms, value = TRUE)
  if (length(hit) == 0) stop("No column matched: ", label)
  if (length(hit) > 1) message("Multiple matches for ", label, " — using first:\n  ", hit[1])
  hit[1]
}
make_comp <- function(df, denom_col, parts_named, add_remainder = TRUE, remainder_name = "Other") {
  denom <- df[[denom_col]]
  parts <- as.data.frame(lapply(parts_named, function(cc) df[[cc]]))
  names(parts) <- names(parts_named)
  
  if (add_remainder) {
    rem <- denom - rowSums(parts, na.rm = TRUE)
    parts[[remainder_name]] <- pmax(rem, 0) 
  }
  
  props <- sweep(as.matrix(parts), 1, denom, "/")
  props[!is.finite(props)] <- 0
  as.data.frame(props)
}
clr_comp <- function(props_df, zero_method = "CZM") {
  X <- as.matrix(props_df)
  Xrep <- cmultRepl(X, 
                    z.delete = FALSE,
                    method = zero_method, 
                    output = "prop")
  as.data.frame(clr(Xrep))
}
marker_logratio <- function(marker_count, parent_count, eps = 0.5) {
  p <- (marker_count + eps) / (parent_count + 2*eps)
  log(p) - log(1 - p)  # logit
}

flow_results_filtered <- flow_results_filtered %>%
  mutate(across(-1, as.numeric))
flow_results_filtered <- as.data.table(flow_results_filtered)

nms <- colnames(flow_results_filtered)

# Define populations and parents ---
cd45_live <- pick_col(nms, "/yZombie-/CD45 \\+/Live Scatter \\| Count$", "CD45+ Live Scatter")
cd45_tot  <- pick_col(nms, "/yZombie-/CD45 \\+ \\| Count$", "CD45+ total")

lymph <- pick_col(nms, "/CD45 \\+/Scatter Lymphs \\| Count$", "Scatter Lymphs")
mono  <- pick_col(nms, "/CD45 \\+/Scatter Mono \\| Count$", "Scatter Mono")

cd3   <- pick_col(nms, "/Scatter Lymphs/CD3\\+ \\| Count$", "CD3+")
cd19  <- pick_col(nms, "/Scatter Lymphs/CD19\\+ \\| Count$", "CD19+")
nonTB <- pick_col(nms, "/Scatter Lymphs/NonT&B \\| Count$", "NonT&B")

cd4   <- pick_col(nms, "/Scatter Lymphs/CD3\\+/CD4 \\+ \\| Count$", "CD4+")
treg  <- pick_col(nms, "/Scatter Lymphs/CD3\\+/CD4 \\+/Treg \\| Count$", "Treg")
cd8   <- pick_col(nms, "/Scatter Lymphs/CD3\\+/CD8\\+ \\| Count$", "CD8+")

nk    <- pick_col(nms, "/Scatter Lymphs/NonT&B/NK cells \\| Count$", "NK cells")

# NK subsets
nk_16_56_pos <- pick_col(nms, "/NK cells/CD16\\+CD56\\+ \\| Count$", "NK CD16+CD56+")
nk_16_pos    <- pick_col(nms, "/NK cells/CD16\\+CD56- \\| Count$",  "NK CD16+CD56-")
nk_56_pos    <- pick_col(nms, "/NK cells/CD56\\+CD16- \\| Count$",  "NK CD56+CD16-")

# Mono subsets
m_cd14_16_pos <- pick_col(nms, "/Scatter Mono/CD14\\+CD16\\+ \\| Count$", "Mono CD14+CD16+")
m_cd14_pos    <- pick_col(nms, "/Scatter Mono/CD14\\+CD16- \\| Count$",  "Mono CD14+CD16-")
m_cd16_pos    <- pick_col(nms, "/Scatter Mono/CD14-CD16\\+ \\| Count$",  "Mono CD14-CD16+")

# CD4 memory quadrants (CD45RA x CCR7)
cd4_q1_mem <- pick_col(nms, "/CD4 \\+/Q1: CD45RA .*CD197 .*\\| Count$", "CD4 Q1 mem")
cd4_q2_mem <- pick_col(nms, "/CD4 \\+/Q2: CD45RA .*CD197 .*\\| Count$", "CD4 Q2 mem")
cd4_q3_mem <- pick_col(nms, "/CD4 \\+/Q3: CD45RA .*CD197 .*\\| Count$", "CD4 Q3 mem")
cd4_q4_mem <- pick_col(nms, "/CD4 \\+/Q4: CD45RA .*CD197 .*\\| Count$", "CD4 Q4 mem")

# CD8 memory quadrants
cd8_q1_mem <- pick_col(nms, "/CD8\\+/Q1: CD45RA .*CD197 .*\\| Count$", "CD8 Q1 mem")
cd8_q2_mem <- pick_col(nms, "/CD8\\+/Q2: CD45RA .*CD197 .*\\| Count$", "CD8 Q2 mem")
cd8_q3_mem <- pick_col(nms, "/CD8\\+/Q3: CD45RA .*CD197 .*\\| Count$", "CD8 Q3 mem")
cd8_q4_mem <- pick_col(nms, "/CD8\\+/Q4: CD45RA .*CD197 .*\\| Count$", "CD8 Q4 mem")

#CD4 activated quadrants 
cd4_38neg_drpos = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q1: CD38 PCPC55- , HLADR BV480+ | Count"  
cd4_38pos_drpos = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q2: CD38 PCPC55+ , HLADR BV480+ | Count"     
cd4_38pos_drneg = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q3: CD38 PCPC55+ , HLADR BV480- | Count"     
cd4_38neg_drneg = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Q4: CD38 PCPC55- , HLADR BV480- | Count"     

#CD8 activated quadrants 
cd8_38neg_drpos = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q1: CD38 PCPC55- , HLADR BV480+ | Count"  
cd8_38pos_drpos = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q2: CD38 PCPC55+ , HLADR BV480+ | Count"     
cd8_38pos_drneg = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q3: CD38 PCPC55+ , HLADR BV480- | Count"     
cd8_38neg_drneg = "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/Q4: CD38 PCPC55- , HLADR BV480- | Count"     

#CD63 and CD68 monos 
cd68_cd163_pos <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Mono/CD68+CD163+ | Count"
cd68_cd163_neg <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Mono/CD68+CD163- | Count"

# Compositional analysis by major populations ----
# 2) Lymph lineages
props_lymph <- make_comp(
  flow_results_filtered,
  denom_col = lymph,
  parts_named = c(CD3 = cd3, CD19 = cd19, NonTB = nonTB),
  add_remainder = FALSE
)

# 3) T lineage within CD3
props_T <- make_comp(
  flow_results_filtered,
  denom_col = cd3,
  parts_named = c(CD4 = cd4, CD8 = cd8),
  add_remainder = TRUE, remainder_name = "Other_T"
)

# 4) Treg within CD4
props_Treg <- make_comp(
  flow_results_filtered,
  denom_col = cd4,
  parts_named = c(Treg = treg),
  add_remainder = TRUE, remainder_name = "NonTreg"
)

# 5) CD4 memory quadrants
props_cd4_mem <- make_comp(
  flow_results_filtered,
  denom_col = cd4,
  parts_named = c(
    Q1_CD45RAneg_CCR7pos = cd4_q1_mem,
    Q2_CD45RApos_CCR7pos = cd4_q2_mem,
    Q3_CD45RApos_CCR7neg = cd4_q3_mem,
    Q4_CD45RAneg_CCR7neg = cd4_q4_mem
  ),
  add_remainder = FALSE
)

# 6) CD8 memory quadrants
props_cd8_mem <- make_comp(
  flow_results_filtered,
  denom_col = cd8,
  parts_named = c(
    Q1_CD45RAneg_CCR7pos = cd8_q1_mem,
    Q2_CD45RApos_CCR7pos = cd8_q2_mem,
    Q3_CD45RApos_CCR7neg = cd8_q3_mem,
    Q4_CD45RAneg_CCR7neg = cd8_q4_mem
  ),
  add_remainder = FALSE
)

# 7) NK vs Other NonT&B
props_nonTB <- make_comp(
  flow_results_filtered,
  denom_col = nonTB,
  parts_named = c(NK = nk),
  add_remainder = TRUE, remainder_name = "Other_NonTB"
)

# 8) NK subsets (+ Other NK)
props_nk <- make_comp(
  flow_results_filtered,
  denom_col = nk,
  parts_named = c(
    NK_CD16pos_CD56pos = nk_16_56_pos,
    NK_CD16pos_CD56neg = nk_16_pos,
    NK_CD16neg_CD56pos = nk_56_pos
  ),
  add_remainder = TRUE, remainder_name = "Other_NK"
)

# 9) Mono CD14/CD16 subsets (+ Other mono)
props_mono <- make_comp(
  flow_results_filtered,
  denom_col = mono,
  parts_named = c(
    Mono_CD14pos_CD16neg = m_cd14_pos,
    Mono_CD14pos_CD16pos = m_cd14_16_pos,
    Mono_CD14neg_CD16pos = m_cd16_pos
  ),
  add_remainder = TRUE, remainder_name = "Other_Mono"
)

#Mono CD68 and CD163 subsets 
props_mono_mac_cd163_cd68 <- make_comp(
  flow_results_filtered,
  denom_col = mono,
  parts_named = c(
    CD68pos_CD163pos = cd68_cd163_pos,
    CD68pos_CD163neg = cd68_cd163_neg
  ),
  add_remainder = TRUE, remainder_name = "Other_Mono"
)

#MDSC gate
# denominator 
cd11b_drneg <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Live Scatter/CD11b+DR- | Count"

# subsets 
gmdsc <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Live Scatter/CD11b+DR-/gMDSC | Count"
mmdsc <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Live Scatter/CD11b+DR-/mMDSC | Count"

props_mdsc <- make_comp(
  flow_results_filtered,
  denom_col = cd11b_drneg,
  parts_named = c(
    gMDSC = gmdsc,
    mMDSC = mmdsc
  ),
  add_remainder = TRUE, remainder_name = "Other_Mono"
)

clr_lymph     <- clr_comp(props_lymph)
clr_lymph <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_lymph))

clr_T         <- clr_comp(props_T)
clr_T <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_T))

clr_Treg      <- clr_comp(props_Treg)
clr_Treg <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_Treg))

clr_cd4_mem   <- clr_comp(props_cd4_mem)
clr_cd4_mem <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_cd4_mem))
clr_cd4_mem = rownames_to_column(clr_cd4_mem, var = "Patient")

clr_cd8_mem   <- clr_comp(props_cd8_mem)
clr_cd8_mem <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_cd8_mem))
clr_cd8_mem = rownames_to_column(clr_cd8_mem, var = "Patient")

clr_nonTB     <- compositions::clr(props_nonTB)
clr_nonTB <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_nonTB))

clr_nk        <- clr_comp(props_nk)
clr_nk <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_nk))

clr_mono      <- clr_comp(props_mono)
clr_mono <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_mono))

clr_mono_cd163_68      <- compositions::clr(props_mono_mac_cd163_cd68)
clr_mono_cd163_68 <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_mono_cd163_68))

clr_mdscs      <- clr_comp(props_mdsc)
clr_mdscs <- cbind(Patient = flow_results_filtered$Patient, as.data.frame(clr_mdscs))

# marker-by-marker log-ratio for overlapping CD4/CD8 markers (not compositional) ----
logit_from_counts <- function(x, n, eps = 0.5) {
  p <- (x + eps) / (n + 2*eps)
  log(p / (1 - p))
}
make_nonoverlap_marker_df <- function(df,
                                      patient_col = "Patient",
                                      parent_col,
                                      markers_named,
                                      eps = 0.5,
                                      prefix = NULL) {
  stopifnot(patient_col %in% names(df))
  stopifnot(parent_col  %in% names(df))
  stopifnot(all(unname(markers_named) %in% names(df)))
  
  Patient <- df[[patient_col]]
  n <- as.numeric(df[[parent_col]])
  
  out <- tibble(
    Patient = Patient,
    !!paste0("n_", ifelse(is.null(prefix), "parent", prefix)) := n
  )
  
  for (nm in names(markers_named)) {
    col <- markers_named[[nm]]
    x <- as.numeric(df[[col]])
    
    frac_raw <- ifelse(n > 0, x / n, NA_real_)
    logit_v  <- logit_from_counts(x, n, eps = eps)
    
    out[[paste0(nm, "_count")]]    <- x
    out[[paste0(nm, "_frac_raw")]] <- frac_raw
    out[[paste0(nm, "_logit")]]    <- logit_v
  }
  
  out
}

cd4_parent <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 + | Count"
pd1_cd4    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD279 + | Count"
cd103_cd4    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD103+ | Count"
ctla4_cd4    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD152+ | Count"
cd39_cd4    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD39+ | Count"
tim3_cd4    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/CD366+ | Count"

cd8_parent <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+ | Count"
pd1_cd8    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD279 + | Count"
cd103_cd8    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD103+ | Count"
ctla4_cd8    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD152+ | Count"
cd39_cd8    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD39+ | Count"
tim3_cd8    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD8+/CD366+ | Count"


tregs_parent <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg | Count"
pd1_tregs    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD279 + | Count"
cd103_tregs    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD103+ | Count"
ctla4_tregs    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD152+ | Count"
cd39_tregs  <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD39+ | Count"
tim3_tregs    <- "PeacoQCGoodEvents/Time/FSC Singlets/SSC Singlets/yZombie-/CD45 +/Scatter Lymphs/CD3+/CD4 +/Treg/CD366+ | Count"


cd4_markers <- c(
  PD1   = pd1_cd4,
  CTLA4 = ctla4_cd4,
  CD103 = cd103_cd4,
  CD39  = cd39_cd4,
  TIM3  = tim3_cd4
)

df_cd4_markers <- make_nonoverlap_marker_df(
  df = flow_results_filtered,
  patient_col = "Patient",
  parent_col  = cd4_parent,
  markers_named = cd4_markers,
  eps = 0.5,
  prefix = "CD4"
)

cd8_markers <- c(
  PD1   = pd1_cd8,
  CTLA4 = ctla4_cd8,
  CD103 = cd103_cd8,
  CD39  = cd39_cd8,
  TIM3  = tim3_cd8
)

df_cd8_markers <- make_nonoverlap_marker_df(
  df = flow_results_filtered,
  patient_col = "Patient",
  parent_col  = cd8_parent,
  markers_named = cd8_markers,
  eps = 0.5,
  prefix = "CD8"
)

tregs_markers <- c(
  PD1   = pd1_tregs,
  CTLA4 = ctla4_tregs,
  CD103 = cd103_tregs,
  CD39  = cd39_tregs,
  TIM3  = tim3_tregs
)

df_tregs_markers <- make_nonoverlap_marker_df(
  df = flow_results_filtered,
  patient_col = "Patient",
  parent_col  = tregs_parent,
  markers_named = tregs_markers,
  eps = 0.5,
  prefix = "Tregs"
)

# Prep generic metadata ----
baseline_meta = SIC_and_SE_assignments

baseline_meta = baseline_meta %>% 
  left_join(pre_ecotype_calls_wide[, c("patient", "Assigned_Max_Ecotype")], 
            by = c("Patient" = "patient"))

baseline_meta$TreatmentArm = baseline_meta$`Treatment Arm`
baseline_meta$SE_pre = baseline_meta$`Sarcoma Ecotype Assignment`
baseline_meta$SIC_pre = baseline_meta$`Sarcoma Immune Class`
baseline_meta$Grade = baseline_meta$`Tumor Grade`
baseline_meta_filtered = baseline_meta[baseline_meta$Patient %in% flow_results_filtered$Patient, ]

###### Figure 6A - CD8+ T Cells Memory Subsets  ###### ----- 
clr_cd8_mem = column_to_rownames(clr_cd8_mem, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(clr_cd8_mem)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-4, -2, 0, 2, 4), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(clr_cd8_mem), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "% of total Lymphocytes",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)

###### Figure 6B - CD4+ T Cells Memory Subsets  ###### ----- 
clr_cd4_mem = column_to_rownames(clr_cd4_mem, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(clr_cd4_mem)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-10, -5, 0, 5, 10), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(clr_cd4_mem), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "CLR",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)

###### Figure 6C - Dot Plot Experimental vs Control Abundance Diff. CD4+ and CD8+ T Memory subsets ###### -----
run_wilcox_arm <- function(df_long, subset_name, filter_expr) {
  dat <- df_long %>%
    filter(!!rlang::enquo(filter_expr)) %>%
    filter(!is.na(CLR), !is.na(TreatmentArm)) %>%
    mutate(TreatmentArm = factor(TreatmentArm))  
  
  res <- dat %>%
    group_by(Group, Feature) %>%
    summarize(
      n = n(),
      n_arm = n_distinct(TreatmentArm),
      effect = mean(CLR[TreatmentArm == levels(TreatmentArm)[2]], na.rm = TRUE) -
        mean(CLR[TreatmentArm == levels(TreatmentArm)[1]], na.rm = TRUE),
      p = ifelse(n_arm == 2,
                 wilcox.test(CLR ~ TreatmentArm, exact = FALSE)$p.value,
                 NA_real_),
      .groups = "drop"
    ) %>%
    mutate(
      subset = subset_name
    ) %>%
    group_by(Group) %>%
    mutate(p_adj = p.adjust(p, method = "BH"), 
           neglog10_fdr = -log10(p_adj))
  
  res
}

clr_list <- list(
  Lymph = clr_lymph,
  Tcells = clr_T,
  Treg = clr_Treg,
  CD4_mem = clr_cd4_mem,
  CD8_mem = clr_cd8_mem,
  NonTB = clr_nonTB,      
  NK = clr_nk,
  Mono = clr_mono,
  Mono_CD68CD163 = clr_mono_cd163_68,
  MDSC = clr_mdscs
)

# 2) Long format
clr_long <- bind_rows(lapply(names(clr_list), function(g) {
  clr_list[[g]] %>%
    pivot_longer(cols = -Patient, names_to = "Feature", values_to = "CLR") %>%
    mutate(Group = g)
}))

meta_flow <- baseline_meta %>%
  distinct(Patient, .keep_all = TRUE) %>%
  transmute(
    Patient,
    TreatmentArm,
    SE_pre,
    SIC_pre,
    Histology,
    Grade,
    Assigned_Max_Ecotype
  )

clr_long <- clr_long %>%
  left_join(meta_flow, by = "Patient")

# Define your subsets (edit labels/values to match your metadata)
res_all   <- run_wilcox_arm(clr_long, "All", TRUE)

res_se1   <- run_wilcox_arm(clr_long %>% dplyr::filter(SE_pre != "NA"), 
                            "SE1", SE_pre == "SE1")
res_nonse <- run_wilcox_arm(clr_long %>% dplyr::filter(SE_pre != "NA"), 
                            "Non-SE1", SE_pre != "SE1")

res_se1_max   <- run_wilcox_arm(clr_long %>% dplyr::filter(SE_pre != "NA"), 
                                "SE1", SE_pre == "SE1")
res_nonse_max <- run_wilcox_arm(clr_long %>% dplyr::filter(SE_pre != "NA"), 
                                "Non-SE1", SE_pre != "SE1")

res_sice  <- run_wilcox_arm(clr_long %>% dplyr::filter(SIC_pre != "NA"), 
                            "SIC E", SIC_pre == "E")
res_nons  <- run_wilcox_arm(clr_long %>% dplyr::filter(SIC_pre != "NA"), 
                            "Non-SIC E", SIC_pre != "E")

wilcox_res <- bind_rows(res_all, res_se1, res_nonse, res_sice, res_nons) %>%
  filter(!is.na(p_adj)) 

sig_markers <- wilcox_res %>%
  filter(p < 0.05) %>%        
  pull(Feature) %>%
  unique()

plot_df <- wilcox_res %>%
  filter(Feature %in% sig_markers) %>%
  mutate(sig = p_adj < 0.05) %>% 
  filter(Group %in% c("CD4_mem", "CD8_mem", "Treg"))

plot_df_cd4_cd8 = plot_df %>% filter(Group %in% c("CD4_mem", 
                                                  "CD8_mem")) 

cd4_cd8_mem_dotplot = ggplot(data = plot_df_cd4_cd8, aes(x = subset, y = Feature)) +
  # non-significant
  geom_point(
    data = dplyr::filter(plot_df_cd4_cd8, !sig),
    aes(size = neglog10_fdr, fill = effect),
    alpha = 0.5, stroke = 0, shape = 21
  ) +
  # significant
  geom_point(
    data = dplyr::filter(plot_df_cd4_cd8, sig),
    aes(size = neglog10_fdr, fill = effect),
    alpha = 1, stroke = 1, shape = 21, colour = "black"
  ) +
  facet_wrap(~ Group, ncol = 2, scales = "free_x") +
  scale_size_continuous(
    name = expression(-log[10]("BH-FDR")),
    range = c(1.2, 8)
  ) +
  scale_fill_distiller(
    name = "Mean delta CLR\n(Exp - Ctrl)",
    palette = "RdBu",
    direction = -1
  ) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(
    strip.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey93")
  )

ggsave(cd4_cd8_mem_dotplot,
       unit = "cm",  width = 20, height = 10, 
       filename = '.../Fgiure_6C_cd4_mem_cd8_mem_bubblePlot.pdf')


###### Figure 6D - Heat map CD4+ T cells exhaustion markers ###### -----
df_cd4_markers_cleaned_logit = df_cd4_markers[, c("Patient", names(df_cd4_markers)[grepl("_logit", colnames(df_cd4_markers))])]
# Logit 
df_cd4_markers_cleaned_logit = column_to_rownames(df_cd4_markers_cleaned_logit, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(df_cd4_markers_cleaned_logit)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-4, -2, 0, 2, 4), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(df_cd4_markers_cleaned_logit), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "CLR",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)


###### Figure 6E - Heat map CD8+ T cells exhaustion markers ###### -----
df_cd8_markers_cleaned_logit = df_cd8_markers[, c("Patient", names(df_cd8_markers)[grepl("_logit", colnames(df_cd8_markers))])]
df_cd8_markers_cleaned_logit = df_cd8_markers_cleaned_logit[df_cd8_markers_cleaned_logit$Patient != "005-013", ]
df_cd8_markers_cleaned_logit = column_to_rownames(df_cd8_markers_cleaned_logit, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(df_cd8_markers_cleaned_logit)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-4, -2, 0, 2, 4), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(df_cd8_markers_cleaned_logit), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "CLR",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)


##### Figure 6F and Figure Supplemental 7B - Dot Plot #####
meta_flow <- baseline_meta %>%
  distinct(Patient, .keep_all = TRUE) %>%
  transmute(
    Patient,
    TreatmentArm,
    SE_pre,
    SIC_pre,
    Histology,
    Grade,
    Assigned_Max_Ecotype
  )

df_cd8_markers <- df_cd8_markers %>%
  left_join(meta_flow, by = "Patient")

df_cd4_markers <- df_cd4_markers %>%
  left_join(meta_flow, by = "Patient")

df_tregs_markers <- df_tregs_markers %>%
  left_join(meta_flow, by = "Patient")

prep_arm <- function(x, control_label = NULL, experimental_label = NULL) {
  x <- as.character(x)
  levs <- sort(unique(x[!is.na(x)]))
  if (!is.null(control_label) && !is.null(experimental_label)) {
    stopifnot(all(c(control_label, experimental_label) %in% levs))
    return(factor(x, levels = c(control_label, experimental_label)))
  }
  
  ctrl_guess <- levs[str_detect(tolower(levs), "control|rt only|radiation only|standard")]
  exp_guess  <- levs[str_detect(tolower(levs), "experimental|pembro|combo|pembrolizumab")]
  
  if (length(ctrl_guess) >= 1 && length(exp_guess) >= 1) {
    return(factor(x, levels = c(ctrl_guess[1], exp_guess[1])))
  }
  
  # Fallback: just use first two unique levels
  if (length(levs) < 2) return(factor(x))
  factor(x, levels = levs[1:2])
}

run_wilcox_features <- function(df, feature_cols, subset_label,
                                subset_filter = rep(TRUE, nrow(df)),
                                control_label = NULL, experimental_label = NULL,
                                min_per_arm = 3) {
  
  d <- df %>%
    filter(subset_filter) %>%
    mutate(TreatmentArm = prep_arm(TreatmentArm,
                                   control_label = control_label,
                                   experimental_label = experimental_label))
  
  arm_levels <- levels(d$TreatmentArm)
  if (length(arm_levels) < 2) {
    stop("TreatmentArm has <2 levels after filtering subset: ", subset_label)
  }
  arm_ctrl <- arm_levels[1]
  arm_exp  <- arm_levels[2]
  
  d_long <- d %>%
    dplyr::select(Patient, TreatmentArm, all_of(feature_cols)) %>%
    pivot_longer(cols = all_of(feature_cols),
                 names_to = "feature", values_to = "value") %>%
    mutate(value = as.numeric(value))
  
  res <- d_long %>%
    group_by(feature) %>%
    summarise(
      subset = subset_label,
      arm_control = arm_ctrl,
      arm_experimental = arm_exp,
      n_control = sum(TreatmentArm == arm_ctrl & !is.na(value)),
      n_experimental = sum(TreatmentArm == arm_exp & !is.na(value)),
      mean_control = mean(value[TreatmentArm == arm_ctrl], na.rm = TRUE),
      mean_experimental = mean(value[TreatmentArm == arm_exp], na.rm = TRUE),
      mean_delta = mean_experimental - mean_control,
      p = if (n_control >= min_per_arm && n_experimental >= min_per_arm) {
        suppressWarnings(wilcox.test(value ~ TreatmentArm, exact = FALSE)$p.value)
      } else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           neglog10_fdr = -log10(p_adj))
  
  res
}
run_wilcox_panel <- function(df, feature_cols,
                             control_label = NULL, experimental_label = NULL,
                             min_per_arm = 3) {
  
  bind_rows(
    run_wilcox_features(df, feature_cols, "All",
                        subset_filter = rep(TRUE, nrow(df)),
                        control_label = control_label,
                        experimental_label = experimental_label,
                        min_per_arm = min_per_arm),
    
    run_wilcox_features(df, feature_cols, "SE1",
                        subset_filter = df$SE_pre == "SE1",
                        control_label = control_label,
                        experimental_label = experimental_label,
                        min_per_arm = min_per_arm),
    
    run_wilcox_features(df, feature_cols, "non-SE1",
                        subset_filter = !is.na(df$SE_pre) & df$SE_pre != "SE1",
                        control_label = control_label,
                        experimental_label = experimental_label,
                        min_per_arm = min_per_arm),
    
    run_wilcox_features(df, feature_cols, "SIC E",
                        subset_filter = df$SIC_pre == "E",
                        control_label = control_label,
                        experimental_label = experimental_label,
                        min_per_arm = min_per_arm),
    
    run_wilcox_features(df, feature_cols, "non-SIC E",
                        subset_filter = !is.na(df$SIC_pre) & df$SIC_pre != "E",
                        control_label = control_label,
                        experimental_label = experimental_label,
                        min_per_arm = min_per_arm)
  )
}

cd4_logit_cols <- names(df_cd4_markers)[grepl("logit", names(df_cd4_markers))]
cd4_frac_cols  <- names(df_cd4_markers)[grepl("frac", names(df_cd4_markers))]

cd8_logit_cols <- names(df_cd8_markers)[grepl("logit", names(df_cd8_markers))]
cd8_frac_cols  <- names(df_cd8_markers)[grepl("frac", names(df_cd8_markers))]

treg_logit_cols <- names(df_tregs_markers)[grepl("logit", names(df_tregs_markers))]
treg_frac_cols  <- names(df_tregs_markers)[grepl("frac", names(df_tregs_markers))]

control_label <- "Control"
experimental_label <- "Experimental"

# Run panels
wilcox_cd4_logit <- run_wilcox_panel(df_cd4_markers, cd4_logit_cols,
                                     control_label, experimental_label)

wilcox_cd4_frac  <- run_wilcox_panel(df_cd4_markers, cd4_frac_cols,
                                     control_label, experimental_label)

wilcox_cd8_logit <- run_wilcox_panel(df_cd8_markers, cd8_logit_cols,
                                     control_label, experimental_label)

wilcox_cd8_frac  <- run_wilcox_panel(df_cd8_markers, cd8_frac_cols,
                                     control_label, experimental_label)

wilcox_treg_logit <- run_wilcox_panel(df_tregs_markers, treg_logit_cols,
                                      control_label, experimental_label)

wilcox_all <- bind_rows(
  wilcox_cd4_logit %>% mutate(compartment = "CD4", scale = "logit"),
  wilcox_cd4_frac  %>% mutate(compartment = "CD4", scale = "fraction"),
  wilcox_cd8_logit %>% mutate(compartment = "CD8", scale = "logit"),
  wilcox_cd8_frac  %>% mutate(compartment = "CD8", scale = "fraction")
)

wilcox_logit <- bind_rows(
  wilcox_cd4_logit %>% mutate(compartment = "CD4", scale = "logit"),
  wilcox_cd8_logit %>% mutate(compartment = "CD8", scale = "logit")
)

sig_markers = wilcox_logit %>% 
  filter(p < 0.05) %>% 
  pull(feature)

## Fig. 6F ----
dotplot_cd4_cd8_overlap = wilcox_logit %>% 
  filter(feature %in% sig_markers) %>% 
  ggplot(aes(y = feature, x = subset)) +
  geom_point(
    data = wilcox_logit %>% filter(p_adj >= 0.05),
    aes(size = neglog10_fdr, fill = mean_delta),
    alpha = 0.5, stroke = 0,  shape = 21
  ) + 
  geom_point(
    data = wilcox_logit %>% filter(p_adj < 0.05),
    aes(size = neglog10_fdr, fill = mean_delta),
    alpha = 1, stroke = 1, shape = 21, colour = "black"
  ) +
  facet_wrap(~ compartment, ncol = 2) +
  scale_size_continuous(
    name = expression(-log[10]("BH-FDR")),
    range = c(1.2, 8)
  ) +
  scale_fill_distiller(
    name = ,
    palette = "RdBu",
    direction = -1
  ) +
  theme_pubr(border = T, legend = "right") +
  theme(
    strip.text = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey93")
  )+
  labs(
    x = NULL, y = NULL,
    size = "-log10(FDR)",
    color = "Mean delta logit\n(Exp - Ctrl)"
  )

ggsave(dotplot_cd4_cd8_overlap,
       unit = "cm",  width = 18.5, height = 10, 
       filename = '.../Figure_6F_cd4_cd8_overlap_dotplot.pdf')

## Fig. S7B ----
sig_markers = wilcox_treg_logit %>% 
  filter(p < 0.05) %>% 
  pull(feature)

wilcox_treg_logit$compartment = "Tregs"

dotplot_treg_overlap = wilcox_treg_logit %>% 
  filter(feature %in% sig_markers) %>% 
  ggplot(aes(y = feature, x = subset)) +
  geom_point(
    data = wilcox_treg_logit %>% filter(p_adj >= 0.05),
    aes(size = neglog10_fdr, fill = mean_delta),
    alpha = 0.5, stroke = 0,  shape = 21
  ) + 
  geom_point(
    data = wilcox_treg_logit %>% filter(p_adj < 0.05),
    aes(size = neglog10_fdr, fill = mean_delta),
    alpha = 1, stroke = 1, shape = 21, colour = "black"
  ) +
  facet_wrap(~ compartment, ncol = 2) +
  scale_size_continuous(
    name = expression(-log[10]("BH-FDR")),
    range = c(1.2, 8)
  ) +
  scale_fill_distiller(
    name = ,
    palette = "RdBu",
    direction = -1
  ) +
  theme_pubr(border = T, legend = "right") +
  theme(
    strip.text = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16), 
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey93")
  )+
  labs(
    x = NULL, y = NULL,
    size = "-log10(FDR)",
    color = "Mean delta logit\n(Exp - Ctrl)"
  )

ggsave(dotplot_treg_overlap,
       unit = "cm",  width = 15, height = 10, 
       filename = '.../Figure_7b_tregs_dotplot.pdf')

##### Figure Supplemental 7A ###### -----
#df_cd8_markers_cleaned_raw = df_cd8_markers[, c("Patient", names(df_cd8_markers)[grepl("_frac", colnames(df_cd8_markers))])]
df_tregs_markers_cleaned_logit = df_tregs_markers[, c("Patient", names(df_tregs_markers)[grepl("_logit", colnames(df_tregs_markers))])]

df_tregs_markers_cleaned_logit = df_tregs_markers_cleaned_logit[df_tregs_markers_cleaned_logit$Patient != "005-013", ]
df_tregs_markers_cleaned_logit = column_to_rownames(df_tregs_markers_cleaned_logit, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(df_tregs_markers_cleaned_logit)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-4, -2, 0, 2, 4), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(df_tregs_markers_cleaned_logit), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "CLR",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)

##### Figure Supplemental 7C ###### -----
##### Figure Supplemental 7C ###### -----
clr_mono_cd163_68 = column_to_rownames(clr_mono_cd163_68, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(clr_mono_cd163_68)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-4, -2, 0, 2, 4), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(clr_mono_cd163_68), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "% of total Lymphocytes",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)

##### Figure Supplemental 7D ##### -----
sig_markers <- wilcox_res %>%
  filter(p < 0.05) %>%          
  pull(Feature) %>%
  unique()

plot_df <- wilcox_res %>%
  mutate(sig = p_adj < 0.05) %>% 
  filter(Group %in% c("Mono_CD68CD163"))

plot_df_mono_cd68 = plot_df %>% filter(Group %in% c("Mono_CD68CD163")) 

mono_cd68_dotplot = ggplot(data = plot_df_mono_cd68, aes(x = subset, y = Feature)) +
  geom_point(
    data = dplyr::filter(plot_df_mono_cd68, !sig),
    aes(size = neglog10_fdr, fill = effect),
    alpha = 0.5, stroke = 0, shape = 21
  ) +
  geom_point(
    data = dplyr::filter(plot_df_mono_cd68, sig),
    aes(size = neglog10_fdr, fill = effect),
    alpha = 1, stroke = 1, shape = 21, colour = "black"
  ) +
  facet_wrap(~ Group, ncol = 2, scales = "free_x") +
  scale_size_continuous(
    name = expression(-log[10]("BH-FDR")),
    range = c(1.2, 8)
  ) +
  scale_fill_distiller(
    name = "Mean delta CLR\n(Exp - Ctrl)",
    palette = "RdBu",
    direction = -1
  ) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(
    strip.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey93")
  )

ggsave(mono_cd68_dotplot,
       unit = "cm",  width = 18.5, height = 10, 
       filename = '.../Fig_S7D_mono_cd68_dotplot.pdf')

##### Figure Supplemental 7E ###### -----
clr_mono = column_to_rownames(clr_mono, var = "Patient")
anno_df = baseline_meta_filtered %>% 
  dplyr::select(Patient, TreatmentArm, Grade, Histology) 
rownames(baseline_meta_filtered) = baseline_meta_filtered$Patient
order = rownames(clr_mono)
anno_df = anno_df[match(order, anno_df$Patient), ]
anno_df$Patient = NULL

col_annotation <- HeatmapAnnotation(
  df = anno_df,
  annotation_height = unit(1, "cm"),
  which = "column", 
  col = list(Grade = c("Grade 2" = "#d682ff", "Grade 3" = "#0095ff"), 
             TreatmentArm = c("Control" = "#61a5c2", "Experimental" = "#9e2a2b"), 
             Histology = c("UPS" = "#7fc97f", "LPS" = "#ffff99", 
                           "Myxoma" = "grey")))

col_custom <- colorRamp2(c(-4, -2, 0, 2, 4), 
                         c("#4500ADFF","#6B58EFFF","#f7f7f7", "#FF5A5AFF", "#D60C00FF"))

ht = Heatmap(t(clr_mono), 
             col = col_custom, 
             column_split = as.factor(anno_df$TreatmentArm),
             show_column_names = F, 
             border = T, 
             heatmap_height = unit(6, "cm"), 
             heatmap_width = unit(12, "cm"),
             cluster_column_slices = F, 
             name = "% of total Lymphocytes",
             top_annotation = col_annotation)

ht <- draw(ht, heatmap_legend_side = "left", 
           annotation_legend_side = "left", 
           merge_legend = TRUE)

##### Figure Supplemental 7F ###### -----
sig_markers <- wilcox_res %>%
  filter(p < 0.05) %>%         
  pull(Feature) %>%
  unique()

plot_df <- wilcox_res %>%
  filter(Feature %in% sig_markers) %>%
  mutate(sig = p_adj < 0.05) %>% 
  filter(Group %in% c("Mono"))

plot_df_mono_cd14 = plot_df %>% filter(Group %in% c("Mono")) 

mono_cd14_dotplot = ggplot(data = plot_df_mono_cd14, aes(x = subset, y = Feature)) +
  geom_point(
    data = dplyr::filter(plot_df_mono_cd14, !sig),
    aes(size = neglog10_fdr, fill = effect),
    alpha = 0.5, stroke = 0, shape = 21
  ) +
  geom_point(
    data = dplyr::filter(plot_df_mono_cd14, sig),
    aes(size = neglog10_fdr, fill = effect),
    alpha = 1, stroke = 1, shape = 21, colour = "black"
  ) +
  facet_wrap(~ Group, ncol = 2, scales = "free_x") +
  scale_size_continuous(
    name = expression(-log[10]("BH-FDR")),
    range = c(1.2, 8)
  ) +
  scale_fill_distiller(
    name = "Mean delta CLR\n(Exp - Ctrl)",
    palette = "RdBu",
    direction = -1
  ) +
  theme_pubr(border = TRUE, legend = "right") +
  theme(
    strip.text = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.35, colour = "grey88"),
    panel.grid.minor = element_line(linewidth = 0.20, colour = "grey93")
  )

ggsave(mono_cd14_dotplot,
       unit = "cm",  width = 20, height = 10, 
       filename = '.../Figure_S7F_mono_cd14_dotplot.pdf')
