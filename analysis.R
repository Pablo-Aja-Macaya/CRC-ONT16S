#######################
# ---- Libraries ---- #
#######################
suppressMessages({
  # Libraries
  library(ggplot2)
  library(patchwork)
  # library(plotly)
  library(vegan)
  library(devtools)
  library(data.table)
  library(RColorBrewer)
  library(egg)
  library(randomcoloR)
  library(ggvenn)
  library(glue)
  library(tidyr)
  library(purrr)
  library(plyr)
  library(dplyr)
  library(phyloseq)
  library(ggalluvial)
  library(stringr)
  library(gtools)
  library(ggsignif)
  library(ggpubr)
  library(rstatix)
  library(ggmagnify)
  library(ggtext)
  library(ggrepel)
})

# Note: ANCOM-BC is loaded in its own section, this is because it loads tidytree and 
# creates an unavoidable warning/message each time a function in the the phyloseq package is used,
# which when outputting to the report makes it difficult to read

# General function
map_signif_level_func = function(x){
  # c("***"=0.001, "**"=0.01, "*"=0.05)
  if (x<=0.001){
    return("***")
  } else if (x<=0.01) {
    return("**")
  } else if (x<=0.05) {
    return("*")
  } else {
    return(NA)
  }
}

# Options
options(dplyr.summarise.inform = FALSE)
options(getClass.msg=FALSE)

# Set seed
set.seed(123)


###################
# ---- Input ---- #
###################
# Set up
# feature_table_file <- opt$features
# contamination_file <- opt$contaminants
# metadata_file <- opt$metadata
# taxmat_file <- opt$tax
# tree_file <- opt$tree
# variables_of_interest <- gsub("-", ".", strsplit(opt$voi, split=";")[[1]])
# output_folder <- opt$output_folder
# biosample_col <- opt$biosample_col
# apply_controls <- FALSE
# 
BASE_DIR = "/home/usuario/Proyectos/Results/ONT16S"
feature_table_file <- as.character(glue("/home/usuario/Proyectos/Results/ONT16S/pipeline/15_biom_conversion/feature_table.biom.json"))
metadata_file <- as.character(glue("/home/usuario/Proyectos/Results/ONT16S/pipeline/15_biom_conversion/metadata_table.tsv"))
taxmat_file <- as.character(glue("/home/usuario/Proyectos/Results/ONT16S/pipeline/15_biom_conversion/taxonomy_table.tsv"))
contamination_file <- as.character("/home/usuario/Proyectos/Results/ONT16S/data/contaminants.tsv")
# tree_file <- as.character(glue("{BASE_DIR}/output/qiime2/export_results/tree.nwk"))
output_folder <- as.character(glue("{BASE_DIR}/r_tests"))
biosample_col <- "biosample_name"
variables_of_interest <- c("platform", "database", "basecalling_model_name", "sample.type")
apply_controls <- FALSE

######################
# ---- Features ---- #
######################
# Transform biome data to JSON.biome
# biom convert -i in.biom -o out.biom.json --table-type="OTU table" --to-json

# Read
otu_df <- import_biom(feature_table_file)

# Remove columns/samples with all 0s (no data)
condition <- colSums(otu_df) != 0
otu_df_empty_cols <- colnames(otu_df)[!condition]
if (length(otu_df_empty_cols)>0){
  print(glue("WARNING: There are {length(otu_df_empty_cols)} empty columns in ASV dataframe (no data, this can happen in controls) and they will be removed: '{paste(otu_df_empty_cols, collapse='; ')}'"))
}
otu_df <- otu_df[, condition]


######################
# ---- Metadata ---- #
######################
# Read metadata
metadata <- read.delim(metadata_file, comment.char = "#")

# Check if there are columns not in otu_df besides the otu_df columns that were empty
condition <- metadata[,'sample.id'] %in% colnames(otu_df)
sample.ids.not.in.otus <- metadata[!condition,]$sample.id
if (length(setdiff(otu_df_empty_cols, sample.ids.not.in.otus))>0){
  warning("There are samples in metadata that do not exist in ASV dataframe. This can happen if they were removed in the previous step.")
  print(setdiff(otu_df_empty_cols, sample.ids.not.in.otus))
}
metadata <- metadata[condition,] # eliminar muestras si no existen en otu_df
rownames(metadata) <- metadata[,'sample.id']

# Create column for grouping by variables_of_interest
if ("grouping_var" %in% names(metadata)){
  stop("Metadata already has grouping_var as column")
} else {
  metadata$grouping_var <- do.call(paste, c(metadata[variables_of_interest], sep=" - "))
}

# Make categorical columns
for (i in variables_of_interest){
  metadata[[i]] <- factor(metadata[[i]])
}

######################
# ---- Taxonomy ---- #
######################
# Taxonomy must be split into columns
taxmat <- read.delim(taxmat_file, sep='\t')
rownames(taxmat) <- taxmat$feature_id
taxmat$feature_id <- NULL

# Propagate NAs
keep_unidentified <- function(tx_table){
  # Keep the ones that have one of c(NA, "", " ", "\t") in target_level column
  # and then find the last level that is known (can be Genus, Family...)
  # Assign to the target_level level the last known level + x__XXX NA
  # Family (last known)   Species
  # f__Lachnospiraceae    f__Lachnospiraceae NA
  elements_to_target <- colnames(tx_table)
  tx_table <- as.data.frame(tx_table)
  previous_level <- elements_to_target[1]
  for (level in elements_to_target[2:length(elements_to_target)]){
    empty_values <- tx_table[[level]] %in% c(NA, "", " ", "\t")
    tx_table[[level]][empty_values] <- paste(tx_table[[previous_level]][empty_values], "NA", sep="_")
    tx_table[[level]][empty_values] <- gsub("(_NA[> ]*)*", "\\1", tx_table[[level]][empty_values]) # Replace repetitions of NA by one NA
    
    previous_level <- level
  }
  
  return(tx_table)
}
taxmat <- as.matrix(keep_unidentified(taxmat))
colnames(taxmat) <- stringr::str_replace(colnames(taxmat), "^\\w{1}", toupper)


############################
# ---- Clean controls ---- #
############################

substract_control <- function(features, metadata, importance_multiplier=1, control_col="origin", control_applies_to_indicator="control.applies.to.origin"){
  # Substract ASVs that appear in controls per sequencing-run
  # - features: ASV table
  # - metadata: dataframe with metadata
  # - importance_multiplier: increase this to a big number to fully remove ASVs that appear in controls (like 999999)
  # cat("\nUsing controls for cleaning...\n")
  
  # Function to calculate the sum of each ASV
  sum_rows <- function(input){
    # if input is a dataframe do rowSums to get a vector
    # if it is already a vector return it
    if (is.null(ncol(input))){
      return (input)
    } else {
      return (rowSums(input))
    }
  }
  
  # Initialize clean dataframe with non-control ids
  clean.otu.df <- features[, metadata[metadata[control_col]!="control",]$sample.id]
  # clean.otu.df <- features
  
  # Split "control.applies.to.origin" column by ";"
  # The resulting elements will be searched in "origin" column 
  # and the control will be applied to corresponding samples
  metadata$control.applies.to_split <- strsplit(metadata[,control_applies_to_indicator], ";")
  
  # For each control ID
  for (control_id in metadata[metadata[control_col]=="control",]$sample.id ){
    # Get the metadata of the control ID
    r <- subset(metadata, sample.id==control_id)
    
    # Log
    control_id_seq_run = r$sequencing.run
    # print(glue("\n==== Using sample.id '{control_id}' as control in sequencing.run='{control_id_seq_run}' ===="))
    
    # Check that if element "all" is present no other elements exist
    # and that an element is not repeated with unique
    applies_to_elements <- unique(r$control.applies.to_split[[1]]) # To access split elements: metadata$control.applies.to_split[65][[1]][1]
    if ("all" %in% applies_to_elements & length(applies_to_elements)>1){
      stop("Element 'all' was specified, but extra elements are present. Please choose either 'all' or the other elements")
    }
    
    # For each origin to which the control applies to
    for (applies_to in applies_to_elements){ 
      # print(glue("\n~~~ Applying control to: '{applies_to}' ~~~"))
      
      tmp <- metadata[metadata[control_col]!="control" & metadata["sequencing.run"]==control_id_seq_run,]
      
      # If the control applies to all samples
      if (applies_to == "all"){
        tmp <- tmp
        # If it only applies to a certain origin type which IS in metadata$origin
      } else if (applies_to %in% metadata[control_col]){
        tmp <- tmp[tmp[control_col]==applies_to]
        # If it is not "all" or the element doesnt exist in metadata${control_col}, raise error 
      } else {
        warning(glue("Origin '{applies_to}' is not in metadata${control_col} column, please fix"))
      }
      
      if (length(tmp$sample.id) >= 1){
        # print(glue("Applying control.id={control_id} from sequencing.run={control_id_seq_run} to element={applies_to} in {length(tmp$sample.id)} samples"), "\n")
        clean.otu.df[,tmp$sample.id] <- clean.otu.df[,tmp$sample.id] - sum_rows(features[,c(control_id)])*importance_multiplier      
      } else {
        stop(glue("WARNING: Control control.id={control_id} from sequencing.run={control_id_seq_run} cant be applied to element={applies_to}, no samples match ({length(tmp$sample.id)} samples)\nThis in unexpected, every control should apply to at least one sample! Check metadata?"))
      }
    }
  }
  
  # Negative numbers to 0
  clean.otu.df[clean.otu.df<0] <- 0
  
  return(clean.otu.df)
}


if (apply_controls==TRUE){
  clean.otu.df <- substract_control(otu_df, metadata)
} else {
  clean.otu.df <- otu_df
}

####################################
# ---- Create phyloseq object ---- #
####################################
OTU = otu_table(clean.otu.df, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
SAM = sample_data(metadata)

physeq = phyloseq(OTU, TAX, SAM)
physeq

physeq = subset_samples(physeq, origin=="faeces" & sample.type!="control" & exclusion!="excluded")
sample_data(physeq)[c("region_16s","sample.type.label","database.label")] <- (data.frame(sample_data(physeq)) %>% dplyr::mutate(
  region_16s = case_when(
    platform == "Illumina" ~ "Illumina-V3V4",
    platform == "Oxford_Nanopore" ~ "ONT-V1V9",
    TRUE ~ ""
  ),
  sample.type.label = case_when(
    sample.type == "crc" ~ "Cancer",
    sample.type == "non-crc" ~ "Control",
    TRUE ~ ""
  ),
  database.label = case_when(
    database == "silva" ~ "SILVA", # SILVA
    database == "default" ~ "Default", # "rrnDB+NCBI" "EmuDefault"
    TRUE ~ ""
  )
))[c("region_16s","sample.type.label","database.label")]


# Fast, HAC, SUP



##################################################
# ---- Remove typical contamination by name ---- #
##################################################
pop_taxa = function(ps, bad_asvs){
  all_asvs = taxa_names(ps)
  keep_these <- all_asvs[!(all_asvs %in% bad_asvs)]
  return(prune_taxa(keep_these, ps))
}
# Only bacteria
bad_asvs <- taxa_names(subset_taxa(tax_table(physeq), Domain!="Bacteria"))
physeq <- pop_taxa(physeq, bad_asvs)

# Remove typical contamination
typical_contamination <- read.csv(contamination_file, sep='\t')

cont_genus <- unique(c(
  # Base genus
  typical_contamination$Genus, 
  paste("g__", typical_contamination$Genus, sep=""), 
  # Genus with NA
  paste(paste("g__", typical_contamination$Genus, sep=""), "NA", sep= "_"),
  paste(paste("g__", typical_contamination$Genus, sep=""), "NA", sep= " "),
  paste(typical_contamination$Genus, "NA", sep= "_"),
  paste(typical_contamination$Genus, "NA", sep= " ")
))
physeq <- subset_taxa(physeq, !Genus %in% cont_genus)

# Manual names
temp <- c("Dermacoccaceae",
          "Mitochondria",
          "Chloroplast")
if (any(temp %in% data.frame(tax_table(physeq))$Family)){
  bad_asvs <- taxa_names(
    subset_taxa(
      tax_table(physeq),
      Family %in% temp #|
        # Genus %in% c(
        #   "g__Burkholderia-Caballeronia-Paraburkholderia",
        #   "g__Methylobacterium-Methylorubrum",
        #   "g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
        #   "g__Acinetobacter"
        # )
    )
  )
  physeq <- pop_taxa(physeq, bad_asvs)
}

# Keep ASVs with more than one count
physeq <- prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq

# # # Remove samples that could be wrong
# # physeq <- subset_samples(physeq,
# #                          !sample.id %in% c(
# #                            "1A"
# #                          )
# # )
# 

# # Truncate features per sample with less than N% relative abundance to 0
# temp <- otu_table(physeq)
# temp[apply(temp, 2, function(col){100*col/sum(col)<0.01})] <- 0
# otu_table(physeq) <- temp

################################################
# ---- Calculate rarefied phyloseq object ---- #
################################################
rar_physeq <- rarefy_even_depth(physeq, sample.size = 15000,
                                replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
rar_physeq


####################
# ---- Export ---- #
####################
save.physeq <- function(physeq_obj, output_file, add_taxonomy=TRUE){
  # Keep taxa with at least one count (not all zeroes)
  physeq_obj <- prune_taxa(taxa_sums(physeq_obj) >= 1, physeq_obj)
  
  # Merge taxonomy with ASVs
  tax_cols <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  merged.taxa.otus <- merge(otu_table(physeq_obj), tax_table(physeq_obj)[,tax_cols], by.x = 0, by.y = 0, all.x = TRUE, all.y = TRUE)
  names(merged.taxa.otus)[1] <- "feature_id"
  
  # Merge tax into one column
  merged.taxa.otus$taxonomy <- do.call(paste, c(merged.taxa.otus[tax_cols], sep=";"))
  for (c in tax_cols) {merged.taxa.otus[c] <- NULL}
  
  # Make feature_id into row index
  rownames(merged.taxa.otus) <- merged.taxa.otus$feature_id
  
  # Count total
  # merged.taxa.otus$TOTAL <- rowSums(merged.taxa.otus[,-which(colnames(merged.taxa.otus) %in% c("taxonomy","feature_id"))])
  
  # Reorder columns
  merged.taxa.otus <- merged.taxa.otus[,c("feature_id","taxonomy", 
                                          colnames(merged.taxa.otus)[!colnames(merged.taxa.otus) %in% c("taxonomy","feature_id")]
  )
  ]
  # merged.taxa.otus <- merged.taxa.otus[order(-merged.taxa.otus$TOTAL),]
  if (add_taxonomy==FALSE){
    merged.taxa.otus$taxonomy <- NULL
  }
  
  # Save
  write.table(merged.taxa.otus, file=output_file, sep='\t', row.names=F)
}
temp <- subset_samples(physeq, platform=="Oxford_Nanopore" & basecalling_model_name == "sup")
temp <- prune_taxa(taxa_sums(temp) >= 1, temp)
save.physeq(temp, glue("{output_folder}/clean_feature_table.tsv"), add_taxonomy=TRUE)
write.table(data.frame(sample_data(temp)), glue("{output_folder}/clean_feature_table_metadata.tsv"), sep="\t", row.names=F) 



# save.physeq(physeq, glue("{output_folder}/clean_feature_table.tsv"), add_taxonomy=FALSE)
# save.physeq(rar_physeq, glue("{output_folder}/rarified_feature_table.tsv"), add_taxonomy=FALSE)



#############################
# ---- Alpha-diversity ---- #
#############################
plot_alpha_diversity <- function(physeq_obj, measures, x, color_col){
  # ---- Plot alpha diversity ----
  alpha_div_plot <- plot_richness(physeq_obj, x=x, color=color_col, measures=measures) + 
    geom_boxplot() + theme_bw() +
    # facet_grid(variable~origin, scale="free") +
    # scale_x_discrete(limits = c("saliva__non-crc", "saliva__crc", "subgingival-fluid__crc", "faeces__non-crc", "faeces__crc", "normal-mucosa__crc", "adenocarcinoma__crc")) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 90, hjust=0.95, vjust=0.2),
          axis.text = element_text(size = 12))
  
  # ---- Wilcoxon-test ----
  richness <- estimate_richness(physeq_obj, measures=measures)
  rownames(richness) <- sample_data(physeq_obj)$sample.id
  # rownames(richness) <- gsub("\\.", '-', rownames(richness))
  richness <- merge(richness, data.frame(sample_data(physeq_obj)), by=0)
  richness$group <- as.factor(richness$group)
  
  # Run test for each combination of groups
  wilcox_list <- list()
  n = 1
  for (m in measures){
    for (i in sort(unique(richness$group))){
      for (j in sort(unique(richness$group))){
        res <- NULL
        # As long as the group is not the same (error if it is)
        if (i!=j){
          # Run the test between groups and print if significative
          res <- wilcox.test(as.formula(glue("{m} ~ group")), 
                             data=subset(richness, 
                                         group %in% c(i, j)
                             )
          )
          wilcox_list[[n]] <- c(i, j, m, res$p.value)
        } else {
          wilcox_list[[n]] <- c(i, j, m, NA)
        }
        n=n+1
      }
    }  
  }
  
  # Format Wilcox results
  wilcox_df <- as.data.frame(do.call(rbind, wilcox_list))
  names(wilcox_df) <- c("C1","C2","measure","pval")
  wilcox_df$pval <- as.numeric(wilcox_df$pval)
  
  # Plot each wilcoxon test
  significance_p_list <- list()
  for (m in measures){
    df <- subset(wilcox_df, measure==m)[c("C1","C2","pval")]
    
    # Turn to wide and block lower triangle
    df <- as.data.frame(df %>% tidyr::pivot_wider(names_from = C1, values_from = pval))
    df[-1][lower.tri(df[-1])] <- NA
    
    # Turn to long
    df <- df %>% pivot_longer(!C2, names_to = "C1", values_to = "pval")
    
    p <- ggplot(df, aes(C1, C2, fill=pval)) + 
      geom_tile() +
      scale_fill_gradient(low="darkgreen", high="white", na.value = "white", limits=c(0,0.1)) +
      # scale_y_discrete(position = "right") +
      theme_bw() + 
      ylab("") + ggtitle(glue("{m}")) + xlab("") + 
      theme(axis.text.x=element_text(angle = 60, hjust=0.95, vjust=0.95), legend.key.size = unit(0.5,"line"))
    significance_p_list[[m]] <- p
  }
  
  return(list(alpha_div_plot=alpha_div_plot, significance_p_list=significance_p_list, richness=richness[,c("sample.id",measures)]))
  
}

measures <- c("ACE","Shannon","Simpson","Observed","Chao1")
ps <- tax_glom(subset_samples(physeq), "Species")
sample_data(ps)$group <- (data.frame(sample_data(ps)) %>% unite("temp", all_of(variables_of_interest), sep=" - "))$temp

res = plot_alpha_diversity(ps, measures, "group", "sample.type")

temp <- res$richness %>% 
  tidyr::pivot_longer(!sample.id, names_to = "alphadiversity_measure", values_to = "alphadiversity_value")
temp <- merge(temp, data.frame(sample_data(ps)), on="sample.id")

model_alphadiversity_p <- ggplot(
    subset(temp, alphadiversity_measure %in% c("Observed","Shannon","Simpson") & platform=="Oxford_Nanopore") %>%
      dplyr::mutate(x_axis=glue("{sample.type.label}\n({database.label})")), 
    aes(x=x_axis, alphadiversity_value, fill=basecalling_model_name)
  ) +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_point(aes(color=basecalling_model_name), position=position_jitterdodge(jitter.width=0.2), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
  facet_wrap(~alphadiversity_measure, scales="free") +
  theme_bw() +
  xlab(NULL) +
  ylab("Alpha-diversity") +
  # ggtitle("") +
  guides(fill=guide_legend(title="Basecalling model")) +
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10, face="bold"),
    axis.text.x=element_text(angle = 0, hjust=0.5, vjust=0.5, size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold")
  ) + 
  ggpubr::geom_pwc(
    aes(group = basecalling_model_name), step.increase=0.05, vjust = 0.5, tip.length = 0.01,
    method = "wilcox_test", label = "{p.adj.signif}",
    p.adjust.method = "holm", hide.ns = TRUE
  ) +
  ggpubr::geom_pwc(
    aes(group = x_axis), step.increase=0.05, vjust = 0.5, size = 0.75, tip.length = 0,
    method = "wilcox_test", label = "{p.adj.signif}",
    p.adjust.method = "holm", hide.ns = TRUE,
    bracket.nudge.y = 0.2
  )

model_alphadiversity_p
ggsave(glue("{output_folder}/alphadiversity_ont.png"), plot=model_alphadiversity_p, dpi="retina", unit="px", height=2000, width=3000)

platform_alphadiversity_p <- ggplot(
    subset(temp, alphadiversity_measure %in% c("Observed","Shannon","Simpson") & database == "silva" & (basecalling_model_name == "sup" | platform=="Illumina")), 
    aes(x=sample.type.label, y=alphadiversity_value, fill=region_16s)
  ) +
  geom_boxplot(alpha=0.5, outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.1), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
  facet_wrap(~alphadiversity_measure, scales="free", nrow=1) +
  ggpubr::geom_pwc(
    aes(group = region_16s), step.increase=0.05, vjust = 0.5, tip.length = 0.01,
    method = "wilcox_test", label = "{p.adj.signif}",
    p.adjust.method = "holm", hide.ns = TRUE
  ) +
  theme_bw() +
  xlab(NULL) +
  ylab("Alpha-diversity") +
  # ggtitle("") +
  guides(fill=guide_legend(title="Method   ")) +
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10, face="bold"),
    axis.text.x=element_text(angle = 0, hjust=0.5, vjust=0.5, size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold")
  )
platform_alphadiversity_p
ggsave(glue("{output_folder}/alphadiversity_platforms.png"), plot=platform_alphadiversity_p, dpi="retina", unit="px", height=2000, width=2500)



############################
# ---- Beta-diversity ---- #
############################

# -- Beta-diveristy comparisons for ONT basecalling models (silva + default) --
rar_ps = tax_glom(subset_samples(rar_physeq, platform == "Oxford_Nanopore"), "Species")
dist <- ordinate(rar_ps, "MDS", "jsd", formula=as.formula(glue("~ {paste(variables_of_interest, collapse = ' + ')}")))
sample_data(rar_ps)$temp <- paste(sample_data(rar_ps)$biosample_name, sample_data(rar_ps)$database.label, sep="_")
sample_data(rar_ps)$model_database <- paste(sample_data(rar_ps)$basecalling_model_name, sample_data(rar_ps)$database.label, sep=" - ")

basecalling_model_silvaanddefault_jsd <- plot_ordination(rar_ps, dist,
                                         type="sites", color="model_database") + 
  geom_path(aes(group=temp), color="black", alpha=0.7, linetype=1, linewidth=0.5) +
  guides(fill=guide_legend(title="Basecalling model")) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=10),
    plot.title=element_text(face="bold", size=17)
  )
basecalling_model_silvaanddefault_jsd

# -- Beta-diveristy comparisons for ONT basecalling models (silva) --
rar_ps = tax_glom(subset_samples(rar_physeq, database == "silva" & platform == "Oxford_Nanopore"), "Species")
dist <- ordinate(rar_ps, "MDS", "jsd", formula=as.formula(glue("~ {paste(variables_of_interest, collapse = ' + ')}")))

basecalling_model_silva_jsd <- plot_ordination(rar_ps, dist,
                                         type="sites", color="basecalling_model_name") + 
  geom_path(aes(group=biosample_name), color="black", alpha=0.7, linetype=1, linewidth=0.5) +
  guides(fill=guide_legend(title="Basecalling model")) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=10),
    plot.title=element_text(face="bold", size=17)
  )
basecalling_model_silva_jsd

# -- Beta-diversity comparison for Illumina vs Nanopore (genus) --
rar_ps = tax_glom(subset_samples(rar_physeq, database == "silva" & (basecalling_model_name == "sup" | platform == "Illumina")), "Genus")
dist <- ordinate(rar_ps, "MDS", "jsd", formula=as.formula(glue("~ {paste(variables_of_interest, collapse = ' + ')}")))

platform_genus_silva_jsd <- plot_ordination(rar_ps, dist,
                                         type="sites", color="region_16s") + 
  stat_ellipse() +
  guides(fill=guide_legend(title="Method")) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=10),
    plot.title=element_text(face="bold", size=17)
  )
platform_genus_silva_jsd

# -- Beta-diversity comparison for Illumina vs Nanopore (species) --
rar_ps = tax_glom(subset_samples(rar_physeq, database == "silva" & (basecalling_model_name == "sup" | platform == "Illumina")), "Species")
dist <- ordinate(rar_ps, "MDS", "jsd", formula=as.formula(glue("~ {paste(variables_of_interest, collapse = ' + ')}")))

platform_species_silva_jsd <- plot_ordination(rar_ps, dist,
                                      type="sites", color="region_16s") + 
  stat_ellipse() +
  guides(fill=guide_legend(title="Method")) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10, face="bold"),
    axis.text.x=element_text(size=10),
    plot.title=element_text(face="bold", size=17)
  )
platform_species_silva_jsd


#####################################
# ---- PERMANOVA using adonis2 ---- #
#####################################
# Filter
rar_ps = tax_glom(subset_samples(rar_physeq), "Genus")

# Matriz de distancias
rar_physeq_dist <- phyloseq::distance(rar_ps, "jsd")
sampledf <- data.frame(sample_data(rar_ps))

# Adonis test
permanova <- adonis2(
  as.formula(paste("rar_physeq_dist ~ ", paste(variables_of_interest, collapse="+"), sep = "")), # rar_physeq_dist ~ origin + sample.type,
  data = sampledf, permutations = 9999, by="onedf",
  parallel = 60)

temp <- data.frame(permanova)
write.table(data.frame(variable = row.names(temp), temp), glue("{output_folder}/permanova.tsv"), sep="\t", row.names=F) 

###########################
# ---- CLR abundance ---- #
###########################


# Functions
plot_boxes <- function(data, x, y, facet="target_level", fill=NULL, title=NULL, boxplot.outlier.shape = 19){
  ggplot(data, aes_string(x=x, y=y, fill=fill)) + 
    facet_wrap(as.formula(glue("~{facet}")), scales="free") +
    geom_boxplot(alpha=0.5, outlier.shape=boxplot.outlier.shape) + 
    # geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, alpha=0.5) +
    theme_bw() +
    xlab(NULL) +
    ylab("Abundance (CLR)") +
    ggtitle(title) +
    guides(fill=guide_legend(title="")) +
    theme(
      axis.text=element_text(size=10, color="black"), 
      axis.title=element_text(size=10, face="bold"),
      axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
      legend.text=element_text(size=10),
      axis.text.x=element_text(angle = 0, hjust=0.5, vjust=0.5, size=10),
      plot.title = element_text(size = 15, color = 'black', face="bold"),
      strip.background = element_rect(fill="white", color="black"),
      strip.text = element_text(colour = 'black', face="bold"),
      strip.text.x = element_text(colour = 'black', face="bold.italic"),
    )
}

like.with.vector.filter <- function(data, searches, column_to_search){
  data <- data[
    rowSums(# Get count of times each pattern appears
      sapply( # Search pattern in column
        sapply(searches, function(x){glue(".*{x}.*")}), # create vector of taxa like ".*Fusobacterium.*" #ORAL_GENERA,
        like, vector = data[[column_to_search]]) 
    )>=1, # If a pattern appears, select row
  ]
  return(data)
}

get.clr.abund.at.zero <- function(rel_abund_data, clr_abund_data){
  # Merge CLR and Relative abundance long objects and get the 
  # point at which clr abundance is mostly 0% relative abundance
  
  # Rename columns
  rel_abund_data <- rel_abund_data %>% dplyr::rename(relative_abundance=Abundance)
  clr_abund_data <- clr_abund_data %>% dplyr::rename(clr_abundance=Abundance)
  
  # Merge
  merge_df <- merge(
    rel_abund_data[c("sample.id","grouping_var","OTU","target_level","relative_abundance")], 
    clr_abund_data[c("sample.id","grouping_var","OTU","target_level","clr_abundance")], 
    by=c("sample.id","grouping_var","OTU","target_level")
  )
  
  # Format
  merge_df$relative_abundance <- merge_df$relative_abundance*100
  merge_df$relative_is_zero <- merge_df$relative_abundance==0
  
  # Plot correspondence
  p <- ggplot(data=merge_df, aes_string(x="relative_abundance", y="clr_abundance", color="relative_is_zero")) +
    geom_point() +
    ggtitle("CLR abundance vs Relative abundance") + 
    geom_hline(yintercept = max(subset(merge_df, relative_abundance==0)$clr_abundance), alpha=0.2) +
    theme_bw()
  # ggplotly(p)
  
  # Get point at which clr abundance is mostly 0% relative abundance
  clr.abund.at.zero <- max(subset(merge_df, relative_abundance==0)$clr_abundance) 
  
  return(clr.abund.at.zero)
}

get.long.clr.and.rel <- function(physeq_obj, selected_taxa_level){
  # Normalize
  physeq_obj <- tax_glom(physeq_obj, selected_taxa_level)
  physeq_obj <- prune_taxa(taxa_sums(physeq_obj) >= 1, physeq_obj)
  
  clr_ps <- microbiome::transform(physeq_obj, transform="clr")
  rel_abund_ps <- microbiome::transform(physeq_obj, transform="compositional")
  
  # Agglomerate and turn into long
  clr_ps_long <- clr_ps %>%
    tax_glom(taxrank = selected_taxa_level, bad_empty = c()) %>%   # Aglomerate to taxnomic range
    psmelt()
  clr_ps_long$target_level <- clr_ps_long[[selected_taxa_level]]
  
  rel_abund_ps_long <- rel_abund_ps %>%
    tax_glom(taxrank = selected_taxa_level, bad_empty = c()) %>%   # Aglomerate to taxnomic range
    psmelt()   
  rel_abund_ps_long$target_level <- rel_abund_ps_long[[selected_taxa_level]]
  
  # Format
  clr_ps_long[["target_level"]] <- gsub(".*__", '', clr_ps_long[["target_level"]])
  clr_ps_long[["target_level"]] <- gsub("_", ' ', clr_ps_long[["target_level"]])
  rel_abund_ps_long[["target_level"]] <- gsub(".*__", '', rel_abund_ps_long[["target_level"]])
  rel_abund_ps_long[["target_level"]] <- gsub("_", ' ', rel_abund_ps_long[["target_level"]])
  
  return(list(clr=clr_ps_long, rel_abund=rel_abund_ps_long))
}



# ---- CLR at certain level ----
# For model sup, database silva at Species level
l <- get.long.clr.and.rel(
  subset_samples(physeq, database == "silva" & (basecalling_model_name=="sup" )), #| platform=="Illumina"
  selected_taxa_level <- "Species"
)
clr_ps_long_sp <- l$clr
rel_abund_ps_long_sp <- l$rel_abund

# For ONT sup with default database at species level
l <- get.long.clr.and.rel(
  subset_samples(physeq, database == "default" & (basecalling_model_name=="sup")),
  selected_taxa_level <- "Species"
)
clr_ps_long_default_sp <- l$clr
rel_abund_ps_long_default_sp <- l$rel_abund

# For all models from ONT, database silva at Species level
l <- get.long.clr.and.rel(
  subset_samples(physeq, database == "silva" & platform=="Oxford_Nanopore"),
  selected_taxa_level <- "Species"
)
clr_ps_long_ont_models_sp <- l$clr
rel_abund_ps_long_ont_models_sp <- l$rel_abund

# For all models from ONT, database silva at genus level
l <- get.long.clr.and.rel(
  subset_samples(physeq, database == "silva" & platform=="Oxford_Nanopore"),
  selected_taxa_level <- "Genus"
)
clr_ps_long_ont_models_genus <- l$clr
rel_abund_ps_long_ont_models_genus <- l$rel_abund

# For ONT vs Illumina, database silva at the genus level
l <- get.long.clr.and.rel(
  subset_samples(physeq, database == "silva" & (basecalling_model_name=="sup"| platform=="Illumina")),
  selected_taxa_level <- "Genus"
)
clr_ps_long_genus <- l$clr
rel_abund_ps_long_genus <- l$rel_abund

# For Illumina (silva) and ONT sup (with both databases) at the genus level
# both databases are used at the same time because we only want relative abundance, not CLR
l <- get.long.clr.and.rel(
  subset_samples(physeq, database %in% c("silva","default") & (basecalling_model_name=="sup"| platform=="Illumina")),
  selected_taxa_level <- "Genus"
)
rel_abund_ps_long_databases_genus <- l$rel_abund

# For comparing both platforms and both databases with CLR
temp <- list(
  get.long.clr.and.rel(
    subset_samples(physeq, database == "silva" & basecalling_model_name=="sup"),
    selected_taxa_level <- "Genus"
  ),
  get.long.clr.and.rel(
    subset_samples(physeq, database == "default" & basecalling_model_name=="sup"),
    selected_taxa_level <- "Genus"
  ),
  get.long.clr.and.rel(
    subset_samples(physeq, database == "silva" & platform=="Illumina"),
    selected_taxa_level <- "Genus"
  )
)
clr_ps_long_platform <- dplyr::bind_rows(temp[[1]]$clr, temp[[2]]$clr, temp[[3]]$clr)
rel_abund_ps_long_platform <- dplyr::bind_rows(temp[[1]]$rel_abund, temp[[2]]$rel_abund, temp[[3]]$rel_abund)




# ---- Plot specific genera, Illumina vs ONT ----
clr.abund.at.zero <- -0.4

clr_ps_long_platform <- clr_ps_long_platform %>%
  dplyr::mutate(
    temp_col = factor(
      glue("{region_16s}\n{sample.type.label}\n{database.label}"), 
      levels=c(
        "Illumina-V3V4\nCancer\nSILVA", "ONT-V1V9\nCancer\nSILVA", "ONT-V1V9\nCancer\nDefault", 
        "Illumina-V3V4\nControl\nSILVA", "ONT-V1V9\nControl\nSILVA", "ONT-V1V9\nControl\nDefault"
      )
    )
  )
temp <- subset(clr_ps_long_platform, target_level %in% c("Parvimonas", "Fusobacterium", "Peptostreptococcus"))
clr_genus_platform <- plot_boxes(
  temp,
  x="temp_col", y="Abundance", facet="target_level", fill="sample.type.label",
  boxplot.outlier.shape = NA # outliers drawn by next geom_point jitter
) +
  # geom_path(aes(group=biosample_name), color="black", alpha=0.6, linetype=2, linewidth=0.2) +
  geom_point(position=position_jitterdodge(jitter.width=0.5), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
  geom_hline(yintercept=clr.abund.at.zero, color="black", alpha=0.5) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
    ymax = clr.abund.at.zero,
    fill = "dodgerblue4", alpha = 0.1, color = NA
  ) +
  geom_signif(
    comparisons = list(
      c("ONT-V1V9\nCancer\nDefault", "ONT-V1V9\nControl\nDefault"),
      c("ONT-V1V9\nCancer\nSILVA", "ONT-V1V9\nControl\nSILVA"),
      c("Illumina-V3V4\nCancer\nSILVA",  "Illumina-V3V4\nControl\nSILVA")
    ),
    map_signif_level=map_signif_level_func,
    step_increase=0.05,
    textsize = 3,
    tip_length = 0.01,
    margin_top = 0.1,
    vjust = 0.5
  )
clr_genus_platform



# ---- Plot specific species, only ONT with both databases ----
# Note: Two different objects, one with database=silva and another with database=default
# are calculated separatedly, as the CLR normalization seems to be influenced by the amount of 
# feature_ids with values equal to 0, and because databases were different the taxonomy ids were also different

clr.abund.at.zero <- mean(
  get.clr.abund.at.zero(rel_abund_ps_long_sp, clr_ps_long_sp),
  get.clr.abund.at.zero(rel_abund_ps_long_default_sp, clr_ps_long_default_sp)
)

temp <- subset(
  dplyr::bind_rows(clr_ps_long_sp, clr_ps_long_default_sp) %>%
    dplyr::mutate(x_label = glue("{sample.type.label}\n{database.label}")) %>%
    dplyr::mutate(x_label = factor(x_label, levels=c("Cancer\nDefault", "Control\nDefault", "Cancer\nSILVA","Control\nSILVA")))  %>%
    dplyr::mutate(
      target_level = dplyr::case_when(    
        target_level == "Clostridium sensu stricto 1 perfringens" ~ "Clostridium perfringens",
        target_level == "Clostridium sensu stricto 1 saudiense" ~ "Clostridium saudiense",
        TRUE ~ target_level
      )
    ),
  target_level %in% c(
    # Species
    "Parvimonas micra",
    "Fusobacterium nucleatum", "Fusobacterium gonidiaformans", "Fusobacterium necrophorum",
    "Bacteroides fragilis",
    "Peptostreptococcus stomatis",
    "Peptostreptococcus anaerobius", "Sutterella wadsworthensis",
    "Paraprevotella clara", "Dialister pneumosintes",
    "Clostridium perfringens",
    "Gemella morbillorum",
    "Clostridium saudiense",
    "Ruficoccus amylovorans",
    # Probable good species
    "Agathobaculum butyriciproducens",
    "Romboutsia ilealis",
    "Anaerostipes rhamnosivorans",
    "Anaerocolumna cellulosilytica"
  )
)
clr_sp_ont_both_p <- plot_boxes(
  temp,
  x="x_label", y="Abundance", facet="target_level", fill="x_label", boxplot.outlier.shape=NA # outliers are drawn by geom_point with jitter
) +
  facet_wrap(~target_level, scales="free_y", nrow=3) +
  geom_point(position=position_jitterdodge(jitter.width=1), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
  geom_hline(yintercept=clr.abund.at.zero, color="black", alpha=0.5) +
  scale_fill_brewer(palette = "Dark2") +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
    ymax = clr.abund.at.zero,
    fill = "dodgerblue4", alpha = 0.1, color = NA
  ) + 
  ggpubr::geom_pwc(
    aes(group = x_label), step.increase=0.05, vjust = 0.5, tip.length = 0.01,
    method = "wilcox_test", label = "{p.adj.signif}",
    p.adjust.method = "holm", hide.ns = TRUE
  )
  # geom_signif(
  #   comparisons = list(
  #     c("Cancer\nDefault","Control\nDefault"),
  #     c("Cancer\nSILVA","Control\nSILVA")
  #   ),
  #   map_signif_level=map_signif_level_func,#map_signif_level_func,
  #   step_increase=0, #0.05,
  #   textsize = 3,
  #   tip_length = 0.01,
  #   margin_top = 0.1,
  #   vjust = 0.5
  # )
clr_sp_ont_both_p



# ---- Plot all Interesting stuff ----
# GENERA <- c(
#   "Fusobacterium", "Bacteroides", "Parvimonas",
#   "Porphyromonas", "Peptostreptococcus", "Faecalibacterium",
#   "Blautia", "Prevotella", "Sutterella", "Roumboutsia", "Longibaculum",
#   "Limosilactobacillus", "Gemella", "Eisenbergiella", "Coraliomargarita",
#   "Anaerostipes", "Anaerocolumna", "Agathobaculum", "Subdoligranulum",
#   "Ruthenibacterium", "Butyricicoccus", "Anaerotaenia", "Anaerocolumna",
#   "Dialister", "Clostridium", "Coprobacter", "Klebsiella", "Paludicola", "Lactobacillus", "Paraprevotella",
#   "Ruficoccus", "Raoultella", "Romboutsia"
# )
# SPECIES <- c(
#   # Species
#   "Parvimonas micra", "Parvimonas NA",
#   "Fusobacterium gonidiaformans", "Fusobacterium necrophorum",
#   "Fusobacterium mortiferum", "Fusobacterium nucleatum", "Fusobacterium varium",
#   "Fusobacterium necrogenes", "Fusobacterium periodonticum", "Fusobacterium ulcerans",
#   "Fusobacterium equinum",
#   "Bacteroides fragilis",
#   "Peptostreptococcus NA", "Peptostreptococcus stomatis",
#   "Peptostreptococcus anaerobius", "Sutterella wadsworthensis",
#   "Paraprevotella clara", "Dialister pneumosintes",
#   "Clostridium sensu stricto 1 perfringens",
#   "Clostridium sensu stricto 1 saudiense",
#   "Coprobacter secundus", "Gemella morbillorum"
# )

# rel_abund_ps_long_sp
# clr_ps_long_sp

wilcoxon.res <- clr_ps_long_default_sp %>% 
  dplyr::group_by(platform, basecalling_model_name, database, target_level) %>%
  rstatix::wilcox_test(Abundance ~ sample.type.label) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  # do(
  #   w = wilcox.test(
  #     x=subset(., sample.type.label=="Cancer")$Abundance,
  #     y=subset(., sample.type.label=="Control")$Abundance
  #   )
  # ) %>% 
  # dplyr::summarise(platform, basecalling_model_name, database, target_level, Wilcox = w$p.value) %>%
  dplyr::filter(p <= 0.1)

clr_ps_long_subset <- like.with.vector.filter(clr_ps_long_default_sp, c(unique(wilcoxon.res$target_level)), "target_level")
rel_abund_ps_long_subset <- like.with.vector.filter(rel_abund_ps_long_default_sp, c(unique(wilcoxon.res$target_level)), "target_level")



# ---- Plot into n split columns ----
# Select taxa
target.taxa <- sort(unique(clr_ps_long_subset$target_level))

# Check if any of them dont exist in data
target.taxa[!target.taxa %in% clr_ps_long_subset$target_level]

# Get point under which CLR abundance is mostly equivalent to 0% relative abundance
clr.abund.at.zero <- get.clr.abund.at.zero(rel_abund_ps_long_subset, clr_ps_long_subset)

# Plot
targets.groups <- split(target.taxa, ceiling(seq_along(target.taxa)/40))
p_list <- list()
n = 1
for (g in targets.groups){
  p <- plot_boxes(
      subset(
        clr_ps_long_subset,target_level %in% g),
      x="sample.type.label", y="Abundance", facet="Species", fill="sample.type.label", boxplot.outlier.shape = NA # outliers are placed by next geom_point jitter
    ) +
    geom_point(position=position_jitterdodge(jitter.width=0.5), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
    # facet_grid(target_level~platform, scales="free") +
    # geom_vline(xintercept = c(4.5, 8.5, 12.5), alpha=0.6, linetype="dotted") +
    geom_hline(yintercept=clr.abund.at.zero, color="black", alpha=0.5) +
    annotate(
      "rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
      ymax = clr.abund.at.zero,
      fill = "dodgerblue4", alpha = 0.1, color = NA
    ) + geom_signif(
      comparisons = list(
        c("Control","Cancer")
      ),
      map_signif_level=map_signif_level_func,
      step_increase=0.05,
      textsize = 3,
      tip_length = 0.01,
      margin_top = 0.1,
      vjust = 0.5
    )
  # p
  ggsave(
    glue("{output_folder}/all_default_p{n}.png"), 
    plot=p, dpi="retina", units="px", height=4000, width=6000, limitsize = FALSE
  )
  
  p_list[[n]] <- p
  n = n+1
}
# 
# p <- ggpubr::ggarrange(
#   p_list[[1]],
#   p_list[[2]],
#   p_list[[3]],
#   p_list[[4]],
#   p_list[[5]],
#   p_list[[6]],
#   # p_list[[7]],
#   # p_list[[8]],
#   # p_list[[9]],
#   # p_list[[10]],
#   # p_list[[11]],
#   # p_list[[12]],
#   # p_list[[13]],
#   # p_list[[14]],
#   # p_list[[15]],
#   # p_list[[16]],
#   # p_list[[17]],
#   # p_list[[18]],
#   # p_list[[19]],
#   ncol = length(p_list),
#   common.legend = TRUE, legend = "bottom"
#   # labels = c("A", "B"), font.label = list(size = 20)
# )
# ggsave(
#   glue("{output_folder}/test.png"), #new_CLR_{selected_taxa_level}_interesting.png
#   plot=p, dpi="retina", units="px", height=4000*3, width=6000*7, limitsize = FALSE
# )


####################################
# ---- Differential abundance ---- #
####################################

# ---- Functions ----
library(ANCOMBC)
run_ancombc_1vs1 <- function(physeq_obj, target_level, target_variables, prv_cut=0.3, qval_min=0.1){
  # IMPORTANT: ANCOM-BC v2.0.1 was used. Its output format varies depending on version
  if (target_level=="ASV"){
    ancom_target <- NULL
  } else {
    ancom_target <- target_level
  }

  # Run ANCOM-BC
  out <- ancombc(data=physeq_obj, tax_level=ancom_target, prv_cut=prv_cut,
                 formula = paste(target_variables, collapse="+"), n_cl=60,
                 p_adj_method="holm", alpha = 0.05)
  # Get results
  res = out$res

  # To long
  lfc_long <- res$lfc[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "lfc")
  diff_abn_long <- res$diff_abn[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "diff_abn")
  se_long <- res$se[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "se")
  qval_long <- res$q_val[,c(-2)] %>% tidyr::pivot_longer(cols = starts_with(target_variables), names_to = "var", values_to = "qval")

  # Merge long
  merged <- Reduce(function(x, y) merge(x, y, by=c("taxon","var")),
                   list(lfc_long, diff_abn_long, se_long, qval_long))

  # Get taxons where any of the tests achived a certain qval
  selected_taxons <- res$q_val[rowSums(res$q_val[,c(-1,-2),drop=F]<=qval_min)!=0,]$taxon
  merged <- subset(merged, taxon %in% selected_taxons)
  merged

  # Change column name
  names(merged)[1] <- "taxon_id"

  # Add names if ASV was selected as target_level
  if (target_level=="ASV"){
    merged$taxon_id <- data.frame(tax_table(physeq_obj))[merged$taxon_id,][,"Species"]
  }

  # Format taxon_id to remove things like "g__"
  # merged$taxon_id <- gsub(".*__", '', merged$taxon_id)

  # Format results for plotting
  df_fig <- merged %>%
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = ifelse(lfc > 0, "Positive LFC", "Negative LFC"))
  df_fig$direction <- factor(df_fig$direction, levels = c("Positive LFC", "Negative LFC"))

  # Significance label
  # - Corrected p-value (qval) under 0.001 will have "***"
  # - Corrected p-value (qval) under 0.01 will have "**"
  # - Corrected p-value (qval) under 0.05 will have "*"
  df_fig$qval.txt <- cut(df_fig$qval, c(-Inf,0.001,0.01,0.05,Inf), labels=c('***','**','*','-'))

  # Orientation for significance level
  df_fig$orientation <- 0
  df_fig$orientation[df_fig$lfc > 0] <- 1
  df_fig$orientation[df_fig$lfc < 0] <- -1

  # Add target level
  df_fig$target_level <- target_level

  return(df_fig)
}

# Function to get the proper alpha value for each significance level
# (non-significant entries will be grayed out)
get_scale_alpha_values <- function(vec){
  l <- length(unique(vec))
  alpha_vec <- rep(1, l)
  if ("-" %in% vec){
    alpha_vec[length(alpha_vec)] <- 0.4
  }
  return(alpha_vec)
}

plot.ancombc <- function(data, fill="direction"){
  p <- ggplot(
    data = subset(
      data %>% 
        dplyr::arrange(lfc) %>% 
        dplyr::mutate(taxon_id=factor(taxon_id, levels=unique(.$taxon_id)))
      # dplyr::filter(grepl('Parvimonas|Fusobacterium|Peptostreptococcus', taxon_id))
    ), 
    aes_string(x = "taxon_id", y = "lfc", fill=fill, alpha="qval.txt")
  ) +
    scale_alpha_manual(values=get_scale_alpha_values(data$qval.txt)) +
    geom_bar(stat = "identity", color="black", linewidth=0.3,
             position = position_dodge2(width = 0.9, preserve = "single")) +
    geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width=0.5,
                  position = position_dodge2(width = 0.9, preserve = "single"),
                  color = "gray36") +
    geom_text(aes(label = qval.txt, y=lfc+orientation*(se+0.22)), vjust = 0.7, color = "gray36",
              position = position_dodge2(width = 0.9, preserve = "single")) +
    # geom_vline(xintercept = (1:length(data$taxon_id)-1)+0.5, alpha=0.5, linetype="dotted") +
    geom_hline(yintercept = 0, alpha=0.5) +
    facet_wrap(~target_level, scale="free_y", ncol=1) +
    labs(x = NULL, y = "Log fold change", title = "") +
    coord_flip() +
    guides(alpha = "none") +
    theme_bw() +
    theme(
      axis.text=element_text(size=10, color="black"), 
      axis.text.x=element_text(angle = 0, hjust=0.95, vjust=0.95, size=10),
      axis.text.y=element_text(size=10, face="italic"),
      axis.title=element_text(size=10, face = "bold"),
      legend.text=element_text(size=8),
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill="white", color="black"),
      strip.text = element_text(colour = 'black', face="bold"),
      strip.text.y = element_text(colour = 'black', face="bold.italic")
    )
  return(p)
} 

ancombc_prv_cut <- 0.1

# ~~~~~~~~~~~~~~~~~~~~~~~~ #
# ---- ONT default db ---- #
# ~~~~~~~~~~~~~~~~~~~~~~~~ #

# ---- Group by and dereplicate ----
# Filter phyloseq object
derep_ps <- subset_samples(physeq, sample.type %in% c("crc","non-crc") & database=="default" & (basecalling_model_name == "sup"))
derep_ps <- prune_taxa(taxa_sums(derep_ps) >= 1, derep_ps)
sample_data(derep_ps)$sample.type.label <- factor(sample_data(derep_ps)$sample.type.label, levels=c("Control", "Cancer"))

# ---- Running ANCOM-BC ----
# Run on each level
proc_genus_ancombc <- NULL
proc_species_ancombc <- NULL

# proc_genus_ancombc <- run_ancombc_1vs1(derep_ps, "Genus", c("sample.type.label"), prv_cut=ancombc_prv_cut)
proc_species_ancombc <- run_ancombc_1vs1(derep_ps, "Species", c("sample.type.label"), prv_cut=ancombc_prv_cut)

# Concatenate results
ont.default.df_fig <- rbind(proc_genus_ancombc, proc_species_ancombc)
ont.default.df_fig$target_level <- factor(ont.default.df_fig$target_level,      # Reordering group factor levels
                              levels = c("Phylum","Class","Order","Family","Genus","Species","ASV")
)

# ~~~~~~~~~~~~~~~~~~~~~~ #
# ---- ONT silva db ---- #
# ~~~~~~~~~~~~~~~~~~~~~~ #

# ---- Group by and dereplicate ----
# Filter phyloseq object
derep_ps <- subset_samples(physeq, sample.type %in% c("crc","non-crc") & database=="silva" & (basecalling_model_name == "sup"))
derep_ps <- prune_taxa(taxa_sums(derep_ps) >= 1, derep_ps)
sample_data(derep_ps)$sample.type.label <- factor(sample_data(derep_ps)$sample.type.label, levels=c("Control", "Cancer"))

# ---- Running ANCOM-BC ----
# Run on each level
proc_genus_ancombc <- NULL
proc_species_ancombc <- NULL

# proc_genus_ancombc <- run_ancombc_1vs1(derep_ps, "Genus", c("sample.type.label"), prv_cut=ancombc_prv_cut)
proc_species_ancombc <- run_ancombc_1vs1(derep_ps, "Species", c("sample.type.label"), prv_cut=ancombc_prv_cut)

# Concatenate results
ont.silva.df_fig <- rbind(proc_genus_ancombc, proc_species_ancombc)
ont.silva.df_fig$target_level <- factor(ont.silva.df_fig$target_level,      # Reordering group factor levels
                              levels = c("Phylum","Class","Order","Family","Genus","Species","ASV")
)




