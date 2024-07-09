

# ---- Read quality stats plot ----
# Read
sup <- data.table::fread(file="/home/usuario/Proyectos/Results/ONT16S/pipeline/09_qualities/sup.tsv.gz")
hac <- data.table::fread(file="/home/usuario/Proyectos/Results/ONT16S/pipeline/09_qualities/hac.tsv.gz")
fast <- data.table::fread(file="/home/usuario/Proyectos/Results/ONT16S/pipeline/09_qualities/fast.tsv.gz")

# Rename
names(sup) <- c("length","avg_quality")
names(hac) <- c("length","avg_quality")
names(fast) <- c("length","avg_quality")

# Filter
sup <- sup %>% dplyr::filter(avg_quality>=12 & length>=1200 & length<=1900)
hac <- hac %>% dplyr::filter(avg_quality>=12 & length>=1200 & length<=1900)
fast <- fast %>% dplyr::filter(avg_quality>=7 & length>=1200 & length<=1900)

# Add basecalling model
sup$basecalling_model_name <- "sup"
hac$basecalling_model_name <- "hac"
fast$basecalling_model_name <- "fast"

# Merge
temp <- rbind(sup, hac, fast)
temp.stats <- temp %>% 
  dplyr::group_by(basecalling_model_name) %>% 
  dplyr::summarise(
    median_length=median(length), 
    median_avg_quality=round(median(avg_quality)), 
    mean_avg_quality=round(mean(avg_quality)), 
    reads=n(),
    reads_m=scales::label_number(accuracy=0.1, scale_cut=scales::cut_short_scale())(n())
  )

# Sample n rows to reduce size of table
temp <- dplyr::sample_n(temp, dim(temp)[1]/10, replace = TRUE)

# Create annotations
annotation <-apply(temp.stats, 1, function(x){glue("- {x['basecalling_model_name']}: Q{x['median_avg_quality']}, {x['median_length']}bp, {x['reads_m']}\n")})
annotation <- paste(
  annotation,
  collapse="\n"
)



# Plot
read_quality_stats_p <- ggplot(temp, aes(x=avg_quality, fill=basecalling_model_name)) +
  # geom_histogram(aes(y=..density..), bins = 100, alpha=0.5, position="identity") +
  geom_density(alpha=0.5) +
  annotate("label", x = 25, y = 0.35, hjust=0, size=4,
           label.padding=unit(0.7, "lines"), label.r = unit(0.4, "lines"),
           label = annotation
  ) +
  theme_bw() +
  xlab("Quality") +
  ylab("Density") +
  ggtitle(NULL) +
  guides(fill=guide_legend(title="Basecalling model")) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) + 
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    # axis.text.x=element_text(angle = 0, hjust=1, vjust=0, size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold"),
    strip.text.y = element_text(colour = 'black', face="bold.italic"),
    # panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )
read_quality_stats_p <- read_quality_stats_p + geom_magnify(
  from = c(xmin = 25, xmax = max(temp$avg_quality), ymin = 0, ymax = 0.005), 
  to = c(xmin = 27, xmax = 37, ymin = 0.05, ymax = 0.25),
  corners = 0.1, axes = "xy"
)
# ggsave(glue("{output_folder}/reads_quality.png"), plot=read_quality_stats_p, dpi="retina", units="px", width=3500, height=2500)
rm(temp, sup, hac, fast) # remove large objects


###################################



# ---- Figure 1 ----

# -- CLR Correspondence between platforms --
temp_cols <- c("Domain", "Phylum", 
               "Class", "Order", "Family", "Genus", "Species","sample.type.label","platform")

x <- clr_ps_long_ont_models_sp %>% 
  dplyr::group_by(platform, basecalling_model_name, across(temp_cols)) %>% 
  dplyr::summarise(mean_abundance=mean(Abundance))
clr.abund.at.zero <- get.clr.abund.at.zero(rel_abund_ps_long_ont_models_sp, clr_ps_long_ont_models_sp)

x <- x[c("sample.type.label","platform","basecalling_model_name","Family","Genus", "Species","mean_abundance")]
x$concat_taxa <- paste(x$Genus, x$Species, sep=";")

# x <- x %>% tidyr::pivot_wider(names_from = basecalling_model_name, values_from = mean_abundance)


x <- dplyr::bind_rows(
  subset(x, basecalling_model_name %in% c("sup","hac")) %>% 
    tidyr::pivot_wider(names_from = basecalling_model_name, values_from = mean_abundance) %>%
    dplyr::rename(x=sup, y=hac) %>%
    mutate(comp = "sup vs. hac"),
  subset(x, basecalling_model_name %in% c("sup","fast")) %>% 
    tidyr::pivot_wider(names_from = basecalling_model_name, values_from = mean_abundance) %>%
    dplyr::rename(x=sup, y=fast) %>%
    mutate(comp = "sup vs. fast"),
  subset(x, basecalling_model_name %in% c("fast","hac")) %>% 
    tidyr::pivot_wider(names_from = basecalling_model_name, values_from = mean_abundance) %>%
    dplyr::rename(x=hac, y=fast) %>%
    mutate(comp = "hac vs. fast")
)



x$comp <- factor(x$comp, levels=c("sup vs. hac", "hac vs. fast", "sup vs. fast"))


# differing.taxa <- unique((x %>% 
#          group_by(sample.type.label, platform, Family, Genus, Species, concat_taxa, comp) %>% 
#          dplyr::summarise(diff=x-y) %>% 
#          dplyr::filter(diff < -0.8 | diff > 0.8))$Species)

x <- x %>% dplyr::mutate(diff=abs(x-y)) 

clr_model_correspondence_p <- ggplot(x, aes(x=x, y=y, color=comp, fill=comp)) +
  geom_smooth(aes(fill=comp), method = "lm", se = TRUE, alpha=0.08, linetype=3, show.legend=FALSE) +
  geom_path(aes(group=concat_taxa), color="black", alpha=0.6, linetype=2, linewidth=0.3) + 
  geom_point(aes(alpha=diff), size = 2) + 
  geom_point(shape=1, size = 2, colour = "black", stroke=0.33) +
  stat_cor(
    aes(label =  paste(..rr.label.., ..p.label.., sep = "~`,`~")), # paste(..rr.label.., ..p.label.., sep = "~`,`~")
    method = "spearman", cor.coef.name="Spearman: R", p.accuracy = 0.001, r.accuracy = 0.01,
    label.x=-0.5, label.y=7, show.legend = FALSE
  ) +
  # scale_colour_gradient2(low = "white", mid = "darkgray", high = "tomato4", midpoint = 0.05) +
  # geom_label_repel(
  #   aes(label=Species, fill=NULL),
  #   subset(x, diff > 0.8 ), #& comp=="sup vs. fast"
  #   color="black", segment.colour="gray50", show.legend = FALSE, nudge_x=0.5, nudge_y=2
  # ) +
  facet_wrap(~sample.type.label) +
  facet_grid(sample.type.label~comp) +
  # guides(color=guide_legend(title="Comparison (D1 vs. D2)          ")) +
  xlab("CLR Abundance (D1)") +
  ylab("CLR Abundance (D2)") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  # scale_x_continuous(expand = c(0,0)) +
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold"),
    strip.text.y = element_text(colour = 'black', face="bold"),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )
clr_model_correspondence_p




xxx <- dplyr::bind_rows(
  subset(clr_ps_long_ont_models_genus[c("biosample_name","basecalling_model_name","Abundance","Genus")], basecalling_model_name %in% c("sup","hac")) %>% 
    tidyr::pivot_wider(names_from = basecalling_model_name, values_from = Abundance) %>%
    dplyr::rename(x=sup, y=hac) %>%
    mutate(comp = "sup vs. hac"),
  subset(clr_ps_long_ont_models_genus[c("biosample_name","basecalling_model_name","Abundance","Genus")], basecalling_model_name %in% c("sup","fast")) %>% 
    tidyr::pivot_wider(names_from = basecalling_model_name, values_from = Abundance) %>%
    dplyr::rename(x=sup, y=fast) %>%
    mutate(comp = "sup vs. fast"),
  subset(clr_ps_long_ont_models_genus[c("biosample_name","basecalling_model_name","Abundance","Genus")], basecalling_model_name %in% c("hac","fast")) %>% 
    tidyr::pivot_wider(names_from = basecalling_model_name, values_from = Abundance) %>%
    dplyr::rename(x=hac, y=fast) %>%
    mutate(comp = "hac vs. fast"),
)
xxx <- xxx %>% dplyr::mutate(diff=abs(x-y)) 


ggplot(xxx, aes(x=x, y=y, color=comp, fill=comp)) +
  geom_smooth(aes(fill=comp), method = "lm", se = TRUE, alpha=0.08, linetype=3, show.legend=FALSE) +
  # geom_path(aes(group=Genus), color="black", alpha=0.6, linetype=2, linewidth=0.3) + 
  geom_point(size = 2) + 
  geom_point(shape=1, size = 2, colour = "black", stroke=0.33) +
  # stat_cor(
  #   aes(label =  paste(..rr.label.., ..p.label.., sep = "~`,`~")), # paste(..rr.label.., ..p.label.., sep = "~`,`~")
  #   method = "spearman", cor.coef.name="Spearman: R", p.accuracy = 0.001, r.accuracy = 0.01,
  #   label.x=-0.5, label.y=7, show.legend = FALSE
  # ) +
  # scale_colour_gradient2(low = "white", mid = "darkgray", high = "tomato4", midpoint = 0.05) +
  geom_label_repel(
    aes(label=Genus, fill=NULL),
    subset(xxx, diff > 2 & biosample_name %in% xxx$biosample_name[1:30] & comp == "sup vs. hac"), #
    color="black"#, segment.colour="gray50", show.legend = FALSE, nudge_x=0.5, nudge_y=2
  ) +
  # facet_wrap(~biosample_name) +
  # facet_grid(biosample_name~comp) +
  # guides(color=guide_legend(title="Comparison (D1 vs. D2)          ")) +
  xlab("CLR Abundance (D1)") +
  ylab("CLR Abundance (D2)") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  # scale_x_continuous(expand = c(0,0)) +
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold"),
    strip.text.y = element_text(colour = 'black', face="bold"),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )

# -- Plot --
p <- ggpubr::ggarrange(
  ggpubr::ggarrange(
    read_quality_stats_p + theme(legend.position="none"), 
    basecalling_model_silva_jsd + theme(legend.position="none"), 
    widths=c(0.6,0.4),
    ncol = 2, labels = c("A", "B")#, common.legend=TRUE, legend="bottom"
  ),
  ggpubr::ggarrange(
    model_alphadiversity_p + theme(legend.position="bottom") + guides(fill=guide_legend(title="Basecalling model     ")),
    # emu_reads_comparison_p + theme(legend.position="bottom"),
    # clr_model_correspondence_p,
    ncol=1,
    labels = c("C")#, common.legend=TRUE, legend="bottom"
  ),
  nrow=2, #heights=c(0.25,0.75),
  font.label = list(size = 30), common.legend=TRUE, legend="bottom"
) + ggpubr::bgcolor("white") + ggpubr::border("white")
p
ggsave(glue("{output_folder}/possible_fig_1.png"), plot=p, dpi="retina", units="px", width=4200, height=3500)


# ---- Figure 2 ----
# -- CLR Correspondence between platforms --
temp_cols <- c("Domain", "Phylum", 
               "Class", "Order", "Family", "Genus","sample.type.label")

x <- clr_ps_long_genus %>% 
  dplyr::group_by(platform, basecalling_model_name, across(temp_cols)) %>% 
  dplyr::summarise(mean_abundance=mean(Abundance))
clr.abund.at.zero <- get.clr.abund.at.zero(rel_abund_ps_long_genus, clr_ps_long_genus)

x <- merge(
  subset(x[c(temp_cols,"platform","basecalling_model_name","mean_abundance")], platform=="Oxford_Nanopore" & basecalling_model_name=="sup"),
  subset(x[c(temp_cols,"platform","basecalling_model_name","mean_abundance")], platform=="Illumina"),
  by=c(temp_cols)
) %>% 
  dplyr::rename(mean_abundance_ont=mean_abundance.x, mean_abundance_illumina=mean_abundance.y)

x <- x[c("sample.type.label","Family","Genus","mean_abundance_ont","mean_abundance_illumina")]
x$concat_taxa <- paste(x$Family, x$Genus, sep=";")

thresh <- clr.abund.at.zero
x <- x %>% dplyr::mutate(
  appears = case_when(
    mean_abundance_ont > thresh & mean_abundance_illumina > thresh ~ "Both",
    mean_abundance_ont > thresh & mean_abundance_illumina <= thresh ~ "Only ONT",
    mean_abundance_ont <= thresh & mean_abundance_illumina > thresh ~ "Only Illumina",
    mean_abundance_ont <= thresh & mean_abundance_illumina <= thresh ~ "None",
    TRUE ~ ""
  )
)
x <- subset(x, Family!="uncultured" & !Family %like% "Unassigned")

# Check unequal taxonomy
x_subset <- subset(x, (appears %like% "Only" | appears=="None") & sample.type.label=="Cancer")
for (i in unique(x_subset$Genus)){
  if (!endsWith(i, " NA") & dim(subset(x_subset, Genus %like% glue("{i} NA")))[1]>=1){
    # print("##################")
    print(i)
    # x <- subset(x, !Genus %like% i) ### WARNING
    # print(subset(x_subset, Genus %like% i))
  }
}

# x <- subset(x, sample.type.label=="Cancer")
clr_technology_correspondence_p <- ggplot(x, aes(x=mean_abundance_ont, y=mean_abundance_illumina, color=appears, fill=appears)) +
  geom_point(alpha=0.5, size = 2, color="black", shape = 21, stroke = 0.5) + 
  geom_smooth(aes(fill=appears), data=subset(x, appears=="Both"), method = "lm", se = TRUE, alpha=0.08, linetype=3, show.legend=FALSE) +
  geom_label_repel(
    aes(label=Genus),
    data=subset(x, sample.type.label=="Control" & Genus %in% c("Fusobacterium","Parvimonas","Peptostreptococcus")),
    segment.colour="gray50", show.legend = FALSE, nudge_x=2, nudge_y=6.5, fontface = "bold.italic", color="black", size = 3
  ) +
  geom_label_repel(
    aes(label=Genus),
    data=subset(x, sample.type.label=="Cancer" & Genus %in% c("Fusobacterium","Parvimonas","Peptostreptococcus")),
    segment.colour="gray50", show.legend = FALSE, nudge_x=1, nudge_y=5, fontface = "bold.italic", color="black", size = 3
  ) +
  stat_cor(
    aes(label =  paste(..rr.label.., ..p.label.., sep = "~`,`~")), # paste(..rr.label.., ..p.label.., sep = "~`,`~")
    data=subset(x, appears=="Both"),
    method = "pearson", cor.coef.name="Pearson: R", p.accuracy = 0.001, r.accuracy = 0.01,
    label.x=3, label.y=7.5, show.legend = FALSE
  ) +
  facet_wrap(~sample.type.label) +
  guides(fill=guide_legend(title="Taxon appears\non average on")) +
  xlab("CLR Abundance (ONT)") +
  ylab("CLR Abundance (Illumina)") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme(
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold"),
    strip.text.y = element_text(colour = 'black', face="bold"),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )
clr_technology_correspondence_p



p <- ggpubr::ggarrange(
  ggpubr::ggarrange(
    platform_alphadiversity_p + 
      theme(legend.position="none"),
      # scale_fill_brewer(palette="Dark2"),
    platform_genus_silva_jsd + 
      # geom_path(aes(group=biosample_name), color="black", alpha=0.6, linetype=2, linewidth=0.1) + 
      theme(legend.position="none"),
      # scale_color_brewer(palette="Dark2"), 
    platform_species_silva_jsd + 
      # geom_path(aes(group=biosample_name), color="black", alpha=0.6, linetype=2, linewidth=0.1) + 
      theme(legend.position="none"),
      # scale_color_brewer(palette="Dark2"), 
    ncol = 3, labels = c("A", "B", "C"), widths=c(0.45,0.27,0.27), common.legend=TRUE, legend="right"#, hjust = -0.2
  ),
  ggpubr::ggarrange(
    clr_technology_correspondence_p,
    labels = c("D"), legend="right"#, hjust = -0.05
  ),
  nrow=2,
  font.label = list(size = 30)
) + ggpubr::bgcolor("white") + ggpubr::border("white")   
p
ggsave(glue("{output_folder}/possible_fig_2.png"), plot=p, dpi="retina", units="px", width=4200, height=2800)





# ---- Figure 3 ----
# p <- ggpubr::ggarrange(
#   clr_genus_platform,# + scale_fill_brewer(palette="Set1"),
#   clr_sp_ont + facet_wrap(~target_level, scales="free", ncol=4),# + scale_fill_brewer(palette="Set1"),
#   nrow=2,
#   heights = c(0.4,0.6),
#   labels=c("A","B"),
#   font.label = list(size = 20),
#   common.legend=TRUE,legend="bottom"
# ) + ggpubr::bgcolor("white") + ggpubr::border("white")   
# p
# ggsave(glue("{output_folder}/possible_fig_3.png"), plot=p, dpi="retina", units="px", width=4100, height=4000)
# 

##### WARNNING TO DO: DONT ALLOW MULTIPLE REPS OF A SAMPLE PER TECHNOLOGY AND DATABASE

temp <- rel_abund_ps_long_databases_genus[c("sample.id","biosample_name","sample.type.label","platform","database.label","region_16s","Genus","Abundance")] %>%
  dplyr::filter(Genus %in% c("Parvimonas", "Fusobacterium", "Peptostreptococcus")) %>%
  dplyr::filter(sample.id %in% sample_data(rar_physeq)$sample.id) %>%
  dplyr::group_by(biosample_name, sample.type.label, region_16s, database.label, Genus) %>%
  dplyr::summarise(sum_abundance=sum(100*Abundance), sample_count=n()) %>%
  tidyr::pivot_wider(names_from = c(region_16s, database.label), values_from = sum_abundance) %>% 
  dplyr::arrange(biosample_name, Genus) %>%
  dplyr::filter(!is.na(`ONT-V1V9_SILVA`) & !is.na(`Illumina-V3V4_SILVA`) & !is.na(`ONT-V1V9_Default`)) # if one of the values is NA that sample is not available for both technologies, remove
temp

# -- Percentage of samples with taxon present per group --
pct.per.class <- temp %>%
  dplyr::group_by(Genus, sample.type.label) %>%
  dplyr::summarise(
    "Illumina-V3V4 SILVA" = round(100 * sum(`Illumina-V3V4_SILVA`!=0) / n(), 2), 
    "ONT-V1V9 SILVA" = round(100 * sum(`ONT-V1V9_SILVA`!=0) / n(), 2),
    "ONT-V1V9 Default" = round(100 * sum(`ONT-V1V9_Default`!=0) / n(), 2),
  )

# pct.per.class.table.p <- ggtexttable(
#     pct.per.class %>% 
#       dplyr::rename("Sample type" = sample.type.label), 
#     rows = NULL, theme = ttheme("blank")
#   ) %>%
#   tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)



pct.per.class.p <- ggplot(
  pct.per.class %>%   
    tidyr::pivot_longer(
      cols=c("Illumina-V3V4 SILVA","ONT-V1V9 SILVA","ONT-V1V9 Default"), 
      names_to = "platform_database", values_to = "presence_pct"
    ) %>% 
    dplyr::mutate(platform_database=gsub(" ","\n",platform_database)), 
  aes(x=platform_database, y=presence_pct, fill=sample.type.label)
) + 
  geom_bar(stat="identity", position=position_dodge(0.9), alpha=0.6, width = 0.8, color="black") +
  geom_text(aes(label = paste(presence_pct, "%", sep="")), hjust = -0.1, size=3, position = position_dodge(0.9)) +
  facet_wrap(~Genus, ncol=1, strip.position="right") +
  ylab("% of samples with taxon") +
  xlab(NULL) +
  guides(fill=guide_legend(title=NULL)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,55)) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold"),
    strip.text.y = element_text(colour = 'black', face="bold.italic"),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )
pct.per.class.p

# -- Relative abundance per taxon and sample --
temp <- temp %>% 
  dplyr::mutate(non_zero_presence_count=sum(c(`Illumina-V3V4_SILVA`!=0, `ONT-V1V9_Default`!=0, `ONT-V1V9_SILVA`!=0)), na.rm = TRUE) %>%
  tidyr::pivot_longer(cols=c("Illumina-V3V4_SILVA","ONT-V1V9_Default","ONT-V1V9_SILVA"), names_to = "platform_database", values_to = "abundance") %>%
  dplyr::arrange(desc(abundance)) %>% 
  dplyr::mutate(biosample_name=factor(biosample_name, levels=unique(.$biosample_name)))  %>%
  dplyr::na_if(0)
temp

sample.matrix.p <- ggplot(
    temp %>% dplyr::mutate(platform_database=gsub("_","\n",platform_database)), 
    aes(x=platform_database, y=biosample_name, fill=abundance)
  ) + 
  geom_tile(show.legend = FALSE) +
  scale_fill_distiller(palette = "YlOrRd", na.value="white", direction=1) +
  geom_text(aes(label = round(abundance, 4)), size=2.5) +
  facet_grid(sample.type.label~Genus, scales="free", space = "free") +
  xlab(NULL) +
  ylab("Samples") +
  scale_x_discrete(expand = c(0,0)) +
  theme_bw() +
  theme(
    # axis.text.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"),
    strip.text = element_text(colour = 'black', face="bold"),
    strip.text.x = element_text(colour = 'black', face="bold.italic"),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )
sample.matrix.p


# -- Figure 4 --
p <- ggpubr::ggarrange(
  clr_genus_platform + facet_wrap(~target_level, scales="free_y", nrow=3,  strip.position="right"),
  pct.per.class.p + theme(legend.position="bottom"),
  ncol=2,
  labels=c("A","B"),
  font.label = list(size = 20),
  common.legend=TRUE,legend="bottom"
) + ggpubr::bgcolor("white") + ggpubr::border("white")   
p

# p <- ggpubr::ggarrange(
#   ggpubr::ggarrange(
#     pct.per.class.p + theme(legend.position="bottom"),
#     nrow=2, heights=c(0.6, 0.4)
#   ),
#   clr_genus_platform,
#   ncol=2,
#   # widths = c(0.45,0.55),
#   labels=c("A","B"),
#   font.label = list(size = 20),
#   common.legend=TRUE,legend="bottom"
# ) + ggpubr::bgcolor("white") + ggpubr::border("white")   
# p
ggsave(glue("{output_folder}/possible_fig_3.png"), plot=p, dpi="retina", units="px", width=4500, height=4000)


# -- Figure 4 --
p <- ggpubr::ggarrange(
  ggpubr::ggarrange(
    plot.ancombc(
      subset(ont.silva.df_fig, target_level=="Species") %>% 
        dplyr::mutate(  
          direction = case_when(
            direction == "Positive LFC" ~ "More in Cancer  ",
            direction == "Negative LFC" ~ "More in Control  ",
            TRUE ~ ""
          ),
          taxon_id = gsub(" NA$", " sp.", taxon_id)
        )
    ) + 
      theme(
        strip.background = element_blank(), 
        strip.text.x = element_blank(), 
        legend.text=element_text(size=12), 
        strip.text=element_text(size=10),
        plot.margin = margin(0,0.25,0,2, "cm")
    ) +
      guides(fill=guide_legend(title=NULL)),
    plot.ancombc(
      subset(ont.default.df_fig, target_level=="Species") %>% 
        dplyr::mutate(  
          direction = case_when(
            direction == "Positive LFC" ~ "More in Cancer  ",
            direction == "Negative LFC" ~ "More in Control  ",
            TRUE ~ ""
          )
        )
    ) + 
      theme(
        strip.background = element_blank(), 
        strip.text.x = element_blank(), 
        legend.text=element_text(size=12), 
        strip.text=element_text(size=10),
        plot.margin = margin(0,0.25,0,2, "cm")
    ) +
      guides(fill=guide_legend(title=NULL)),
    labels=c("A","B"),
    font.label = list(size = 20),
    common.legend=TRUE, legend="bottom"
  ),
  ggpubr::ggarrange(
    clr_sp_ont_both_p + 
      theme(legend.text=element_text(size=12), 
            strip.text=element_text(size=10)),
    font.label = list(size = 20),
    labels=c("C"), legend="bottom"
  ),
  nrow=2, heights=c(0.3,0.7)
) + ggpubr::bgcolor("white") + ggpubr::border("white")



ggsave(glue("{output_folder}/possible_fig_4.png"), plot=p,
       dpi="retina", unit="px", height=6000, width=6000)




#####################################################
#####################################################



################################################################
# -- Get proportion of feature counts not matching patterns -- #
################################################################
selected_taxa_level = "Species"
temp <- physeq %>%
  tax_glom(taxrank = selected_taxa_level, bad_empty = c()) %>%   # Aglomerate to taxnomic range
  psmelt()
temp$target_level <- temp[[selected_taxa_level]]


pattern <- " NA$|uncultured|unidentified|unclassified|human gut|^bacterium | sp.$|gut metagenome|metagenome"
species_level_idenfitication_percentages_p <- temp %>%
  # Establish taxa not matching pattern
  dplyr::mutate(
    taxa_is_identified = !grepl(pattern, target_level)
  ) %>%
    # Group by and sum counts according to if taxa is identified
    dplyr::group_by(biosample_name, platform, region_16s, basecalling_model_name, database, database.label, sample.type.label, taxa_is_identified) %>%
    dplyr::summarise(sum_abundance=sum(Abundance)) %>%
    pivot_wider(names_from = taxa_is_identified, values_from = sum_abundance, names_prefix="taxa_is_identified_") %>%
    dplyr::mutate(
      identified_pct = 100 * taxa_is_identified_TRUE/(taxa_is_identified_FALSE+taxa_is_identified_TRUE),
      unidentified_pct = 100 * taxa_is_identified_FALSE/(taxa_is_identified_FALSE+taxa_is_identified_TRUE)
    ) %>%
    pivot_longer(cols=c(identified_pct, unidentified_pct), names_to = "type", values_to = "pct") %>%
  ##############
  # -- Plot -- #
  ##############
  dplyr::filter(database=="silva") %>%
    dplyr::filter(type=="identified_pct") %>%
    dplyr::mutate(label=glue("{region_16s}\n{database.label}\n{basecalling_model_name}")) %>%
    ggplot(
      aes(x=label, y=pct, fill=label)
    ) + 
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width=0.25), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
    facet_wrap(~sample.type.label) + 
    ggpubr::geom_pwc(
      aes(group = label), step.increase=0.05, vjust = 0.5, tip.length = 0.01,
      method = "wilcox_test", label = "{p.adj.signif}",
      p.adjust.method = "holm", hide.ns = TRUE
    ) +
    ylab(glue("Feature counts percentage properly identified at {selected_taxa_level} level")) +
    xlab(NULL) +
    guides(fill=guide_legend(title=NULL)) +
    # scale_y_continuous(expand = c(0,0), limits=c(0,110)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text=element_text(size=10, color="black"), 
      axis.title=element_text(size=10, face="bold"),
      axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
      legend.title=element_text(size=10, face="bold"),
      legend.text=element_text(size=10),
      plot.title = element_text(size = 15, color = 'black', face="bold"),
      strip.background = element_rect(fill="white", color="black"), # #3263b0
      strip.text = element_text(colour = 'black', face="bold"),
      panel.background = element_blank(),
      axis.ticks.x=element_blank()
    )




##############################################
# ---- Compare Emu reads classification ---- #
##############################################

get.matched.pct.per.sample <- function(filepath, selected_taxa_level){
  data <- read.delim(filepath) %>%
    dplyr::rename(level_0=glue("{selected_taxa_level}_0"), level_1=glue("{selected_taxa_level}_1")) %>%
    dplyr::group_by(comparison_id, level_0, level_1) %>%
    dplyr::summarise(count=sum(read_id)) %>%
    dplyr::mutate(
      matches=case_when(
        level_0 == level_1 ~ "match",
        TRUE ~ "nomatch"
      )
    ) %>%
    dplyr::group_by(comparison_id, matches) %>%
    dplyr::summarise(count=sum(count)) %>%
    tidyr::pivot_wider(names_from="matches", values_from="count", values_fill=0) %>%
    dplyr::summarise(non_matching_reads_pct=100*nomatch/(match+nomatch)) %>%
    dplyr::mutate(level=selected_taxa_level) 
  
  return(data)
}


temp_folder = "/home/usuario/Proyectos/Results/ONT16S/pipeline/10.2_compare_emu_reads"
data <- dplyr::bind_rows(
  # hac vs sup silva
  get.matched.pct.per.sample(glue("{temp_folder}/hac_sup_silva.tsv"), "family") %>% dplyr::mutate(comp="hac_vs_sup (silva)"),
  get.matched.pct.per.sample(glue("{temp_folder}/hac_sup_silva.tsv"), "genus") %>% dplyr::mutate(comp="hac_vs_sup (silva)"),
  get.matched.pct.per.sample(glue("{temp_folder}/hac_sup_silva.tsv"), "species") %>% dplyr::mutate(comp="hac_vs_sup (silva)"),
  # fast vs sup silva
  get.matched.pct.per.sample(glue("{temp_folder}/fast_sup_silva.tsv"), "family") %>% dplyr::mutate(comp="fast_vs_sup (silva)"),
  get.matched.pct.per.sample(glue("{temp_folder}/fast_sup_silva.tsv"), "genus") %>% dplyr::mutate(comp="fast_vs_sup (silva)"),
  get.matched.pct.per.sample(glue("{temp_folder}/fast_sup_silva.tsv"), "species") %>% dplyr::mutate(comp="fast_vs_sup (silva)"),
  # hac vs sup default
  get.matched.pct.per.sample(glue("{temp_folder}/hac_sup_default.tsv"), "family") %>% dplyr::mutate(comp="hac_vs_sup (default)"),
  get.matched.pct.per.sample(glue("{temp_folder}/hac_sup_default.tsv"), "genus") %>% dplyr::mutate(comp="hac_vs_sup (default)"),
  get.matched.pct.per.sample(glue("{temp_folder}/hac_sup_default.tsv"), "species") %>% dplyr::mutate(comp="hac_vs_sup (default)"),
  # fast vs sup default
  get.matched.pct.per.sample(glue("{temp_folder}/fast_sup_default.tsv"), "family") %>% dplyr::mutate(comp="fast_vs_sup (default)"),
  get.matched.pct.per.sample(glue("{temp_folder}/fast_sup_default.tsv"), "genus") %>% dplyr::mutate(comp="fast_vs_sup (default)"),
  get.matched.pct.per.sample(glue("{temp_folder}/fast_sup_default.tsv"), "species") %>% dplyr::mutate(comp="fast_vs_sup (default)")
)

data %>%
  dplyr::group_by(comp, level) %>%
  rstatix::get_summary_stats()

stat.test <- data %>% 
  dplyr::group_by(comp) %>%
  rstatix::wilcox_test(non_matching_reads_pct ~ level) %>%
  rstatix::adjust_pvalue(method = "holm") %>%
  # rstatix::add_xy_position(x = "comp", group = "level") %>%
  dplyr::filter(p.adj<=0.1) 

emu_reads_comparison_p <- ggplot(data, aes(x=comp, y=non_matching_reads_pct, fill=level)) +
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.1), color="black", shape = 21, stroke = 0.5, show.legend=FALSE, size=1, alpha=0.7) +
  ggpubr::geom_pwc(
    aes(group = level), step.increase=0.05, vjust = 0.5, tip.length = 0.01,
    method = "wilcox_test", label = "{p.adj.signif}",
    p.adjust.method = "holm", hide.ns = TRUE
  ) +
  ggpubr::geom_pwc(
    aes(group = comp), step.increase=0.05, vjust = 0.5, tip.length = 0,
    method = "wilcox_test", label = "{p.adj.signif}",
    p.adjust.method = "holm", hide.ns = TRUE,
    bracket.nudge.y = 0.2, size=0.75
  ) +
  scale_fill_brewer(palette="Dark2") +
  ylab("Percentage of reads identified differently") +
  xlab(NULL) +
  guides(fill=guide_legend(title="Taxa level     ")) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text=element_text(size=10, color="black"), 
    axis.title=element_text(size=10, face="bold"),
    axis.title.y=element_text(size=10, face="bold", margin = margin(t = 0, r = 5, b = 0, l = 0)),
    legend.title=element_text(size=10, face="bold"),
    legend.text=element_text(size=10),
    plot.title = element_text(size = 15, color = 'black', face="bold"),
    strip.background = element_rect(fill="white", color="black"), # #3263b0
    strip.text = element_text(colour = 'black', face="bold"),
    panel.background = element_blank(),
    axis.ticks.x=element_blank()
  )


























