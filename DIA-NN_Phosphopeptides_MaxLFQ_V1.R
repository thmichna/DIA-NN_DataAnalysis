
library(tidyverse)
library(data.table)


##### LOAD DATA ################################################################

raw_data = read_delim("output.tsv", delim="\t")


# add raw file mapping from .csv file
rawfile_mapping = read_delim("rawfile_mapping.csv", delim=",")
raw = left_join(raw_data, rawfile_mapping, by = c("Run"="raw_file") )




##### OVERVIEW OF PEPTIDE COUNTS ###############################################

pep_summary = raw %>%
  group_by(identifier) %>%
  filter( !duplicated(Modified.Sequence) ) %>%
  summarize(
    pep_all = n(),
    pep_phos = sum( str_detect(Modified.Sequence, fixed("UniMod:21")), na.rm=T  ),
    enrichment = pep_phos / pep_all
  )


ggplot(pep_summary) +
  geom_col( aes(x=identifier, y=pep_all, fill="All Peptides"), fill="#335584", color="black", size=0.3, width=0.75 ) +
  geom_col( aes(x=identifier, y=pep_phos, fill="Phosphopeptides"), fill="darkgreen", color="black", size=0.3, width=0.75 ) +
  geom_text( aes(x=identifier, y=pep_all, label=pep_all), angle=90, size=3, hjust=-0.5 ) +
  geom_text( aes(x=identifier, y=pep_phos, label=pep_phos), angle=90, size=3, hjust=-0.1 ) +
  geom_text( aes(x=identifier, y=55000, label=round(enrichment, 2)), angle=90, size=3, hjust=0.5, color="red" ) +
  geom_text( aes(x=identifier[1], y=52500, label="enrichment"), angle=90, size=3, hjust=1, vjust=0.5, color="red" ) +
  scale_y_continuous(breaks=seq(0,100000,5000), limits=c(0,60000)) +
  theme_minimal() +
  theme( axis.text.x = element_text(angle=90, vjust=0.5, size=10),
         axis.text.y = element_text(size=10))



ggsave("01_peptide_counts.png", width = 12)




##### SUMMARIZE Modified.Sequence ##############################################

raw_summary = raw %>% 
                  filter( str_detect(Modified.Sequence, fixed("UniMod:21")) ) %>%
                  group_by( identifier, Modified.Sequence ) %>%
                  summarize( Site.confidence = max(PTM.Site.Confidence, na.rm=T),
                             CScore = max(CScore, na.rm=T),
                             Max.Int.Norm = max(Precursor.Normalised, na.rm=T),
                             prot.names =  paste(unique(Protein.Names), collapse = "_") ) %>%
                  ungroup() 



#### WIDE FORMAT ##############################################################

res_wide = raw_summary %>%
  group_by(Modified.Sequence) %>%
  mutate( Site.confidence = max(Site.confidence, na.rm=T),
          CScore = max(CScore, na.rm=T),
          prot.names =  paste(unique(prot.names), collapse = "_") ) %>%
  spread(key=identifier, value=Max.Int.Norm)

write.csv(res_wide, "PrecInt_PPep_WideFormat.csv")



#### MAKE SAMPLE CORRELATION ###################################################

sample_correlation = tibble(
  "file1" = NULL,
  "file2" = NULL,
  "pearson" = NULL,
  "spearman" = NULL
)

for (f1 in unique(raw_summary$identifier) ) {
  
  for (f2 in unique(raw_summary$identifier) ) {
    
  f1_int = res_wide %>% ungroup() %>% select(f1) %>% unlist() %>% as.numeric()
  f1_int[is.na(f1_int)] = 0
  
  f2_int = res_wide %>% ungroup() %>% select(f2) %>% unlist() %>% as.numeric()
  f2_int[is.na(f2_int)] = 0
  
  pearson_corr = cor(f1_int, f2_int, method="pearson")
  spearman_corr = cor(f1_int, f2_int, method="spearman")
  
  current_correlations = tibble(
    "file1" = f1,
    "file2" = f2,
    "pearson" = pearson_corr,
    "spearman" = spearman_corr
  )
  
  sample_correlation = bind_rows(sample_correlation, current_correlations)
  
  }
  
}

#  "heatmap"
ggplot(sample_correlation) +
  scale_fill_gradient2(low = "darkred", mid="yellow", high = "darkgreen", midpoint=0.5) +
  geom_raster( aes(x=file1, y=file2, fill=pearson)  ) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
  
ggsave("02_phosphopeptides_correlation.png", width=10, height=8 )


