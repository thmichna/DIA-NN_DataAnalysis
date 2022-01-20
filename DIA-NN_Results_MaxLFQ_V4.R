
library(tidyverse)
library(data.table)


##### LOAD DATA ################################################################

raw_data = read_delim("output.tsv", delim="\t")


# add raw file mapping from .csv file
rawfile_mapping = read_delim("rawfile_mapping.csv", delim=",")
raw = left_join(raw_data, rawfile_mapping, by = c("Run"="raw_file") )




##### SUMMARIZE Protein.Group ##################################################

prot_count = raw %>% group_by(Protein.Group, identifier) %>% summarize( count=n() ) %>% group_by(identifier) %>% summarize(count=n())
unique_protIDs = raw %>% select(Protein.Group) %>% filter(!duplicated(Protein.Group))

write.csv(prot_count, "prot_count.csv")
write_csv(unique_protIDs, "unique_protIDs.csv")

ggplot(prot_count) +
  geom_col( aes(x=identifier, y=count), fill="#335584" ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  scale_y_continuous(limits=c(0,6000), breaks=seq(0,10000,500)) +
  geom_text(aes( x=identifier, y=count-500, label = count ), angle=90, color="white", size=3.5 )

ggsave("01_protein_counts.png", units="cm", dpi = 300, width = 24, height = 10)



##### SUMMARIZE Protein.Group ##################################################

raw_summary = raw %>% 
                  group_by( identifier, Protein.Group) %>%
                  summarize( no_pep = n(),
                             proteotypic_pep = sum(Proteotypic) ,
                             shared_pep = sum(!Proteotypic) ,
                             PG.Q.val = paste(unique(PG.Q.Value), collapse = "_"),
                             PG.MaxLFQ = paste(unique(Genes.MaxLFQ), collapse = "_"),
                             prot.names =  paste(unique(Protein.Names), collapse = "_"),
                             cond = paste(unique(condition), collapse="_"),
                             submission_no = paste(unique(submission_no), collapse="_"),
                             run = paste(unique(run), collapse="_") ) %>%
                  filter( str_detect(prot.names, "MOUSE"), no_pep>=2 ) %>% # optional - maybe include filter for no_pep >= 2
                  ungroup() %>%
                  group_by( Protein.Group ) %>%
                  mutate( max_pep = max(no_pep, na.rm = T),
                          max_proteotypic_pep = max(proteotypic_pep, na.rm = T),
                          max_shared_pep = max(shared_pep, na.rm = T),
                          max_PG.Q.val = max(PG.Q.val, na.rm = T),
                          prot.names = paste(unique(prot.names), collapse = "_"))



#### WIDE FORMAT ##############################################################

res_wide = raw_summary %>%
  select(identifier, Protein.Group, prot.names, max_PG.Q.val,  max_pep, max_proteotypic_pep, max_shared_pep,  PG.MaxLFQ) %>%
  spread(key=identifier, value=PG.MaxLFQ)

write.csv(res_wide, "MaxLFQ_Intensities.csv")



#### MAKE SAMPLE CORRELATION ##############################################################

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
  scale_fill_gradient2(low = "darkred", mid="yellow", high = "darkgreen", midpoint=0.75 ) +
  geom_raster( aes(x=file1, y=file2, fill=pearson)  ) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5))
  
ggsave("02_proteome_correlation.png", width=10, height=8 )

# for better vizualization for saving
pearson_wide = sample_correlation %>%
                  select("file1","file2","pearson") %>%
                  pivot_wider(names_from="file2", values_from="pearson")
