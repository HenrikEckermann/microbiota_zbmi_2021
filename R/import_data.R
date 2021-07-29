library(tidyverse)
library(phyloseq)
library(microbiome)

# load helper functions
source(here::here("R/mb_helper.R"))

# 16S data was analyzed by Yang O. using ngtax 2 with SILVA DB 132
map_file <- here::here("data/2019_0054 to 2019_0062 20190074 75 82 mapping file.csv")
biom_file <- here::here("data/Galaxy224-[NG-Tax__BIBO_1m-10y_(20200116)].biom1")
tree_file <- read_tree(here::here("data/Galaxy226-[Tree_files__NG-Tax__BIBO_1m-10y_(20200116)]/all_otus.tree"))

# merge files
bibo_temp <- read_phyloseq(otu.file = biom_file, taxonomy.file = NULL, metadata.file = map_file, type = "biom")

# merge with tree file 
bibo <- merge_phyloseq(bibo_temp, tree_file) 
bibo_genus <- tax_glom(bibo, taxrank = "Genus")
bibo_genus_rel <- microbiome::transform(bibo_genus, transform = "compositional")

# write helper function to deduce sampling time point from sample name
deduce_time_var <- function(df) {
  df %>% mutate(time = ifelse(startsWith(Sample, "a"), "1m", 
                  ifelse(startsWith(Sample, "b"), "3m", 
                  ifelse(startsWith(Sample, "c"), "4m",
                  ifelse(startsWith(Sample, "d"), "6y",
                  ifelse(startsWith(Sample, "e"), "10y", NA))))))
}

# get clr abundances for firm and bact:
bibo_phy <- tax_glom(bibo, taxrank = "Phylum")
bibo_phy_clr <- microbiome::transform(bibo_phy, "clr")
fb_clr <- bibo_phy_clr %>% 
  otu_to_df() %>%
  select(Sample = sample_id, Firmicutes = "15796961", Bacteroidetes = "157969124") %>%
  filter(!grepl("^L.+", Sample)) %>% # disregard MOCK samples)
  deduce_time_var()



# create data objects according to description in methods section:

# list of SCFA
genus_scfa <- c("g__Bacteroides", "g__Prevotella", "g__Prevotella_2", 
                "g__Prevotella_7", "g__Prevotella_9", "g__Alistipes", "g__[Eubacterium]_hallii_group",
                "g__[Eubacterium]_rectale_group", "g__Roseburia", "g__Anaerostipes", 
                "g__Coprococcus_1", "g__Coprococcus_2", "g__Coprococcus_3", 
                "g__Blautia", "g__Faecalibacterium", "g__Subdoligranulum", "g__Dialister", "g__Phascolarctobacterium",
                "g__Akkermansia", "g__Holdemanella", "g__Bifidobacterium")
     

# we need both clr and relative abundance for analysis and plotting
bibo_rel <- microbiome::transform(bibo, "compositional")
df_rel <- psmelt(bibo_rel) %>% 
  filter(!grepl("^L.+", Sample)) %>% # disregard MOCK samples)
  deduce_time_var() %>%
  mutate(id = str_replace( Sample, "\\w", ""))

bibo_clr <- microbiome::transform(bibo, "clr")
# a df that contains clr values of the scfa genera:
df_clr <- psmelt(bibo_clr) %>% 
  filter(!grepl("^L.+", Sample)) %>% # disregard MOCK samples)
  deduce_time_var() %>% 
  filter(Genus %in% genus_scfa) %>%
  group_by(Sample, Genus) %>%
  summarise(abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = abundance) %>% 
  mutate(id = str_replace(Sample, "\\w", "")) %>%   
  mutate(time = ifelse(startsWith(Sample, "a"), "1", 
                ifelse(startsWith(Sample, "b"), "2", 
                ifelse(startsWith(Sample, "c"), "3",
                ifelse(startsWith(Sample, "d"), "4",
                ifelse(startsWith(Sample, "e"), "5", NA))))),
         time = as.integer(time)) 


# in this section I create the scfa vector:
# to deduce the scfa count, we first need all genera counts
df_count <- psmelt(bibo) %>% 
  filter(!grepl("^L.+", Sample)) %>% # disregard MOCK samples)
  deduce_time_var() %>%
  mutate(
    # dont aggreg different families/orders of undefined genera/families
    Genus = ifelse(Genus %in% c(     
      "g__uncultured_organism", 
      "g__", 
      "g__uncultured"
    ), Family, Genus),
    Genus = ifelse(Genus == "g__<empty>", Order, Genus),
    Genus = ifelse(Genus == "f__", Order, Genus)
  )


# the sum up sfca counts 
otu_temp <- df_count %>% 
  group_by(Sample, Genus) %>%
  summarise(abundance_temp = sum(Abundance), .groups = "drop") %>%
  mutate(genus = ifelse(Genus %in% genus_scfa, "scfa", Genus)) %>% 
  group_by(Sample, genus) %>%
  summarise(abundance = sum(abundance_temp), .groups = "drop") %>%
  pivot_wider(names_from = "Sample", values_from = "abundance") %>%
  column_to_rownames("genus") %>%
  as.matrix()

tax_temp <- tax_table(bibo) %>% 
  as.data.frame() %>%
  as_tibble(rownames = NA) %>%
  mutate(
    Genus = ifelse(Genus %in% c(
      "g__uncultured_organism", 
      "g__", 
      "g__uncultured"
    ), Family, Genus),
    Genus = ifelse(Genus == "g__<empty>", Order, Genus),
    Genus = ifelse(Genus == "f__", Order, Genus)
  ) %>%
  group_by(Genus) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(genus = Genus) %>%
  column_to_rownames("genus") %>%
  filter(Genus %in% rownames(otu_temp)) %>%
  select(-n) %>%
  add_row(Genus = "scfa")
rownames(tax_temp) <- c(rownames(tax_temp)[-length(rownames(tax_temp))], "scfa")

pseq_manip <- phyloseq(
    otu_table(otu_temp, taxa_are_rows = TRUE), 
    tax_table(as.matrix(tax_temp)))

pseq_manip_clr <- microbiome::transform(pseq_manip, "clr")
df_scfa_clr <- psmelt(pseq_manip_clr) %>% 
  filter(!grepl("^L.+", Sample)) %>% # disregard MOCK samples)
  deduce_time_var() %>% 
  filter(Genus == "scfa") %>%
  select(Sample, scfa_clr = Abundance, time) %>%
  mutate(id = str_replace( Sample, "\\w", ""))



# relative abundance for random forest models
otu_rf <- otu_to_df(bibo_genus_rel)
otu_rf <- otu_rf %>% 
  filter(!grepl("^L.+", sample_id)) %>%
  select(Sample = sample_id, everything()) %>%
  deduce_time_var() 


# F/B ratio 
df_ratio <- df_rel %>% group_by(Sample, Phylum) %>%
  summarise(abundance = sum(Abundance), .groups = "drop") %>%
  filter(Phylum %in% c("p__Firmicutes", "p__Bacteroidetes")) %>%
  pivot_wider(names_from = "Phylum", values_from = "abundance") %>%
  mutate(fb_ratio = p__Firmicutes / p__Bacteroidetes)  %>% 
  mutate(id = str_replace(Sample, "\\w", "")) %>%   
  mutate(time = ifelse(startsWith(Sample, "a"), "1", 
                ifelse(startsWith(Sample, "b"), "2", 
                ifelse(startsWith(Sample, "c"), "3",
                ifelse(startsWith(Sample, "d"), "4",
                ifelse(startsWith(Sample, "e"), "5", NA))))),
         time = as.integer(time)) %>%
  select(-contains("p__"))


df_ratio_pseudo <- df_rel %>% group_by(Sample, Phylum) %>%
  summarise(abundance = sum(Abundance), .groups = "drop") %>%
  filter(Phylum %in% c("p__Firmicutes", "p__Bacteroidetes")) %>%
  pivot_wider(names_from = "Phylum", values_from = "abundance") %>%
  mutate(
    p__Bacteroidetes = ifelse(p__Bacteroidetes == 0, 0.0001, p__Bacteroidetes),
    fb_ratio = p__Firmicutes / p__Bacteroidetes)  %>% 
    mutate(id = str_replace(Sample, "\\w", "")) %>%   
    mutate(time = ifelse(startsWith(Sample, "a"), "1", 
                  ifelse(startsWith(Sample, "b"), "2", 
                  ifelse(startsWith(Sample, "c"), "3",
                  ifelse(startsWith(Sample, "d"), "4",
                  ifelse(startsWith(Sample, "e"), "5", NA))))),
           time = as.integer(time)) %>%
    select(-contains("p__"))

           

# SCFA relab for each sample
df_scfa <- df_rel %>% 
  filter(Genus %in% genus_scfa) %>%
  group_by(Sample) %>%
  summarise(scfa = sum(Abundance), .groups = "drop") %>%
  deduce_time_var() %>%
  mutate(id = str_replace( Sample, "\\w", ""))




# create df that we use for the brms model per remote
ids <- unique(df_ratio_pseudo$id)  
temp1 <- df_ratio_pseudo %>% 
    full_join(
      select(df_scfa_clr, -time, -id), 
      by = "Sample"
    ) %>%
    full_join(
      select(fb_clr, Firmicutes, Bacteroidetes, Sample),
      by = "Sample"
    ) %>%
    filter(id %in% ids) %>%
    select(subject = id, everything())  

# bmi z scores
temp2 <- readxl::read_xlsx(here::here("data/tm_df_2.xlsx")) %>%
  rename(time = TIME) %>%
  filter(subject %in% ids) %>%
  mutate(subject = as.character(subject))


df <- full_join(temp1, temp2, by = c("subject", "time")) %>% 
  select(-f_subject, -FB_ratio_Tminus, scfa_Tminus) %>%
  select(
    subject,
    time, 
    fb_ratio, 
    scfa_tminus = scfa_clr, 
    bmi, 
    bmi_tminus = bmi_Tminus, 
    bw_grams,
    Firmicutes,
    Bacteroidetes
  ) %>%
  mutate(
    l_fb_ratio = ifelse(fb_ratio != 0, log(fb_ratio), log(0.000001)), 
    s_fb_ratio = scale(fb_ratio)[, 1], 
    sl_fb_ratio = scale(l_fb_ratio)[, 1],
    s_bw_grams = scale(bw_grams)[, 1]
  ) %>%
  arrange(subject, time)
df <- select(df, id = subject , everything())


# the following objects can be used for analyses or visualization
save(
  df,
  df_clr, 
  otu_rf,
  df_rel,
  df_scfa_clr,
  genus_scfa,
  fb_clr,
  file = here::here("rdata/mb_import_H.Rds")
)

