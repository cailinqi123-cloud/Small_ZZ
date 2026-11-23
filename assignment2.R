##***************************
## Assignment 1
##
## Ali Al-Reda El-Zein 
##
## 2026-11-02
##
## Assignment 2
##
## Modified by Luis Villasana 
## 
##***************************

## clear environment  
rm(list = ls())

### ADDING PACKAGES USED 
## _ Packages used -------
library(tidyverse)
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())


## Import data -----
# trichinella_spiralis
df_trich_s <- read_csv("../Small_ZZ/data/trichinella_spiralis_te.csv")
# trichinella zimbabwensis
trich_zb <- read_csv("../Small_ZZ/data/trich_z.csv")
# trichinella nativa 
df_te_nat <-  read.csv("../Small_ZZ/data/GCA_001447565.2 .csv")


## Filtering data for DNA transposons -----
# T. zimbabwensis 
trich_zb$matching_class |> unique()

non_te_z <- c("Simple_repeat", "Low_complexity", "Satellite", "snRNA", "rRNA", "Unknown")

dna_te_z <- c("DNA/PIF-Harbinger", "DNA/TcMar-Tc1", "DNA/MULE-MuDR", 
              "DNA/TcMar-Tc4", "DNA/TcMar-Pogo", "DNA/CMC-Chapaev", 
              "DNA/TcMar-Mariner", "DNA/PiggyBac","RC/Helitron")

trich_z <- trich_zb |> 
  filter(matching_class %notin% non_te_z) |>
  filter(matching_class %in% dna_te_z) 

trich_z$species <-  c("Trichinella zimbabwensis") # to create a new column with the species information "T. zimbabwensis"

# T. spiralis
df_trich_s$matching_class |> unique()

non_te_s <-  c("simple_repeat", "Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")

dna_te_s <- c("DNA/TcMar-Tc4", "DNA/TcMar-Tc1", "DNA/hAT", 
              "DNA/PIF-Harbinger", "DNA/TcMar-Mariner", "RC/Helitron", "DNA/hAT-Ac", "DNA/MULE-MuDR", "DNA/TcMar-Pogo", "DNA/CMC-Chapaev")

trich_s <- df_trich_s |> 
  filter(matching_class %notin% non_te_s) |> 
  filter(matching_class %in% dna_te_s) 

trich_s$species  <-  c("Trichinella spiralis") # to create a new column with the species information "T. spiralis"

# T. nativa
df_te_nat$matching_class |> unique()

ls_no_te_nat <-  c("Simple_repeat", "simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")

dna_te_nat <- c("DNA/TcMar-ISRm11", "DNA/hAT-Ac", "DNA/TcMar-Tc1", "DNA/hAT-Tip100", "DNA/TcMar-Tc4", "DNA/TcMar-Mariner", "DNA/MULE-MuDR", "DNA/PIF-Harbinger", "DNA/hAT", "DNA/TcMar-Pogo", "DNA/TcMar-Tc2")

te_nat <- df_te_nat |> 
  filter(matching_class %notin% ls_no_te_nat) |>
  filter(matching_class %in% dna_te_nat) 

te_nat$species  <-  c("Trichinella nativa") # to create a new column with the species information "T. nativa"


## identifying and filtering the common DNA transposons between the three species -----
trich_s$matching_class |> unique()

trich_z$matching_class |> unique()

te_nat$matching_class |> unique()

# 1ST SUBSTANTIVE EDIT: automatically find the DNA transposon classes shared by all species ----
dte_common <- Reduce(
  intersect,                            
  list(
    unique(trich_s$matching_class),     
    unique(trich_z$matching_class),     
    unique(te_nat$matching_class)      
  )
)

# using dte_common to filter all the unique te in each species
# Trichinella spiralis
ts_common <- trich_s |> 
  group_by(matching_class) |> 
  select(matching_class, qry_width, perc_insert, species) |> 
  filter(matching_class %in% dte_common) |> 
  distinct(perc_insert, .keep_all = TRUE) |> 
  distinct(qry_width, .keep_all = TRUE)

#Trichinella zimbabwensis
tz_common <- trich_z |> 
  group_by(matching_class) |> 
  select(matching_class, qry_width, perc_insert, species) |> 
  filter(matching_class %in% dte_common) |> 
  distinct(perc_insert, .keep_all = TRUE) |> 
  distinct(qry_width, .keep_all = TRUE)

# Trichinella nativa
tn_common <- te_nat |> 
  group_by(matching_class) |> 
  select(matching_class, qry_width, perc_insert, species) |> 
  filter(matching_class %in% dte_common) |> 
  distinct(perc_insert, .keep_all = TRUE) |> 
  distinct(qry_width, .keep_all = TRUE)


## filtering the top 5 per class largest DNA transposons in the three species -----
# Trichinella spiralis
ts_5 <- ts_common|> 
  top_n(5, qry_width)

# Trichinella zimbabwensis
tz_5 <- tz_common |> 
  top_n(5, qry_width)

# Trichinella nativa
tn_5 <- tn_common |> 
  top_n(5, qry_width)

# to combine the previous tables into a single one for easier posterior use
combined_5 <- rbind(ts_5, tz_5, tn_5)


## filtering the largest transposons for each species -----
# Trichinella spiralis
ts_largest <- ts_common |> 
  slice(which.max(qry_width))  

# Trichinella zimbabwensis
tz_largest <- tz_common |> 
  slice(which.max(qry_width))

# Trichinella nativa
tn_largest <- tn_common |> 
  slice(which.max(qry_width))

## creating a combined dataframe for the largest transposons in all species 
combined_largest <- rbind(tz_largest, ts_largest, tn_largest)

## creating plots of the largest te per species and the largest 5 of each species in a direct comparison
#creating T. zimbabwensis plot
p_tz <- ggplot(
  tz_largest,
  aes(x = qry_width, y=perc_insert, colour = matching_class)
) +
  geom_point(size = 5) +
  labs(title = "How Sequence Length Affects Percent Insertion in T. zimbabwensis for Each DNA Transposon Type.",
       subtitle = "Using the Largest Sequence per Transposon Type (In Base Pairs)",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type")

print(p_tz)

ggsave("tz.png", path = "../Small_ZZ/output/")

# creating T. spiralis plot
p_ts <- ggplot(
  ts_largest,
  aes(x = qry_width, y=perc_insert, colour = matching_class)
) +
  geom_point(size = 5)+
  labs(title = "How Sequence Length Affects Percent Insertion in T. spiralis for Each DNA Transposon Type.",
       subtitle = "Using the Largest Sequence per Transposon Type (In Base Pairs)",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type")

print(p_ts)   

ggsave("ts.png", path = "../Small_ZZ/output")

# creating T. nativa plot 
p_tn <- ggplot(
  tn_largest,
  aes(x = qry_width, y=perc_insert, colour = matching_class)
) +
  geom_point(size = 5)+
  labs(title = "How Sequence Length Affects Percent Insertion in T. nativa for Each DNA Transposon Type.",
       subtitle = "Using the Largest Sequence per Transposon Type (In Base Pairs)",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type")

print(p_tn)   

ggsave("tn.png", path = "../Small_ZZ/output")


# creating combined plot for top 5 largest sequences per DNA transposon
combined_top5 <- ggplot(
  combined_5,
  aes(x =qry_width, y =perc_insert, 
      colour = matching_class, shape = combined_5$species)) +
  geom_point(size = 3) + 
  geom_line(size = 0.5) +
  xlim(0,750) +
  ylim(0,15) +
  labs(title = "The Effects of Sequence Length on Percent Insertion in T. spiralis, T.zimbabwensis & T. nativa for Each DNA Transposon Type.",
       subtitle = "Using the Top 5 Largest DNA Transposon Sequences in Base Pairs.",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type",
       shape = "Species")


print(combined_top5)

ggsave("combined_top5.png", path = "../Small_ZZ/output")


## creating combined plot for largest transposon length per DNA transposon-----
plot_combined_largest <- ggplot() +
  geom_point(data = combined_largest, 
             aes(x =qry_width, y =perc_insert, 
                 colour = matching_class, shape = combined_largest$species,
                 size = 5)) +
  guides(size = "none") +
  labs(title = "The Effects of Sequence Length on Percent Insertion in T. spiralis, T.zimbabwensis & T. nativa for Each DNA Transposon Type.",
       subtitle = "Using the Largest DNA Transposon Sequences in Base Pairs.",
       x = "Sequence Length (bp)",
       y = "Percent Insertion",
       colour = "DNA Transposon Type",
       shape = "Species")

print(plot_combined_largest)

ggsave("combined_largest.png", path = "../Small_ZZ/output")


## 2ND SUBSTANTIVE EDIT: Statistics to analyze the trend of sequence length in relation to percent insertion ----

# spearman tests to quantify the length vs. percent insertion trend

# Trichinella spiralis
count(ts_largest, species)  
cor.test(ts_largest$qry_width, ts_largest$perc_insert, method = "spearman", exact = FALSE)
## ==========> S = 10, p-value = 0.1108 rho = 0.7142857 

# Trichinella zimbabwensis
count(tz_largest, species) 
cor.test(tz_largest$qry_width, tz_largest$perc_insert, method = "spearman", exact = FALSE)
## ==========> S = 11.664, p-value = 0.1481 rho = 0.6667367 

# Trichinella nativa
count(tn_largest, species) 
cor.test(tn_largest$qry_width, tn_largest$perc_insert, method = "spearman", exact = FALSE)
## ==========> S = 60, p-value = 0.1108 rho = -0.7142857  


## 3RD SUBSTANTIVE EDIT: Using te of all sizes, not just the largest, to compare the previous results with a new plot and spearman test ----
# Creating a table that has the information of all the te of the three species with only the pertinent information
all_common <- bind_rows(
  trich_s |> select(matching_class, qry_width, perc_insert, species),
  trich_z |> select(matching_class, qry_width, perc_insert, species),
  te_nat  |> select(matching_class, qry_width, perc_insert, species)
) |> 
  filter(matching_class %in% dte_common)

# Making sure that there is enough information 
count(all_common, species)  


## creating combined plot for all the transposon length per DNA transposon
p_cor <- ggplot(all_common, aes(x = qry_width, y = perc_insert)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ species) +
  labs(title = "Sequence Length vs. Percent Insertion in Shared DNA Transposons",
       x = "Sequence Length (bp)", y = "Percent Insertion")

print(p_cor)

ggsave("length_vs_insertion_all.png", plot = p_cor, path = "../Small_ZZ/output/")


# obtaining statistical results using all the data
all_common |>
  group_by(species) |>
  summarise(rho = cor(qry_width, perc_insert, method = "spearman"),
            p = cor.test(qry_width, perc_insert, method = "spearman", exact = FALSE)$p.value)
#==========> Trichinella nativa. rho= 0.239,  p= 0.0180
#==========> Trichinella spiralis rho= 0.193 p= 0.0302
#==========> Trichinella zimbabwensis rho= 0.135 p= 0.132 
