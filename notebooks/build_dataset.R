#!/usr/bin/env Rscript

library(tidyverse)

samples <- read_tsv("input_data/samples.tsv.gz")
proteins <- read_tsv("input_data/proteins.tsv.gz")

samples <- samples %>%
           filter(
               `Characteristics[sex]` == "female",
               `Characteristics[disease]` %in% c("normal", "systemic lupus erythematosus")
           ) %>%
           extract(`Source Name`, "id", "(\\d+)", convert = TRUE) %>%
           transmute(
               id,
               country = recode(`Characteristics[cohort]`,
                                "United Kingdom" = "UK",
                                "United States of America" = "USA"),
               ethnicity = recode(`Characteristics[ethnic group]`,
                                  "African American" = "Afro-American",
                                  "UK Afro-Caribbean" = "Afro-Caribbean"),
               case = `Characteristics[disease]` == "systemic lupus erythematosus"
           ) %>%
           arrange(id)

protein_ids <- proteins$`Reporter Identifier`
proteins <- t(as.matrix(proteins %>% select(-`Reporter Identifier`)))
colnames(proteins) <- protein_ids

proteins <- as_tibble(proteins, rownames = "id") %>%
            extract(id, "id", "(\\d+)", convert = TRUE) %>%
            select_if(colSums(is.na(.)) == 0)

samples <- samples %>%
           inner_join(proteins, by = "id") %>%
           arrange(id)

write_rds(samples, "datasets/sle_proteomics.rds", compress = "gz")

