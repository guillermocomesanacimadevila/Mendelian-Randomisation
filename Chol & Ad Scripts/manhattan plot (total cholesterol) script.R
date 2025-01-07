install.packages("tidyverse")
install.packages("qqman")

library(tidyverse)
library(qqman)

setwd("/Users/guillermocomesanacimadevila/Desktop/mendelian randomisation/Data")

tc_path = "/Users/guillermocomesanacimadevila/Desktop/mendelian randomisation/Data/Willer2013tc.chrall.CPRA_b37.tsv"

coltypes = cols(
  ID = col_character(),
  CHROM = col_double(),
  POS = col_double(),
  REF = col_character(),
  ALT = col_character(),
  AF = col_double(),
  TRAIT = col_character(),
  BETA = col_double(),
  SE = col_double(),
  Z = col_double(),
  P = col_double(),
  N = col_double(),
  OR = col_double(),
  OR_L95 = col_double(),
  OR_U95 = col_double(),
  DIR = col_character(),
  G1000_ID = col_character(),
  G1000_VARIANT = col_character(),
  DBSNP_ID = col_character(),
  DBSNP_VARIANT = col_character(),
  OLD_ID = col_character(),
  OLD_VARIANT = col_character()
)

tc_data <- read_tsv(tc_path, comment = "##", col_types = coltypes, 
                  col_select = c(DBSNP_ID, CHROM, POS, REF, ALT, AF, BETA, SE, Z, P, N, TRAIT))

# filter for p < 0.05
filtered_data <- filter(tc_data, P < 0.05 & P > 1e-100)
manhattan_data <- filtered_data %>%
  rename(
    chr = CHROM,
    bp = POS,
    p = P,
    snp = DBSNP_ID
  )

custom_colors <- c("skyblue", "tomato1")
manhattan(
  manhattan_data,
  chr = "chr",
  bp = "bp",
  p = "p",
  snp = "snp",
  col = custom_colors,
  genomewideline = -log10(5e-8),   
  suggestiveline = -log10(1e-5), 
  main = "Manhattan Plot"
)