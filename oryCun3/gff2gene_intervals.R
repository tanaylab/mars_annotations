library(tidyverse)
library(tgutil)

gff <- fread("lifted.gff3", skip=3)

gene_intervals <- fread("lifted.gff3", skip=3) %>% select(chrom=V1, start=V4, end=V5, strand=V7, gene_name=V2) %>% mutate_at(vars(start, end), as.numeric) %>% mutate(strand = ifelse(strand == "+", +1, -1))

fwrite(gene_intervals, "gene_intervals_oryCun3.txt.gz", sep="\t", quote=FALSE)

ercc_intervs <- fread('../oryCun2/ERCC_intervals.tsv.gz') %>% as_tibble()

fwrite(bind_rows(gene_intervals, ercc_intervs), "gene_intervals_oryCun3_with_ERCC.txt", sep="\t", quote=FALSE)

