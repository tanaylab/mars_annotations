library(tgutil)
library(tidyverse)

umis_dir <- "/net/mraid14/export/tgdata/db/tgdb/mars_runs/stelzer_gene_annot/plates/umis/"
annots_file <- "/net/mraid14/export/tgdata/db/tgdb/mars_runs/stelzer_gene_annot/work/210426_03/scdb_oryCun2/annotations/gene_intervals_oryCun2.txt"


# +
# library(metacell)
# dataset_table_fn <- "/net/mraid14/export/tgdata/users/atanay/proj/rabemb/config/samp_key_190421.txt"
# scdb_init("scrna_db", force_reinit=T)
# mcell_import_multi_mars(mat_nm = "embs",
#                      dataset_table_fn = dataset_table_fn,
#                      base_dir = umis_dir,
#                      patch_cell_name = TRUE,
#                      force=TRUE)
# -


umi_files <- list.files(umis_dir, full.names=TRUE)

doMC::registerDoMC(20)
marginals <- plyr::llply(umi_files, function(x) tgutil::fread_rownames(x) %>% as.data.frame() %>% filter(!is.na(rowname)) %>% column_to_rownames("rowname") %>% as.matrix() %>% rowSums(), .parallel=TRUE) 

mars <- do.call('cbind', marginals) %>% rowSums() %>% enframe("gene_name", "tot_umis") %fcache_df% "tot_umis_gene.tsv" %>% as_tibble()

gene_types <- fread("gene_types.tsv") %>% 
    # keep a single organism per gene
    mutate(organism = factor(organism, levels=c("orycun2", "mm9", "hg19"))) %>% 
    arrange(gene_name, organism) %>% 
    distinct(gene_name, .keep_all=TRUE) %>% 
    as_tibble()

# The policy is: 
# 1. If there is orycun refseq - use it
# 2. If there is "other" gene from mouse or human annotation - take it (all the genes separated by ";", first mouse and then human)
# 3. If there is a single "other" name (even from other organisms) - take it. 

mars_types <- mars %>% 
    distinct(gene_name)  %>%     
    # save the orginial long name in "gene_name_long" column
    mutate(gene_name_long = gene_name) %>%
    # suggest an alternative name by removing RIKEN,LOC,ORF etc. genes    
    mutate(gene_name = gsub(".*\\d+.*\\d+RIK;*","",gsub("LOC\\d+;*","",gene_name))) %>% 
    mutate(gene_name = gsub("C\\d+ORF\\d+;*","",gsub("^\\d+.*;","",gene_name))) %>%  
    separate_rows("gene_name", sep=";") %>% 
    mutate(gene_name = ifelse(gene_name == "", NA, gene_name)) %>% 
    # add organism annotation
    left_join(gene_types, by = "gene_name")

# Add a genome suffix for the 'other' genes
mars_types <- mars_types %>%     
    mutate(gene_name = case_when(
        organism == "mm9" & !is.na(gene_name) ~ paste0(gene_name, "_mm9"),
        organism == "hg19" & !is.na(gene_name) ~ paste0(gene_name, "_hg19"),
        is.na(organism) & !is.na(gene_name) ~ paste0(gene_name, "_other"),
        TRUE ~ gene_name # for orycun genes do not add a suffix
    )) 

mars_types <- mars_types %>%             
    arrange(gene_name_long, organism) %>% 
    group_by(gene_name_long) %>% 
    # for each long name - extract the number of refseq and "other" genes
    summarise(        
        refseq = paste(sort(unique(gene_name[type == "refseq" & !is.na(type)])), collapse=";"), 
        other =  paste(sort(unique(gene_name[type == "other" & !is.na(type)])), collapse=";"),

        # if we have an organism (note that names are sorted by the organism factor) - take all the genes from it (separated by '"'', can be either orycun refseq or other mm9 and hg19, but only from a single organism)        
        new_name = ifelse(!is.na(organism[1]), paste(gene_name[organism == organism[1] & !is.na(organism)], collapse=";"), NA),

        # if there are no "other" genes - use the concatenation of the refseq genes
        new_name = ifelse(is.na(new_name) & other == "", gene_name_long, new_name)
    ) %>%     
    mutate(new_name = ifelse(new_name == '', NA, new_name))

# if we have a single "other" name - use it. 
mars_types <- mars_types %>%         
    mutate(new_name = ifelse(is.na(new_name) & !grepl(";", other) & !is.na(other), other, new_name)) 

# save the data frame
mars_types <- mars_types %>%     
    rename(gene_name=gene_name_long) %>%     
    right_join(mars) %fcache_df% "tot_umis_gene_types.tsv"

trans_tab <- mars_types %>% 
    mutate(more_than_one_gene = grepl(";", new_name), other_organism = grepl("_mm9", new_name) | grepl("_hg19", new_name)) %>% 
    filter(tot_umis > 10 & (is.na(new_name) | (more_than_one_gene & other_organism))) %>% 
    select(-more_than_one_gene, -other_organism) %>% 
    arrange(desc(tot_umis)) %fcache_df% "missing_gene_annots.tsv"

# At this point we went manually and curated the translation table. We now merge the output with a suffix of "_m"
trans_table_manual <- readxl::read_xlsx("missing_gene_annots.xlsx")

# Use the following heuristic - take the names with the least number characters excluding names with parenthesis
trans_table_curated <- trans_table_manual %>% 
    separate_rows(new_name, sep=";")  %>% 
    mutate(bad = grepl("\\(", new_name)) %>% 
    mutate(numbers = map(str_extract_all(gsub("_hg19", "", new_name) %>% gsub("_mm9", "", .), "(\\d+)"), str_length) %>% map_dbl(max)) %>% 
    arrange(gene_name, as.numeric(bad), numbers)  %>% 
    group_by(gene_name) %>% 
    summarise(ADD = new_name[1])

writexl::write_xlsx(trans_table_curated, "missing_gene_annots_curated.xlsx")

mars_types <- mars_types %>% 
    left_join(trans_table_curated %>% select(gene_name, ADD) %>% mutate(ADD = paste0(ADD, "_m"))) %>% 
    mutate(new_name = ifelse(!is.na(ADD), ADD, new_name)) %>% 
    select(-ADD)

# remove not-curated genes that are not sufficiently expressed
mars_types <- mars_types %>% filter(!(tot_umis <= 10 & is.na(new_name))) %fcache_df% "tot_umis_gene_types_curated.tsv"

gene_annots <- fread(annots_file) %>% as_tibble()

# we remove genes that were not expressed (less than 11 umis) and didn't come from refseq
new_gene_annots <- gene_annots %>% 
    left_join(mars_types %>% select(gene_name, new_name, tot_umis), by = "gene_name") %>% 
    mutate(other_organism = grepl("_mm9", new_name) | grepl("_hg19", new_name) | grepl("_other", new_name)) %>% 
    filter(tot_umis > 10 | !other_organism) %>% 
    select(chrom, start, end, strand, gene_name = new_name) %>% 
    arrange(chrom, start, end, strand) %>% 
    mutate(start = as.numeric(start), end = as.numeric(end), strand = as.numeric(strand)) %>% 
    filter(!is.na(start), !is.na(gene_name), gene_name != "gene_name")    

fwrite(new_gene_annots, "curated_gene_annots.tsv", sep="\t")

ercc_intervs <- fread('ERCC_intervals.tsv.gz') %>% as_tibble()

fwrite(bind_rows(new_gene_annots, ercc_intervs), "gene_intervals_oryCun2_with_ERCC_curated.txt", sep="\t", quote=FALSE)

