#!/usr/bin/env Rscript
library(dplyr, warn.conflicts = FALSE)

.misha_root <- "/net/mraid14/export/tgdata/db/tgdb/orycun2/trackdb"


########################################################################
main <- function(argv) {
    if (length(argv) != 2) {
        cat("Usage: gene_intervals.R  <refseq.txt.gz>  <output.tsv>\n", file = stderr())
        return(2)
    }

    misha::gdb.init(.misha_root)

    gene_tab <- read_refseq(argv[1])
    gene_tab <- merge_genes(gene_tab)

    data.table::fwrite(gene_tab, argv[2], sep = "\t")

    return(0)
}


########################################################################
read_refseq <- function(fname) {
    cnames <- c(
        "bin",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "score",
        "name2",
        "cdsStartStat",
        "cdsEndStat",
        "exonFrames"
    )

    refseq <- data.table::fread(fname, header = FALSE, col.names = cnames, data.table = FALSE)

    refseq <- refseq %>% filter(grepl("^chr[0-9X]+$", chrom))
    refseq <- refseq %>% mutate(strand = case_when(strand == "+" ~ +1L, strand == "-" ~ -1L))
    refseq <- refseq %>%
              mutate(exonStarts = stringr::str_remove(exonStarts, ",$"),
                     exonEnds = stringr::str_remove(exonEnds, ",$")) %>%
              tidyr::separate_rows(exonStarts, exonEnds, sep = ",") %>%
              mutate(exonStarts = as.integer(exonStarts), exonEnds = as.integer(exonEnds))
    refseq <- refseq %>%
              select(chrom, start = exonStarts, end = exonEnds, strand, gene_name = name2)

    return(refseq)
}


########################################################################
merge_genes <- function(gene_tab) {
    pairs <- bind_rows(
        gene_pairs(gene_tab %>% filter(strand == +1)),
        gene_pairs(gene_tab %>% filter(strand == -1))
    )
    pairs <- pairs %>%
             distinct()
    pairs <- igraph::graph_from_edgelist(as.matrix(pairs), directed = FALSE)
    pairs <- igraph::components(pairs)
    pairs <- tibble::tibble(gene_name = names(pairs$membership), component = pairs$membership)
    pairs <- pairs %>%
             group_by(component) %>%
             arrange(gene_name) %>%
             mutate(merged_name = paste0(gene_name, collapse = ";")) %>%
             ungroup() %>%
             select(-component)

    stopifnot(nrow(pairs %>% group_by(gene_name) %>% filter(n() > 1)) == 0)

    gene_tab <- gene_tab %>%
                left_join(pairs, by = "gene_name") %>%
                mutate(gene_name = if_else(is.na(merged_name), gene_name, merged_name)) %>%
                select(-merged_name)

    gene_tab <- gene_tab %>%
                group_by(gene_name, strand) %>%
                do(merge_intervals(.)) %>%
                ungroup() %>%
                select(chrom, start, end, strand, gene_name)

    tab_order <- gene_tab %>%
                 mutate(chrom = factor(chrom, levels = levels(misha::gintervals.all()$chrom), ordered = TRUE)) %>%
                 group_by(gene_name) %>%
                 summarize(chrom = min(chrom), start = min(start), .groups = "drop") %>%
                 arrange(chrom, start)
    stopifnot(all(!is.na(tab_order$chrom)))

    gene_tab <- tab_order %>%
                distinct(gene_name) %>%
                left_join(gene_tab, by = "gene_name") %>%
                select(chrom, start, end, strand, gene_name)

    stopifnot(!has_overlap(gene_tab %>% filter(strand == +1)))
    stopifnot(!has_overlap(gene_tab %>% filter(strand == -1)))

    return(gene_tab)
}


########################################################################
gene_pairs <- function(gene_tab) {
    canonic_intervals <- misha::gintervals.canonic(gene_tab, unify_touching_intervals = TRUE)
    gene_tab <- gene_tab %>% mutate(canonic = attr(canonic_intervals, "mapping"))

    pairs <- gene_tab %>%
             distinct(gene_name, canonic)
    pairs <- pairs %>%
             rename(gene1 = gene_name) %>%
             left_join(pairs %>% rename(gene2 = gene_name), by = "canonic") %>%
             distinct(gene1, gene2) %>%
             filter(gene1 < gene2)

    return(pairs)
}


########################################################################
merge_intervals <- function(gene_tab) {
    return(misha::gintervals.canonic(gene_tab, unify_touching_intervals = TRUE))
}


########################################################################
has_overlap <- function(gene_tab) {
    canonic_intervals <- misha::gintervals.canonic(gene_tab, unify_touching_intervals = TRUE)
    return(nrow(canonic_intervals) != nrow(gene_tab))
}


########################################################################
if (sys.nframe() == 0) {
    rc <- main(commandArgs(trailingOnly = TRUE))
    if (is.null(rc)) {
        rc <- 0
    }
    if (!is.numeric(rc)) {
        cat(rc, "\n", file = stderr(), sep = "")
        rc <- 1
    }
    rc <- floor(rc)
    q(save = "no", status = rc)
}