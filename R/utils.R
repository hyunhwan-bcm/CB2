#' Repeat a row
#'
#' @param x A vector to repeat
#' @param n Number of repetitions
#' @return A matrix with `n` copies of `x`
#' @export
rep_row <- function(x, n) matrix(rep(x, each = n), nrow = n)

#' A function to normalize sgRNA read counts.
#' 
#' @param sgcount The input table contains read counts of sgRNAs for each sample
#'
#' A function to calculate the CPM (Counts Per Million) (required)
#' 
#' @return a normalized CPM table will be returned
#' 
#' @examples
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' get_CPM(Evers_CRISPRn_RT112$count)
#' 
#' @export
get_CPM <- function(sgcount) {
    cols <- which(sapply(sgcount, class) == "numeric")
    nmat <- rep_row(colSums(sgcount[,cols]), nrow(sgcount[,cols]))
    sgcount[,cols] <- sgcount[,cols]/nmat * 10^6
    sgcount
}

#' A function to plot the first two principal components of samples.
#'
#' This function will perform a principal component analysis, and it returns a ggplot object of the PCA plot.
#' 
#' @param sgcount The input matrix contains read counts of sgRNAs for each sample.
#' @param df_design The table contains a study design.
#'
#' @importFrom magrittr %>%
#' @return A ggplot2 object contains a PCA plot for the input.
#' 
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' plot_PCA(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design)
#'  
#' @export
plot_PCA <- function(sgcount, df_design) {
    cols <- which(sapply(sgcount, class) == "numeric")
    pca_obj <- sgcount[,cols] %>% t %>% prcomp
    importance <- summary(pca_obj)$importance
    prop_pc1 <- importance[2,1]
    prop_pc2 <- importance[2,2]
    pca_obj$x %>% as.data.frame() %>% 
        tibble::rownames_to_column("sample_name") %>% 
        dplyr::left_join(df_design, by = "sample_name") %>% 
        ggplot2::ggplot(ggplot2::aes_string(x = "PC1", y = "PC2")) + 
        ggplot2::geom_point(ggplot2::aes_string(color = "group"), size = 2) + 
        ggplot2::geom_text(ggplot2::aes_string(label = "sample_name")) +
        ggplot2::xlab(sprintf("PC1 (%.2f%%)", prop_pc1*100)) +
        ggplot2::ylab(sprintf("PC2 (%.2f%%)", prop_pc2*100))
}

#' A function to show a heatmap sgRNA-level corrleations of the NGS samples.
#' @param sgcount The input matrix contains read counts of sgRNAs for each sample.
#' @param df_design The table contains a study design.
#' @param cor_method A string parameter of the correlation measure. One of the three - "pearson", "kendall", or "spearman" will be the string. 
#'
#' @importFrom magrittr %>%
#' @importFrom pheatmap pheatmap
#' @importFrom stats cor
#' @return A pheatmap object contains the correlation heatmap
#' 
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' plot_corr_heatmap(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design)
#' @export
plot_corr_heatmap <- function(sgcount, df_design, cor_method = "pearson") {
    cols <- which(sapply(sgcount, class) == "numeric") 
    sgcount[,cols] %>% cor(method = cor_method) %>% 
        pheatmap::pheatmap(display_numbers = T, 
                           number_format = "%.2f", 
                           annotation_col = df_design %>%
                               tibble::column_to_rownames("sample_name") %>%
                               dplyr::select(dplyr::all_of("group")))    
}

#' A function to calculate the mappabilities of each NGS sample.
#' 
#' @param count_obj A list object is created by `run_sgrna_quant`.
#' @param df_design The table contains a study design.
#' @importFrom magrittr %>%
#' @examples
#' library(CB2)
#' library(magrittr)
#' library(tibble)
#' library(dplyr)
#' library(glue)
#' FASTA <- system.file("extdata", "toydata", "small_sample.fasta", package = "CB2")
#' ex_path <- system.file("extdata", "toydata", package = "CB2")
#' 
#' df_design <- tribble(
#'   ~group, ~sample_name,
#'   "Base", "Base1",  
#'   "Base", "Base2", 
#'   "High", "High1",
#'   "High", "High2") %>% 
#'     mutate(fastq_path = glue("{ex_path}/{sample_name}.fastq"))
#' 
#' cb2_count <- run_sgrna_quant(FASTA, df_design)
#' calc_mappability(cb2_count, df_design)
#'
#' @export
calc_mappability <- function(count_obj, df_design) {
    csum <- count_obj$count %>% colSums()
    mp <- csum/count_obj$total * 100
    df_design %>%
        dplyr::mutate(total_reads = count_obj$total) %>%
        dplyr::mutate(mapped_reads = csum) %>%
        dplyr::mutate(mappability = mp) %>%
        dplyr::select(-dplyr::any_of("fastq_path"))
}

#' A function to join a count table and a design table.
#'
#' @param sgcount The input matrix contains read counts of sgRNAs for each sample.
#' @param df_design The table contains a study design.
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @return A tall-thin and combined table of the sgRNA read counts and study design will be returned.
#'
#' @examples 
#' library(CB2)
#' data(Evers_CRISPRn_RT112) 
#' head(join_count_and_design(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design))
#' 
#' @export
join_count_and_design <- function(sgcount, df_design) {
    cols <- colnames(sgcount)
    
    if(all(sapply(sgcount, class) == "numeric")) {
        sgcount %>% as.data.frame(stringsAsFactors=F) %>%
            tibble::rownames_to_column("sgRNA") %>%
            tidyr::pivot_longer(cols = dplyr::all_of(cols),
                                names_to = "sample_name",
                                values_to = "count") %>%
            dplyr::left_join(df_design, by = "sample_name")
    } else {
        cols <- colnames(sgcount)[sapply(sgcount, class) == "numeric"]
        sgcount %>% as.data.frame(stringsAsFactors=F) %>%
            tidyr::pivot_longer(cols = dplyr::all_of(cols),
                                names_to = "sample_name",
                                values_to = "count") %>%
            dplyr::left_join(df_design, by = "sample_name")
    }
}

#' A function to plot read count distribution.
#' 
#' @param sgcount The input matrix contains read counts of sgRNAs for each sample.
#' @param df_design The table contains a study design.
#' @param add_dots The function will display dots of sgRNA counts if it is set to `TRUE`.
#' @importFrom magrittr %>%
#' @return A ggplot2 object contains a read count distribution plot for `sgcount`.
#' 
#' @examples
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' cpm <- get_CPM(Evers_CRISPRn_RT112$count)
#' plot_count_distribution(cpm, Evers_CRISPRn_RT112$design)
#' 
#' @export
plot_count_distribution <- function(sgcount, df_design, add_dots = FALSE) {
    p <- join_count_and_design(sgcount, df_design) %>%
        dplyr::mutate(count = log2(1+count)) %>% 
        ggplot2::ggplot(ggplot2::aes_string(y="count", x="sample_name")) + 
        ggplot2::geom_violin(ggplot2::aes_string(fill = "group")) + 
        ggplot2::ylab("log2(1+count)")
    
    if( add_dots == T ) {
        p <- p + ggplot2::geom_jitter(width=0.1, alpha=0.5)   
    }
    
    p
}

#' A function to visualize dot plots for a gene.
#' @param sgcount The input matrix contains read counts of sgRNAs for each sample.
#' @param df_design The table contains a study design.
#' @param gene The gene to be shown.
#' @param ge_id A name of the column contains gene names.
#' @param sg_id A name of the column contains sgRNA IDs.
#' 
#' @importFrom stats prcomp as.formula
#' @return A ggplot2 object contains dot plots of sgRNA read counts for a gene.
#' 
#' @examples 
#' library(CB2)
#' data(Evers_CRISPRn_RT112)
#' plot_dotplot(get_CPM(Evers_CRISPRn_RT112$count), Evers_CRISPRn_RT112$design, "RPS7")
#' 
#' @export 
plot_dotplot <- function(sgcount, df_design, gene, ge_id = NULL, sg_id = NULL) {
    if(all(sapply(sgcount, class) == "numeric")) {
        if(sum(stringr::str_detect(rownames(sgcount), glue::glue("^{gene}")))==0) {
            stop(glue::glue("{gene} is not in sgcount."))
        }
        
        join_count_and_design(sgcount, df_design) %>%
            dplyr::filter(stringr::str_detect(.data$sgRNA, glue::glue("^{gene}"))) %>% 
            ggplot2::ggplot(ggplot2::aes_string(x = "group", y = "count")) + 
            ggplot2::geom_dotplot(ggplot2::aes_string(fill = "group", color = "group"), binaxis = "y", stackdir = "center", stackratio = 1.5, dotsize = 1.2) + 
            ggplot2::facet_wrap(~sgRNA, scales = "free_y") + ggplot2::ggtitle(gene)
    } else {
        if(is.null(ge_id) || is.null(sg_id)) {
            stop("ge_id and sg_id should not be null.")
        }
        if(!(ge_id %in% colnames(sgcount))) {
            stop("ge_id is not found in sgcount.")
        } 
        if(!(sg_id %in% colnames(sgcount))) {
            stop("sg_id is not found in sgcount.")
        } 
        
        df <- join_count_and_design(sgcount, df_design)
        df <- df[df[,ge_id]==gene,]
        if(nrow(df) == 0) {
            stop(glue::glue("{gene} is not in sgcount."))
        }
        df %>% ggplot2::ggplot(ggplot2::aes_string(x = "group", y = "count")) + 
            ggplot2::geom_dotplot(ggplot2::aes_string(fill = "group", color = "group"), binaxis = "y", stackdir = "center", stackratio = 1.5, dotsize = 1.2) + 
            ggplot2::facet_wrap(as.formula(paste("~", sg_id)), scales = "free_y") + ggplot2::ggtitle(gene)
    }
    
}

