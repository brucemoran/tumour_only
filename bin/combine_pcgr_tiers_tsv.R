#! R

#' Take PCGR and CPSR TSVs from current dir and combine per sample
#' naming comes from input files, each single file generated by PCGR
#'
#' @param path character, path to multiple *tiers.tsv files (default: getwd())
#' @param pattern character, pattern to match files (default: "tiers.tsv$")
#' @param delim character, delimiter to split filenames, first element is name for output (default: "\\.")
#' @param out_tag character, name of output file (default: NULL)

combine_tsvs <- function(path = getwd(),
                         pattern = "tiers.tsv$",
                         delim = "\\.",
                         out_tag){

  paths_vec <- dir(path = path,
                   pattern = pattern,
                   full.names = TRUE)
  files_vec <- dir(path = path,
                   pattern = pattern,
                   full.names = FALSE)

  paste(strsplit(path_vec[1], delim)[[1]][-1], collapse = gsub("\\\\", "", delim))

  ##find CPSR, PCGR files
  cpsr_ind <- pcgr_ind <- c()

  for(x in 1:length(files_vec)){
    ifelse(!is.na(match("cpsr", strsplit(files_vec[x], delim)[[1]][2])),
           cpsr_ind <- c(cpsr_ind, x),
           pcgr_ind <- c(pcgr_ind, x))
  }

  pcgr_list <- process_tsvs(path_vec = paths_vec[pcgr_ind],
                            delim = delim,
                            pcgr = TRUE)
  cpsr_list <- process_tsvs(path_vec = paths_vec[cpsr_ind],
                            delim = delim,
                            pcgr = FALSE)


  pcgr_tb <- dplyr::left_join(pcgr_list[[1]], pcgr_list[[2]], by = "GENOMIC_CHANGE")
  if(dim(pcgr_tb)[1] > 0){
    readr::write_tsv(x = pcgr_tb,
                     file = paste0(out_tag, ".pcgr.combined.tsv"))

  } else {
    print("No results for PCGR found...")
  }

  cpsr_tb <- dplyr::left_join(cpsr_list[[1]], cpsr_list[[2]], by = "GENOMIC_CHANGE")
  if(dim(cpsr_tb)[1] > 0){
    readr::write_tsv(x = cpsr_tb,
                     file = paste0(out_tag, ".cpsr.combined.tsv"))
  } else {
    print("No results for CPSR found...")
  }

}

process_tsvs <- function(path_vec, delim, pcgr){

  tsv_gc_list <- lapply(path_vec, function(f){
    nm <- strsplit(x = f,
                   split = delim)[[1]][1]
    tf <- readr::read_tsv(f)
    gc <- tf$GENOMIC_CHANGE
    list("tsv" = tf, "GENOMIC_CHANGE" = gc)
  })

  ##create named list of tsv input
  name_func <- function(f){
    unlist(lapply(f, function(ff){
      strsplit(x = rev(strsplit(ff, "/")[[1]][-1]),
               split = delim)[[1]][1]
    }))}
  names(tsv_gc_list) <- name_func(path_vec)

  ##take all GENOMIC_CHANGE as unique ID, tabulate, then per-sample
  all_gc_tab <- table(unlist(lapply(tsv_gc_list, function(f){
    f[[2]]
  })))
  all_gc_tb <- tibble::tibble(GENOMIC_CHANGE = names(all_gc_tab),
                              n_samples = c(all_gc_tab))

  all_gc_smp <- purrr::reduce(lapply(seq_along(tsv_gc_list), function(f){
    tibble::tibble(GENOMIC_CHANGE = tsv_gc_list[[f]][[2]],
                   sampleID = names(tsv_gc_list)[f])
  }),
  dplyr::full_join,
  by = "GENOMIC_CHANGE")

  for(x in 1:dim(all_gc_smp)[1]){
    all_gc_smp$sampleID.x[x] <- gsub(",NA|NA,", "",
                                     paste(all_gc_smp[x, 2:dim(all_gc_smp)[2]],
                                           collapse = ","))
  }

  all_gc_tbsmp <- dplyr::left_join(dplyr::select(.data = all_gc_smp,
                                                 GENOMIC_CHANGE,
                                                 sampleID = sampleID.x),
                                   all_gc_tb)

  all_tsv_tb <- tibble::as_tibble(
    data.table::rbindlist(
      lapply(tsv_gc_list, function(f){
        return(f[[1]])
      })))

  if(pcgr == TRUE){
    all_tsv_tbu <- unique(dplyr::select(.data = all_tsv_tb,
                                        SYMBOL,
                                        GENE_NAME,
                                        HGVSp,
                                        TIER,
                                        CLINVAR,
                                        CLINVAR_CLNSIG,
                                        DP_TUMOR,
                                        AF_TUMOR,
                                        tidyselect::everything(),
                                        -VCF_SAMPLE_ID))
  } else {
    all_tsv_tbu <- unique(dplyr::select(.data = all_tsv_tb,
                                        SYMBOL,
                                        GENE_NAME,
                                        HGVSp,
                                        CONSEQUENCE,
                                        DBSNP,
                                        CLINVAR_CLASSIFICATION,
                                        CLINVAR_PHENOTYPE,
                                        tidyselect::everything(),
                                        -VCF_SAMPLE_ID,
                                        -GENOTYPE))

  }
  return(list(tbsmp = all_gc_tbsmp, tbu = all_tsv_tbu))
}
