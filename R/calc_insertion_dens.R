#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom rlang .data

################################################################################
############################## format transposon data ##########################
################################################################################

calculate_insertion_density <- function(file_path){
  nas_gone <- utils::read.csv(file_path, header = T) %>%
    stats::na.omit()

  gene_level_data <- nas_gone[!duplicated(nas_gone), ] %>%
    dplyr::group_by(rlang::.data$gene) %>%
    dplyr::summarise(n = n(), length = first(rlang::.data$length)) %>%
    dplyr::transmute(gene = rlang::.data$gene,
                     insdens = (n / length))
  return(gene_level_data)
}
