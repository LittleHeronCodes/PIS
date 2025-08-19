#' Get KEGG pathways as a named list (pathway name/ID as name, gene vector as entry)
#'
#' @param species Character. KEGG species code (e.g., 'hsa' for human)
#' @return Named list: names are pathway IDs, entries are gene vectors
#' @import KEGGREST
#' @export
get_pathways_KEGG <- function(species = "hsa") {
  kegg_pathways <- KEGGREST::keggList("pathway", species)
  kegg_list <- list()
  for (pid in names(kegg_pathways)) {
    genes <- KEGGREST::keggGet(pid)[[1]]$GENE
    if (is.null(genes)) next
    gene_ids <- genes[seq(1, length(genes), 2)]
    kegg_list[[pid]] <- gene_ids
  }
  kegg_list
}

#' Get Reactome pathways as a named list (pathway name/ID as name, gene vector as entry)
#'
#' @param species Character. Reactome species code (e.g., 'Homo sapiens')
#' @return Named list: names are pathway IDs, entries are gene vectors
#' @import httr jsonlite
#' @export
get_pathways_Reactome <- function(species = "hsa") {
  reactome_url <- paste0("https://reactome.org/ContentService/data/pathways/low/", toupper(species))
  reactome_res <- httr::GET(reactome_url)
  reactome_list <- list()
  if (httr::status_code(reactome_res) == 200) {
    reactome_data <- jsonlite::fromJSON(httr::content(reactome_res, as = "text"))
    for (pw in reactome_data) {
      pw_id <- pw$dbId
      gene_url <- paste0("https://reactome.org/ContentService/data/pathway/", pw_id, "/participants")
      gene_res <- httr::GET(gene_url)
      if (httr::status_code(gene_res) == 200) {
        gene_data <- jsonlite::fromJSON(httr::content(gene_res, as = "text"))
        gene_symbols <- unique(unlist(lapply(gene_data, function(x) x$geneName)))
        if (length(gene_symbols) > 0) {
          reactome_list[[as.character(pw_id)]] <- gene_symbols
        }
      }
    }
  }
  reactome_list
}

#' Get GO terms as a named list (GO ID as name, gene vector as entry)
#'
#' @param species Character. KEGG species code (e.g., 'hsa' for human)
#' @param go_ontology Character. GO ontology: 'BP', 'MF', or 'CC'.
#' @return Named list: names are GO IDs, entries are gene vectors
#' @import biomaRt
#' @export
get_pathways_GO <- function(species = "hsa", go_ontology = "BP") {
  mart <- biomaRt::useMart("ensembl", dataset = paste0(species, "_gene_ensembl"))
  go_terms <- biomaRt::getBM(
		attributes = c("go_id", "external_gene_name", "namespace_1003"),
    filters = "namespace_1003", values = go_ontology, mart = mart
	)
  go_split <- split(go_terms, go_terms$go_id)
  go_list <- lapply(go_split, function(df) unique(df$external_gene_name))
  go_list
}

#' Get all pathway lists (KEGG, Reactome, GO) as a named list by ontology
#'
#' @param species Character. KEGG species code (e.g., 'hsa' for human)
#' @param go_ontology Character. GO ontology: 'BP', 'MF', or 'CC'.
#' @return Named list with ontologies as names and pathway lists as entries
#' @export
get_pathway_list <- function(species = "hsa", go_ontology = "BP") {
  list(
    KEGG = get_pathways_KEGG(species),
    Reactome = get_pathways_Reactome(species),
    GO = get_pathways_GO(species, go_ontology)
  )
}
