getID <- function(symbol,DB=Symbol2ID){
  error_symbol = c("ABP6","LARS","STSL2","ANX10","CARS","WARS")
  right_symbol = c("FABP6","LARS1","ABCG5","ANXA10","CARS1","WARS1")
  symbol[symbol %in% error_symbol] = right_symbol[error_symbol %in% symbol]
  DB[symbol]
}

#' @title Return meta Data
#' Get HCC risk-score related metaData.
#' @examples 
#' \dontrun{
#' meta()
#' metaDf = meta()
#' }
#' @export
meta <- function(){
  metaData
}


#' @title score
#' Accept a dataframe of HCC RNA data, and return 48 risk-scores in a new dataframe.
#' @section Waining:
#' The missing gene expression will be fill with 0.
#' @param data RNA expression dataframe, samples as in column and gene Ensembl ID as row name.
#' @examples 
#' \dontrun{
#' df = read.delim("GEOdata.txt")
#' score_df = score(df)
#' }
#' @export
score <- function(data){
  S = as.data.frame(matrix(NA,48,ncol(data)))
  colnames(S) = colnames(data)
  rownames(S) = paste0('score',1:48)

  # score 
  for (i in 1:48) {
    G = as.matrix(data[getID(Genes[[i]]),]) 
    G[is.na(G)] = 0
    S[i,] = Weight[[i]] %*% G
    if (i >= 44) {
      S[i,] = exp(S[i,]) 
    }
  }
  return(S)
}

