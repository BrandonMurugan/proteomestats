#' ID Converter Function
#'
#' This function takes an annotation package, and conbversion types as input to create an ID converter function
#' This code was modified from the SetRank package (Cedric Simillion (2016). SetRank: Advanced Gene Set Enrichment Analysis. R package version 1.1.0. https://CRAN.R-project.org/package=SetRank)
#'
#'
#' @param annotationPackageName Annotation package name. Format = String.
#' @param from Input ID type. See keytypes(<annotation package name>) for annotation package-specific ID types
#' @param to Output ID type. See keytypes(<annotation package name>) for annotation package-specific ID types
#'
#'
#' @author Brandon D. Murugan

createIDconverter2 <- function (annotationPackageName, from, to)
{
  require(annotationPackageName, character.only = TRUE)
  keySet = keys(eval(parse(text = annotationPackageName)),
                from)
  conversionTable = AnnotationDbi::select(eval(parse(text = annotationPackageName)),
                           keys = keySet, keytype = from, columns = to)
  conversion = as.list(by(conversionTable, as.factor(conversionTable[[from]]),
                          function(t) unique(t[[to]]), simplify = FALSE))
  outputFunction = function(x, na.rm = TRUE, drop.ambiguous = FALSE) {
    knownX = unique(x)
    outputList = as.list(rep(NA, length(knownX)))
    names(outputList) = knownX
    knownX = intersect(knownX, names(conversion))
    if (drop.ambiguous) {
      ambiguous = sapply(conversion, function(x) length(x) >
                           1)
      knownX = setdiff(knownX, names(conversion[ambiguous]))
    }
    outputList[knownX] = conversion[knownX]
    output = unlist(outputList[x], use.names = FALSE)
    if (na.rm) {
      output = output[!is.na(output)]
    }
    output
  }
  return(outputFunction)
}
