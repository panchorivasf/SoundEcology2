#' Calculate Diversity Indices
#'
#' This function from the `vegan` package calculates various diversity indices (Shannon, Simpson, and Inverse Simpson) for a given community data matrix.
#' Reference: Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin
#' P, O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour
#'  M, Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De
#'  Caceres M, Durand S, Evangelista H, FitzJohn R, Friendly M,
#'  Furneaux B, Hannigan G, Hill M, Lahti L, McGlinn D, Ouellette
#'  M, Ribeiro Cunha E, Smith T, Stier A, Ter Braak C, Weedon J
#'  (2022). _vegan: Community Ecology Package_. R package version
#'  2.6-4, <https://CRAN.R-project.org/package=vegan>.
#'
#' @param x A numeric matrix or data frame where rows typically represent sites or samples and columns represent species or categories.
#' @param index A character string specifying the diversity index to compute. Options are `"shannon"` (default), `"simpson"`, or `"invsimpson"`.
#' @param groups A vector or factor specifying groupings for pooling sites or units. If provided, diversity indices are calculated by these groups.
#' @param MARGIN An integer indicating whether to compute indices by rows (`MARGIN = 1`, default) or by columns (`MARGIN = 2`).
#' @param base The logarithm base used for calculating the Shannon index. Default is the natural logarithm (base `exp(1)`).
#'
#' @return A numeric vector of diversity index values for each site, sample, or group. Returns \code{NA} for any row or column containing all zero values.
#' @export
#'
#' @examples
#' # Example data: a matrix with species counts for different sites
#' data <- matrix(c(10, 20, 10, 15, 0, 5, 15, 25), nrow = 2)
#' colnames(data) <- c("Species1", "Species2", "Species3", "Species4")
#' rownames(data) <- c("Site1", "Site2")
#'
#' # Shannon diversity index
#' diversity(data, index = "shannon")
#'
#' # Simpson diversity index
#' diversity(data, index = "simpson")
#'
#' # Inverse Simpson diversity index with grouping
#' groups <- c("Group1", "Group1")
#' diversity(data, index = "invsimpson", groups = groups)
#'
#' # Using MARGIN to calculate diversity by columns
#' diversity(data, index = "shannon", MARGIN = 2)
diversity <- function(x, 
                      index = "shannon", 
                      groups, 
                      MARGIN = 1, 
                      base = exp(1)) {
  x <- drop(as.matrix(x))
  if (!is.numeric(x))
    stop("input data must be numeric")
  
  if (any(x < 0, na.rm = TRUE))
    stop("input data must be non-negative")
  
  ## sum communities for groups
  if (!missing(groups)) {
    if (MARGIN == 2)
      x <- t(x)
    
    if (length(groups) == 1) # total for all SU
      groups <- rep(groups, NROW(x))
    
    x <- aggregate(x, list(groups), sum) # pool SUs by groups
    
    rownames(x) <- x[,1]
    x <- x[,-1, drop=FALSE]
    
    if (MARGIN == 2)
      x <- t(x)
  }
  
  INDICES <- c("shannon", "simpson", "invsimpson")
  
  index <- match.arg(index, INDICES)
  
  if (length(dim(x)) > 1) {
    total <- apply(x, MARGIN, sum)
    x <- sweep(x, MARGIN, total, "/")
  } else {
    x <- x/(total <- sum(x))
  }
  
  if (index == "shannon")
    x <- -x * log(x, base)
  else
    x <- x*x
  if (length(dim(x)) > 1)
    H <- apply(x, MARGIN, sum, na.rm = TRUE)
  else
    H <- sum(x, na.rm = TRUE)
  if (index == "simpson")
    H <- 1 - H
  else if (index == "invsimpson")
    H <- 1/H
  ## check NA in data
  if (any(NAS <- is.na(total)))
    H[NAS] <- NA
  H
}