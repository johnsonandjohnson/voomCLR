#' @include utils.R


#' @title Plot distribution of coefficients across features.
#' @name plotBeta
#' @param fit The result obtained after running \code{lmFit}.
#' @examples
#' #' @examples
#' library(limma)
#' set.seed(495212344)
#' n <- 40 # sample size
#' P <- 10 # number of cell types
#' mu0 <- rnbinom(n = P, size = 1 / 2, mu = 400)
#' mu0 # absolute counts in group 0
#' beta <- rlnorm(n = P, meanlog = 0, sdlog = 2) * # these are log-fold-changes
#'   rbinom(n = P, size = 1, prob = .15) *
#'   sample(c(-1, 1), size = P, replace = TRUE) # fold change on log scale
#' mu1 <- exp(beta) * mu0 # because we want log(mu2/mu1) = beta
#' relAbundances <- data.frame(
#'   g0 = mu0 / sum(mu0),
#'   g1 = mu1 / sum(mu1)
#' ) # relative abundance of absolute count
#' # relative abundance information (observed data in typical experiment)
#' Y0 <- rmultinom(n = 10, size = 1e4, prob = relAbundances$g0)
#' Y1 <- rmultinom(n = 10, size = 1e4, prob = relAbundances$g1)
#' Y <- cbind(Y0, Y1)
#' rownames(Y) <- paste0("celltype", 1:10)
#' colnames(Y) <- paste0("sample", 1:20)
#' group <- factor(rep(0:1, each = 10))
#' design <- model.matrix(~group)
#' v <- voomCLR(
#'   counts = Y,
#'   design = design,
#'   lib.size = NULL,
#'   plot = TRUE
#' )
#' fit <- lmFit(v, design)
#' plotBeta(fit)
#' @export
plotBeta <- function(fit) {
  for (cc in seq_len(ncol(fit$coefficients))) {
    curBeta <- fit$coefficients[, cc]
    plot(density(curBeta), main = colnames(fit$coefficients)[cc])
    abline(
      v = .getMode(beta = curBeta, n = nrow(fit$design)),
      col = "orange", lwd = 2
    )
  }
}
