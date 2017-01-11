#'  Form Gaussian Matrix
#'
#'  In probability theory, Gaussian distribution is also called as normal distribution
#'  It is a continuous probability distribution used to represent real-valued random variables
#' The elements in the random matrix are drawn from N(0,1/k),
#' where k value calculated based on JL - Lemma.
#'
#'
#' @param n_rows - number of rows in the sample
#' @param n_cols - number of columns in the sample
#' @param JLT - Boolean to set JL transform (TRUE or FALSE)
#' @param eps - error tolerance level with default value 0.1
#'
#' @details
#' The function uses pnorm() function from stats package to generate random matrix with mean is zero and
#' standard deviation is 1/sqrt(k),  where k is the minimum number of dimension return from the find_dim_JL() function
#'
#' @export
#'
#' @examples
#' # Load Library
#' library(RandPro)
#'
#' # Without JLT
#' mat <- form_gauss_matrix(600,1000,FALSE)
#'
#' # With JLT of eps = 0.5
#' mat <- form_gauss_matrix(300,100000,TRUE,0.5)
#'
#' # With JLT of default eps value = 0.1
#' mat <- form_gauss_matrix(300,100000,TRUE)

#' @keywords RandomProjection, Johnson-Lindenstrauss_Lemma, Dimension_Reduction, Gaussian-Distribution
#'
#' @return Random Dense Matrix
#'
#' @author Aghila G
#' @author Siddharth R
#'
#' @importFrom stats rnorm
#'
#' @references [1] N.I.R. Ailon and B.Chazelle, "The Fast Johnson Lindenstrauss Transform and Approximate Nearest Neighbors(2009)"
#'
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

form_gauss_matrix <- function(n_rows,n_cols,JLT,eps=0.1)
{
  if(JLT == TRUE)
  {
    if(missing(eps))
    {
      message("Function uses deafault value 0.1 for epsilon")
    }
    n_comp = n_cols
    n_cols=find_dim_JL(sample =  n_cols,epsilon = eps)
  }

  if(n_rows <= 0.0)
  {
    stop("The Number of rows should be greater than zero, got",dQuote(n_rows))
  }

  if(n_cols <= 0.0)
  {
    stop("The Number of columns should be greater than zero, got",dQuote(n_cols))
  }
  gauss_matrix <- matrix(rnorm(n_rows*n_cols,mean=0,sd=(1/sqrt(n_cols))),n_rows,n_cols)
  return(gauss_matrix)
}
