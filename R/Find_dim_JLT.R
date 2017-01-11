#'  Johnson - Lindenstrauss Function
#'
#' Johnson Lindenstrauss Transform[JLT] is the heart of random projection.
#' The lemma states that a small set of points in a high-dimensional space can be embedded into
#' a space of much lower dimension in such a way that distances between the points are nearly preserved.
#' The lemma has used in dimensionality reduction, compressed sensing, manifold learning and graph embedding.
#' \eqn{(1 - epsilon) ||x - y||^2 < ||RP(x) - RP(y)||^2 < (1 + epsilon) ||x - y||^2}
#' where x and y are number of rows and columns respectively
#'
#' @details
#' The function find_dim_JL() is used to find the minimum dimension required to project the data from high dimensional
#' space to low dimensional space. The number of sample and error tolerant level was passed as an input argument
#' to the function find_dim_JL() .
#' It will return the minimal size of the random subspace to guarantee a bounded distortion introduced by
#' the random projection.
#'
#' @param sample - number of samples
#' @param epsilon - error tolerance level with default value 0.1
#'
#' @export
#'
#' @examples
#' #load library
#' library(RandPro)
#'
#' #Calculate minimum dimension using eps =0.5 for 1000000 sample
#' y <- find_dim_JL(1000000,0.5)
#'
#' #Calculating minimum dimension using different epsilon value for 1000000 sample
#' d <-  c(0.5,0.1)
#' x<- find_dim_JL(103260,d)
#'
#' @keywords Random projection , Johnson - Lindenstrauss Lemma, Dimension Reduction
#'
#' @return minimum number of dimension required to maintain the pai wise distance with the controlled amount of error
#'
#' @author Aghila G
#' @author Siddharth R
#'
#' @references [1] William B.Johnson, Joram Lindenstrauss, "Extension of Lipschitz mappings into a Hilbert space (1984)"
#' @references [2] Sanjoy Dasgupta , Anupam Gupta "An elementary proof of a theorem of Johnson and Lindenstrauss (2003)"
#'
#'
#'
#' @seealso \href{http://cseweb.ucsd.edu/~dasgupta/papers/jl.pdf}{Johnson-Lindenstrauss Elementary Proof}
#'
#'

#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

find_dim_JL <- function(sample,epsilon=0.1)
{
  if(nargs()==1)
  {
    message("Function uses deafault value 0.1 for epsilon")
  }
  sam <- sample
  eps <- epsilon
  for(i in seq(eps))
  {
    try(
    if((eps[i]>1.0)|(eps[i]<0.0))
    {
      stop("The JL bound for epsilon is [0.0,1.0], got",dQuote(eps))
    }
    else if(sam[i]<=0.0)
    {
      stop("The JL bound for number of samples is greater than zero, got",dQuote(sam))
    }
    else
    {
      denominator = ( (eps ^ 2) / 2 - (eps ^ 3) / 3 )
      return ( floor((4 * log(sample) / denominator)))
    }
    )
  }
}
