#'  Forms the Projection Matrix
#'
#'  The projection function is used to generate the random projection matrix. It will form either dense or spare
#'  projection matrix. The Package supports 4 projection functions namely gaussian, probability, achlioptas and li.
#'  The number of rows and columns of the input sample is passed with the boolean value JLT. If JLT is set to TRUE,
#'  the dimension of the input data is reduced to the value returned by dimension() method.
#'  For Dense Matrix -  "gaussian" method
#'  For Sparse Matrix - "probability, achlioptas and li" method.

#'
#' @param rows - number of rows
#' @param cols - number of columns
#' @param JLT - Boolean to set JL transform (TRUE or FALSE)
#' @param eps - error tolerance level with default value 0.1
#' @param projection - projection function with default value "gaussian"
#'
#' @details
#' The 4 projection functions are
#'
#' 1."gaussian" - The default projection function is "gaussian". In probability theory, Gaussian distribution is also
#'  called as normal distribution. It is a continuous probability distribution used to represent real-valued random
#'  variables. The elements in the random matrix are drawn from N(0,1/k), N is a Natural number and
#'  k value calculated based on JL - Lemma using dimension() function.
#'
#' 2. "probability" - In this method, the matrix was generated using the equal probability
#' distribution with the elements [-1, 1].
#'
#' 3. "achlioptas" - Achlioptas matrix is easy to generate and also the 2/3rd of the matrix was filled
#' with zero which makes it as more sparse and cut-off the 2/3rd computation.
#'
#' 4. "li" - This method generalizes the achlioptas method and generate very sparse random matrix
#'  to improve the computational speed up of random projection.
#'
#' When comparing to gaussian function,  the other projection functions creates sparse matrix by filling with
#' zero's or one's to reduce the computation even more.
#'
#' @export
#'
#' @examples
#' # Load Library
#' library(RandPro)
#'
#' # Default Gaussian projection matrix without JL transform
#' mat <- form_matrix(600,1000,FALSE)
#'
#' # Default Gaussian projection matrix with JL transform of 50% Error tolerance
#' mat <- form_matrix(300,100000,TRUE,0.5)
#'
#' # Projection matrix with probability distribution of 50% Error tolerance
#' mat <- form_matrix(250,1000000,TRUE,0.5,"probability")
#'
#' # Projection matrix with li distribution of 50% Error tolerance
#' mat <- form_matrix(250,1000000,TRUE,0.5,"li")
#'
#' # Projection matrix with achlioptas distribution of 50% Error tolerance
#' mat <- form_matrix(250,1000000,TRUE,0.5,"achlioptas")
#'
#'
#' @keywords RandomProjection, Johnson-Lindenstrauss_Lemma, Dimension_Reduction, Gaussian-Distribution,
#' sparse_matrix, Achlioptas, Li,Probability distribution
#'
#' @return Projection Matrix
#'
#' @author Aghila G
#' @author Siddharth R
#'
#' @importFrom stats rnorm runif
#'
#' @references [1] N.I.R. Ailon and B.Chazelle, "The Fast Johnson Lindenstrauss Transform and Approximate Nearest Neighbors(2009)"
#' @references [2] Ping Li, Trevor J. Hastie, and Kenneth W. Church,  "Very sparse random projections(2006)".
#' @references [3] D. Achlioptas, "Database-friendly random projections(2002)"
#'
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

form_matrix <- function(rows,cols,JLT,eps=0.1,projection="gaussian")
{
  if(rows <= 0.0)
  {
    stop("The Number of rows should be greater than zero, got",dQuote(rows))
  }

  if(cols <= 0.0)
  {
    stop("The Number of columns should be greater than zero, got",dQuote(cols))
  }
  if(JLT == TRUE)
  {
    if(missing(eps))
    {
      message("Function uses deafault value 0.1 for epsilon")
    }
    cols=dimension(sample = cols,epsilon = eps)
  }
  if(missing(projection) || projection=="gaussian")
  {
  message("Function uses Gaussian Projection function")
  random_matrix <- matrix(rnorm(rows*cols,mean=0,sd=(1/sqrt(cols))),rows,cols)
  return(random_matrix)
  }
  if(projection=="achlioptas")
  {
    message("Function uses Achlioptas's Projection function")
    random_matrix <- floor(runif(cols*rows,1,7))
    sqr_3 <- sqrt(3);
    random_matrix[random_matrix==1] <- sqr_3;
    random_matrix[random_matrix==6] <- -sqr_3;
    random_matrix[random_matrix==2 | random_matrix==3 | random_matrix==4 | random_matrix==5] <- 0;
    random_matrix<- matrix(random_matrix,rows);
    return(random_matrix)
  }
  if(projection=="li")
  {
    message("Function uses Li's Projection function")
    s <- ceiling(sqrt(cols))
    pro <- c((1/(2*s)),(1-(1/s)),(1/(2*s)))
    random_matrix <- matrix(sample(c(-1,0,1), size=rows*cols, replace=TRUE, prob=pro), nrow=rows)
    return(random_matrix)
  }
  if(projection=="probability")
  {
    message("Function uses Probability Projection function")
    pro <- c(0.5,0.5)
    random_matrix <- matrix(sample(c(-1,1), size=rows*cols, replace=TRUE, prob =pro), nrow=rows)
    return(random_matrix)
  }
}

