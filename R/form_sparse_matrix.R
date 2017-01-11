#'  Form Sparse Random Matrix
#'
#' The purpose of sparse matrix is to reduce the computational complexity even more when comparing to gaussian matrix.
#' The sparse matrix can be generated using a simpler approach by filling with zero's or one's to reduce the computation.
#' This function supports 3 methods that are used to generate the sparse matrix.
#'
#' @details
#'
#' The 3 types of methods used in the function are
#'
#' 1. "probability" - In this method, the sparse matrix was generated using the equal probability
#' distribution with the elements [-1, 1].
#'
#' 2. "achlioptas" - Achlioptas matrix is easy to generate and also the 2/3rd of the matrix was filled
#' with zero which makes it as more sparse and cut-off the 2/3rd computation.
#'
#' 3. "li" - This method generalizes the achlioptas method and generate very sparse random matrix
#'  to improve the computational speed up of random projection.
#'
#' @param n_rows - number of rows in the sample
#' @param n_cols - number of columns in the sample
#' @param JLT - Boolean to set JL transform (TRUE or FALSE)
#' @param eps - error tolerance level
#' @param method - method to generate projection matrix
#'
#' @export
#'
#' @examples
#' #load library
#' library(RandPro)
#'
#' # Matrix generated with deafult method "achlioptas" without JLT
#' mat <- form_sparse_matrix(n_rows = 300,n_cols = 1000,JLT = FALSE)
#'
#' #create matrix with JLT and default method
#' mat <- form_sparse_matrix(n_rows = 250,n_cols = 1000000,JLT = TRUE,eps = 0.5)
#'
#' #create matrix with JLT and probability method
#' mat <- form_sparse_matrix(250,1000000,JLT = TRUE,eps = 0.5,method = "probability")
#'
#' #create matrix with JLT and li method
#' mat <- form_sparse_matrix(n_rows = 250,n_cols = 1000000,JLT = TRUE,eps = 0.5,method = "li")
#'
#' @keywords RandomProjection, Johnson-Lindenstrauss_Lemma, Dimension_Reduction, sparse_matrix, Achlioptas, Li, Probability distribution
#'
#' @return Random Sparse Matrix
#'
#' @author Aghila G
#' @author Siddharth R
#'
#' @importFrom  stats runif
#'
#' @references [1] Ping Li, Trevor J. Hastie, and Kenneth W. Church,  "Very sparse random projections(2006)".
#' @references [2] D. Achlioptas, "Database-friendly random projections(2002)"
#'
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

form_sparse_matrix <- function(n_rows,n_cols,JLT,eps=0.1,method="achlioptas")
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
  if(missing(method))
  {
    message("Projection Matrix was formed using default method Achlioptas")
  }
  if(method == "achlioptas")
  {
    sparse_matrix <- floor(runif(n_cols*n_rows,1,7))
    sqr_3 <- sqrt(3);
    sparse_matrix[sparse_matrix==1] <- sqr_3;
    sparse_matrix[sparse_matrix==6] <- -sqr_3;
    sparse_matrix[sparse_matrix==2 | sparse_matrix==3 | sparse_matrix==4 | sparse_matrix==5] <- 0;
    sparse_matrix<- matrix(sparse_matrix,n_rows);
    return(sparse_matrix)
  }
  if(method == "probability")
  {
    pro <- c(0.5,0.5)
    sparse_matrix <- matrix(sample(c(-1,1), size=n_rows*n_cols, replace=TRUE, prob =pro), nrow=n_rows)
    return(sparse_matrix)
  }
  if(method == "li")
  {
    s <- ceiling(sqrt(n_cols))
    pro <- c((1/(2*s)),(1-(1/s)),(1/(2*s)))
    sparse_matrix <- matrix(sample(c(-1,0,1), size=n_rows*n_cols, replace=TRUE, prob=pro), nrow=n_rows)
    return(sparse_matrix)
  }
}
