#'  Classification Function
#'
#'  The classify() function allows the user to combine the task of random projection based dimension reduction
#'  and classification within a single function. The dimension of the training data and test data was reduced
#'  by the  value returned from the dimension() method. Then the projection matrix was generated using
#'  form_matrix() function based on the input paramater "projection".Then the training data and test data was
#'  projected into the low dimensional space by multiplying with the projection matrix. At last the reduced matrix
#'  was given to the classifier. The confusion matrix is the output of the classifier where we can calculate
#'  the performance of the classifier.
#'
#'
#' @param train_data - Training data of either matrix or data frame
#' @param test_data - Test data of either matrix or data frame
#' @param train_label - Training label of either vector or data frame
#' @param test_label - Test label of either vector or data frame
#' @param eps - Epsilon with default 0.1
#' @param projection - projection function with default "gaussian"
#' @param classifier - classifier with default "knn"
#'
#' @details
#'
#'  The parameters train_data,test_data,train_label and test_label are mandatory arguments. The eps is the error
#'  tolerance paramater. The value of eps must be \eqn{0.0<eps<1.0}. The default value of eps is 0.1 that means 10 percentage of
#'  error is acceptable during projection. The supported projection functions are gaussian, probability, li, and
#'  achlioptas.The default projection method is "gaussian". The complete detail of the projection function is given in
#'  form_matrix() function. The final argument "classifier" in the function defines the classifier to train the model.
#'  The supported classifier for classification task are
#'
#'  "knn" - k-nearest neighbor classification
#'
#'  "svmlinear" - Support Vector Machine
#'
#'  "nb" - Naive Bayes Classifier
#'
#' @importFrom stats predict
#' @import caret
#' @import e1071
#'
#' @export
#'
#' @examples
#' # Load Library
#' library(RandPro)
#'
#' #Load Iris Data
#' data("iris")
#'
#' #Split the data into training set and test set of 75:25 ratio.
#' set.seed(101)
#' sample <- sample.int(n = nrow(iris), size = floor(.75*nrow(iris)), replace = FALSE)
#' trainn <- iris[sample, ]
#' testt  <- iris[-sample,]
#'
#' #Extract the train label and test label
#' trainl <- trainn$Species
#' testl <- testt$Species
#' typeof(trainl)
#'
#' #Remove the label from training set and test set
#' trainn <- trainn[,1:4]
#' testt <- testt[,1:4]
#'
#' #classify the Iris data with default K-NN Classifier.
#' res <- classify(trainn,testt,trainl,testl)
#' res
#'
#' @keywords classification feature_extraction k-nn svm confusion_matrix 
#'
#' @return Confusion Matrix
#'
#' @author Aghila G
#' @author Siddharth R
#'
#'
#' @references [1] Cannings, T. I. and Samworth, R. J. "Random projection ensemble classification(2015)".
#' @references [2] Ella Bingham and Heikki Mannila, "Random projection in dimensionality reduction: Applications to image and text data(2001)".
#'
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

classify <- function(train_data,test_data,train_label,test_label,eps=0.1,projection="gaussian",classifier="knn")
{
  train <-train_data
  test <-test_data
  j<- ncol(train)
  class <- train_label
  if(is.data.frame(train))
  {
    train <- as.matrix(train)
    colnames(train) <- NULL
    rownames(train) <- NULL
  }
  if(is.data.frame(test))
  {
    test <- as.matrix(test)
    colnames(test) <- NULL
    rownames(test) <- NULL
  }
  if((nrow(train)||nrow(test)||ncol(train)||ncol(test))==0)
  {
    stop("The Number of rows/columns should be greater than zero in training set and test set")
  }
  if(ncol(train)!=ncol(test))
  {
    stop("The Number of columns in training set and test set must be same")
  }
  if(missing(eps))
  {
      message("Function uses default value 0.5 for epsilon")
  }
    k=dimension(j,eps)

    rmat <- form_matrix(j,k,FALSE,projection = projection)
  ftrain <- train %*% rmat
  ftest <- test %*% rmat

  ftrain <- as.data.frame(cbind(ftrain,class))
  ftest <- as.data.frame(ftest)
  colnames(ftrain)[ncol(ftrain)]<-"class"

  ftrain[["class"]] <- factor(class)

  if(missing(classifier))
  {
    message("Function uses default K-NN classifier")
  }
    fit <- train(class~.,ftrain,method=classifier)
    pred <- predict(fit,ftest)
    print(confusionMatrix(pred,test_label))
}






