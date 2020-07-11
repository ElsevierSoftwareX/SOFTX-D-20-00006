
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RandPro

<!-- badges: start -->

<!-- badges: end -->

The goal of RandPro is to provide the platform for R users to perform
random projection (RP) based feature extraction with ease and simple
functions. Random Projection is one of the projection based feature
extraction technique where the projection matrix is filled with
independent and identically distributed (i.i.d) random values. Even
though the projection is random, RP nearly preserves all the pairwise
distance between any two samples in the projected low dimensional
subspace with the controlled amount of error. RP even outperforms
principal component analysis (PCA) in very high dimensional data
analysis and streaming applications.

## Functions in RandPro

**dimension()** - determines the number of dimension required to project
the data in low dimensional space using Johnson-Lindenstrauss Lemma.

**form\_matrix()** - creates a projection matrix filled with four
suitable distributions that are primarily used in the research
community.

**classify()** - performs the classification with random projection
based feature extraction. It uses classifiers in the caret package.

## Installation

You can install the released version of RandPro from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RandPro")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cran/RandPro")
```

## Example - IRIS dataset

This is a basic example which shows you how to perform random projection
based feature extraction with classification

``` r
# Load Library
library(RandPro)
```

    ## Loading required package: caret

    ## Loading required package: lattice

    ## Loading required package: ggplot2

``` r
#Load Iris Data
data("iris")
#Split the data into training set and test set of 75:25 ratio.
set.seed(101)
sample <- sample.int(n = nrow(iris), size = floor(.75*nrow(iris)), replace = FALSE)
trainn <- iris[sample, ]
testt <- iris[-sample,]
#Extract the train label and test label
trainl <- trainn$Species
testl <- testt$Species
typeof(trainl)
#Remove the label from training set and test set
trainn <- trainn[,1:4]
testt <- testt[,1:4]
#classify the Iris data with default K-NN Classifier.
res <- classify(trainn,testt,trainl,testl)
```

    ## Function uses default value 0.5 for epsilon

    ## Function uses Gaussian Projection function

    ## Function uses default K-NN classifier

``` r
res
```

## Additional Information

The detailed explanation of the RandPro package and its functions are
available in this manual - [Click
here](https://CRAN.R-project.org/package=RandPro)
