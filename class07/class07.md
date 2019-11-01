Class 7 R functions and packages
================
Fabian Lim
10/23/2019

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](class07_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

\#\#Revisit our functions from last day

``` r
source("http://tinyurl.com/rescale-R")
```

## Let’s try our rescale function from last day

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
rescale(c(3, 10, NA, 7))
```

    ## [1] 0.0000000 1.0000000        NA 0.5714286

``` r
# Lets define an example x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
sum(is.na(x))
```

    ## [1] 2

``` r
# Our working snippet
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
both_na <- function(x,y) {
sum(is.na(x) & is.na(y))
}
```

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
```

``` r
both_na(x, y1)
```

    ## [1] 2

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na3 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be vectors of the same length")
 }

 na.in.both <- ( is.na(x) & is.na(y) )
 na.number <- sum(na.in.both)
 na.which <- which(na.in.both)
 message("Found ", na.number, " NA's at position(s):",
 paste(na.which, collapse=", ") )

 return( list(number=na.number, which=na.which) )
}
```

``` r
grade <- function(student1, student2) {
# student 1
student1<- c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
# assign value 0 to NA
student2[is.na(student2)] <- 0
# sort scores in decreasing order
sort(student1, decreasing=TRUE)
sort(student2, decreasing=TRUE)
# derive mean for top 7 scores
mean(student1[1:7])
mean(student2[1:7])
}
```

``` r
grade(student1, student2)
```

    ## [1] 79.57143
