# Network-inference-LASSO-stabsel-extensions
**R functions to infer gene co-expression networks using LASSO with stability selection, with experimental extensions.**

This is an R implementation of LASSO with stability selection, originally introduced by Meinshausen & Bühlmann (2010). It is solving a feature selection regression task, i.e.: what is the set of predictor variables that best predict a response vector? In addition, I implement an extension of stability selection for arbitrary sub-divisions during sub-sampling (Beinrucker et al., 2016), as well as an experimental combination of stability selection with the distance precision matrix method by Ghanbari et al. (2016).

I make use of the R packages *glmnet* (LASSO), *corpcor*, and *foreach*/*doParallel* for parallalization. 

### Example application
A common task in computational biology is inferring which genes are acting together in response to certain stimuli (construct a co-expression network). The input data are measurements of gene expression in different conditions. The naive approach would be to calculate pairwise correlations between all genes, and draw network edges between genes that are highly correlated or anti-correlated. However, this will typically result in a very dense network with many redundant edges, and other problems (e.g. potentially losing local structure).

Many more sophisticated methods exist (partial correlation, multiple regression methods, random forests, ...), one of them being regularized multiple regression using LASSO. This can be combined with the stability selection approach by Meinshausen & Bühlmann to:
* obtain significance scores
* avoid having to explicitly select a regularization parameter

Briefly: for each gene, we perform multiple regression with L1 regularization (LASSO) to infer which subset of other genes (features) can best explain the observed behavior of that gene. We repeat this for a number of iterations, using only a randomly chosen 1/2 of observations on each iteration. We note down the top predictors (genes) on each iteration, by order of their appearance in the LASSO regularization path. If a predictor is consistently among the top predictors, regardless of sub-sample of observations, we regard it as reliable ("stable") and draw an edge to our predicted gene in the network.

### Extension: arbitrary sub-divisions
The functions implemented here support arbitrary sub-divisions of observations (instead of the original 1/2 split), as introduced by Beinrucker et al. (2016). This is useful to reduce complexity in cases with very many observations and may have other interesting properties.
* function for single response: *rs.stabsel*
* parallelized function for pairwise associations in a matrix of predictors: *rs.stabsel.matrix*

### Extension: distance partial correlations
One down-side of LASSO is that it is an inherently linear method, whereas gene co-expression patterns may be non-linear. 
Introduced by Ghanbari et al. (2016), the distance partial correlation method uses the distance covariance metric by Szekely et al. (2013) to solve this problem: observations are first mapped into a high-dimensional space where they behave linearly. Then, the fast approximation of the precision matrix (Schafer et al., 2005) is used to calculate partial correlations on the extended representation.
This method is implemented here in the function *rs.distancePartialCorrelation*.

In addition, I have implemented an experimental method (*rs.stabsel.dcor*), using the stability selection procedure that employs the full partial correlations (as in Schafer and Strimmer) instead of LASSO, with or without mapping into the extended space (N2). Note that the number of observations in the extended space is the square of the original number, so this can get expensive.
(Note: this was not terribly beneficial in my limited set of tests, but it's included for curiosity)

### References
Meinshausen, N. and Bühlmann, P., 2010. Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), pp.417-473.

Beinrucker, A., Dogan, Ü. and Blanchard, G., 2016. Extensions of stability selection using subsamples of observations and covariates. Statistics and Computing, 26(5), pp.1059-1077.

Schaefer, J., Opgen-Rhein, R., Zuber, V., Ahdesmäki, M., Duarte Silva, A.P. and Strimmer, K., 2013. corpcor: Efficient estimation of covariance and (partial) correlation. R package version, 1(6).

Ghanbari, M., Lasserre, J. and Vingron, M., 2016. The Distance Precision Matrix: computing networks from nonlinear relationships. arXiv preprint arXiv:1605.03378.

Friedman, J., Hastie, T. and Tibshirani, R., 2009. glmnet: Lasso and elastic-net regularized generalized linear models. R package version, 1(4).
