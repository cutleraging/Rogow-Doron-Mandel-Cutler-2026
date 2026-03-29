## BIOINFORMATICS ==============================================================

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")

## core packages
BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "biomaRt",
                       "GSVA", "GSEABase", "GSVAdata", "Biobase", 
                       "clusterProfiler", "SummarizedExperiment", 
                       "ComplexHeatmap", "EnhancedVolcano", "PCAtools",
                       "org.Mm.eg.db", "org.Hs.eg.db"))

## GenomicFeatures: functionality for transcript database or TxDb objects
## AnnotationDbi: querying SQLite-based annotation data packages
## biomaRt: Interface to BioMart databases (i.e. Ensembl)

install.packages("msigdbr")

## msigdbr: Molecular Signatures Database (MSigDB) gene sets for GSEA


## BULK TRANSCRIPTOMICS =============================================================

## diann: process and analyse mass spectrometry-based single cell proteomics data

BiocManager::install(c("limma", "DESeq2", "genefilter", "zFPKM"))

## limma: Linear Models for Microarray Data


## CLEANING DATA ===============================================================

install.packages(c("janitor","outliers","missForest","frequency","Amelia",
                   "diffobj","mice","VIM","mi",
                   "wrangle"))

## janitor: Simple tools for data cleaning
## outliers: Detect outliers in data sets
## missForest: Impute missing values using the random forest algorithm
## frequency: Easy frequency tables from data sets
## Amelia: Multiple imputation of missing data
## diffobj: Compute and visualize differences of R objects
## mice: Multivariate imputation by chained equations
## VIM: Visualization and imputation of missing values
## mi: Tools for handling missing data in R
## wrangle: A collection of tools for data manipulation


## DATA TYPES AND FORMATS ======================================================

install.packages(c("stringr","lubridate","glue",
                   "scales","hablar","readr","readxl","haven", "arrow"))

## stringr: String manipulation package
## lubridate: Work with dates and times in R
## glue: Glue strings together in a flexible way
## scales: Graphical scales map data to aesthetics in plots
## hablar: Convert natural language to logical expressions in R
## readr: Read rectangular text data from file or string
## readxl: Read Excel files (.xls and .xlsx) into R
## haven: Import and export 'SPSS', 'Stata' and 'SAS' files


## WRANGLING, SUBSETTING, PLOTTING =============================================

install.packages(c("tidyverse", "data.table", "pheatmap", "viridis", 
                   "RColorBrewer"))

## tidyverse: A collection of packages for data manipulation and visualization
## - ggplot2: A system for creating graphics in R
## - dplyr: A grammar of data manipulation
## - tidyr: Tidy messy data
## - readr: A fast and friendly way to read rectangular data
## - purrr: Functional programming tools
## - tibble: A modern re-imagining of the data.frame
## - stringr: A consistent, simple and easy-to-use set of wrappers around the 
##   stringr functions of stringi package
## - forcats: Tools for working with categorical variables (factors)
##
## data.table: Fast and efficient data manipulation using data.table syntax


## STATISTICAL TESTS ===========================================================

install.packages(c("stats","ggpubr","lme4","MASS","car"))


## PARALLEL COMPUTING ==========================================================

install.packages(c("parallel", "foreach", "doParallel"))

## parallel: Base R package for parallel computing
## foreach: Provides a simple and flexible way to write parallel code
## doParallel: Provides a parallel backend for the foreach package


## SAMPLING ====================================================================

# install.packages(c("sampling","icarus","sampler","SamplingStrata",
#                    "survey","laeken","stratification","simPop"))


## MULTIVARIATE ANALYSIS =======================================================

# install.packages(c("CCA","CCP","MASS","gvlma","smacof",
#                    "MVN","rpca","EFA.MRFA","MFAg","MVar","fabMix",
#                    "fad","spBFA","mnlfa","GFA","lmds","SPCALDA",
#                    "semds", "vcd", "vcdExtra"))


## CLASSIFICATION AND CLUSTERING ===============================================

# install.packages(c("fpc","cluster","treeClust","e1071","NbClust","skmeans",
#                    "kml","protoclust","pvclust","genie", "tclust",
#                    "ClusterR","dbscan","CEC","GMCM","EMCluster","randomLCA",
#                    "MOCCA","factoextra",'poLCA', 'pdfCluster','flexclust',
#                    "teigen"))

## fpc: Flexible Procedures for Clustering
## cluster: Cluster Analysis Extended Rousseeuw et al.
## treeClust: Cluster analysis with trees
## e1071: Misc Functions of the Department of Statistics (e1071), TU Wien
## NbClust: Determining the Best Number of Clusters in a Data Set
## skmeans: Spherical k-Means Clustering
## kml: K-Means for Longitudinal Data
## compHclust: Complementary Hierarchical Clustering Analysis
## protoclust: Hierarchical Clustering with Prototypes
## pvclust: Hierarchical Clustering with P-Values via Multiscale Bootstrap Resampling
## genie: genie: A New, Fast, and Outlier Resistant Hierarchical Clustering Algorithm
## tclust: Robust Trimmed Clustering
## ClusterR: Gaussian Mixture Models, K-Means, Mini-Batch-Kmeans, K-Medoids and Affinity Propagation Clustering
## dbscan: Density-Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms
## CEC: Cross-Entropy Clustering
## GMCM: Generalized Maximum Contrast Method
## EMCluster: EM Algorithm for Model-Based Clustering of Finite Mixture Gaussian Distribution
## randomLCA: Random Effects Latent Class Analysis
## MOCCA: Multi-Objective Optimization for Collecting Cluster Analysis Results
## factoextra: Extract and Visualize the Results of Multivariate Data Analyses
## poLCA: Polytomous variable Latent Class Analysis
## pdfCluster: Cluster analysis via nonparametric density estimation
## flexclust: Flexible Cluster Algorithms
## EMMIXskew: The EM Algorithm and Skew-Elliptical Distributions
## teigen: Model-Based Clustering and Dimension Reduction for the Wrapped Normal Distribution


## TIME SERIES =================================================================

# install.packages(c("zoo","xts","timeSeries","tsModel",
#                    "TSA","fma","fpp2","fpp3","TSdist","TSclust","feasts",
#                    "MTS","sazedR","kza","fable","forecast","tseries",
#                    "nnfor","quantmod",'meboot','rugarch','betategarch','GAS'))

## ts: Time Series analysis and computation
## zoo: S3 Infrastructure for Regular and Irregular Time Series
## xts: Extensible time series
## timeSeries: Rmetrics - Financial Time Series Objects
## tsModel: Time Series Modeling for Air Pollution and Health
## TSMining: Mining Univariate and Multivariate Motifs in Time-Series Data
## TSA: Time Series Analysis
## fma: Data Sets from Forecasting: Methods and Applications
## fpp2: Data for "Forecasting: Principles and Practice" (2nd Edition)
## fpp3: Data for "Forecasting: Principles and Practice" (3rd Edition)
## tsfa: Time Series Factor Analysis
## TSdist: Distance Measures for Time Series Data
## TSclust: Time Series Clustering Utilities
## feasts: Feature Extraction And Statistics for Time Series
## MTS: All-Purpose Toolkit for Analyzing Multivariate Time Series (MTS) and Estimating Multivariate Volatility Models
## dse: Dynamic Systems Estimation (Time Series Package)
## sazedR: Statistical Analysis of Ziggurat Method (R)
## kza: Kolmogorov-Zurbenko Adaptive Filters
## fable: Forecasting Models for Tidy Time Series
## forecast: Forecasting functions for time series and linear models
## tseries: Time Series Analysis and Computational Finance
## nnfor: Time Series Forecasting with Neural Networks
## quantmod: Quantitative Financial Modelling Framework
## meboot: Maximum Entropy Bootstrap for Time Series
## rugarch: Univariate GARCH models
## betategarch: Beta-t-EGARCH models
## GAS: Generalized Autoregressive Score Models


## FUNCTIONAL DATA ==========================================================

# install.packages(c("fda","FDboost","fds","ftsa","fdasrvf","refund","fdapace",
#                    "fdatest","fdaMixed","goffda","mlr","fda.usc"))

## fda It is a basic reference to work in R with functional data
## FDboost: Boosting Functional Regression Models
## fds: Functional Data Sets
## ftsa: Functional Time Series Analysis
## fdasrvf: Elastic Functional Data Analysis
## refund: Regression with Functional Data
## fdapace: Functional Data Analysis and Empirical Dynamics
## StatFda: exploratory analysis and functional regression models
## tidyfun: makes data wrangling and exploratory analysis of functional data easier
## fdatest: Interval Testing Procedure for Functional Data
## fdakma: Functional Data Analysis: K-Mean Alignment
## fdaMixed: Functional data analysis in a mixed model framework
## goffda: Goodness-of-Fit Tests for Functional Data


## BUILDING AND VALIDATING ML ==================================================

# install.packages(c("tree", "e1071","crossval","caret","rpart","bcv",
#                    "klaR","cvAUC","CVThresh",
#                    "cvTools","cvms","blockCV"))

## tree: Decision tree models
## e1071: Misc Functions of the Department of Statistics, Probability Theory Group (Formerly: E1071), TU Wien
## crossval: Generic functions for cross validation
## caret: Classification and Regression Training package
## rpart: Recursive Partitioning and Regression Trees
## bcv: Cross-Validation for the SVD (Bi-Cross-Validation)
## klaR: Classification and visualization
## EnsembleCV: Extensible package for Cross-validation-Based Integration of Base Learners
## gencve: General Cross Validation Engine
## cvAUC: Cross-validated area under the curve (AUC) for the ROC curve
## CVThresh: Estimating the optimal threshold for a prediction model with cross-validation
## cvTools: Tools for cross-validation in R
## dcv: Cross-Validation for Discriminant Analysis
## cvms: Cross-Validation for Model Selection
## blockCV: Block-wise Cross-Validation for Covariate-Adjusted Models



## RANDOM FORESTS ==============================================================

# install.packages(c("randomForest","grf","ipred","party","randomForestSRC",
#                    "grf","BART","Boruta","LTRCtrees","REEMtree",
#                    "binomialRF","superml"))

## randomForest: Classification and regression based on a forest of trees using random inputs
## grf: Generalized random forest for classification, regression and survival analysis
## ipred: Improved predictive models using random forests
## party: A laboratory for recursive partytioning
## randomForestSRC: Random Forests for survival, regression and classification (RF-SRC)
## grf: Generalized random forest for classification, regression and survival analysis
## BART: Bayesian Additive Regression Trees
## Boruta: Wrapper algorithm for all-relevant feature selection
## LTRCtrees: Survival trees to model left-truncated and right-censored data
## REEMtree: Regression trees with random effects for longitudinal (panel) data
## refr: An implementation of the Randomized Ensembles Framework
## binomialRF: Binomial random forest for imbalanced data
## superml: Build machine learning models as easily as a formula


## MODEL INTERPRETATION =====================================================

# install.packages(c("lime","localModel","iml","EIX","flashlight",
#                    "interpret","outliertree","breakDown"))

## Lime: Machine learning explanation and model interpretation package
## localModel: Fit local regression models to non-linear data
## iml: Interpretable Machine Learning package
## EIX: Explainability of Importance in Random Forest
## flashlight: Model Agnostic Feature Selection
## interpret: Fit interpretable models and explain blackbox models
## outliertree: Identify and visualize outliers in regression trees
## breakDown: Model-agnostic methods for decomposition of prediction results


## DOCUMENTS CREATION ==========================================================

install.packages(c("devtools","usethis","roxygen2","knitr",
                   "rmarkdown","flexdashboard","shiny",
                   "xtable","httr","profvis","officedown"))

## devtools: Tools to make developing R packages easier
## usethis: Automate package and project setup tasks
## roxygen2: In-line documentation for R
## knitr: A general-purpose package for dynamic report generation in R
## rmarkdown: Dynamic documents for R
## flexdashboard: Easy interactive dashboards for R
## Shiny: Web application framework for R
## xtable: Export tables to LaTeX or HTML
## httr: Tools for working with URLs and HTTP
## profvis: Interactive visualizations for profiling R code
## officedown: Convert R Markdown to Office documents (Word, PowerPoint and Excel)


## BAYESIAN INFERENCE =======================================================

# install.packages(c("rstan","brms","rstanarm","tidybayes","bayestestR",
#                    "bayesrules","bayesplot","loo"))


## rstan: Provides R interface to the Stan probabilistic programming language for Bayesian inference.
## brms: Implements Bayesian regression models using Stan for fitting, estimation, and posterior prediction.
## rstanarm: Implements Bayesian regression models using Stan for fitting, estimation, and posterior prediction, but with a simplified R syntax.
## tidybayes: Provides tidy data structures and visualization tools for Bayesian models using ggplot2.
## bayestestR: Provides various tools for Bayesian model checking, hypothesis testing, and model comparison.
## bayesrules: Implements Bayesian networks for modeling probabilistic relationships between variables.
## bayesplot: Provides plotting functions for visualizing Bayesian inference results.
## loo: Implements the leave-one-out cross-validation method for comparing models and assessing model fit.


## TIDY MODELLING ===========================================================

# install.packages("tidymodels")

## rsamples: split the data into training and testing sets (as well as cross validation sets)
## recipes: prepare the data with preprocessing (assign variables and preprocessing steps)
## parsnip: specify and fit the data to a model 
## yardstick and tune: evaluate model performance by metrics
## workflows: combining recipe and parsnip objects into a workflow (this makes it easier to keep track of what you have done and it makes it easier to modify specific steps)
## tune and dials: model optimization (hyperparameters)
## broom: make the output from fitting a model easier to read
## baguette: speeding up bagging pipelines
## butcher: dealing with pipelines that create model objects that take up too much memory
## discrim: more model options for classification
## embed: extra preprocessing options for categorical predictors
## corrr: more options for looking at correlation matrices
## rules: more model options for prediction rule ensembles
## tidypredict: running predictions inside SQL databases
## modeldb: working within SQL databases and it allows for dplyr and tidyeval use within a database
## tidyposterior: compares models using resampling statistics


## GAMs ========================================================================

# install.packages(c("gam","mgcv","gamlss","mboost","gss","scam","gamm4"))

## gam: Generalized Additive Models for regression with the original backfitting approach *
## mgcv: Mixed GAM computation vehicle with GCV/AIC/REML smoothness estimation 
## vgam: Vector Generalized Linear and Additive Models *
##     - pospois and posnegbinom functions in the VGAM package fit truncated poisson and negbin models
## gamlss: Generalized Additive Models for location, scale and shape **
## mboost: Model-based boosting
## gss: General smoothing splines
## scam: Shape constrained additive models
## gamm4: Generalized Additive Mixed Models using 'mgcv' and 'lme4'
## bayesx: Bayesian inference in structured additive regression models ***
## Methods discussed here are in R recommended package mgcv

## *No smoothing parameter selection
## **Limited smoothing parameter selection
## ***see also mgcv::jagam


## NETWORK MODELS ==============================================================

# install.packages(c("igraph", "igraphdata", "blockmodels",
#                    "latentnet", "huge", "covglasso"))


## igraph: Network analysis and visualization package
## igraphdata: Example datasets for the 'igraph' package
## networkdata: Example datasets for network analysis
## blockmodels: Fit and analyze block models for networks
## latentnet: Latent position and cluster models for network analysis
## huge: High-dimensional undirected graph estimation
## covglasso: Gaussian graphical models with the lasso penalty


## NLP ========================================================================

# install.packages(c("tm", "quanteda", "tidytext", "wordcloud", "topicmodels",
#                    "SentimentAnalysis", "text2vec", "syuzhet", "openNLP",
#                    "openNLPdata", "NLP", "RWeka", "tau", "koRpus", "stringr",
#                    "tokenizers", "textmineR", "rvest", "hunspell", "lexicon"))


## tm: Text mining tools for managing, processing, and analyzing text data
## quanteda: Quantitative text analysis methods and utilities
## tidytext: Text mining using tidy data principles
## wordcloud: Creating word cloud visualizations
## topicmodels: Fitting and analyzing topic models for text data
## SentimentAnalysis: Analyzing sentiment in text data
## text2vec: Efficient text vectorization and topic modeling
## syuzhet: Sentiment and emotion analysis using narrative arcs
## openNLP: Natural language processing with the Apache OpenNLP tools
## openNLPdata: Model files required for the 'openNLP' library
## NLP: Basic classes and methods for natural language processing
## RWeka: R interface to the Weka machine learning toolkit, including text classification algorithms
## tau: Text analysis utilities
## koRpus: Analyzing and processing text data, including readability measures
## stringr: Simple and consistent manipulation of strings
## tokenizers: Various text tokenization methods
## textmineR: Text mining and topic modeling
## rvest: Web scraping, including text extraction from HTML
## hunspell: Spell checking and stemming based on the Hunspell library
## lexicon: Collection of lexicons and dictionaries for text analysis

## Record the state in renv.lock
renv::snapshot()

sessionInfo()
