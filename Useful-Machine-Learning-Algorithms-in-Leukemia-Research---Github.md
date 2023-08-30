Useful Machine Learning Algorithms in Leukemia Research
================
Erica Laidler
2022-05-20

**Introduction**

One of the most powerful methods for scientists to better understand,
prevent, and predict the development of cancer is to assess the
relationships between different types of cancers and their associated
risk factors. These may include old age, obesity, alcohol use, smoking,
personal or family history of cancer, past viral infections, certain
chemicals, and exposure to radiation, depending on the type of cancer
(ASCO). It can be helpful to determine the diseases which individuals
are most likely to develop so that they can make informed decisions and
take precautionary measures. This paper will focus on developing models
which can help predict whether a patient is likely to develop certain
types of leukemia, based on a large sampling of their genes. I will also
present and assess a variety of methods which could be useful for
isolating a smaller subset of highly relevant genes in the context of
leukemia or other diseases.

When trying to determine the risk factors for a particular form of
cancer, it is useful for researchers to analyze data from previous
cohorts of patients. This analysis is quite challenging when we are
focusing on factors on the biochemical level. In particular, it is
difficult to summarize and derive general patterns about genes because
there are thousands of them, they interact with one another, and they do
not always function the same way in everyone. For instance, the same
gene might perform slightly different function or have different effects
on the body for different people. Thus, a high rate of expression for a
particular gene might not produce the same effect in two different
people. This being said, it is possible to uncover general trends. If we
have data to show that a huge percentage of people who shared certain
similarities in their genetic makeup ended up devleping acute myeloid
leukemia, this may suggest that future patients with a similar
combination of genes may be at high risk. Though there is almost never
certainty in these cases, it is often possible to construct models which
can help guide our predictions about patients’ level of risk.

When determining patterns in cancer risk, there are different ways of
assessing genetic makeup. Observing genetic mutations is one common way.
In fact, leukemia itself occurs as a result of mutations in the DNA of
bone marrow cells (Medical News Today). A mutation is a “change that
occurs in our DNA sequence, either due to mistakes when the DNA is
copied or as the result of environmental factors,” either before birth
or during one’s lifetime (Your Genome). Many of the cancer-causing risk
factors we mentioned earlier, like cigarette smoke or excessive exposure
to UV light, are carcinogenic because they can lead to mutations in DNA
which put individuals at a higher risk of developing cancer. Some
established examples of cancer-causing mutations are changes in certain
genes such as FLT3, c-KIT, and RAS, which are common in acute-myeloid
leukemia cells (ARUP Consult). It is possible to test for these
mutations by taking a sample of saliva, blood, or tissue (Family
Doctor.org).

However, mutations are not the only way to describe genetic makeup. Some
researchers have taken a more holistic approach when assessing the types
of genetic sequences which are associated with a higher risk of cancer,
and developing models to explain them. While an awareness of key
mutations can be a simple and extremely useful and effective method for
assessing cancer risk, these mutations only describe the deviations of a
few genes. Other methods, which look for patterns in gene expression,
protein-to-protein interactions, and regulatory-sequence information of
brain genes, can help predict using many more genes. (Asif). For
instance, a model which makes a classification prediction for a
patient’s likelihood of developing a particular disease, based on the
data from hundreds of past patients and thousands of genes, is likely to
be more informative than a prediction based on testing for the mutations
of just a few genes.

This paper focuses on assessing the use of different machine learning
methods in the analysis of the genetic component behind the development
of five different types of leukemia. I will look at the gene expression
for over 22,000 genes to develop a model which can help predict a
patient’s risk of developing each of the five subtypes. I will also
explain how these tools may be useful in other cancer-related analysis,
including the determination of important genes in the early stages of
researching a disease. Leukemia is an extremely deadly form of cancer,
with an approximate 23,660 deaths expencted to be attributed to it in
2021 (LLS). The more that can be understood about the genetic components
behind it, the better we can assess risk as well as develop better
treatments and prevention strategies.

**Data**

The original data set comes from the an “extensively curated microarray
database” from the Structural Bioinformatics and Computational Biology
Lab. This lab collects gene expression data on a variety of cancers, and
publishes it for public analysis and use. The data set I chose contains
the gene expression data for 22,283 genes, collected from 64 different
patients diagnosed with one of the five subtypes of leukemia. Thus it
contains 64 rows to represent the 64 patients, 22,283 columns to
represent the genes, and 1 column called ‘type’ which contains the
leukemia subypte diagnosis for a particular patient.

The first step is to download and import the data. This presents some
challenges due to the extreme size of the data set. First I download the
data frame via the application Numbers. However, I find that I have to
download it as the tranpose of the original, because the application
only offers a limited number of (1000) rows. This is not sufficient for
the 22, 000+ rows in my data set. Downloading the transpose, however, is
a functional solution. I then import the transposed data set into R, as
shown below.

``` r
leukemia_og <- read.csv("leukemia.csv")
```

Next, it is important to tidy the data sets, and put them in a structure
that will prove useful later. I build three slightly different data
sets.

Data Frame 1 is the original data set, called ‘leukemia’. It contains 64
rows and 22,284 columns, including the 22,283 gene columns and 1 type
column.

Data Frame 2 is identical to the first, except that it does not contain
the ‘type’ column and all of the entries of the data frame are converted
to numeric. This will be useful during the clustering procedure.

Data Frame 3 is the transpose of Data Frame 2, with 22,283 rows and 64
columns. It will be useful in exploring methods for the isolation of
important genes related to leukemia.

Below I create Data Frame 1, and output the first six lines.

``` r
#Create Data Frame 1: leukemia:

#Previously, in order to import the data set into R, I used the transpose of the original data. 
#Now, to return the data to its original form, I again take the transpose, 
#and name the resultingn data frame 'leukemia_t':

library(sjmisc)
leukemia_t <- leukemia_og %>% rotate_df()

#Now that the data is transposed back to its original form, there is still a bit more
#data tidying to do. Here, the first row does not actually contain data - instead, it #simply contains the column names. So I edit the data set to make sure its columns #names are labeled correctly. 

#1) Name the columxns based on the entries of the first row.
colnames(leukemia_t) <- leukemia_t[1,]

#2) Delete the first row, which is now unnecessary.
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
    ## ✓ tibble  3.1.5     ✓ dplyr   1.0.7
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.0.2     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x tibble::add_case()  masks sjmisc::add_case()
    ## x dplyr::filter()     masks stats::filter()
    ## x purrr::is_empty()   masks sjmisc::is_empty()
    ## x dplyr::lag()        masks stats::lag()
    ## x tidyr::replace_na() masks sjmisc::replace_na()

``` r
leukemia <- leukemia_t %>%
  filter(!row_number() %in% 1)

#Now, the data frame is sufficiently tidied. It contains 64 rows (observations) and 
#22,284 columns (genes)
```

Below, I create Data Frames 2 and 3.

``` r
#Create Data Frame 2 (leukemia_wo2):
#Contains alterations which will be useful for certain methods, like unlabeled #clustering
#Eliminate the first column of the data frame because it contains the categorical 
#data/response variable
leukemia_wo <- subset(leukemia, select = -c(type) )

#Change all of the values in the columns from character to numeric (doubles):
leukemia_wo2 <- mutate_all(leukemia_wo, function(x) as.numeric(as.character(x)))

#Finished Data Frame 2: 64 rows (observations) and 22,283 columns (genes)
```

``` r
#Data Frame 3: leukemia_tr_wo2:

#Now I will use the initial 'leukemia_og' data frame I imported (the one that was the
#transposed version of the original due to the limitations of the Numbers app) to 
#produce a data frame with 22,283 rows (genes) and 64 columns (patients). This is in 
#contrast to data frames 1 and 2, which had 64 rows and 22,284 columns.
#Begin with the initially imported leukemia_og, in which 22,284 rows = genes, and 64 
#columns = patients. Rename it as leukemia_tr, as a reminder that the initial import 

#data frame is the transposed version of the data frame found on the website.
leukemia_tr <- leukemia_og

#Remove the first row, as this contains the 'type' variable.
library(dplyr)
leukemia_tr_wo <- leukemia_tr %>%
  filter(!row_number() %in% 1)

#Make it so that the names of the rows are equal to the values of the first column 

#(gene name), and then delete the first column.
rownames(leukemia_tr_wo) <- leukemia_tr_wo$samples
leukemia_tr_wo <- subset(leukemia_tr_wo, select = c(-samples))

#Convert the entries to numeric.
leukemia_tr_wo2 <- mutate_all(leukemia_tr_wo, function(x) as.numeric(as.character(x)))

#Finished data frame: 22,283 rows; 64 columns
```

Now that I have tidied the data and developed relevant data frames, I
can continually refer to them during the analysis portion of the
project.

**Methods**

The research goal is to assess different techniques in the analysis of
gene expression data for multiple applications, such as (1) Acting as a
predictive model for the identification/classification of patients,
based on their gene expression, into low-risk or high-risk groups for
the development of different subtypes of cancers, and (2) Isolating a
subset of genes which are highly influential in the context of a
particular form of cancer.

Past research has used a variety of strategies, like support vector
machines, clustering, AI and neural networks, random forest, and random
walk with restart (RWR). In this project, I use clustering, because it
is a powerful method for effectively reducing the dimensions of the
data, and logistic regression.

There are two main stages of the plan. First, I cluster the data, and
then I use the clustering information to inform a logistic regression
model.

Though the original data set contains the patient leukemia diagnosis,
the clustering method I use is unsupervised, meaning that it does not
depend on labels. The labels in this case refer to the patients’
leukemia diagnosis. Thus, in the clustering methods, I use a version of
the data set which does not contain information about the patients’
leuekemia subtype. Data Frame 2 is useful here because it does not
contain the column ‘type’.

In theory, since the clustering algorithms do not know the patient
subtype, they will be organizing patients based purely on similarities
in their gene expression. I hypothesize that patients sorted into the
same cluster will be more likely to have the same leukemia diagnosis. If
this is true, it suggests that the clustering method is able to
correctly discern patterns in gene expression which are nearly
impossible to synthesize without algorithmic methods. It also suggests
that these methods may be used in future to sort new patients, who have
not yet been diagnosed, into risk groups. New patients who are sorted
into a cluster that is associated with a particular subtype may be more
likely to develop this subtype than any of the others.

There are many different types of clustering. It is not always easy to
“properly tune” different parameters of the clustering procedure (Biomed
Central). This includes determining the clustering method as well as the
distance which will be used to measure between data points. In this
project, I use hierarchical clustering with complete linkage and
multiple distance measures, k-means clustering with a chosen value k =
5, and k-means clustering based on principle components.

**Hierarchical Clustering**

Within hierarchical clustering, there are different types of linkage and
distance measures. I choose to use complete linkage because it tends to
perform well when there are well-defined clusters of data, which is
plausibly the case, as there are five concrete subtypes.

Choosing the distance measure well is also important. The distance
chosen has the potential to strongly impact the clustering. Deciding
which distance measure to use tells the algorithm how the similarity
between two elements is quantified, thus affecting how the clusters will
be chosen. In a study published in BMC Bioinformatics, a team of
researchers analyzed how different distances and clustering methods
interact with “regards to their ability to cluster gene expression.”
They use 15 different distance metrics. They noted that “the selection
of an appropriate distance measure can make the difference between
meaningful and poor clustering outcomes, even for a suitable clustering
method,” emphasizing that the choice of optimal distance measure depends
on the context (BMC).

The best choice of distance measure in any given scenario is not always
obvious, but there are some guidelines. A blog on DisplayR notes that
the choice should be based on “theoretical concerns from the domain of
study” (DisplayR) In other words, the distance metric should be able to
define similarity in a way that makes sense in the context of the data.
“For example, if clustering crime sites in a city, city block distance
may be appropriate. Or, better yet, the time taken to travel between
each location.” The authors note that Euclidean distance tends to be
reasonable when there is no obvious justification for any other
particular metric, as it is the way that we measure distance in the
physical world (DisplayR).

The first distance metric I try is Euclidean distance. When used, this
metric chooses the clusters in a way such that it minimizes the
following quantity:

$$\\frac{1}{(C_k)} \\sum\_{i, i’ \\in C_k}^{} \\sum\_{j=1}^{p} (x\_{ij} - x\_{i’j})^2$$

where the absolute value of *C*<sub>*k*</sub> refers to the number of
observations in the 4th cluster. *x*<sub>*i**j*</sub> and
*x*<sub>*i*′*j*</sub> are two vectors of length n. It is used in most
cases, when we are finding the distance between data rows that are
numeric (ints or floats).

Below I perform hierarhical clustering with complete linkage and
Euclidean distance, using scaled data.

``` r
#Hierarchical Clustering using Complete Linkage and Euclidean Distance 

library(dplyr)

#Use the data frame with 64 rows (patients) and 22,283 columns (genes). Create a new 
#data frame, leukemia_wo2_scaled. The scaled version has columns of data frame which 
#contain the original values subtracted by the column mean (center = TRUE), and 
#divided by the column standard deviations (scale = TRUE).  
leukemia_wo2_scaled <- scale(leukemia_wo2, center=TRUE, scale=TRUE)

#Create the Euclidean distance matrix.
leuk_wo2_eucdis <- dist(leukemia_wo2_scaled, method= "euclidean")

#hierarchical clustering using complete linkage
set.seed(1) 

leukemia.hclust = hclust(leuk_wo2_eucdis, method= "complete")

plot(leukemia.hclust, labels=leukemia$type, main='Hierarchical Clustering of Leukemia with Euclidean Distance')
```

![](Useful-Machine-Learning-Algorithms-in-Leukemia-Research---Github_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The dendogram above depicts the output of hierarchical clustering. Each
node represents a patient, and its label (AML, PB, PBSC_CD34,
Bone_Marrow, or Bone_Marrow_CD34) represents the patient subtype of
leukemia. Hierarchical clustering does a reasonable good job of
clustering the patients based on their subtype. However, it is not
perfect, as seen in the far right end, which depicts some patients with
AML (acute myeloid leukemia) in a cluster filled with patients with PB
(peripheral blood leukemia).

In hierarchical clustering, there is no fixed number of clusters; it is
possible to cut the tree at different heights in order to sever them
into clusters of varying sizes. Below, I cut the tree at a height which
creates a preferred number of groups (k = 5 groups).

``` r
clusters <- cutree(leukemia.hclust, k=5)
table(clusters)
```

    ## clusters
    ##  1  2  3  4  5 
    ## 42 11  9  1  1

``` r
table(leukemia$type)
```

    ## 
    ##              AML      Bone_Marrow Bone_Marrow_CD34               PB 
    ##               26               10                8               10 
    ##        PBSC_CD34 
    ##               10

The clusters, which are named arbitrarily, are of varying sizes. The
first cluster contains the vast majority of the patients (42 of them),
while the second contains 11, the third contains 9, and the last two
clusters contain one patient each. Below the clustering output, I print
the real subtypes for the patients. As shown, while most have AML
leukemia (26 patients), the rest of the subtypes are represented by a
similar number of patients. None of the subtypes have less than 8
patients. Upon further inspection of the original dataset, the
hierarhical method clusters many patients with different subtypes into
cluster 1. Thus, hierarchical clustering does not seem to be a very good
option.

I also tried the Pearson correlation measure as the distance measure.
The output was more difficult to interpret and it also performed very
flawed clustering.

**K-Means Clustering**

Now I try the k-means clustering method with chosen value of k = 5
clusters.

``` r
set.seed(16)
cluster_k <- kmeans(leukemia_wo2, 5)
```

The next step is to assess how well k-means clustering worked.

``` r
#Create a list which represents the true patient subtype. This will next be used to 
#compare the real subtype to the cluster which the patients were sorted in.

#BMCD34 = 1, BM = 2, AML = 3, PB = 4, PBSC_CD34
nrow(leukemia)
```

    ## [1] 64

``` r
types_list <- c()
for (i in 1:nrow(leukemia)){
  if (leukemia$type[i] == "Bone_Marrow_CD34"){
    types_list <- append(types_list, 1)
  }
  else if (leukemia$type[i] == "Bone_Marrow"){
    types_list <- append(types_list, 2)
  }
  else if (leukemia$type[i] == "AML"){
    types_list <- append(types_list, 3)
  }
  else if (leukemia$type[i] == "PB"){
    types_list <- append(types_list, 4)
  }
  else if (leukemia$type[i] == "PBSC_CD34"){
    types_list <- append(types_list, 5)
  }
}
types_list
```

    ##  [1] 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
    ## [39] 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5

``` r
table(cluster_k$cluster, types_list, dnn=c("Cluster","Class"))
```

    ##        Class
    ## Cluster  1  2  3  4  5
    ##       1  0  0 16  0  0
    ##       2  0  0  3  0 10
    ##       3  0  0  5  0  0
    ##       4  8  0  0  0  0
    ##       5  0 10  2 10  0

``` r
ggsave("table_now.png", device="png", width=7, height=2, units="in",dpi=300)
```

Above I compare the number of patients in each cluster and class to
assess whether the k-means algorithm sorted patients into clusters based
on their subtype reasonably well. Here, the cluster represents the
k-means cluster the patient was sorted into by the algorithm, and the
class represents their leukemia subtype.

It is important to note that the cluster number is arbitrary; hence it
is not important whether cluster 2 is equivalent to class 2. Overall,
the k-means algorithm performed relatively well. For instance, every
patient in class 1 (diagnosed with Bone Marrow CD34 leukemia) was sorted
into cluster 4. Every patient in class 2 (bone marrow leukemia) was
sorted into cluster 5. Though the clustering algorithm does not produce
a perfect one-to-one sorting, there are certainly a lot of positive
signs to suggest that being sorted in a particular cluster is associated
with having a certain subtype of leukemia. This makes it likely that it
will be useful to use the cluster as a variable in a logistic regression
model for classification of future patients.

Note that k-means clustering based on two principle components of the
data significantly increased the training error, such that I decided it
was not worth any advantages and was not the optimal option.

**Logistic Regression**

Logistic regression is a very useful, robust classification model. I
implement a logistic model which will theoretically have the capacity to
classify new patients into one of the five leukemia subtypes. The model
has one categorical variable, cluster, which contains the value of the
cluster which the patient was sorted into via k-means clustering. It
will be a multinomial logistic function because the dependent variable
has five levels (5 subtypes of leukemia).

One advantage of logistic regression is that it makes no assumptions
about the distributions of classes for its variables. This is useful
because since the clustering algorithm determined the classes, it is
difficult to determine their exact distribution. One potential drawback
is that logistic regression should not be used if the number of
observations is less than the number of features. However, this is not
the case here because there are 64 patients (observations) and only one
5-level variable (cluster).

``` r
#class (dependent variable)
class <- as.factor(types_list)

#cluster (independent variable)
cluster <- as.factor(cluster_k$cluster)
require(nnet)
```

    ## Loading required package: nnet

``` r
#Multinomial Logistic Model
multinom_model <- multinom(class ~ cluster)
```

    ## # weights:  30 (20 variable)
    ## initial  value 103.004026 
    ## iter  10 value 27.784784
    ## iter  20 value 27.588002
    ## iter  30 value 27.587592
    ## final  value 27.587592 
    ## converged

``` r
summary(multinom_model)
```

    ## Call:
    ## multinom(formula = class ~ cluster)
    ## 
    ## Coefficients:
    ##   (Intercept)  cluster2  cluster3  cluster4  cluster5
    ## 2  -0.9004305  3.505044  1.602115 -24.19701 27.262212
    ## 3  29.3878985 -9.448305 -4.119276 -52.89487 -4.635559
    ## 4  -0.9005490  3.505018  1.602112 -24.19719 27.262332
    ## 5  -1.5634752 22.707050  1.605320 -23.92393  6.734998
    ## 
    ## Std. Errors:
    ##   (Intercept)     cluster2     cluster3     cluster4     cluster5
    ## 2   0.1776998 6.757492e-09 6.834599e-12 3.061303e-12 1.776998e-01
    ## 3   0.2350200 1.917734e-01 1.623586e-11 1.756279e-11 2.762406e-01
    ## 4   0.1776998 6.756514e-09 6.833803e-12 3.060405e-12 1.776998e-01
    ## 5   0.1917734 1.917734e-01 1.793930e-12 2.207384e-12 4.769593e-10
    ## 
    ## Residual Deviance: 55.17518 
    ## AIC: 95.17518

``` r
set.seed(232)
#split data into train and test 
library(caret)
```

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

``` r
#create a dataframe containing the x and y variable as columns
clus_class <- data.frame(cluster, class)
index <- createDataPartition(clus_class[,2], p = .70, list = FALSE)
#create test and train data subsets
test <- clus_class[-index,]
train <- clus_class[index,]

train$ClassPredicted <- predict(multinom_model, newdata = train, "class")
# Create classification table
table_mn <- table(train[,2], train$ClassPredicted)
table_mn
```

    ##    
    ##      1  2  3  4  5
    ##   1  6  0  0  0  0
    ##   2  0  5  0  2  0
    ##   3  0  1 15  0  3
    ##   4  0  4  0  3  0
    ##   5  0  0  0  0  7

``` r
# To get test success rate, divide the sum of the diagonals (number of observations 
#for which class = predicted class) by the total number of observations 
round((sum(diag(table_mn))/sum(table_mn))*100,2)
```

    ## [1] 78.26

**Analysis**

Overall, the use of a logistic regression model with a k-means
clustering variable had a very effective, 78.26 percent training success
rate. In the output above, I print the summary of the multinomial
function, and a table to compare the patient’s leukemia subtype to their
predicted subtype under the model.

The summary of the multinomial logistic function is slightly difficult
to interpret, because cluster A is not always positively associated with
class A, as we would expect. This is likely because the intercepts
associated with each cluster vary.

Nevertheless, the function effectively performs its predictive goal,
with a training success rate of 78.26 percent. Note that in the
comparison table, the class numbers are not arbitrary and match up. The
table shows that the majority of the patients with a particular subtype
of leukemia are assigned the correct predicted class, and vice-versa. In
other words, the model’s predictions match reality quite well. For
instance, the model correctly classifies all of the patients with
leukemia subtype \#1 to be in group 1. Conversely, all of the patients
in the predicted group 1 have the leukemia subtype 1.

**Further Notes**

While the majority of this project focused on the use of machine
learning algorithms in the classification of patients based on cancer
subtype, some of the methods are useful for other purposes in cancer
research. For instance, one of the important goals which emerge when
researching a particular disease or form of cancer is to determine a
reasonable subset of genes which may be highly influential in mediating
or predicting cancer.

I explored the use of hierarchical clustering with complete linkage and
the Pearson correlation distance measure for identifying highly relevant
genes.

``` r
library(mlbench)
library(corrplot)
```

    ## corrplot 0.92 loaded

``` r
data(leukemia_wo2)
```

    ## Warning in data(leukemia_wo2): data set 'leukemia_wo2' not found

``` r
library(corrplot)
plot1 <- corrplot(cor(leukemia_wo2[,1:20]), 
                  method="square",
                  order="hclust", tl.cex=0.7, cl.cex=0.5, tl.col="black", addrect=2)
```

![](Useful-Machine-Learning-Algorithms-in-Leukemia-Research---Github_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
cor.info <- cor(leukemia_wo2[, 1:50])
sim.by.hclust <- hclust(dist(cor.info), method = "complete")
plot(sim.by.hclust)
```

![](Useful-Machine-Learning-Algorithms-in-Leukemia-Research---Github_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

This method helps us identify genes which are highly correlated, and
sorts them into clusters. Though this does not immediately tell us which
genes are influential in leukemia, it helps reduce the dimension. It may
be possible to determine certain clusters which are particularly
relevant. Above, I include both a correlation plot to depict highly
correlated genes, in which dark colors represent high correlation, as
well as the dendogram output of the hierarchical clustering.

**Conclusion/Extension of Analysis**

Overall, the main goal of this research was to explore machine learning
algorithms which could be useful in the context of cancer research, by
specifically focusing on extensive genetic data from a sample of
leukemia patients. While there are many possible methods for studying
genes, by for instance focusing on genetic mutations or regulatory
sequence interactions, in this project I focused on gene expression
data.

I found that clustering can serve as a very powerful, effective tool for
simplifying data into a smaller number of clusters. There are several
points to consider, such as clustering method, distance measure, and
linkage type. In this project, after trying several methods, I found
that k-means clustering with a chosen value of k=5 worked best at
clustering the data.

Clustering alone does not exactly provide a model for predicting the
leukemia subtype of future patients. However, incorporating the results
of clustering as a variable in a logistic classification model resulted
in a predictive training accuracy of about 78%. Though the test accuracy
would naturally be lower, this decent accuracy rate suggests that it
could be useful for predicting leukemia subtype. I posit that using a
combination of clustering and logistic regression might be useful in the
context of other cancers and diseases besides leukemia as well.

There is certainly room for growth in this model, and certain caveats to
be aware of. In this project, I tried two different types of clustering,
hierarchical clustering and k-means, as well as two distance metrics
(Euclidean and correlation). In future studies, it would be a good idea
to try out a greater number of combinations of clustering parameters.

Also, these methods might not be the most practical due to the fact that
they depend on a very large number of genes. Each of the patients in the
data set had the expression measurements for over 22,000 genes on file.
Doctors would not be likely to have this many gene readings available
for the average patient. That being said, this argument points to the
utility of isolating a smaller subset of highly relevant genes for a
particular disease. In my exploration I found that clustering based on
correlation might be a very useful way of isolating these highly
relevant genes. Once we isolate a subset of genes which have a strong
effect on leukemia (or another disease), we can build logistic models
based on the results of clustering on this smaller subset of genes. It
would be more practical for doctors to order tests on patients for this
smaller number of genes.

In my current model, the logistic function is only able to classify
based on one of the five leukemia subtypes. I would like to add data
such that it would be able to distinguish between leukemia and
non-leukemia patients as well.

Also, in future research, I would like to experiment with adding more
variables to the logistic model. I could add risk factors for leukemia,
like exposure to radiation and certain blood diseases, as well as other
genetic measurements like protein-protein interactions, or a binary
variable for the presence of certain genetic mutations. Considering that
the model already works quite well for patients in the data set, it is
likely that the incorporation of further relevant variables could boost
its function even further.

Overall, research like this is hugely important because a better
understanding of different diseases and forms of cancer is crucial to
improving the correct assessment of risk, the development of better
preventative care, and even the discovery of new treatment options.

**References**

“Acute Myeloid Leukemia Molecular Genetic Testing.” ARUP Consult, ARUP,
<https://arupconsult.com/ati/acute-myeloid-leukemia-molecular-genetic-testing>.

Asif, Muhammad, et al. “Identifying Disease Genes Using Machine Learning
and Gene Functional Similarities, Assessed through Gene Ontology.” PloS
One - PubMed Central , National Library of Medicine, 10 Dec. 2018,
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6287949/>.

Bock, Tim. “What Is Hierarchical Clustering?” DisplayR, 17 Jan. 2022,
<https://www.displayr.com/what-is-hierarchical-clustering/#>:\~:text=Hierarchical%20clustering%20starts%20by%20treating,the%20clusters%20are%20merged%20together.

“Facts and Statistics Overview: General Blood Cancers.” Leukemia &
Lymphoma Society ,
<https://www.lls.org/facts-and-statistics/facts-and-statistics-overview#>:\~:text=Approximately%2023%2C660%20deaths%20(13%2C900%20males,in%20females%20in%20the%20US.

Feltes, B.C.; Chandelier, E.B.; Grisci, B.I.; Dorn, M. CuMiDa: An
Extensively Curated Microarray Database for Benchmarking and Testing of
Machine Learning Approaches in Cancer Research. Journal of Computational
Biology, 2019.

“Genetic Testing: What You Should Know.” Family Doctor.org, 12
Aug. 2020,
<https://familydoctor.org/genetic-testing-what-you-should-know/>.

“Is Leukemia Hereditary? Causes, Risk Factors, and Prevention.” Medical
News Today, Healthline Media,
<https://www.medicalnewstoday.com/articles/325332>.

James, Gareth, Daniela Witten, Trevor Hastie, and Robert Tibshirani.
(2021). Introduction to Statistical Learning with Applications in R, 2nd
ed. Springer-Verlag.

Jaskowiak, Pablo A, et al. “On the Selection of Appropriate Distances
for Gene Expression Data Clustering.” BioMed Central Informatics, BioMed
Central, 24 Jan. 2014,
<https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-S2-S2>.

“Multinomial and Ordinal Logistic Regression in R.” Analytics Vidhya, 5
July 2020,
<https://www.analyticsvidhya.com/blog/2016/02/multinomial-ordinal-logistic-regression/>.

Puthier, Denis. “Distance Metrics and Clustering.” Pedagogix, 21
Nov. 2017,
<http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/distances/distances.html#>:\~:text=The%20Pearson’s%20correlation%20coefficient%20is,microarray%20or%20RNA%2DSeq%20data.

“Understanding Cancer Risk.” Cancer.Net, 2005-2022 American Society of
Clinical Oncology (ASCO), 28 Jan. 2022,
<https://www.cancer.net/navigating-cancer-care/prevention-and-healthy-living/understanding-cancer-risk>.

“What Is a Mutation?” Your Genome Topics, The Public Engagement Team at
the Wellcome Genome Campus, 21 July 2021,
<https://www.yourgenome.org/facts/what-is-a-mutation>.
