# PALMER
<!--
PALMER: A Constrained Biclustering Algorithm to Improve Pathway Annotation Based on the Biomedical Literature Mining
-->

PALMER (a constrained biclustering algorithm to improve **P**athway **A**nnotation based on the biomedical **L**iterature **M**ining) is a constrained biclustering approach that allows to identify indirect relationships among genes based on the text mining of biomedical literature, which allows researchers to utilize prior biological knowledge to guide identification of gene-gene associations.
'palmer' package provides computationally efficient and user friendly interface to fit the PALMER models. 
[The 'palmer' vignette](https://github.com/dongjunchung/chunglab_binary_packages/blob/master/palmer.pdf?raw=true) provides a good start point for the step-by-step data analysis using 'palmer' package.The following help pages provide a good start point for the genetic analysis using the 'GPA' package, including the overview of 'GPA' package and the example command lines:

```
library(palmer)
package?palmer
class?palmer
vignette("palmer")
```

Installation
============ 

The stable versions of 'palmer' package can be obtained from the following URLs:

Package source: [https://github.com/dongjunchung/chunglab_binary_packages/blob/master/palmer_0.1.tar.gz](https://github.com/dongjunchung/chunglab_binary_packages/blob/master/palmer_0.1.tar.gz?raw=true)

Windows binary: [https://github.com/dongjunchung/chunglab_binary_packages/blob/master/palmer_0.1.zip](https://github.com/dongjunchung/chunglab_binary_packages/blob/master/palmer_0.1.zip?raw=true)

Mac OS/X binary: [comming soon](https://?raw=true)

To install the developmental versions of 'palmer' package, it's easiest to use the 'devtools' package.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/palmer", build_vignettes= TRUE)
```

References
==========
Nam JH, Couch D, Silveira W.A, Yu Z and Chung D (2019) ''PALMER: A constrained biclustering Algorithm to improve pathway annotation based on the biomedical literature mining''.


