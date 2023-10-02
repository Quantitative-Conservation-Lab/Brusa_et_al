## Correcting for measurement errors in a long-term aerial survey with auxiliary photographic data  

#### *In Review*  

#### Jamie Brusa, Matthew T. Farr, Joseph Evenson, Bryan Murphie, Thomas A. Cyra, Heather J. Tschaekofske, Kyle A. Spragens, Emily Silverman, Sarah J. Converse  

##### Please contact the first author for questions about the code or data: Jamie Brusa (jlbwcc@gmail.com)
_______________________________________________________________________________________

## Abstract

Long-term, large-scale monitoring of wildlife populations is an integral part of conservation research and management. However, some traditional monitoring protocols lack the information needed to account for sources of measurement error in data analyses. Ignoring measurement error, such as partial availability, imperfect detection, and species misidentification, can lead to mischaracterizations of population states and processes. Accounting for measurement error is key to robust monitoring of populations, which can inform a wide variety of decisions, including harvest, habitat restoration, and determination of the legal status of species. We undertook an effort to retroactively minimize bias in a large-scale, long-term monitoring program for marine birds in the Salish Sea, Washington, USA by conducting an auxiliary study to jointly estimate components of measurement error. We built a novel model in a Bayesian framework that simultaneously harnessed human observer and photographic data types to produce estimates necessary to correct for the effects of partial availability, imperfect detection, and species misidentification. Across all 31 species identified in photographs, both of the two participating observers had instances of undercounting and overcounting birds but tended to undercount (observers undercounted totals across all species on 69.3% – 78.9% of transects). We estimated species-specific correction factors that can be used to correct both historical and future counts from the Salish Sea survey, which has been running since 1992. Our novel modeling framework can be applied in other multi-species monitoring contexts where minimal photographic data can be collected for the purposes of correcting for measurement error in large-scale, long-term datasets.  

### Table of Contents

### [Scripts](./scripts)

+ [Marine birds model code with data.R](<./scripts/Marine birds model code with data.R>) does ...
+ [Marine birds out of sample sim code.R](<./scripts/Marine birds out of sample sim code.R>) does ...
+ [Marine birds simulation model code.R](<./scripts/Marine birds simulation model code.R>) does ...

### [Data](./data)

+ [FF.csv](./data/FF.csv) are forward-facing camera data for the fine-grained (i.e., species-specific) analysis.
+ [Obs1.csv](./data/Obs1.csv) are Observer One's data for the fine-grained (i.e., species-specific) analysis.
+ [Obs2.csv](./data/Obs2.csv) are Observer Two's data for the fine-grained (i.e., species-specific) analysis.
+ [FF_c.csv](./data/FF_c.csv) are forward-facing camera data for the coarse-grained (i.e., species-groups) analysis.
+ [Obs1_c.csv](./data/Obs1_c.csv) are Observer One's data for the the coarse-grained (i.e., species-groups) analysis.
+ [Obs2_c.csv](./data/Obs2_c.csv) are Observer Two's data for the the coarse-grained (i.e., species-groups) analysis.
+ [seat.csv](./data/seat.csv) are the seat position during a flight transect for each observer.

### [Figures](./figures)

+ [Figure 1](./figures/Figure1.png) is a directed acyclic graph showing relationship between in-sample data and out-of-sample correction.
+ [Figure 2](./figures/Figure2.png) is the observer count accuracy per species.
+ [Figure 3](./figures/Figure3.png) is the model performance for count corrections.

### Required Packages and Versions Used
tidyverse
nimble
ggmcmc
coda
MCMCvis
viridis
here

### Details of Article

### How to Use this Repository
The model code can be used to estimate correction factors for count data of marine bird species taken from aerial surveys. We provide code for simulations (as a proof-of-concept) and to use with our data, which are provided in the data folder. The model code in the "Marine birds model code with data.R" script can be adapted to fit data for future surveys in the Salish Sea or data for other studies using a similar protocol.
