Modeling the Effects of Nutrition with Mixed-Effect Bayesian Network
--------------------------------------------------------------------

This work proposes Mixed-Effect Bayesian Network (MEBN) as a method for modeling the effects of nutrition. It allows to identify both typical and personal correlations between nutrients and their bodily responses. Predicting a personal network of nutritional reactions would allow interesting applications at personal diets and in understanding this complex system. Brief theory of MEBN is given followed by implementation in R and Stan. A real life dataset from a nutritional Sysdimet study is then analyzed with the method and the results are visualized with an interactive JavaScript-visualization.

[View HTML-version of the notebook here](http://htmlpreview.github.io/?https://github.com/turkiaj/StanCon2018/blob/master/PersonalEffectsOfNutrition.html)

Contents
--------

Main body of the presentation is found here as fully functional RMarkdown notebook, and also HTML and PDF renderitions of it. The notebook uses my R-package "MEBN" for constructing a Mixed-Effect Bayesian Network that models the effects of nutrition from a real life nutritional dataset. The fully Bayesian estimation of the local mixed-effect models at the network is done with Stan. Visualization of the network is created with JavaScript library sigma.js and some customizations of it.

-   PersonalEffectsOfNutrition.Rmd : RMarkdown notebook
-   PersonalEffectsOfNutrition.html : HTML knit of the notebook with active JavaScript visualization
-   PersonalEffectsOfNutrition.pdf : PDF knit of the notebook
-   population\_graph.htm : HTML document for graph visualization
-   biblio.bib : References for the notebook in BibLatex format
-   data-folder: The dataset that is used in the notebook analysis
-   Data description.csv : CSV file that describes the metadata for analyzed dataset
-   mebn-folder: The R-code for constructing MEBN and stan-model definitions
-   visualization-folder: JavaScript code for graph visualization
-   models-folder: Sampled Stan-models are cached in this folder once the notebook is executed. It is empty by default.

License
-------

Material in this repository is licensed under CC 4.0 <https://creativecommons.org/licenses/by/4.0/>

Requirements
------------

Besides working installation of Stan, following R-packages are required for this notebook to execute correctly.

Use install\_packages.r script to install.

\*\* Session Info \*\*

R version 3.4.4 (2018-03-15) Platform: x86\_64-w64-mingw32/x64 (64-bit) Running under: Windows Server 2008 R2 x64 (build 7601) Service Pack 1

Matrix products: default

locale: \[1\] LC\_COLLATE=Finnish\_Finland.1252 LC\_CTYPE=Finnish\_Finland.1252
\[3\] LC\_MONETARY=Finnish\_Finland.1252 LC\_NUMERIC=C
\[5\] LC\_TIME=Finnish\_Finland.1252

attached base packages: \[1\] stats graphics grDevices utils datasets methods base

loaded via a namespace (and not attached): \[1\] Rcpp\_0.12.16 digest\_0.6.15 rprojroot\_1.3-2 plyr\_1.8.4 grid\_3.4.4
\[6\] gtable\_0.2.0 backports\_1.1.2 magrittr\_1.5 evaluate\_0.10.1 scales\_0.5.0
\[11\] pillar\_1.2.1 ggplot2\_2.2.1 rlang\_0.2.0 stringi\_1.1.7 lazyeval\_0.2.1
\[16\] rmarkdown\_1.9 tools\_3.4.4 stringr\_1.3.0 munsell\_0.4.3 yaml\_2.1.19
\[21\] compiler\_3.4.4 colorspace\_1.3-2 htmltools\_0.3.6 knitr\_1.20 tibble\_1.4.2
