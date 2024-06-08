# VARDetect

This repository is for R package **VARDetect**. 



##  Install the R package 

To install the package from R CRAN(https://CRAN.R-project.org/package=VARDetect):
```
install.packages("VARDetect")
```


You could download the R package [VARDetect_0.1.8.tar.gz](VARDetect_0.1.8.tar.gz) and locally install the package using the tar.gz File:
```
install.packages("VARDetect_0.1.8.tar.gz", repos = NULL, type="source")
```

To install the latest version of the package from GitHub:
```
library(devtools)
devtools::install_github("peiliangbai92/VARDetect")
```

To cite package **VARDetect** in publications use:

  Peiliang Bai, Yue Bai, Abolfazl Safikhani and George Michailidis (2021). VARDetect: Multiple Change Point Detection in Structural VAR Models
  R package version 0.1.8. https://CRAN.R-project.org/package=VARDetect

A BibTeX entry for LaTeX users is
```
  @Manual{,
    title = {VARDetect: Multiple Change Point Detection in Structural VAR Models},
    author = {Peiliang Bai, Yue Bai, Abolfazl Safikhani and George Michailidis},
    year = {2021},
    note = {R package version 0.1.8},
    url = {https://CRAN.R-project.org/package=VARDetect},
  }
 ```


To apply the dynamical programming method introduced by [Wang et al. (2019)](https://arxiv.org/abs/1909.06359), one can find the functions in the file ["VARDetect_DP.R"](https://github.com/peiliangbai92/VARDetect/blob/main/VARDetect_DP.R). 
