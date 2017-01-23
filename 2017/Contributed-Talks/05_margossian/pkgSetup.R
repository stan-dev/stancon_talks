## Package Setup for Stan Con. (Margossian & Gillespie)
pkgs <- c('gglplot2', 'rstan', 'plyr', 'dplyr', 'tidyr', 'tibble', 'DBI', 'StanHeaders',
          'metrumrg', 'tibble')

user <- Sys.info()["user"]
parentScriptDir <- getwd()  ## You may need to mod this to be in the top level of scriptDir
if(!exists("pkgDir")) pkgDir <- file.path(parentScriptDir, "pkg")
libDir <- file.path(parentScriptDir, "lib")
dir.create(libDir)
dir.create(pkgDir)

.libPaths(libDir)
mycran <- paste("file://", pkgDir, sep="")
library(tools)
if(file.exists(file.path(pkgDir,"PACKAGES"))){ 
  available <- available.packages(contriburl=mycran)[,"Package"]
}else{
  available <- NULL
}


## Only authors can install from CRAN and write_PACKAGES
##fromCRAN <- user %in% author
fromCRAN <- TRUE
if(fromCRAN){
    newpkgs <- setdiff(pkgs, available)
    write_PACKAGES(pkgDir)
   if("audited" %in% newpkgs){
      download.file("https://metrumrg-soft.s3.amazonaws.com/audited/audited_1.9.tar.gz",destfile=file.path(pkgDir,"audited_1.9.tar.gz"))
      write_PACKAGES(pkgDir)
    }
    if("fork" %in% newpkgs){
      download.file("https://metrumrg-soft.s3.amazonaws.com/fork/fork_1.2.4.tar.gz",destfile=file.path(pkgDir,"fork_1.2.4.tar.gz"))
      write_PACKAGES(pkgDir)
    }
    if("review" %in% newpkgs){
      download.file("https://metrumrg-soft.s3.amazonaws.com/review/review_2.5.tar.gz",destfile=file.path(pkgDir,"review_2.5.tar.gz"))
      write_PACKAGES(pkgDir)
    }
    if(length(newpkgs)>0){
      install.packages(newpkgs,
                       lib=libDir,
                       contriburl=c(mycran,
                         contrib.url("http://R-Forge.R-project.org","source"),
                         contrib.url("http://cran.rstudio.com","source")),
                       destdir=pkgDir,
                       type="source",
                       INSTALL_opts="--no-multiarch")
              write_PACKAGES(pkgDir)
            }
    ## If multiple authors qcing each other, a package could be available but uninstalled.  Install from local.
    uninstalled <- setdiff(pkgs, installed.packages(libDir))
    if(length(uninstalled)>0){
      install.packages(uninstalled,
                       lib = libDir,
                       contriburl = mycran,
                       type = "source",
                       INSTALL_opts="--no-multiarch")
    }    
  }
if(!fromCRAN){
  installed <- row.names(installed.packages(libDir))
  newpkgs <- setdiff(pkgs, installed)
  if(length(newpkgs)>0){
    if("metrumrg" %in% newpkgs){
                                        # XML frequently fails to build from source.  Use the binary and place it in lib.
      if(!"XML" %in% rownames(installed.packages())) install.packages("XML",lib=libDir,repos="https://cran.rstudio.com/",INSTALL_opts="--no-multiarch",destdir=pkgDir)
      install.packages('metrumrg', repos=c('http://R-Forge.R-project.org','https://cran.rstudio.org'))
    }
    
    install.packages(newpkgs,
                     lib = libDir,
                     contriburl = mycran,
                     type = "source",
                     INSTALL_opts="--no-multiarch")
  }
}