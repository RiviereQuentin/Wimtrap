# Wimtrap

Wimtrap: Integrative tools to predict the location of transcription factor binding sites

## Installation

Wimtrap is an R package that requires the last version of R (R 4.0.4), BiocManager and remotes to be installed. 

*Important: the installation of Wimtrap might take up to 1h if all the packages on which it depends need also to be installed.*

In R, type the following lines:
```
    if(!require("remotes", quietly = TRUE)){  
        install.packages("remotes")
        }
    if(!require("BiocManager", quietly = TRUE)){  
        install.packages("BiocManager")
        }
 ```
  
Then, you can enter:
```
  options(repos = BiocManager::repositories())
  getOption("repos")
  BiocManager::install("RiviereQuentin/Wimtrap",                     
    dependencies = TRUE,                     
    build_vignettes = TRUE,
    force = TRUE)
````

If an error occurs, it might be because the version of R, BiocManager and/or remotes is not updated. 

In some cases, we also noticed that it might be necessary to install beforehand the dependence e1071:

```
  if(!require("e1071", quietly = TRUE)){  
    install.packages("e1071")
    }
  BiocManager::install("RiviereQuentin/Wimtrap",                     
    dependencies = TRUE,                     
    build_vignettes = TRUE)
```
On linux, it might be also necessary to install as a prerequisite some software. On ubuntu, this might be achieved by entering the following in the terminal:

```
sudo apt install libcurl4-gnutls-dev icu-devtools libicu-dev libxml2-dev bzip2-doc libbz2-dev liblzma-dev
```
## Predictions of transcription factor binding sites in *Arabidopsis thaliana* or *Solanum lycopersicum*

Pre-built models are available for Arabidopsis and the tomato. 

Predictions can be made taking into consideration chromatin state features related to different conditions (for **Arabidopsis**: whole seedlings, seedling roots, non-hair part of seedling roots, flowers in stages 4-5, seed coats, heat-shocked seedlings, dark-grown seedlings, dark-grown seedlings exposed to 30min or 3h of light and dark-grown seedlings exposed to a long-day cycle; for the **tomato**: immature and ripening fruits).

To predict the binding sites of "AT2G46830" in flowers of Arabidopsis or of "Solyc00g0224680.1" in immature fruits of tomato, type:

```
    library(Wimtrap)
    #Predictions of the binding sites of "AT2G46830" in flowers of Arabidopsis
    CCA1predictions.flowers <- carepat(organism = "Arabidopsis thaliana",
                                       condition = "flowers",
                                       TFnames = "AT2G46830")
    #Predictions of the binding sites of "Solyc00g024680.1" in immature fruits of tomato
    DOF24predictions.immature <- carepat(organism = "Solanum lycopersicum",
                                         condition = "immature_fruits",
                                         TFnames = "Solyc00g024680.1")
```

## Build and apply a predictive model

Extensive details about the functionalities offered by Wimtrap are provide in the [user guide](https://htmlpreview.github.io/?https://github.com/RiviereQuentin/Wimtrap/blob/main/vignettes/Wimtrap.html) and in the manual pages of the functions, which can be accessed by entering in R:

```
    library(Wimtrap)
    ?importGenomicData()
    ?getTFBSdata()
    ?buildTFBSmodel()
    ?predictTFBS()
    ?plotPredictions()
```
