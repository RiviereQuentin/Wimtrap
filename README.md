# Wimtrap

Wimtrap: Integrative tools to predict the location of transcription factor binding sites

## Installation

Wimtrap is an R package that requires BiocManager to be installed. 

In R, type the following lines:
```
    if(!require("BiocManager", quietly = TRUE)){  
        install.packages("BiocManager")
        }
 ```
  
Then, you can enter:
```
  BiocManager::install("RiviereQuentin/Wimtrap",                     
  dependencies = TRUE,                     
  build_vignettes = TRUE)
````

If an error occurs, it might be fixed in some cases by installing beforehand the dependence e1071:

```
  if(!require("e1071", quietly = TRUE)){  
    install.packages("e1071")
    }
  BiocManager::install("RiviereQuentin/Wimtrap",                     
  dependencies = TRUE,                     
  build_vignettes = TRUE)
```
## Predictions of transcription factor binding sites in *Arabidopsis thaliana* or *Solanum lycopersicum*

Pre-built models are available for Arabidopsis and the tomato. Predictions can be made taking into consideration chromatin state features related to different conditions ().

## Build a predictive model

Extensive details about the functionalities offered by Wimtrap are provide in the [user guide](http://lpgmp.ulb.be/wp-content/uploads/2021/02/Wimtrap.pdf)
