# Wimtrap

Wimtrap: Integrative tools to predict the location of transcription factor binding sites

## Installation

Installing Wimtrap requires at first the installation of BiocManager, if it has not been done yet:
    if(!require("BiocManager", quietly = TRUE)){  
        install.packages("BiocManager")
        }
  
The package can then be installed by typing the following:

  BiocManager::install("RiviereQuentin/Wimtrap",                     
  dependencies = TRUE,                     
  build_vignettes = TRUE)

If an error occurs, it might be fixed in some cases by installing beforehand the dependence e1071:

  if(!require("e1071", quietly = TRUE)){  
    install.packages("e1071")
    }
  BiocManager::install("RiviereQuentin/Wimtrap",                     
  dependencies = TRUE,                     
  build_vignettes = TRUE)
  
To take advantage of the package, please refer to the user guide available on our website: https://lpgmp.ulb.be/ressources/
