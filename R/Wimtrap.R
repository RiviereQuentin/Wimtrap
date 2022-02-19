#' Imports genomic data into the R session
#' @export
#' @description
#' Imports genomic data that will allow to define contextual features around matches with the primary motifs of
#' transcription factors. Data about gene structures might be optionally automatically downloaded from Biomart.
#' @param organism Binomial name of the organism. Can be set to \code{NULL} if you provide
#' the location of the transcription start sites (TSS), transcription termination site (TTS) and structures
#' of the protein-coding genes of the organism (see the arguments  \code{tss}, \code{tts} and \code{genomic_data}).
#' @param genomic_data A named character vector defining the local paths to BED files describing genomic features.
#' The vector has to be named according to the features described by the files indicated. All the data related to
#' the chromatin state have to be specific of the samed condition. The properties of the BED files are the following,
#' depending on the type of feature: if 'numeric', the score field of the file is fulfilled or empty - if empty, the
#' score will automatically be set to '1'; if 'categorical', the score field of the file is empty while its name field
#' is fulfilled with the name of the categories of features. If you want to input the location of gene structures, name
#' the file paths with exactly the following names: '\code{ProximalPromoter}', '\code{Promoter}', '\code{X5UTR}', '\code{X3UTR}',
#' '\code{CDS}', '\code{Intron}', '\code{Downstream}'. If you use these names, the potential binding sites will be annotated only
#' with the gene structure that their centers overlap. If you don not use these names, the data related to gene structure will
#' be extracted in the same way than for the others (average of the signal on windows of 3 different lengths centered on the potential binding sites).
#' 'X5UTR' stands for 5'Untranslated Region, 'X3UTR' for '3'Untranslated Region', 'CDS' for coding sequence, 'Donwstream'
#' for the regions downstream from the transcription termination site.
#' @param biomart Logical. Should be automatically downloaded through biomart the location of the transcription
#' start sites (TSS), transcription termination site (TTS) and structures of the protein-coding genes of the organism?
#' Default is \code{TRUE}.
#' @param tss \code{NULL} (by default) or local path to a BED file defining the transcription stat site (TSS), name and
#' orientation of each protein-coding transcript of the organism. The default value allows to download automatically
#' these informations from biomart.
#' @param tts \code{NULL} (by default) or local path to a BED file defining the transcription termination site (TTS), name
#' and orientation of each protein-coding transcript of the organism. The default value allows to download automatically
#' these informations from biomart.
#' @param promoter_length An integer setting the length of the promoters. By default, the promoter is defined as the region spanning the 2000bp
#' upstream of the transcription start site (TSS).
#' @param downstream_length An integer setting the length of the downstream regions. By default, the downstream region is defined as the region
#' spanning the 1000bp downstream of the transcription termnination site (TTS).
#' @param proximal_length An integer setting the length of the proximal promoters. By default, the proximal promoter is defined as the region
#' spanning the 500bp upstream of the transcription start site (TSS).
#' @return A named list of \code{GRanges} objects. Each component of the list is named according to the nature of the genomic
#' data. The last two components describe respectively the position of the transcritpion termination site (TTS) and the
#' transcription start site (TSS) of each protein-coding transcripts of the organism considered.
#' @examples
#'## Without automatic download from Biomart of data related to gene structure
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("extdata/x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("extdata/cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("extdata/intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("extdata/x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex,
#'                                               tss = system.file("extdata/tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("extdata/tts_example.bed", package = "Wimtrap"))
#'
#' ##With automatic download from Biomart of data related to gene structure
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                      DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                      DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(organism = "Arabidopsis thaliana",
#'                                               biomart = TRUE,
#'                                               genomic_data = genomic_data.ex)
importGenomicData <- function(organism = NULL,
                              genomic_data,
                              biomart = TRUE,
                              tss = NULL,
                              tts = NULL,
                              promoter_length = 2000,
                              downstream_length = 1000,
                              proximal_length = 500)
{
  #@ Defining all the parameters required by the function
  #@ _inputGenomicData()_
  if (biomart == TRUE) {
    if (length(organism)==1){
      proximalpromoter <- "biomart"
      promoter <- "biomart"
      x5utr <- "biomart"
      cds <- "biomart"
      intron <- "biomart"
      x3utr <- "biomart"
      downstream <- "biomart"
    } else {
      if(is.na(genomic_data["ProximalPromoter"])) {proximalpromoter <- NULL} else {proximalpromoter <- genomic_data["ProximalPromoter"]}
      if(is.na(genomic_data["Promoter"])) {promoter <- NULL} else {promoter <- genomic_data["Promoter"]}
      if(is.na(genomic_data["X5UTR"])) {x5utr <- NULL} else {x5utr <- genomic_data["X5UTR"]}
      if(is.na(genomic_data["CDS"])) {cds <- NULL} else {cds <- genomic_data["CDS"]}
      if(is.na(genomic_data["Intron"])) {intron <- NULL} else {intron <- genomic_data["Intron"]}
      if(is.na(genomic_data["Downstream"])) {downstream <- NULL} else {downstream <- genomic_data["Downstream"]}
      if(is.na(genomic_data["X3UTR"])) {x3utr <- NULL} else {x3utr <- genomic_data["X3UTR"]}
      }
  } else {
    if(is.na(genomic_data["ProximalPromoter"])) {proximalpromoter <- NULL} else {proximalpromoter <- genomic_data["ProximalPromoter"]}
    if(is.na(genomic_data["Promoter"])) {promoter <- NULL} else {promoter <- genomic_data["Promoter"]}
    if(is.na(genomic_data["X5UTR"])) {x5utr <- NULL} else {x5utr <- genomic_data["X5UTR"]}
    if(is.na(genomic_data["CDS"])) {cds <- NULL} else {cds <- genomic_data["CDS"]}
    if(is.na(genomic_data["Intron"])) {intron <- NULL} else {intron <- genomic_data["Intron"]}
    if(is.na(genomic_data["Downstream"])) {downstream <- NULL} else {downstream <- genomic_data["Downstream"]}
    if(is.na(genomic_data["X3UTR"])) {x3utr <- NULL} else {x3utr <- genomic_data["X3UTR"]}
  }

  if (length(tss) == 0){
    tss <- "biomart"
    if (length(organism)==0){
      tss <- NULL
    }
  }
  if (length(tts) == 0){
    tts <- "biomart"
    if (length(organism)==0){
      tts <- NULL
    }
  }
  DNAregionsName = "Motifs"

  #@ The location of each structural genomic feature
  #@ as well as the location of the TSS of protein-coding
  #@ genes are input separetely but are from now
  #@ put together. The _annotateOccurences_ function
  #@ will then determine which ones have to be ignored,
  #@ downloaded or imported from source files.
  structuralFeatures <- list(proximalpromoter,promoter, x5utr, cds,
                             intron, x3utr, downstream, tts, tss)
  names(structuralFeatures) <-
    c("ProximalPromoter",
      "Promoter",
      "X5UTR",
      "CDS",
      "Intron",
      "X3UTR",
      "Downstream",
      "TTS",
      "TSS")

  #@ Additional features can be eiher ignored,
  #@ imported from source file(s) or input as a list of
  #@ GRanges objects.
  if (is.null(genomic_data)) {
    directInput_featureGRanges <- NULL
    feature_paths <- NULL
  } else {
    if (class(genomic_data) == "character") {
      feature_paths <- genomic_data[!(names(genomic_data) %in%  c("ProximalPromoter",
                                                                  "Promoter",
                                                                  "X5UTR",
                                                                  "CDS",
                                                                  "Intron",
                                                                  "X3UTR",
                                                                  "Downstream",
                                                                  "TTS",
                                                                  "TSS"))]
      directInput_featureGRanges <- NULL
    } else {
      feature_paths <- NULL
      directInput_featureGRanges <- genomic_datagenomic_data[!(names(genomic_data) %in%  c("ProximalPromoter",
                                                                                           "Promoter",
                                                                                           "X5UTR",
                                                                                           "CDS",
                                                                                           "Intron",
                                                                                           "X3UTR",
                                                                                           "Downstream",
                                                                                           "TTS",
                                                                                           "TSS"))]
    }
  }

  ####################################################################

  ####################################################################
  ### Carry on all the preliminary steps to the modeling: import the
  ### data, identify the potential binding sites by pattern-matching,
  ### annotate them with the considered features and label them as
  ### TRUE or FALSE according to the ChIP-seq data.
  ### Return as well the objects storing the sources data and options
  dna_regions_name = DNAregionsName

  message("Importing genomic features...")

  #@ The first part of the code is mainly format ckecking
  #@ Globally, data can be input from objects in the R environment
  #@ (cf. "directInput") or from a source file (cf. "paths").
  #@ In certain cases data can be as well downloaded

  #@ Treatment of additional features
  if (length(directInput_featureGRanges) == 0) {
    if (!missing(feature_paths) | !(is.null(feature_paths))) {
      FeatureRanges_list <- importFeatures(feature_paths)
      if (!(is.null(names(feature_paths)))){
        names(FeatureRanges_list) <- names(feature_paths)
      }
    } else {

    }
  } else {
    if (is.list(directInput_featureGRanges)) {
      FeatureRanges_list <- directInput_featureGRanges
    } else {
      FeatureRanges_list <- list(directInput_featureGRanges)
    }
  }

  #@ Treatment of the structural genomic features

  #@ Identify the location of which genonmic
  #@ structure(s) has to be downloaded (OPTIONAL)

  UseOfBiomart <- rep(0, length(structuralFeatures))
  for (i in seq_along(UseOfBiomart)) {
    considered <- structuralFeatures[[i]][1]
    if (is.character(considered)) {
      UseOfBiomart[i] <- (considered == "biomart")
    } else {

    }
  }
  if (length(UseOfBiomart[UseOfBiomart == 1]) > 0) {
    whichToBiomart <- which(UseOfBiomart == 1)
    UseOfBiomart <- TRUE
  } else {
    UseOfBiomart <- FALSE
  }

  if (UseOfBiomart & (length(organism)==1)) {
    message("Downloading structural genomic ranges...")
    #@ Load the biomart database
    Ensembl <-
      loadEnsembl(organism, NULL, NULL)
    #@ Query the database
    StructureBiomart <-
      getStructure(Ensembl, promoter_length, downstream_length, proximal_length)
  } else {

  }

  #@ List the location of each genomic structures in "StructuralFeatures"
  #@ from the specified source (GRanges, source file, biomart)
  if (length(structuralFeatures) > 0) {
    for (structure in seq_along(structuralFeatures)) {
      considered <- structuralFeatures[[structure]]
      if (length(considered) > 1) {
        #@ The input is a GRanges objectS
      } else {
        if (length(considered) == 1 & is.character(considered)) {
          if (considered == "biomart") {
            structuralFeatures[[structure]] <- StructureBiomart[[structure]]
          } else {
              structuralFeatures[[structure]] <- rtracklayer::import(structuralFeatures[[structure]])
          }
        } else {

        }
      }
    }
  }

  #@ Change the name of "structuralFeatures" to "StructureTracks"
  if (biomart & length(organism)==1){
    StructureTracks <- StructureBiomart
  } else {
    StructureTracks <- structuralFeatures
    if (!all(unlist(lapply(StructureTracks, is.null)))){
      StructureTracks <- StructureTracks[!unlist(lapply(StructureTracks, is.null))]
      GenomeInfoDb::seqlevels(StructureTracks$TSS) <- getRiddChr(GenomeInfoDb::seqlevels(StructureTracks$TSS))

      if (is.null(StructureTracks$Promoter)){
        PartA <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="+"])
        GenomicRanges::start(PartA) <- GenomicRanges::start(PartA) - promoter_length
        PartB <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="-"])
        GenomicRanges::end(PartB) <- GenomicRanges::start(PartB) + promoter_length
        StructureTracks$Promoter <- c(PartA, PartB)
      }
      if (is.null(StructureTracks$Downstream)){
        PartA <- GenomicRanges::granges(StructureTracks$TTS[GenomicRanges::strand(StructureTracks$TTS)=="+"])
        GenomicRanges::end(PartA) <- GenomicRanges::end(PartA) + downstream_length
        PartB <- GenomicRanges::granges(StructureTracks$TTS[GenomicRanges::strand(StructureTracks$TTS)=="-"])
        GenomicRanges::start(PartB) <- GenomicRanges::start(PartB) - downstream_length
        StructureTracks$Downstream <- c(PartA, PartB)
      }
      if (is.null(StructureTracks$ProximalPromoter)){
        PartA <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="+"])
        GenomicRanges::start(PartA) <- GenomicRanges::start(PartA) - proximal_length
        PartB <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="-"])
        GenomicRanges::end(PartB) <- GenomicRanges::start(PartB) + proximal_length
        StructureTracks$ProximalPromoter <- c(PartA, PartB)
      }
      fulfilled <- lapply(StructureTracks, is.null)
      StructureTracks <- StructureTracks[unlist(fulfilled)==FALSE]
      StructureTracks <- StructureTracks[c(which(!(names(StructureTracks) %in% c("TTS", "TSS"))),
                                           which((names(StructureTracks) %in% c("TTS", "TSS"))))]
    } else {
      StructureTracks <- NULL
    }
  }


  sourceData <- c(
    FeatureRanges_list,
    StructureTracks
  )

  for (feature in 1:length(sourceData)){
    sourceData[[feature]] <- sourceData[[feature]][GenomicRanges::order(sourceData[[feature]]),]
  }

  return(sourceData)
}


#' Performs pattern-matching and extraction of genomic features at match location
#' @export
#' @importFrom rtracklayer mcols import
#' @importFrom GenomicRanges strand makeGRangesFromDataFrame GRanges reduce as.data.frame start end resize findOverlaps nearest distanceToNearest
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings readDNAStringSet width Views matchPWM reverseComplement letterFrequency
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom utils read.csv
#' @importFrom biomartr getENSEMBLGENOMESInfo getENSEMBLInfo getGenome getKingdoms
#' @importFrom readr read_tsv cols col_character col_integer col_date write_tsv
#' @importFrom stringr str_to_lower str_replace_all str_to_upper str_sub str_replace_all str_detect
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom utils download.file
#' @importFrom dplyr bind_rows mutate filter setdiff
#' @importFrom curl curl_fetch_memory
#' @importFrom tibble as_tibble
#' @importFrom biomaRt listEnsembl listEnsemblGenomes useEnsemblGenomes searchDatasets getBM
#' @importFrom methods as
#' @importFrom IRanges IRanges resize
#' @importFrom readr write_tsv
#' @importFrom data.table fread data.table
#' @importFrom graphics hist
#' @importFrom stats runif median IQR
#' @importFrom universalmotif read_cisbp read_homer read_jaspar read_meme read_transfac
#' @importFrom httr GET content_type stop_for_status content
#' @description
#' `getTFBSdata` writes files encoding for datasets characterizing the genomic context
#' around motif occurences along the genome for the considered transcription factors
#' (training and/or studied TFs).
#' @param  pfm Path to a file including the position frequency or weight matrices (PFMs or PWMs) of the motifs recognized
#' by the considered transcription factors (training and/or studied TFs). This file can be in different formats, determined based
#' on the file extension: raw pfm (".pfm"), jaspar (".jaspar"), meme (".meme"), transfac (".transfac"), homer (".motif") or
#' cis-bp (".txt").\code{pfm} can be set to \code{NULL} (default value) if you provide the results of pattern-matching
#' obtained from an external source (see the argument \code{matches}).
#' @param TFnames names of the considered transcription factors among those described in the \code{pfm} file.
#' @param organism Binomial name of the organism. Can be set to \code{NULL} if you provide the genome sequence
#' (see the argument \code{genome_sequence})
#' @param genome_sequence "getGenome" (by default) or local path to a FASTA file encoding the genomic sequence
#'  of the organism. The default value allows the automatic download of the genomic sequence (when \code{organism} is input)
#'  from ENSEMBL or ENSEMBL GENOMES.
#' @param imported_genomic_data An object output by [importGenomicData()] and that includes data related to the chromatin
#' state that are specific to the training or studied condition.
#' @param matches \code{NULL} (by default) or a named list of \code{GRanges} objects. Each \code{GRanges} object is related to a
#' given transcritpion factor and defines the location along the genome of the matches with the primary motif of the latter.
#' The \code{GRanges} objects contain also a metadata column named '\code{matchLogPval}'that gives the p-value of the matches.
#' The list input through \code{matches} has to be named according to the names of the transcription factors considered.
#' These names have to be consistent with those provided through the \code{ChIP-peaks} argument. The default value allows to
#' perform the pattern-matching analysis with the function encoded by the Wimtrap package.
#' @param strand_as_feature A logical. Should be considered as feature the orientation of the matches in relation to the
#' direction of transcription of the closest transcript? Default is \code{FALSE}.
#' @param pval_threshold P-value threshold to identify the matches with the primary motif of the transcription factors.
#' Default is set to 0.001.
#' @param short_window An integer (20 by default). Sets the length of the short-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param medium_window An integer (400 by default). Sets the length of the medium-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param long_window An integer (1000 by default). Sets the length of the long-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @seealso [importGenomicData()] for importing genomic data and [buildTFBSmodel()] to train a predictive model of
#' transcription factor binding sites.
#' @return A vector indicating the local paths to the tab-delimited files in which are written the results of pattern-matching
#' and genomic feature extraction for each of the transcription factors considered. The 5 first fields of these files describe
#' the location of the potential binding sites identified by pattern-matching. The following fields contain the
#' raw score and/or p-value of the matches, and the the genomic features extracted at location of the matches
#' on short-, medium- and long-ranges-centered windows the label ('1' = "positive" = "ChIP-validated in
#' the considered condition" or '0' = "negative") for the the training TFsr.
#' @examples
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("extdata/x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("extdata/cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("extdata/intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("extdata/x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex,
#'                                               tss = system.file("extdata/tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("extdata/tts_example.bed", package = "Wimtrap"))
#' TFBSdata.ex <- getTFBSdata(pfm = system.file("extdata/pfm_example.pfm", package = "Wimtrap"),
#'                            TFnames = c("PIF3", "TOC1"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("extdata/genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
getTFBSdata <- function(pfm = NULL,
                        TFnames = NULL,
                        organism = NULL,
                        genome_sequence = "getGenome",
                        imported_genomic_data,
                        matches = NULL,
                        strand_as_feature = FALSE,
                        pval_threshold = 0.001,
                        short_window  = 20,
                        medium_window = 400,
                        long_window = 1000)
{
  #@ Import the source data into the R session
  sourceData <- inputSourceData(NULL,
                                pfm,
                                organism,
                                genome_sequence,
                                matches,
                                strand_as_feature,
                                TFnames)
  #@ Build the dataset
  DataSet <- getCandidatesRegions(directInput_matches = sourceData$matches,
                                  Pwm = sourceData$pwm,
                                  StructureTracks = imported_genomic_data[names(imported_genomic_data) %in% c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                                                                                           "X3UTR", "Downstream", "TTS", "TSS")],

                                  FeatureRanges_list = imported_genomic_data[!(names(imported_genomic_data) %in% c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                                                                                            "X3UTR", "Downstream", "TTS", "TSS"))],
                                  genome = sourceData$genome,
                                  TFNames = names(sourceData$pwm),
                                  pval_threshold,
                                  ChIP_regions = sourceData$ChIP,
                                  short_window,
                                  medium_window,
                                  long_window)
  return(DataSet)
}

#' Builds a predictive model of transcription factor binding site (TFBS) location
#' @export
#' @importFrom pROC roc auc
#' @importFrom stats model.matrix cor predict
#' @importFrom caret findCorrelation confusionMatrix
#' @importFrom xgboost xgb.DMatrix xgb.cv xgb.train xgb.importance xgb.plot.importance
#' @description
#' `buildTFBSmodel` learns (and optionally validates) a predictive model of transcription factor(TF) binding sites based
#' on the integration of genomic features extracted at the location of potential binding
#' sites identified by a prior analysis (such as pattern-matching). This function implements
#' an algorithm of supervised machine learning, the extreme gradient boosting (XGBoost), which
#' will be fed by a dataset of potential binding sites of training TFs either labelled as 'positive' or 'negative'
#' (i.e. ChIP-seq 'validated' or 'not-validated' in a given condition).
#' @param TFBSdata A named character vector as output by the [getTFBSdata()] function, defining the local paths
#' to files encoding for the results of pattern-matching and geonmic feature extraction for the training TFs and/or
#' studied TFs.
#' @param  ChIPpeaks A named character vector defining the local paths to BED files encoding the
#' location of ChIP-peaks. The vector is named according to the training transcription factors
#' that are described by the files indicated. **Caution**: the names of the \code{ChIPpeaks} have to find
#' a match with those of \code{TFBSdata}.
#' @param ChIPpeaks_length An integer setting a fixed length for the ChIP-peaks, that are defined as the intervals of
#' \code{ChIPpeaks_length} bp that are centered on the regions encoded in the \code{ChIPpeaks} files. Default velue = 400.
#' @param TFs_validation \code{NULL} (by default) or a character vector of names of training TFs. If \code{NULL}, the model performances will be
#' assessed from a validation dataset obtained by sampling randomly 20% of the potential binding sites of all the training TFs.
#' Otherwise the model will be validated with the data related to the training TFs named in the \code{TFs_validation} parameter and
#' trained using the remaining training TFs.
#' @param model_assessment A logical. If \code{TRUE} (by default), 1) the imortance of features is represented,
#' 2) the confusion matrix is printed and 3) the ROC and AUC are plotted.
#' @param xgb_modeling A logical. If \code{TRUE} (by default), the step of modeling will be achieved and a model will be
#' output. If \code{FALSE}, only the steps of balancing, labelling and splitting between a training and a validation
#' datasets will be carried on and the output will be a list of two `data.table` corresponding to the training and validation
#' datasets.
#' @return An object of \code{xgb.Booster} class if \code{xgb_modeling = TRUE}. A list of two `data.table` corresponding
#' each to the training and validation datasets otherwise.
#' @seealso [getTFBSdata()] for obtaining the datasets and [predictTFBS()] to predict transcription factor
#' binding site location
#' @examples
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("extdata/x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("extdata/cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("extdata/intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("extdata/x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex,
#'                                               tss = system.file("extdata/tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("extdata/tts_example.bed", package = "Wimtrap"))
#' TFBSdata.ex <- getTFBSdata(pfm = system.file("extdata/pfm_example.pfm", package = "Wimtrap"),
#'                            TFnames = c("PIF3", "TOC1"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("extdata/genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
#' TFBSmodel.ex <- buildTFBSmodel(TFBSdata = TFBSdata.ex,
#'                                ChIPpeaks = c(PIF3 = system.file("extdata/PIF3_example.bed", package = "Wimtrap"),
#'                                              TOC1 = system.file("extdata/TOC1_example.bed", package = "Wimtrap")),
#'                                TFs_validation = "PIF3")

buildTFBSmodel <- function(TFBSdata,
                           ChIPpeaks,
                           ChIPpeaks_length = 400,
                           TFs_validation = NULL,
                           model_assessment = TRUE,
                           xgb_modeling = TRUE){
  ChIP_regions <- listChIPRegions(ChIPpeaks, NULL, ChIPpeaks_length)
  DataSet <- data.frame()
  if (length(names(ChIPpeaks))==0 &length(TFBSdata) == 1) {names(ChIPpeaks) == names(TFBSdata)}
  for (trainingTF in names(ChIPpeaks)){
     considered <- data.table::fread(TFBSdata[trainingTF],
                                     stringsAsFactors = TRUE)
     considered$TF <- trainingTF
     considered <- GenomicRanges::makeGRangesFromDataFrame(considered,
                                                           keep.extra.columns = TRUE)
     validated_TFBS <- GenomicRanges::findOverlaps(considered, ChIP_regions[[trainingTF]], select = "all")
     considered <- as.data.frame(considered)
     considered$ChIP.peak <- 0
     considered$ChIP.peak[validated_TFBS@from[!duplicated(validated_TFBS@from)]] <- 1
     NbTrueBs <- nrow(considered[considered$ChIP.peak == 1,])
     DataSet <- rbind(DataSet,
                      considered[considered$ChIP.peak == 1,],
                      considered[sample(which(considered$ChIP.peak == 0), NbTrueBs),])
     rm(considered)
  }
  DataSet <- data.table::as.data.table(DataSet)
  TFBSdata <- DataSet[,-seq(1,5), with = FALSE]
  rm(DataSet)
  #Split the dataset into a training and a validation datasets
  if(is.null(TFs_validation)){
    trainind <- sample(seq(1,nrow(TFBSdata)), as.integer(nrow(TFBSdata)*0.8))
    testind <- seq(1,nrow(TFBSdata))[!(seq(1,nrow(TFBSdata)) %in% trainind)]
  } else {
    trainind <- which(TFBSdata$TF %in% TFs_validation)
    testind <- which(!(TFBSdata$TF) %in% TFs_validation)
  }
  TFBSdata.training <- TFBSdata[trainind,]
  TFBSdata.validation <- TFBSdata[testind,]
  #Pre-processing of the training dataset
  ##Remove columns that do not have to be taken into account
  if (length(grep(pattern = "matchScore", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "matchScore", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "TF", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "TF", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "ClosestTSS", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "ClosestTSS", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "ClosestTTS", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "ClosestTTS", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  ## Remove the infinite p-values associated to P-M,
  ## that occurs when the PWM is not flexible (i.e. is a consensus)
  TFBSdata.training$matchLogPval[which(is.infinite(TFBSdata.training$matchLogPval))] <- max(TFBSdata$matchLogPval)
  ##Create dummy variables
  dummy <- stats::model.matrix(~.+0, data = TFBSdata.training[,-c("ChIP.peak"),with=F])
  TFBSdata.training <- cbind(dummy, TFBSdata.training[,"ChIP.peak"])
  ##Remove highly correlated features
  descrCor <- stats::cor(TFBSdata.training, use = 'complete')
  highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .95)
  filteredDescr <- TFBSdata.training[,colnames(TFBSdata.training)[highlyCorDescr] := NULL]
  NAs <- is.na(as.data.frame(filteredDescr))
  NAs <- apply(NAs, 1, function(x) {if (length(which(x==TRUE)) > 0 ) {return(TRUE)} else {return(FALSE)} })
  train <- filteredDescr[!NAs,]
  #Pre-processing of the validation dataset
  ##Remove columns that do not have to be taken into account
  if (length(grep(pattern = "TF", colnames(TFBSdata.validation))) > 0){
    torm <- grep(pattern = "TF", colnames(TFBSdata.validation))
     TFBSdata.validation <- TFBSdata.validation[, colnames(TFBSdata.validation)[torm] := NULL]
   } else {}
  if (length(grep(pattern = "TSS", colnames(TFBSdata.validation))) > 0){
    torm <- grep(pattern = "TSS", colnames(TFBSdata.validation))
    TFBSdata.validation <- TFBSdata.validation[, colnames(TFBSdata.validation)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "TTS", colnames(TFBSdata.validation))) > 0){
    torm <- grep(pattern = "TTS", colnames(TFBSdata.validation))
    TFBSdata.validation <- TFBSdata.validation[, colnames(TFBSdata.validation)[torm] := NULL]
  } else {}
  ## Remove the infinite p-values associated to P-M,
  ## that occurs when the PWM is not flexible (i.e. is a consensus)
  TFBSdata.validation$matchLogPval[which(is.infinite(TFBSdata.validation$matchLogPval))] <- max(TFBSdata$matchLogPval)
  #Create dummy variables
  dummy <- stats::model.matrix(~.+0, data = TFBSdata.validation[,-c("ChIP.peak"),with=F])
  TFBSdata.validation <- cbind(dummy, TFBSdata.validation[,"ChIP.peak"])
  ##Remove highly correlated features
  descrCor <- stats::cor(TFBSdata.validation, use = 'complete')
  highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .95)
  filteredDescr <- TFBSdata.validation[,colnames(TFBSdata.validation)[highlyCorDescr] := NULL]
  NAs <- is.na(as.data.frame(filteredDescr))
  NAs <- apply(NAs, 1, function(x) {if (length(which(x==TRUE)) > 0 ) {return(TRUE)} else {return(FALSE)} })
  test <- filteredDescr[!NAs,]

  train <- train[,which(colnames(train) %in% colnames(test)), with = FALSE]
  test <- test[,which(colnames(test) %in% colnames(train)), with = FALSE]
  #Build a model by extreme gradient boosting
  labels <- train$ChIP.peak
  ts_label <- test$ChIP.peak

  new_tr <- stats::model.matrix(~.+0,data = train[,-c("ChIP.peak"),with=F])
  new_ts <- stats::model.matrix(~.+0,data = test[,-c("ChIP.peak"),with=F])

  if (xgb_modeling == TRUE) {
    rm(train)
    rm(filteredDescr)
    rm(descrCor)
    rm(highlyCorDescr)
    rm(TFBSdata.validation)
    rm(TFBSdata.training)
    rm(NAs)
    dtrain <- xgboost::xgb.DMatrix(data = new_tr,label = labels)
    dtest <- xgboost::xgb.DMatrix(data = new_ts,label=ts_label)

    params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)
    xgbcv <- xgboost::xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, print.every.n = 10, early.stop.round = 100, maximize = F, metrics = "auc")
    
    xgb1 <- xgboost::xgb.train( params = params, data = dtrain, watchlist = list(train = dtrain, eval = dtest), nrounds = 100, eval_metric = "auc")

    xgbpred <- stats::predict(xgb1,dtest)

    if (model_assessment){
      xgbpred <- ifelse(xgbpred > 0.5,1,0)
      cat("Performances of the model\n")
      print(caret::confusionMatrix(as.factor(xgbpred), as.factor(ts_label)))
      mat <- xgboost::xgb.importance(model = xgb1)
      cat("Features importance")
      print(mat)
      nb_max <- ifelse(nrow(mat) < 35, nrow(mat), 35)
      print(xgboost::xgb.plot.importance(importance_matrix = mat[1:nb_max]))
      print(plot(pROC::roc(ts_label, xgbpred)))
      print(lines(pROC::roc(ts_label, test$matchLogPval), col = "red"))
      print(legend(x = "bottom", legend = c("Model", "Pattern-Matching"), fill = c("black", "red")))
      print(mtext(text = round(pROC::auc(ts_label, xgbpred), 2), side = 2))
    } else {
    }
    save(xgb1, file = "model.RData")
    return(xgb1)
  } else {
    return(list(training.dataset = train, validation.dataset = test))
  }
}


#' @name predictTFBS
#' @title Predicts transcription factor binding sites
#' @export
#' @importFrom rtracklayer mcols import
#' @importFrom GenomicRanges strand makeGRangesFromDataFrame GRanges reduce as.data.frame start end resize findOverlaps nearest distanceToNearest
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings readDNAStringSet width Views matchPWM reverseComplement letterFrequency
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom utils read.csv
#' @importFrom biomartr getENSEMBLGENOMESInfo getENSEMBLInfo getGenome
#' @importFrom biomaRt listEnsembl listEnsemblGenomes useEnsemblGenomes searchDatasets getBM
#' @importFrom methods as
#' @importFrom IRanges IRanges resize
#' @importFrom readr write_tsv
#' @importFrom data.table fread data.table
#' @importFrom graphics hist
#' @importFrom stats runif
#' @description
#' Predicts the binding sites of studied transcription factors in a studied condition
#' by applying a predictive model on potential binding sites
#' located by a prior pattern-matching analysis and annotated with genomic features extracted at their location.
#' @param TFBSmodel A \code{xgb.Booster} object as output by the function [buildTFBSmodel()].
#' @param TFBSdata A named character vector as output by the [getTFBSdata()] function, defining the local paths
#' to files encoding for the results of pattern-matching and genomic feature extraction for the training TFs and/or
#' studied TFs.
#' @param studiedTFs A character vector setting the names of the studied transcription factors (for which the
#' location of the binding sites has to be predicted). This names have to match with one of the names of the
#' \code{TFBSdata} argument.
#' @param show_annotations A logical. Default = \code{FALSE}. If \code{TRUE}, the annotation of the potential binding sites
#' with the genomic features extracted from their genomic regions will be output.
#' @param score_threshold A numeric (comprised between 0 and 1). Sets the minimum prediction score output by the
#' \code{TFBSmodel} to predict a potential binding site as a binding site of a studied transcription factor in the
#' studied condition. Higher the prediction score, higher is the specificity and lower the sensitivity. Default = 0.5.
#' @details  Remark: the features included in the datasets encoded in the files given by \code{TFBSdata} have to correspond
#' to those integrated by the predictive model set by \code{TFBSmodel}.
#' @seealso [importGenomicData()] for importing genomic data, [buildTFBSmodel()] to train a predictive model of
#' transcription factor binding sites, and [plotPredictions()] to vizualize the results for a given potential target gene.
#' @return A \code{data.table} listing the predicted binding sites. The '\code{TF}' column annotates the potential binding sites
#' with their cognate transcription factor. Additionally, the \code{data.table} describes, for the potential
#' binding sites, the chromosomic coordinates, the closest transcript (relatively to the transcript start site) and the prediction score.
#' Optionally, the \code{data.table} might also include the genomic features used to make the predictions.
#' @examples
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("extdata/x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("extdata/cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("extdata/intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("extdata/x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex,
#'                                               tss = system.file("extdata/tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("extdata/tts_example.bed", package = "Wimtrap"))
#' TFBSdata.ex <- getTFBSdata(pfm = system.file("extdata/pfm_example.pfm", package = "Wimtrap"),
#'                            TFnames = c("PIF3", "TOC1"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("extdata/genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
#' TFBSmodel.ex <- buildTFBSmodel(TFBSdata = TFBSdata.ex,
#'                                ChIPpeaks = c(PIF3 = system.file("extdata/PIF3_example.bed", package = "Wimtrap"),
#'                                              TOC1 = system.file("extdata/TOC1_example.bed", package = "Wimtrap")),
#'                                TFs_validation = "PIF3")
#' PIF3BS.predictions <- predictTFBS(TFBSmodel.ex,
#'                                   TFBSdata.ex,
#'                                   studiedTFs = "PIF3")
#' ##To get the transcripts whose expression is potentially regulated by PIF3 do as follows:
#' PIF3_regulated.predictions <- as.character(PIF3BS.predictions$transcript[!duplicated(PIF3BS.predictions)])
#' ###If you want to consider only the gene model,
#' ###then do as follows:
#' PIF3_regulated.predictions <- unlist(strsplit(PIF3_regulated.predictions, "[.]"))[seq(1, 2*length(PIF3_regulated.predictions),2)]
#' PIF3_regulated.predictions <- PIF3_regulated.predictions[!duplicated(PIF3_regulated.predictions)]

predictTFBS <- function(TFBSmodel,
                        TFBSdata,
                        studiedTFs = NULL,
                        show_annotations = FALSE,
                        score_threshold = 0.86){
  results <- list()
  paths <- TFBSdata
  if (length(studiedTFs)==0){studiedTFs <- names(TFBSdata)}
  for (TF in studiedTFs){
    DataSet <- data.table::fread(paths[TF],
                                  stringsAsFactors = TRUE)

    TFBSdata <- DataSet[,-seq(1,5), with = FALSE]
    #Pre-processing of the  dataset
    ##Remove columns that do not have to be taken into account
    if (length(grep(pattern = "matchScore", colnames(TFBSdata))) > 0){
      torm <- grep(pattern = "matchScore", colnames(TFBSdata))
      TFBSdata <- TFBSdata[, colnames(TFBSdata)[torm] := NULL]
    } else {}
    if (length(grep(pattern = "TF", colnames(TFBSdata))) > 0){
      torm <- grep(pattern = "TF", colnames(TFBSdata))
      TFBSdata <- TFBSdata[, colnames(TFBSdata)[torm] := NULL]
    } else {}
    if (length(grep(pattern = "ClosestTSS", colnames(TFBSdata))) > 0){
      torm <- grep(pattern = "ClosestTSS", colnames(TFBSdata))
      TFBSdata <- TFBSdata[, colnames(TFBSdata)[torm] := NULL]
    } else {}
    if (length(grep(pattern = "ClosestTTS", colnames(TFBSdata))) > 0){
      torm <- grep(pattern = "ClosestTTS", colnames(TFBSdata))
      TFBSdata <- TFBSdata[, colnames(TFBSdata)[torm] := NULL]
    } else {}
    ## Remove the infinite p-values associated to P-M,
    ## that occurs when the PWM is not flexible (i.e. is a consensus)
    TFBSdata$matchLogPval[which(is.infinite(TFBSdata$matchLogPval))] <- max(TFBSdata$matchLogPval[which(!is.infinite(TFBSdata$matchLogPval))])
    TFBSdata <- TFBSdata[,TFBSmodel$feature_names, with=FALSE]
    ##Create dummy variables
    TFBSdata <- stats::model.matrix(~.+0, data = TFBSdata)
    NAs <- is.na(as.data.frame(TFBSdata))
    NAs <- apply(NAs, 1, function(x) {if (length(which(x==TRUE)) > 0 ) {return(TRUE)} else {return(FALSE)} })
    TFBSdata <- TFBSdata[!NAs,]
    TFBSdata <- TFBSdata[,colnames(TFBSdata)[colnames(TFBSdata) %in% TFBSmodel$feature_names]]
    TFBSdata <- data.table::as.data.table(TFBSdata)
    TFBSdata <- stats::model.matrix(~.+0,data = TFBSdata)

    TFBSdata <- xgboost::xgb.DMatrix(data = TFBSdata)
    xgbpred <- stats::predict(TFBSmodel,TFBSdata)
    if (show_annotations){
      results[[TF]] <- cbind(DataSet, prediction.score = xgbpred)
      colnames(results[[TF]])[which(colnames(results[[TF]]) == "ClosestTSS")] <- "transcript"
    } else {
      results[[TF]] <- cbind(DataSet[,c(seq(1,5), which(colnames(DataSet)=="ClosestTSS")),with=FALSE], prediction.score = xgbpred)
      colnames(results[[TF]])[which(colnames(results[[TF]]) == "ClosestTSS")] <- "transcript"
  }
  }
  results <- lapply(seq_along(results), function(i){
    x <- cbind(results[[i]], TF = names(results)[i]);
    return(x)
  })
  results <- do.call(rbind, results)
  results <- results[results$prediction.score >= score_threshold,]
  return(results)
}

#' Plots genomic data and predicted binding sites along a gene track
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame mcols strand seqnames findOverlaps sort reduce start end
#' @importFrom Gviz GeneRegionTrack GenomeAxisTrack DataTrack plotTracks
#' @importFrom grDevices rainbow topo.colors
#' @importFrom GenomeInfoDb keepSeqlevels
#' @importFrom grid grid.newpage pushViewport viewport grid.layout popViewport
#' @description
#' Allows the vizualisation of the predicted binding sites on a given gene for the studied condition.
#' The function plots three panels of gene tracks: (i) the different transcript models of the gene of interest;
#' (ii) the predicted binding sites, annotated with their prediction score; and (iii), optionally, the signals
#' of the genomic data that have been used to extract the predictive features at the location of the potential
#' binding sites.
#' @param TFBSprediction A \code{data.table} as output by [predictTFBS()]. This \code{data.table} comprises the columns '\code{seqnames}',
#' '\code{start}', '\code{end}', '\code{width}', '\code{strand}', '\code{transcript}', '\code{prediction.score}', '\code{TF}'
#' and, optionally, the annotations with the predictive features extracted from the genomic data.
#' @param imported_genomic_data A list of \code{GRanges} objects as output by [importGenomicData()]. It includes at least
#' the three \code{GRanges} named '\code{CDS}', '\code{X5UTR}' and '\code{X3UTR}' that give the position of the coding sequence(s), 5'UTR and
#' 3'UTR of the transcript models that are indicated in the '\code{name}' metadata column.
#' @param gene A character. Sets the name of the gene to be plotted. The system of nomenclature corresponds to that
#' used to name the transcript models in the \code{GRanges} objects of \code{imported_genomic_data}.
#' @param TFs \code{NULL} or a vector of characters indicating the transcription factors for which the predicted binding
#' sites have to be plotted. By default, all the studied transcription factors will be considered.
#' @param genomic_data \code{NULL} or a vector of characters allowing to name the genomic data for which the signal has to
#' be plotted. The names that are accepted correspond to those of the GRanges objects listed in \code{imported_genomic_data}.
#' By default, no genomic data is plotted.
#' @param promoter_length A numeric allowing to set the start of the gene track (=lowest start among the transcript
#' models of the considered gene - \code{promoter_length}). Default = 2000.
#' @param downstream_length A numeric allowing to set the end of the gene track (=highest among the transcript
#' models of the considered gene - \code{downstreamr_length}). Default = 2000.
#' @return A list of \code{GenomeGraph} tracks.
#' @seealso [predictTFBS()] to get predictions of transcription factor binding sites.
#' @examples
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("extdata/x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("extdata/cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("extdata/intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("extdata/x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex,
#'                                               tss = system.file("extdata/tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("extdata/tts_example.bed", package = "Wimtrap"))
#' TFBSdata.ex <- getTFBSdata(pfm = system.file("extdata/pfm_example.pfm", package = "Wimtrap"),
#'                            TFnames = c("PIF3", "TOC1"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("extdata/genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
#' TFBSmodel.ex <- buildTFBSmodel(TFBSdata = TFBSdata.ex,
#'                                ChIPpeaks = c(PIF3 = system.file("extdata/PIF3_example.bed", package = "Wimtrap"),
#'                                              TOC1 = system.file("extdata/TOC1_example.bed", package = "Wimtrap")),
#'                                TFs_validation = "PIF3")
#' PIF3BS.predictions <- predictTFBS(TFBSmodel.ex,
#'                                   TFBSdata.ex,
#'                                   studiedTFs = "PIF3")
#' plotPredictions(PIF3BS.predictions,
#'                 imported_genomic_data.ex,
#'                 gene = "AT1G01140",
#'                 genomic_data = "DHS")

plotPredictions <- function(TFBSpredictions,
                            imported_genomic_data,
                            gene,
                            TFs = NULL,
                            genomic_data = NULL,
                            promoter_length = 2000,
                            downstream_length = 2000){
  TFBSpredictions <- TFBSpredictions[grep(pattern = gene, TFBSpredictions$transcript),]
  TFBSpredictions$seqnames <- droplevels(as.factor(TFBSpredictions$seqnames))
  TFBSpredictions <- GenomicRanges::makeGRangesFromDataFrame(TFBSpredictions,
                                                             keep.extra.columns = TRUE)

  structures <- imported_genomic_data[names(imported_genomic_data) %in% c("X5UTR", "CDS", "Intron",
                                                                          "X3UTR","TTS", "TSS")]
  structures <- lapply(seq_along(structures), function(i){x <- as.data.frame(structures[[i]]);
  x <-x[grep(pattern = gene, x$name),];
  x <- cbind(x[,seq(1,5)],
             data.frame(feature=names(structures)[i],
                        gene = gene,
                        exon =x$name,
                        transcript = x$name,
                        symbol = gene))
  return(x)})
  structures <- do.call(rbind, structures)
  models <- structures[structures$feature %in% c("CDS", "X5UTR", "X3UTR"),]
  models$feature <- droplevels(as.factor(models$feature))
  levels(models$feature) <- c("utr5", "cds", "utr3")
  models$seqnames <- droplevels(as.factor(models$seqnames))
  grtrack <- Gviz::GeneRegionTrack(range = GenomicRanges::makeGRangesFromDataFrame(models),
                                   symbol = as.character(models$transcript),
                                   transcript =  as.character(models$transcript),
                                   name="Gene Model", background.title="brown",  showId = TRUE, geneSymbol = TRUE)
  gtrack <- Gviz::GenomeAxisTrack()
  strack <- list()
  if (length(TFs)==0){
    TFs <- as.factor(GenomicRanges::mcols(TFBSpredictions)$TF)
  } else {
    TFs <- TFs[TFs %in% as.factor(GenomicRanges::mcols(TFBSpredictions)$TF)]
  }
  col <- grDevices::rainbow(length(TFs))
  names(col) <- TFs
  GenomicRanges::strand(TFBSpredictions) <- "*"
  strack <- list()
  for (TF in TFs){
    strack[[TF]] <- Gviz::DataTrack(TFBSpredictions[GenomicRanges::mcols(TFBSpredictions)$TF == TF],
                                    name = TF, col = col[TF], background.title = col[TF])
  }
  if (length(genomic_data) > 0){
    col <- grDevices::topo.colors(length(genomic_data))
    names(col) <- genomic_data
    ctrack <- list()
    for (feature in genomic_data){
      extracted <- imported_genomic_data[[feature]][GenomicRanges::seqnames(imported_genomic_data[[feature]])
                                                    == as.character(models$seqnames[1])]
      extracted <- GenomeInfoDb::keepSeqlevels(extracted,
                                               as.character(models$seqnames[1]))
      overlapping <- GenomicRanges::findOverlaps(extracted, GenomicRanges::makeGRangesFromDataFrame(data.frame(seqname = as.character(models$seqnames[1]),
                                                                                                               start = min(models$start)-promoter_length,
                                                                                                               end = max(models$end)+downstream_length,
                                                                                                               strand = "*")))
      extracted <- extracted[overlapping@from[!duplicated(overlapping@from)]]
      if (length(extracted) > 0){
        extracted <-GenomicRanges::sort(extracted, ignore.strand = TRUE)
        GenomicRanges::strand(extracted) <- "*"
        if(is.factor(GenomicRanges::mcols(extracted)[,1]) | is.character(GenomicRanges::mcols(extracted)[,1])){
          GenomicRanges::mcols(extracted)[,1] <- 1
        }
        intermediate <- GenomicRanges::reduce(extracted, ignore.strand=TRUE)
        complements <- intermediate[-length(intermediate)]
        GenomicRanges::mcols(complements)[,1] <- 0
        names(GenomicRanges::mcols(complements)) <- names(GenomicRanges::mcols(extracted))
        for (i in seq_along(complements)){
          GenomicRanges::start(complements[i]) <- GenomicRanges::end(complements[i])
          GenomicRanges::end(complements[i]) <- GenomicRanges::start(intermediate[i+1])
        }
        complements <- c(extracted[1], complements)
        GenomicRanges::end(complements)[1] <- GenomicRanges::start(complements)[1]
        GenomicRanges::start(complements)[1] <- min(models$start)-promoter_length
        complements <- c(complements, extracted[length(complements)])
        GenomicRanges::start(complements)[length(complements)] <- GenomicRanges::end(complements[length(complements)])
        GenomicRanges::end(complements[length(complements)]) <- max(models$end)+downstream_length
        GenomicRanges::mcols(complements)[,1] <- 0
        extracted <- c(extracted, complements)
        extracted <- GenomicRanges::sort(extracted)
        tmp_extracted <- extracted
        GenomicRanges::end(tmp_extracted) <- GenomicRanges::start(tmp_extracted)
        GenomicRanges::start(extracted) <- GenomicRanges::end(extracted)
        extracted <- c(extracted, tmp_extracted)
        extracted <- GenomicRanges::sort(extracted)
        ctrack[[feature]] <- Gviz::DataTrack(extracted, col = col[feature], name = feature)
      } else {
        extracted <- GenomicRanges::makeGRangesFromDataFrame(
          data.frame(seqname = models$seqnames[1],
                     start = c(min(models$start)-promoter_length, max(models$end)+downstream_length),
                     end = c(min(models$start)-promoter_length, max(models$end)+downstream_length),
                     feature = 0
          ), keep.extra.columns = TRUE)
        ctrack[[feature]] <- Gviz::DataTrack(extracted, col = col[feature], name = feature)
      }
    }
    ncols <- 1
    nrows <- 2
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrows, ncols)))
    grid::pushViewport(grid::viewport(layout.pos.col=((1-1)%%ncols)+1, layout.pos.row=(((1)-1)%/%ncols)+1))
    Gviz::plotTracks(c(list(gtrack, grtrack), strack), from = min(models$start)-promoter_length, to = max(models$end)+downstream_length, type = "p", add=TRUE)
    grid::popViewport(1)
    grid::pushViewport(grid::viewport(layout.pos.col=1, layout.pos.row=2))
    Gviz::plotTracks(ctrack, from = min(models$start)-promoter_length, to = max(models$end)+downstream_length, type = "s", add=TRUE)
    grid::popViewport(1)
  } else {
    Gviz::plotTracks(c(list(gtrack, grtrack), strack), from = min(models$start)-promoter_length, to = max(models$end)+downstream_length, type ="p")
  }
}

#' Cis-acting regulatory element predictions for Arabidopsis and tomato
#' @export
#' @importFrom utils download.file unzip
#' @description
#' Predicts the location of transcription factor binding sites (=cis-acting regulatory elements) in various conditions
#' for *Arabidopsis thaliana* and *Solanum lycopersicum*. The function integrates 5 pre-built general models obtained based on a more
#' or less extended set of genomic data and trained from different organisms/conditions. These models almost all integrate the degree of opening
#' of the chromating (DHS: DNAseI hypersensitive sites) and results of digital genomic footprinting (DGF: digital genomic footprints) in
#' the conditions that can be studied using \code{carepat}. These represent genomic data with high potenital of prodectivity (see details).
#' @param organism "Arabidopsis thaliana" or "Solanum lycopersicum"
#' @param condition Character indicating the studied condition. For *Arabidopsis thaliana*: "seedlings", "flowers", "roots",
#' "roots_non_hairs","seed_coats", "seedlings_dark7d", "seedlings_dark7dLD24h","seedlings_dark7dlight3h",
#'  "seedlings_dark7dlight30min" or "seedlings_heatshock"". For *Solanum lycopersicum*: "ripening_fruits" or "immaturefruits".
#' @param TFnames Character vector setting the name(s) of the studied transcription factors. These names have to follow the
#' AGI (Arabidopsis)/Solyc (Tomato) nomenclature to allow the retrieval of the motis from [PlantTFDB](http://planttfdb.gao-lab.org/)
#' database.Otherwise, if you input the motifs from a local file through \code{pfm}, these names have to be among those described in that file.
#' @param  pfm Path to a file including the position frequency or weight matrices (PFMs or PWMs) of the motifs recognized
#' by the considered transcription factors (training and/or studied TFs). This file can be in different formats, determined based
#' on the file extension: raw pfm (".pfm"), jaspar (".jaspar"), meme (".meme"), transfac (".transfac"), homer (".motif") or
#' cis-bp (".txt").\code{pfm} can be set to \code{NULL} (default value) if you provide the results of pattern-matching
#' obtained from an external source (see the argument \code{matches}).
#' @param show_annotations A logical. Default = \code{FALSE}. If \code{TRUE}, the annotation of the potential binding sites
#' with the genomic features extracted from their genomic regions will be output.
#' @param score_threshold A numeric (comprised between 0 and 1). Sets the minimum prediction score output by the
#' \code{TFBSmodel} to predict a potential binding site as a binding site of a studied transcription factor in the
#' studied condition. Higher the prediction score, higher is the specificity and lower the sensitivity. Default = 0.5.
#' @details  The following table details, for each organism-condition that can be studied using \code{carepat}, the model
#' that is considered: from which training organism-condition it has been obtained and which genomic features.
#'
#' |studied  = training organism | studied condition | training condition | genomic features |
#' |:----------------------------|:-----------------:|:------------------:|:-----------------------:|
#' | *Arabidopsis thaliana* | whole seedlings ("seedlings") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + Full layer 5
#' | *Arabidopsis thaliana* | flowers in stages 4-5 ("flowers") | flowers in stages 4-5 ("flowers")| Layers 1, 2, 3, 4 + DHS, Cme
#' | *Arabidopsis thaliana* | seedling roots ("roots") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | non-hair part of seedling roots ("roots_non_hair") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | seed coats, 4 days after anthesis ("seedlings_coats") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | heat-shocked seedlings ("seedlings_heatshock) | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | dark-grown seedlings ("seedlings_dark7d") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | dark-grown seedlings exposed to 30 min of light ("seedlings_dark7d30min") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | dark-grown seedlings exposed to 3h of light ("seedlings_dark7d3h") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Arabidopsis thaliana* | dark-grown seedlings exposed to a long day cycle ("seedlings_dark7dLD24h") | whole seedlings ("seedlings") | Layers 1, 2, 3, 4 + DHS
#' | *Solanum lycopersicum* | ripening fruits ("ripening_fruits") | ripening fruits ("ripening_fruits") | Layers 1, 2, 3, 4 + DHS, Cme, H3K27me3
#' | *Solanum lycopersicum* | immature fruits ("ripening_fruits") | ripening fruits ("ripening_fruits") | Layers 1, 2, 3 + DHS, Cme, H3K27me3
#'
#' The different layers of genomic features are composed of the following:
#' - **Layer 1**: results of pattern-matching (log10 p-value of the score and local density of matches)
#' - **Layer 2**: phastcons-scored conserved elements (for Arabidopsis and the tomato) and conserved non-coding sequences
#' (for Arabiopsis only)
#' - **Layer 3**: position on the gene (promoter, proximal promoter, 5'untranslated region, coding sequence, intron, 3'untranslated
#' region, downstream region, distance to the transcription start and termination site of the gene)
#' - **Layer 4**: local signals of digital footprints
#' - **Layer 5**: local signals of Chromatin state, Cytosine methylation (Cme), Histone 2A.Z positioning (H2AZ), DNA looping (Dloop),
#' Nucleosomes positioning (Nuc), Histone2B monoubiquitination (H2BuB), Monomethylation on lysine 4 of the histone 3 (H3K4me1),
#' Dimethylation on lysine 4 of the histone 3 (H3K4me2), Trimethylation on lysine 4 of the histone 4 (H3K4me3),
#' Dimethylation on lysine 9 of the histone 3 (H3K9me2), Monomethylation on lysine 27 of the histone 3 (H3K27me1),
#' Trimethylation on lysine 27 of the histone 3 (H3K27me3), Trimethylation on lysine 36 of the histone 4 (H3K36me3),
#' Acetylation on lysine 9 of histone 3 (H3K9ac), Acetylation on lysine 14 of histone 3 (H3K14ac), Acetylation on lysine 18 of histone 3 (H3K18ac),
#' Acetylation on lysine 27 of histone 3 (H3K27ac), Acetylation on lysine 56 of histone 3 (H3K56ac),
#' Phosphorylation on tyrosine 3 of histone 3 (H3T3ph), Acetylation on lysine 5 of histone 4 (H4K5ac),
#' Acetylation on lysine 8 of histone 4 (H4K8ac), Acetylation on lysine 12 of histone 4 (H4K12ac),
#' Acetylation on lysine 16 of histone 4 (H4K16ac).
#'
#' The **source** of the data used to train the models and, if applicable, to transfer them the studied conditions are described in the file "Sources.ods"
#' on the "RivereQuentin/carepat" repository.
#' @seealso [plotPredictions()] to vizualize the results for a given potential target gene.
#' @return A \code{data.table} listing the predicted binding sites. The '\code{TF}' column annotates the potential binding sites
#' with their cognate transcription factor. Additionally, the \code{data.table} describes, for the potential
#' binding sites, the chromosomic coordinates, the closest transcript (relatively to the transcript start site) and the prediction score.
#' Optionally, the \code{data.table} might also include the genomic features used to make the predictions.
#' @examples
#'#Predictions of the binding sites of "AT2G46830" in flowers of Arabidopsis
#'CCA1predictions.flowers <- carepat(organism = "Arabidopsis thaliana",
#'                                   condition = "flowers",
#'                                   TFnames = "AT2G46830")
#'#Predictions of the binding sites of "Solyc00g024680.1" in immature fruits of tomato
#'DOF24predictions.immature <- carepat(organism = "Solanum lycopersicum",
#'                                   condition = "immature_fruits",
#'                                   TFnames = "Solyc00g024680.1")
#'

carepat <- function(organism = c("Arabidopsis thaliana", "Solanum lycopersicum"),
                    condition = c("seedlings", "flowers", "roots", "roots_non_hairs",
                                  "seed_coats", "seedlings_dark7d", "seedlings_dark7dLD24h",
                                  "seedlings_dark7dlight3h", "seedlings_dark7dlight30min", "seedlings_heatshock",
                                  "ripening_fruits", "immature_fruits"),
                    TFnames = NULL,
                    pfm = NULL,
                    show_annotations = FALSE,
                    score_threshold = 0.5){

  #Download the required data and models from the github
  #repository RiviereQuentin/caretrap if necessary
  package.dir <- system.file(package = "Wimtrap")
  test.file = paste0(package.dir, "/carepat-main/data/Arabidopsis_thaliana/PWMs_athal.meme")
  if (!file.exists(test.file)){
    message("Downloading data and models - needs to be done once")
    utils::download.file(url = "https://github.com/RiviereQuentin/carepat/archive/main.zip",
                         destfile = paste0(package.dir, "/carepat.zip" ),
                         timeout = 600,
                         quiet = FALSE)
    utils::unzip(zipfile = paste0(package.dir, "/carepat.zip"),
                 exdir = package.dir)
  }
  dir.models <- paste0(package.dir, "/carepat-main/models/")
  if (organism == "Solanum lycopersicum"){
    dir.data <- paste0(package.dir, "/carepat-main/data/Solanum_lycopersicum/")
    if (condition == "ripening_fruits"){
      genomic_data <- c("CE_sl.bed", "CDS_sl.bed", "Intron_sl.bed", "X5UTR_sl.bed",
                               "X3UTR_sl.bed", "ripeningfruits/DHS_sl_ripening.bed",
                               "ripeningfruits/DGF_sl_ripening.bed", "ripeningfruits/H3K27me3_sl_ripening.bed",
                               paste0("ripeningfruits/Methylome_sl_ripening_part", seq(1,14), ".bed"))
      names(genomic_data) <- c("phastcons", "CDS", "Intron", "X5UTR", "X3UTR", "DHS", "DGF", "H3K27me3", paste0("Methylome_", seq(1,14)))
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      TFBSmodel <- paste0(dir.models, "TFBSmodel_sl_ripening_full.RData")
      load(TFBSmodel)
      TFBSmodel <- TFBSmodel.sl.ripening_full
    } else if (condition == "immature_fruits"){
      genomic_data <- c("CE_sl.bed", "CDS_sl.bed", "Intron_sl.bed", "X5UTR_sl.bed",
                               "X3UTR_sl.bed", "immaturefruits/DHS_sl_immature.bed",
                               "immaturefruits/H3K27me3_sl_immature.bed",
                               paste0("ripeningfruits/Methylome_sl_ripening_part", seq(1,14), ".bed"))
      names(genomic_data) <- c("phastcons", "CDS", "Intron", "X5UTR", "X3UTR", "DHS", "H3K27me3", paste0("Methylome_", seq(1,14)))
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      TFBSmodel <- paste0(dir.models, "TFBSmodel_sl_ripening_reduced.RData")
      load(TFBSmodel)
      TFBSmodel <- TFBSmodel.sl.ripening_reduced
    }

    imported_genomic_data <- Wimtrap::importGenomicData(genomic_data = genomic_data,
                                                        tts = paste0(dir.data, "TTS_sl.bed"),
                                                        tss = paste0(dir.data, "TSS_sl.bed"))
    imported_genomic_data$Methylome_1 <- imported_genomic_data[grep(pattern = "Methylome", names(imported_genomic_data))]
    imported_genomic_data$Methylome_1 <- lapply(imported_genomic_data$Methylome_1, function(x){x <- data.table::as.data.table(x); colnames(x)[ncol(x)] <- "CpG"; return(x)})
    imported_genomic_data$Methylome_1 <- do.call(rbind, imported_genomic_data$Methylome_1)
    imported_genomic_data$Methylome_1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data$Methylome_1)
    names(imported_genomic_data)[which(names(imported_genomic_data) == "Methylome_1")] <- "Cme"
    imported_genomic_data <- imported_genomic_data[!(seq(1, length(imported_genomic_data)) %in% grep(pattern = "Methylome", names(imported_genomic_data)))]
    genome_sequence <- Biostrings::readDNAStringSet(paste0(dir.data, c("genome_sl_chr0.fa.gz",
                                                                       "genome_sl_chr1_part1.fa.gz",
                                                                       "genome_sl_chr1_part2.fa.gz",
                                                                        paste0("genome_sl_chr", seq(2,12),".fa.gz"))))
    genome_sequence[2] <- Biostrings::xscat(genome_sequence[2], genome_sequence[3])
    genome_sequence <- genome_sequence[c(1,2,seq(4,14))]
    PWMs.file <- paste0(dir.data, "PWMs_sl.meme")

  } else if (organism == "Arabidopsis thaliana"){
    dir.data <- paste0(package.dir, "/carepat-main/data/Arabidopsis_thaliana/")
    if (condition == "seedlings"){
      genomic_data <- c(DHS = "seedlings/DHS_athal_seedlings.bed",
                        DGF = "seedlings/DGF_athal_seedlings.bed",
                        Dloops = "seedlings/Dloops_athal_seedlings.bed",
                        H2AZ = "seedlings/H2AZ_athal_seedlings.bed",
                        H2BuB = "seedlings/H2BUB_athal_seedlings.bed",
                        H3K4me1 = "seedlings/H3K4me1_athal_seedlings.bed",
                        H3K4me2 = "seedlings/H3K4me2_athal_seedlings.bed",
                        H3K9me2 = "seedlings/H3K9me2_athal_seedlings.bed",
                        H3K14ac = "seedlings/H3K14ac_athal_seedlings.bed",
                        H3K18ac = "seedlings/H3K18ac_athal_seedlings.bed",
                        H3K27ac = "seedlings/H3K27ac_athal_seedlings.bed",
                        H3K27me1 = "seedlings/H3K27me1_athal_seedlings.bed",
                        H3K56ac = "seedlings/H3K56ac_athal_seedlings.bed",
                        H4K5ac = "seedlings/H4K5ac_athal_seedlings.bed",
                        H4K8ac = "seedlings/H4K8ac_athal_seedlings.bed",
                        H4K12ac = "seedlings/H4K12ac_athal_seedlings.bed",
                        H4K16ac = "seedlings/H4K16ac_athal_seedlings.bed",
                        Methylome = "seedlings/Methylome_athal_seedlings.bed",
                        phast1 = "phastcons_athal_part1.bed",
                        phast2 = "phastcons_athal_part2.bed",
                        CDS = "CDS_athal.bed",
                        Intron = "Intron_athal.bed",
                        X3UTR = "X3UTR_athal.bed",
                        X5UTR = "X5UTR_athal.bed",
                        CNS = "CNS_athal.bed")
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      load(file = paste0(dir.models, "TFBSmodel_athal_seedlings_full.RData"))
      TFBSmodel <- TFBSmodel.athal.seedlings_full
    } else if (condition == "flowers"){
      genomic_data <- c(DHS = "flowers/DHS_athal_flowers.bed",
                        DGF = "flowers/DGF_athal_flowers.bed",
                        Methylome = "flowers/Methylome_athal_flowers.bed",
                        phast1 = "phastcons_athal_part1.bed",
                        phast2 = "phastcons_athal_part2.bed",
                        CDS = "CDS_athal.bed",
                        Intron = "Intron_athal.bed",
                        X3UTR = "X3UTR_athal.bed",
                        X5UTR = "X5UTR_athal.bed",
                        CNS = "CNS_athal.bed"
                      )
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      load(file = paste0(dir.models, "TFBSmodel_athal_flowers.RData"))
      TFBSmodel <- TFBSmodel.athal.flowers
    } else if (condition %in% c("roots", "roots_non_hairs", "seed_coats", "seedlings_dark7d", "seedlings_dark7dLD24h",
                                "seedlings_dark7dlight3h", "seedlings_dark7dlight30min", "seedlings_heatshock")){
      load(file = paste0(dir.models, "TFBSmodel_athal_seedlings_reduced.RData"))
      TFBSmodel <- TFBSmodel.athal.seedlings_reduced
      genomic_data <- c(phast1 = "phastcons_athal_part1.bed",
                                         phast2 = "phastcons_athal_part2.bed",
                                         CDS = "CDS_athal.bed",
                                         Intron = "Intron_athal.bed",
                                         X3UTR = "X3UTR_athal.bed",
                                         X5UTR = "X5UTR_athal.bed",
                                         CNS = "CNS_athal.bed")
      if (condition == "roots"){
        genomic_data <- c(genomic_data,
                          c(DGF = "roots/DGF_athal_roots.bed",
                                   DHS = "roots/DHS_athal_roots.bed"))
      } else if (condition == "roots_non_hairs"){
        genomic_data <- c(genomic_data,
                          c(DGF = "roots/DGF_athal_roots_non_hairs.bed",
                                   DHS = "roots/DHS_athal_roots_non_hairs.bed"))
      } else if (condition == "seed_coats"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seed_coats/DGF_athal_seed_coats.bed",
                                   DHS = "seed_coats/DHS_athaliana_seed_coats.bed"))
      } else if (condition == "seedlings_dark7d"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seedlings/DGF_athal_seedlings_dark7d.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7d.bed"))
      } else if (condition == "seedlings_dark7dLD24h"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seedlings/DGF_athal_seedlings_dark7dLD24h.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7dLD24h.bed"))
      } else if (condition == "seedlings_dark7dlight3h"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seedlings/DGF_athal_seedlings_dark7dlight3h.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7dlight3h.bed"))
      } else if (condition == "seedlings_dark7dlight30min"){
        genomic_data <- c(genomic_data,
                         c(DGF = "seedlings/DGF_athal_seedlings_dark7dlight30min.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7dlight30min.bed"))
      } else if (condition == "seedlings_heatshock"){
        genomic_data <- c(genomic_data,
                         c(DGF = "seedlings/DGF_athal_seedlings_heatshock.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_heatshock.bed"))
      }
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp

    }

    imported_genomic_data <- importGenomicData(genomic_data = genomic_data,
                                               tts = paste0(dir.data, "TTS_athal.bed"),
                                               tss = paste0(dir.data, "TSS_athal.bed"))

    imported_genomic_data$phast1 <- imported_genomic_data[grep(pattern = "phast", names(imported_genomic_data))]
    imported_genomic_data$phast1 <- lapply(imported_genomic_data$phast1,
                                           function(x){x <- as.data.frame(x);
                                           colnames(x)[ncol(x)] <- "phastcons";
                                           return(x)})
    imported_genomic_data$phast1 <- do.call(rbind, imported_genomic_data$phast1)
    imported_genomic_data$phast1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data$phast1, keep.extra.columns = TRUE)
    imported_genomic_data <- imported_genomic_data[!(names(imported_genomic_data)=="phast2")]
    names(imported_genomic_data)[which(names(imported_genomic_data)=="phast1")] <- "phastcons"
    genome_sequence <- Biostrings::readDNAStringSet(paste0(dir.data, paste0("genome_athal_chr", seq(1,5), ".fa.gzip")))
    PWMs.file <- paste0(dir.data, "PWMs_athal.meme")

  }

  if(length(pfm) == 1){
    PWMs.file <- pfm
  }
  TFBSdata <- getTFBSdata(pfm = PWMs.file,
                          TFnames = TFnames,
                          genome_sequence = genome_sequence,
                          imported_genomic_data = imported_genomic_data)
  TFBSpredictions <- predictTFBS(TFBSmodel = TFBSmodel,
                                 TFBSdata = TFBSdata,
                                 studiedTFs = TFnames,
                                 score_threshold = score_threshold,
                                 show_annotations = show_annotations)
  return(TFBSpredictions)
}

#' Test the performances of predicting gene targets based on the location of potential TFBS identified by Wimtrap
#' @export
#' @importFrom rtracklayer mcols
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom xgboost xgb.DMatrix
#' @importFrom stats predict
#' @importFrom pROC roc coords
#' @description
#' This function aims at defining the optimal threshold to set on the TFBS prediction score output by Wimtrap in order to infer
#' the potential gene targets of transcription factors. Subsequently, the performances at predicting gene targets are assessed
#' for each transcription factor considered by considering the whole set of potential TFBS on a given chromosome (unbalanced dataset).
#' @param TFBSdata A named character vector as output by the [getTFBSdata()] function, defining the local paths
#' to files encoding for the results of pattern-matching and geonmic feature extraction for the training TFs and/or
#' studied TFs.
#' @param TFBSmodel A \code{xgb.Booster} object as output by the function [buildTFBSmodel()].
#' @param chrTest An integer specifying the number of the chromosome that will be considered to assess
#' the performances at predicting the TF gene targets. Default = 1.
#' @param  ChIPpeaks A named character vector defining the local paths to BED files encoding the
#' location of ChIP-peaks. The vector is named according to the transcription factors that are described
#'  by the files indicated. **Caution**: the names of the \code{ChIPpeaks} have to find a match with those of \code{TFBSdata}.
#'  Default is \code{NULL} and 
#' @param ChIPpeaks_length An integer setting a fixed length for the ChIP-peaks, that are defined as the intervals of
#' \code{ChIPpeaks_length} bp that are centered on the regions encoded in the \code{ChIPpeaks} files. Default value = 400.
#' @param tss A list of \code{GRanges} objects as output by [importGenomicData()] or local path to a BED file defining 
#' the transcription stat site (TSS), name and orientation of each protein-coding transcript of the organism.
#' @param targets A named character vector defining the local paths to text files encoding the
#' manually curated transcriptional targets of each transcription factor. The vector is named according to the transcription factors
#' that are described by the files indicated. **Caution**: the names of the \code{ChIPpeaks} have to find
#' a match with those of \code{TFBSdata}. Default = NULL (Transcriptional targets are predicted from ChIP-peaks)
#' @details  Each gene is at first scored with the highest prediction score among the TFBSs associated with it and predicted by Wimtrap. 
#' Each gene is then labelled as positive or negative. The positive genes are the genes whose the TSS is the closest to an
#' occurrence on a ChIP-peak of the cognate TF-primary motif. This allows to draw a ROC curve based on a balanced dataset obtained from all the
#' chromosomes but one and to identify the best threshold to set on the prediction score in order to predict TF gene targets.
#' Finally, the performances are assessed for each TF based on the whole dataset of predicted TFBS on the left-over chromosome.
#' @seealso [plotPredictions()] to vizualize the results for a given potential target gene.
#' @return A \code{data.frame} that gives, for each TF considered, the performances of prediction of the transcriptional
#' targets encoded on the test chromosome, taking into consideration all the TFBSs predicted by Wimtrap (prediction score >= 0.5) 
#' on that chromosome.
#' Due to the highly imbalanced dataset, the performances are expressed in terms of recall, precision, accuracy and F-score.
#' In addition, in the last column, is presented the optimal threshold obtained when including all the input TFs.
#' @examples
#' genomic_data.ex <- c(CE = system.file("extdata/conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("extdata/DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("extdata/DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("extdata/x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("extdata/cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("extdata/intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("extdata/x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex,
#'                                               tss = system.file("extdata/tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("extdata/tts_example.bed", package = "Wimtrap"))
#' TFBSdata.ex <- getTFBSdata(pfm = system.file("extdata/pfm_example.pfm", package = "Wimtrap"),
#'                            TFnames = c("PIF3", "TOC1"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("extdata/genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
#' TFBSmodel.ex <- buildTFBSmodel(TFBSdata = TFBSdata.ex,
#'                                ChIPpeaks = c(PIF3 = system.file("extdata/PIF3_example.bed", package = "Wimtrap"),
#'                                              TOC1 = system.file("extdata/TOC1_example.bed", package = "Wimtrap")),
#'                                TFs_validation = "PIF3")
#' ##Determine the optimal score threshold
#' targetPerformances <- testTargetPredictions(
#' TFBSdata = TFBSdata.ex["TOC1"],
#' TFBSmodel = TFBSmodel.ex, 
#' tss = imported_genomic_data.ex,
#' ChIPpeaks =  c(TOC1 = system.file("extdata/TOC1_example.bed", package = "Wimtrap")))
#' optimal_threshold <- targetPerformances$threshold[1]
#' PIF3BS.predictions <- predictTFBS(TFBSmodel.ex,
#'                                   TFBSdata.ex,
#'                                   studiedTFs = "PIF3",
#'                                   score_threshold = optimal_threshold)
#' ##To get the transcripts whose expression is potentially regulated by PIF3 do as follows:
#' PIF3_regulated.predictions <- as.character(PIF3BS.predictions$transcript[!duplicated(PIF3BS.predictions)])
#' ###If you want to consider only the gene model,
#' ###then do as follows:
#' PIF3_regulated.predictions <- unlist(strsplit(PIF3_regulated.predictions, "[.]"))[seq(1, 2*length(PIF3_regulated.predictions),2)]
#' PIF3_regulated.predictions <- PIF3_regulated.predictions[!duplicated(PIF3_regulated.predictions)]

testTargetPredictions <- function(
  TFBSdata,
  TFBSmodel,
  chrTest = 1,
  tss,
  ChIPpeaks = NULL,
  ChIPpeaks_length = 400,
  targets = NULL
  )
{
  if (is.character(tss)){
    tss <- Wimtrap::importGenomicData(genomic_data = tss,
                                      tts = tss,
                                      tss = tss)$TSS
  } else {
    tss <- tss$TSS
  }
  allgenes <- unlist(strsplit(x = rtracklayer::mcols(tss)[,1], split = "[.]"))[
    seq(1, length(tss)*2, 2)
    ]
  tss$name <- allgenes
  allgenes <- allgenes[!duplicated(allgenes)]
  tss <- as.data.frame(tss)
  allgenes_Chr5 <- tss$name[tss$seqnames == chrTest]
  allgenes_Chr5 <- allgenes_Chr5[!duplicated(allgenes_Chr5)]
  PSOnBins_all <- list()
  for (TF in names(TFBSdata)) {
    #Obtain the full dataset
    annotations <- read.table(TFBSdata[TF],
                              header = TRUE, sep = "\t")
    annotations <- GenomicRanges::makeGRangesFromDataFrame(annotations,
                                                           keep.extra.columns = TRUE)
    if (is.null(targets)){
      ChIPdata <- read.table(ChIPpeaks[TF], sep = "\t")
      colnames(ChIPdata) <- c("seqnames", "start", "end")
      ChIPdata$strand <- "*"
      ChIPdata$seqnames <- Wimtrap:::getRiddChr(ChIPdata$seqnames)
      ChIPdata$start <- (ChIPdata$start + (ChIPdata$end-ChIPdata$start)/2) - 200
      ChIPdata$end <- ChIPdata$start + 400
      ChIPdata <- GenomicRanges::makeGRangesFromDataFrame(ChIPdata)
      
      overlaps <- GenomicRanges::findOverlaps(query = annotations, 
                                              subject = ChIPdata)
      
      ChIP.peaks <- rep(0, length(annotations))
      ChIP.peaks[overlaps@from[!duplicated(overlaps@from)]] <- 1
      
      #Get the gene targets of the TF
      targets <- annotations$ClosestTSS[ChIP.peaks == 1]
      targets <- unlist(strsplit(x = as.character(targets), split = "[.]"))[seq(1, 2*length(targets),2)]
      targets <- targets[!duplicated(targets)]
    } else {
      targets <- read.table(targets[TF], sep = "\t")
      targets <- as.character(targets[,1])
    }
    annotations <- as.data.frame(annotations)
    
    #Get predictions
    assessed <- annotations[,-seq(1:6)]
    assessed <- assessed[,colnames(assessed) %in% TFBSmodel$feature_names]
    assessed <- as.matrix(assessed)
    assessed <- xgboost::xgb.DMatrix(data = assessed)
    
    predictions <- stats::predict(TFBSmodel, assessed)
    annotations$predictions <- predictions
    annotations <- annotations[predictions >= 0.5,]
    
    kept <- annotations[,c("seqnames", "start", "end", "width", "strand", "predictions", "ClosestTSS")]
    kept$ClosestTSS <- unlist(strsplit(x = as.character(kept$ClosestTSS), split = "[.]"))[seq(1, 2*length(kept$ClosestTSS),2)]
    
    PSOnBins <- data.frame(seqnames = kept$seqnames[!duplicated(kept$ClosestTSS)],
                           GeneID = kept$ClosestTSS[!duplicated(kept$ClosestTSS)],
                           Max = 0)
    rownames(PSOnBins) <- PSOnBins$GeneID
    max_TFBS <- tapply(kept$predictions, kept$ClosestTSS, max, na.rm = TRUE)
    PSOnBins[as.character(names(max_TFBS)), "Max"] <- max_TFBS
    PSOnBins$Target <- 0
    PSOnBins$Target[rownames(PSOnBins) %in% targets] <- 1
    PSOnBins_all[[TF]] <- PSOnBins
  }
  PSOnBins <- do.call(rbind, PSOnBins_all)
  #Select randomly the individuals used for training
  PSOnBins <- PSOnBins[order(PSOnBins$Target),]
  
  #Keep an unbalanced dataset on the Chr5
  
  PSOnBins_test <- PSOnBins[which(PSOnBins$seqnames == chrTest),]
  PSOnBins_train <- PSOnBins[-which(PSOnBins$seqnames == chrTest),] 
  
  
  #Balance the dataset on the other chromosomes
  NbPosNeg <- table(PSOnBins_train$Target)
  PSOnBins <- PSOnBins_train[c(sample.int(NbPosNeg[1], NbPosNeg[2]),
                               (seq((NbPosNeg[1]+1), nrow(PSOnBins_train)))),]
  
  NbToSample <- NbPosNeg[2]
  if (NbToSample >= 10000) {NbToSample <- 10000}
  
  training_ids <- c(sample.int(NbPosNeg[2], NbToSample),
                    sample.int(NbPosNeg[2], NbToSample)+NbPosNeg[2])
  
  
  ##Performances considering Max
  roc_test <- pROC::roc(response = as.integer(PSOnBins_test$Target), predictor = PSOnBins_test$Max)
  threshold <- pROC::coords(roc_test, "best", ret="threshold", transpose = FALSE)
  
  #Assess the performances of the general model on each TF
  
  ##Performances considerind Max, Nb, Mean
  performances <- list()
  for (TF in names(TFBSdata)){
    print(TF)
    PSOnBins <- PSOnBins_all[[TF]]
    PSOnBins_test <- PSOnBins[PSOnBins$seqnames == chrTest,]
    geneswithoutTFBS <- allgenes_Chr5[!(allgenes_Chr5 %in% PSOnBins_test$GeneID)]
    geneswithoutTFBS <- data.frame(seqnames = chrTest,
                                   GeneID = geneswithoutTFBS,
                                   Max = 0,
                                   Target = 0)
    geneswithoutTFBS$Target[geneswithoutTFBS$GeneID %in% targets] <- 1
    PSOnBins_test <- rbind(PSOnBins_test, geneswithoutTFBS)
    targetactual <- as.character(PSOnBins_test[,"Target"])
    targetpredictions <- PSOnBins_test$Max >= as.numeric(threshold2[1])
    targetpredictions <- as.character(as.integer(targetpredictions))
    cm <- table(targetactual, targetpredictions)
    recall <- (cm[2,2])/(cm[2,1]+cm[2,2])
    precision <- (cm[2,2])/(cm[1,2]+cm[2,2])
    f.score <- (precision * recall) / (precision + recall)
    accuracy <- (cm[1,1] + cm[2,2])/length(targetactual)
    performance <- data.frame(recall = recall, precision = precision, 
                              accuracy = accuracy, f.score = f.score, 
                              TF = TF, threshold = threshold)
    performances[[TF]] <- performance
  }
  performances <- as.data.frame(do.call(rbind, performances))
  for (i in c(1:4, 6)){
    performances[,i] <- as.numeric(performances[,i])
  }
  performances$threshold <- round(as.numeric(performances$threshold), 4)
  return(performances)
}
#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Allows to import as R objects all the source data that will be used to carry on the
### pattern-matching analysis and to annotate the PWM-matches with contextual genomic
### features.
#========================================================================================#
inputSourceData <- function(ChIPpeaks,
                            pfm = NULL,
                            organism = NULL,
                            genome_sequence = "getGenome",
                            matches = NULL,
                            strand_as_feature = FALSE,
                            studiedTF = NULL)
{
  DNAregionsName = "Motifs"

  if (is.null(pfm)) {
    pfm_path <- NULL
    directInput_pwmMatrices <- NULL
  } else {
    if (class(pfm) == "character") {
      pfm_path <- pfm
      directInput_pwmMatrices <- NULL
    } else {
      pfm_path <- NULL
      directInput_pwmMatrices <- pfm
    }
  }

  #@ Location of the ChIP-peaks related to the TFs for which
  #@ the PFMs are provided can be input from a BED3 file
  #@ or from a GRanges object
  if (!is.null(ChIPpeaks)){

    if (class(ChIPpeaks) == "character") {
      ChIP_paths <- ChIPpeaks
      directInput_ChIPGRanges <- NULL
    } else {
      ChIP_paths <- NULL
      directInput_ChIPGRanges <- ChIPpeaks
    }
  }

  #@ The genome sequence is not necessary to get
  #@ if the user achieves the patern-matching analysis
  #@ with an external tool.
  #@ Otherwise the genome sequence can be imported either
  #@ from a FASTA file, a DNAStringSet or can be downloaded.
  if (is.null(genome_sequence)) {
    genome_path <- NULL
    directInput_DNAsequences <- NULL
  } else {
    if (class(genome_sequence) != "DNAStringSet") {
      if (genome_sequence == "getGenome") {
        genome_path <- NULL
        directInput_DNAsequences <- NULL
      } else {
        genome_path <- genome_sequence
        directInput_DNAsequences <- NULL
      }
    } else {
      genome_path <- NULL
      directInput_DNAsequences <- genome_sequence
    }
  }

  #@ We have to indicate whether the user has input his
  #@ own matches or not, as a list of GRanges objects.
  #@ The user can optionally indicate the matches scores
  #@ in the first metadata column of the GRanges objects,
  #@ as raw scores or as p-values scores.
  #@ If the pattern-matching analysis is carried on by wimtrap,
  #@ the raw score is reported in the first metadata column
  #@ and the p-value in the second. But only the p-value will
  #@ be considered as predictive feature.
  if (is.null(matches)) {
    directInput_matches <- NULL
    DirectMatches <- FALSE
    NoScore <- FALSE
    IndexMatch <- 2
  } else {
    directInput_matches <- matches
    if (!(is.list(directInput_matches))) {
      directInput_matches <- list(directInput_matches)
    } else {

    }
    DirectMatches <- TRUE
    if (length(rtracklayer::mcols(directInput_matches[[1]])) == 0) {
      NoScore <- TRUE
      IndexMatch <- NULL
    } else {
      NoScore <- FALSE
      IndexMatch <- 1
    }
    for (match in seq_along(directInput_matches)) {
      #@ The code will crash if the matches input by the user are not
      #@ oriented.
      if (length(which(as.logical(GenomicRanges::strand(directInput_matches[[match]]) == "*"))) ==
          length(GenomicRanges::strand(directInput_matches[[match]]) == "*")) {
        if (strand_as_feature == FALSE) {
          #@ Inventing the orientation of the matches won't affect the annotations,
          #@ except if the orientation of the matches to the genes whose the TSS
          #@ is the closest is taken into account.
          GenomicRanges::strand(directInput_matches[[match]]) <-
            S4Vectors::Rle(values = c("+", "-"),
                           length = c(1, length(directInput_matches[[match]]) -
                                        1))
        } else {
          stop("Please provide oriented matches")
        }
      } else {

      }
    }
  }
  ####################################################################

  ####################################################################
  ### Carry on all the preliminary steps to the modeling: import the
  ### data, identify the potential binding sites by pattern-matching,
  ### annotate them with the considered features and label them as
  ### TRUE or FALSE according to the ChIP-seq data.
  ### Return as well the objects storing the sources data and options
  no_score = NoScore
  index_match = IndexMatch
  dna_regions_name = DNAregionsName


  #@ The first part of the code is mainly format ckecking
  #@ Globally, data can be input from objects in the R environment
  #@ (cf. "directInput") or from a source file (cf. "paths").
  #@ In certain cases data can be as well downloaded

  if(!is.null(ChIPpeaks)){
    ChIP_regions <- listChIPRegions(ChIP_paths, directInput_ChIPGRanges, 400)
  } else {
    ChIP_regions <- NULL
  }
  #@ Get the sequence of the chromosomes of the genome either from
  #@ the fasta file whose the path is specified in the genome_path
  #@ parameter, either directly from a database (Ensembl, Ensembl Genomes,
  #@ Refseq or Genebank) through the biomartr package
  if (missing(genome_path) | is.null(genome_path)) {
    if (length(directInput_matches) == 0) {
      if (length(directInput_DNAsequences) > 0) {
        genome <- directInput_DNAsequences
      } else {
        message("Downloading chromosomes DNA sequences from Ensembl...")
        genome <- getChromosomes(organism)
      }
    } else {
      genome <- NULL
    }
  } else {
    genome <- Biostrings::readDNAStringSet(genome_path)
  }

  #@ Avoid any discreapency between sources data
  #@ by deleting the "Chr", "CHR" and "chr" prefixes
  #@ from the chromosomes names
  if (!(is.null(genome))) {
    names(genome) <- getRiddChr(names(genome))
  } else {

  }

  #@ Treatment of the PFM or PWM of the TFs

  #@ If the considered organism is A. thaliana,
  #@ the PWM of the TFs can be optionally automatically
  #@ retrieved from the PlantTFDB database

  #@ The following steps are not necessary if the user
  #@ has performed the pattern-matching analysis with an
  #@ external tool
  if(!is.null(ChIPpeaks)){
    TFs <- names(ChIP_regions)
  } else {
    if(!is.null(studiedTF)){
      TFs <- studiedTF
    } else {
      TFs <- NULL
    }
  }

  if (length(pfm_path) >0){
   Pwm <- getPwm(pfm_path, directInput_pwmMatrices, organism, TFs ,genome)
   Pwm <- Pwm[names(Pwm) %in% TFs]
  } else {
   Pwm <- NULL
  }

  sourceData <- list(
    ChIP = ChIP_regions,
    genome = genome,
    pwm = Pwm,
    matches = directInput_matches
  )

  return(sourceData)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Import information encoded as BED or GFF files as GRanges objects
#========================================================================================#
importFeatures <- function(file_paths)
{
  feature_GRanges_list <- list()
  for (i in seq_along(file_paths)) {
    BEDorGFF <- unlist(strsplit(file_paths[i], "[.]"))[length(strsplit(file_paths[i],
                                                                       "[.]")[[1]])]
    if (BEDorGFF == "bed") {
      BEDorGFF <- 5
    }
    else {
      if (BEDorGFF == "gff") {
        BEDorGFF <- 6
      }
      else {
      }
    }
    if (is.numeric(BEDorGFF)) {
      NumOrCat_Data <- utils::read.csv(file_paths[[i]],
                                       sep = "\t")
      if (BEDorGFF == 5) {
        if (length(NumOrCat_Data) > 4) {
          if (is.factor(NumOrCat_Data[, BEDorGFF])) {
            if (length(NumOrCat_Data) > 5) {
              Strands = NumOrCat_Data[, 6]
            }
            else {
              Strands = rep("*", dim(NumOrCat_Data)[1])
            }
            feature_GRanges <- data.frame(seqnames = NumOrCat_Data[,
                                                                   1], start = NumOrCat_Data[, 2], end = NumOrCat_Data[,
                                                                                                                       3], strand = Strands, score = as.character(NumOrCat_Data[,
                                                                                                                                                                                5]))
          }
          else {
            feature_GRanges <- as.data.frame(rtracklayer::import(file_paths[i]))
          }
        }
        else {
          feature_GRanges <- as.data.frame(rtracklayer::import(file_paths[i]))
        }
      }
      else {
        if (length(NumOrCat_Data) > 5) {
          if (is.factor(NumOrCat_Data[, BEDorGFF])) {
            if (length(NumOrCat_Data) > 6) {
              Strands = NumOrCat_Data[, 7]
            }
            else {
              Strands = rep("*", dim(NumOrCat_Data)[1])
            }
            feature_GRanges <- data.frame(seqnames = NumOrCat_Data[,
                                                                   1], start = NumOrCat_Data[, 4], end = NumOrCat_Data[,
                                                                                                                       5], strand = Strands, score = as.character(NumOrCat_Data[,
                                                                                                                                                                                6]))
          }
          else {
            feature_GRanges <- as.data.frame(rtracklayer::import(file_paths[i]))
          }
        }
        else {
          feature_GRanges <- as.data.frame(rtracklayer::import(file_paths[i]))
        }
      }
    }
    else {
      feature_GRanges <- as.data.frame(rtracklayer::import(file_paths[i]))
    }
    if (is.null(feature_GRanges$score)) {
      feature_GRanges <- data.frame(seqnames = feature_GRanges$seqnames,
                                    start = feature_GRanges$start, end = feature_GRanges$end,
                                    strand = feature_GRanges$strand)
    }
    else {
      feature_GRanges <- data.frame(seqnames = feature_GRanges$seqnames,
                                    start = feature_GRanges$start, end = feature_GRanges$end,
                                    strand = feature_GRanges$strand, score = feature_GRanges$score)
    }
    feature_GRanges$seqnames <- getRiddChr(feature_GRanges$seqnames)
    FeatureName <- strsplit(file_paths[i], "/")[[1]]
    FeatureName <- FeatureName[length(FeatureName)]
    FeatureName <- strsplit(FeatureName, "[.]")[[1]][1]
    FeatureName <- FeatureName[length(FeatureName)]
    if (!(is.null(feature_GRanges$score))) {
      colnames(feature_GRanges)[5] <- FeatureName
    }
    feature_GRanges <- GenomicRanges::makeGRangesFromDataFrame(feature_GRanges,
                                                               keep.extra.columns = TRUE)
    feature_GRanges_list[[i]] <- feature_GRanges
    names(feature_GRanges_list)[[i]] <- FeatureName
  }
  return(feature_GRanges_list)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Lists for each TF as GRanges object the location of their ChIP-peaks,
### which is input by user as GRanges objects or from BED source file(s)
### The regions defining the ChIP-peaks are limited to a fixed lenght
### anchored to the center
#========================================================================================#
listChIPRegions <- function(ChIP_paths, directInput_ChIPGRanges, width){
  #@Import regions from BED files
  ChIP_regions <- lapply(ChIP_paths, read.delim, header = FALSE)
  #Limit the length of the ChIP-peaks to
  #400bp around the center
  for (i in 1:length(ChIP_regions)){
    considered <- ChIP_regions[[i]]
    considered <- data.frame(
      seqnames = considered[,1],
      start = considered[,2],
      end = considered[,3]
    )
    considered$seqnames <- getRiddChr(considered$seqnames)
    considered <- GenomicRanges::makeGRangesFromDataFrame(considered)
    considered <- IRanges::resize(considered, width, fix = "center")
    ChIP_regions[[i]] <- considered
  }
  names(ChIP_regions) <- names(ChIP_paths)
  ChIP_regions <- c(ChIP_regions, directInput_ChIPGRanges)
  return(ChIP_regions)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Returns a list of matrices describing the Pwm of each considered TF.
### The PWM might be directly input as a matrix or imported from a source
### file encoded in jaspar.
### If the organism is Arabidopsis thaliana, the PWM can be automatically
### searched in the PlanTFDB database
#========================================================================================#
getPwm <- function(pfm_path = NULL, directInput_pwmMatrices = NULL, organism, TFs, genome){
  TFNames <- TFs
  if (!is.null(directInput_pwmMatrices)) {
    if (is.list(directInput_pwmMatrices)) {
      directInput_pwmMatrices <- directInput_pwmMatrices
    } else {
      directInput_pwmMatrices <- list(directInput_pwmMatrices)
    }
    input_pwm <- list()
    for (i in seq_along(directInput_pwmMatrices)) {
      input_pwm[[i]] <-
        buildPwm(directInput_pwmMatrices[[i]], genome)
      names(input_pwm)[[i]] <-
        names(directInput_pwmMatrices)[i]
      rownames(input_pwm[[i]]) <- c("A", "C", "G", "T")
    }
  } else {

  }

  if (!is.null(pfm_path)) {
    List_pfm <- readPSFM(pfm_path)
    input_pwm <- list()

    for (i in seq_along(List_pfm)) {
      input_pwm[[i]] <- buildPwm(List_pfm[[i]], genome)
      names(input_pwm)[[i]] <- names(List_pfm)[i]
      rownames(input_pwm[[i]]) <- c("A", "C", "G", "T")
    }
  } else {

  }

  if (length(organism) > 0 && organism == "Arabidopsis thaliana") {
    #@ Identification of the TFs for which no PFM/PWM has been
    #@ specified and thus for which the PWM has to be retrieved
    #@ from the database
    if (!is.null(pfm_path) | !is.null(directInput_pwmMatrices)) {
      TFDBquery <- TFNames[!(TFNames %in% names(input_pwm))]
    } else {
      input_pwm <- NULL
      TFDBquery <- TFNames
    }

    #@ Query the database
    if (length(TFDBquery) == 0) {
      Pwm <- input_pwm
    } else {
      Pwm <- list()
      for (TF in seq_along(TFDBquery)) {
        data("ath_ID_mapping", package = "Wimtrap")
        data("ath_PWM", package = "Wimtrap")
        list_columns <- list()

        for (i in seq_along(ath_ID_mapping)) {
          list_columns[[i]] <- ath_ID_mapping[, i]
        }

        GeneLine <- sapply(list_columns,
                           function(x) {
                             Line = grep(pattern = TFDBquery[TF], as.character(x))
                             return(Line[1])
                           })

        if (all(is.na(GeneLine))) {
          stop(paste0(
            c(
              "The motif of ",
              TFDBquery[TF],
              " is not described by plantTFDB.\nInput it through the \"psfm\" argument."
            )
          ))
        } else {
          GeneID <- ath_ID_mapping[GeneLine[which(!(is.na(GeneLine)))[1]], 1]
          Pwm[[TF]] <- ath_PWM[[TFDBquery[TF]]]
          names(Pwm)[[TF]] <- TFDBquery[TF]
        }

      }

      if (!(is.null(input_pwm))) {
        Pwm <- c(Pwm, input_pwm)
      } else {

      }
    }
  } else {
    Pwm <- input_pwm
  }
  return(Pwm)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Construct Position Pseudo--Weighted Matrices from Position Frequency Matrices
#========================================================================================#
buildPwm <- function (pfm, genome)
{
  NbChromosomes <- length(genome)
  if (Biostrings::width(genome)[1] > 2e+06) {
    Start <- stats::runif(1000, 1, Biostrings::width(genome)[1] -
                            2000)
    SequencesSample <- Biostrings::Views(genome[[1]], start = Start,
                                         end = Start + 2000)
    NucleotideComp <- apply(Biostrings::letterFrequency(SequencesSample,
                                                        c("A", "C", "G", "T")), 2, mean)/2000
  }
  else {
    NucleotideComp <- Biostrings::letterFrequency(genome[[1]],
                                                  c("A", "C", "G", "T"))/length(genome[[1]])
  }
  for (i in seq(2, NbChromosomes)) {
    if (Biostrings::width(genome)[i] > 2e+06) {
      Start <- stats::runif(1000, 1, Biostrings::width(genome)[i] -
                              2000)
      SequencesSample <- Biostrings::Views(genome[[i]],
                                           start = Start, end = Start + 2000)
      NucleotideComp <- (NucleotideComp + apply(Biostrings::letterFrequency(SequencesSample,
                                                                            c("A", "C", "G", "T")), 2, mean)/2000)/2
    }
    else {
      NucleotideComp <- (NucleotideComp + Biostrings::letterFrequency(genome[[i]],
                                                                      c("A", "C", "G", "T"))/length(genome[[i]]))/2
    }
  }
  NucleotideComp <- round(NucleotideComp, 3)
  TotalCounts <- apply(pfm, 2, sum)
  b <- rep(sqrt(TotalCounts), 4)
  pwm <- (pfm + NucleotideComp * b)/(TotalCounts + b)
  return(pwm)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Import as list of matrices Pseudo Frequency Matrices encoded in a raw Jaspar file
#========================================================================================#
readPSFM <- function (file_path)
{
  extension <- unlist(strsplit(file_path, "[.]"))[length(unlist(strsplit(file_path, "[.]")))]
  if (extension == "pfm"){
    list_pfm <- list()
    readRaw <- suppressWarnings(readLines(file(file_path)))
    MotifIDLines <- which(regexpr("^>", readRaw) >= 1)
    for (i in seq_along(MotifIDLines)) {
      MotifID <- unlist(strsplit(readRaw[MotifIDLines[i]],
                                 ">"))[2]
      MotifID <- unlist(strsplit(MotifID, "\\s+"))[[1]]
      A_frequencies <- unlist(strsplit(readRaw[MotifIDLines[i] +
                                                 1], "\\s+"))
      A_frequencies <- suppressWarnings(as.numeric(A_frequencies))
      A_frequencies <- A_frequencies[!(is.na(A_frequencies))]
      Frequencies <- matrix(A_frequencies, nrow = 1)
      for (j in seq(2, 4)) {
        X_frequencies <- unlist(strsplit(readRaw[MotifIDLines[i] +
                                                   j], "\\s+"))
        X_frequencies <- suppressWarnings(as.numeric(X_frequencies))
        X_frequencies <- X_frequencies[!(is.na(X_frequencies))]
        Frequencies <- rbind(Frequencies, X_frequencies)
      }
      rownames(Frequencies) <- NULL
      list_pfm[[i]] <- Frequencies
      rownames(list_pfm[[i]]) <- c("A", "C", "G", "T")
      names(list_pfm)[i] <- MotifID
    }
  } else if (extension == "txt"){
    list_pfm <- universalmotif::read_cisbp(file_path)
    TFnames <- c()
    for (TF in seq_along(list_pfm)){
      TFnames <- c(TFnames, list_pfm[[TF]]@name)
      list_pfm[[TF]] <- list_pfm[[TF]]@motif
      colnames(list_pfm[[TF]]) <- NULL
    }
    names(list_pfm) <- TFnames
  } else if (extension == "motif"){
    list_pfm <- universalmotif::read_homer(file_path)
    TFnames <- c()
    for (TF in seq_along(list_pfm)){
      TFnames <- c(TFnames, list_pfm[[TF]]@name)
      list_pfm[[TF]] <- list_pfm[[TF]]@motif
      colnames(list_pfm[[TF]]) <- NULL
    }
    names(list_pfm) <- TFnames
  } else if (extension == "jaspar"){
    list_pfm <- universalmotif::read_jaspar(file_path)
    TFnames <- c()
    for (TF in seq_along(list_pfm)){
      TFnames <- c(TFnames, list_pfm[[TF]]@name)
      list_pfm[[TF]] <- list_pfm[[TF]]@motif
      colnames(list_pfm[[TF]]) <- NULL
    }
    names(list_pfm) <- TFnames
  } else if (extension == "transfac"){
    list_pfm <- universalmotif::read_transfac(file_path)
    TFnames <- c()
    for (TF in seq_along(list_pfm)){
      TFnames <- c(TFnames, list_pfm[[TF]]@name)
      list_pfm[[TF]] <- list_pfm[[TF]]@motif
      colnames(list_pfm[[TF]]) <- NULL
    }
    names(list_pfm) <- TFnames
  }  else if (extension == "meme"){
    list_pfm <- universalmotif::read_meme(file_path)
    TFnames <- c()
    for (TF in seq_along(list_pfm)){
      TFnames <- c(TFnames, list_pfm[[TF]]@name)
      list_pfm[[TF]] <- list_pfm[[TF]]@motif
      colnames(list_pfm[[TF]]) <- NULL
    }
    names(list_pfm) <- TFnames
  }
  return(list_pfm)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTIONt
#========================================================================================#
### Download the genome sequence from ENSEMBL or ENSEMBL GENOMES
#========================================================================================#
getChromosomes <- function(organism)
{
  file_path <- tryCatch(getsequence(organism = organism), error = function(e) {return(NULL)}, finally = message("Interrogating Ensembl Plants"))
  if(is.null(file_path)) {
    file_path <- biomartr::getGenome( db = "ensemblgenomes",
                          organism,
                          path = file.path("_ncbi_downloads","genomes"))
    Genome <- Biostrings::readDNAStringSet(file_path)
    ChrNames <- GenomeInfoDb::seqlevels(Genome)
    SplitChrNames <- unlist(lapply(as.list(ChrNames), base::strsplit,
                                   split = " "))
    FieldNumber <- which(SplitChrNames=="chromosome")
    ChrNames <- SplitChrNames[FieldNumber+1]
    ChrNames <- getRiddChr(ChrNames)
    Genome <- Genome[1:length(ChrNames)]
    names(Genome) <- ChrNames

  } else {
    Genome <- Genome[1:length(ChrNames)]
    names(Genome) <- ChrNames
    Genome <- Biostrings::readDNAStringSet(file_path)
    ChrNames <- GenomeInfoDb::seqlevels(Genome)
    SplitChrNames <- lapply(as.list(ChrNames), base::strsplit,
                            split = " ")
    SplitChrNames <- unlist(lapply(SplitChrNames, function(names) {
      return(names[[1]][1])}))
    ChrNames <- getRiddChr(SplitChrNames)
    Genome <- Genome[1:length(ChrNames)]
    names(Genome) <- ChrNames

    }
  
  return(Genome)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Remove the 'chr' prefix from chromosome names
#========================================================================================#
getRiddChr <- function (names)
{
  Character <- is.character(names)
  names <- base::tolower(names)
  names <- factor(names)
  regChr <- regexpr("^chr", levels(names))
  if (length(which(regChr > 0))) {
    NewNames <- c()
    for (i in which(regChr > 0)) {
      NewNames <- c(NewNames, unlist(strsplit(levels(names)[[i]],
                                              "chr"))[2])
    }
    levels(names)[which(regChr > 0)] <- NewNames
  }
  else {
  }
  if (Character) {
    names <- as.character(names)
  }
  return(names)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Make a connection to the biomart database (related to the given organism)
#========================================================================================#
loadEnsembl <- function (organism, biomart = NULL, dataset = NULL)
{
  OrganismsList <- data.frame(organism = c(), dataset = c(), mart = c())
  
  biomartsEnsembl <- biomaRt::listEnsembl()$biomart[grep(pattern = "Genes", 
                                                         x = biomaRt::listEnsembl()$version)]
  biomartsEnsemblGenomes <- biomaRt::listEnsemblGenomes()$biomart[grep(pattern = "Genes", 
                                                                       x = biomaRt::listEnsemblGenomes()$version)]
  #Let's test at first whether the species considered is a plant species
  
  biomart_considered <- biomartsEnsemblGenomes[grep(pattern = "plant", x = biomartsEnsemblGenomes)][1]
  ensembl_considered <- biomaRt::useEnsemblGenomes(biomart_considered)
  dataset_searched <- biomaRt::searchDatasets(ensembl_considered, pattern = organism)[1,]
  
  Ensembl <- c()
  
  if (!is.null(dataset_searched)){
    Ensembl <- biomaRt::useEnsemblGenomes(biomart_considered, dataset_searched$dataset)
    return(Ensembl)
  } else {
    #Let's query for databases related to other life kingdoms
    biomartsEnsemblGenomes <- biomartsEnsemblGenomes[-grep(pattern = "plant", x = biomartsEnsemblGenomes)]
    for (biomart_considered in biomartsEnsemblGenomes) {
      ensembl_considered <- biomaRt::useEnsemblGenomes(biomart_considered)
      dataset_searched <- biomaRt::searchDatasets(ensembl_considered, pattern = organism)[1,]
      if (!is.null(dataset_searched)){
        Ensembl <- biomaRt::useEnsemblGenomes(biomart_considered, dataset_searched$dataset)
        return(Ensembl)
      }
    }
  }
  
  if(length(Ensembl)==0){
    #Do we have a Ensembl genome?
    for (biomart_considered in biomartsEnsembl) {
      ensembl_considered <- biomaRt::useEnsembl(biomart_considered)
      dataset_searched <- biomaRt::searchDatasets(ensembl_considered, pattern = organism)[1,]
      if (!is.null(dataset_searched)){
        Ensembl <- biomaRt::useEnsembl(biomart_considered, dataset_searched$dataset)
        return(Ensembl)
      }
    }
  }
  
  if (length(Ensembl) == 0) {
    warning(base::paste0(c("Error: no dataset related to ",
                           organism, " has been found.")))
    stop()
  }
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Retrieve the structure of protein-coding genes from biomart database
#========================================================================================#
getStructure <- function (ensembl, promoter_length = 2000, downstream_length = 1000, proximal_length = 500)
{
  Transcripts <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                               "transcript_start", "transcript_end", "chromosome_name",
                                               "strand"), filters = "transcript_biotype", values = "protein_coding",
                                mart = ensembl)
  Transcripts <- Transcripts[order(Transcripts$chromosome_name, Transcripts$transcript_start),]
  X5UTRs <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                          "5_utr_start", "5_utr_end", "chromosome_name", "strand"),
                           filters = "transcript_biotype", values = "protein_coding",
                           mart = ensembl)
  X3UTRs <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                          "3_utr_start", "3_utr_end", "chromosome_name", "strand"),
                           filters = "transcript_biotype", values = "protein_coding",
                           mart = ensembl)
  CDSs <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                        "genomic_coding_start", "genomic_coding_end", "chromosome_name",
                                        "strand"), filters = "transcript_biotype", values = "protein_coding",
                         mart = ensembl)
  Exons <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                         "exon_chrom_start", "exon_chrom_end", "chromosome_name",
                                         "strand"), filters = "transcript_biotype", values = "protein_coding",
                          mart = ensembl)
  TTS <- Transcripts
  TTS$transcript_start[TTS$strand == 1] <- TTS$transcript_end[TTS$strand == 1]
  TTS$transcript_end[TTS$strand == -1] <- TTS$transcript_start[TTS$strand == -1]
  TSS <- Transcripts
  TSS$transcript_end[TSS$strand == 1] <- TSS$transcript_start[TSS$strand == 1]
  TSS$transcript_start[TSS$strand == -1] <- TSS$transcript_end[TSS$strand == -1]
  Promoters <- Transcripts
  Promoters$transcript_start[Promoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                      1] - promoter_length
  Promoters$transcript_end[Promoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                    1]
  Promoters$transcript_start[Promoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                     -1]
  Promoters$transcript_end[Promoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                   -1] + promoter_length
  colnames(Promoters)[c(2, 3)] <- c("promoter_start", "promoter_end")
  ProximalPromoters <- Transcripts
  ProximalPromoters$transcript_start[ProximalPromoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                      1] - proximal_length
  ProximalPromoters$transcript_end[ProximalPromoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                    1]
  ProximalPromoters$transcript_start[ProximalPromoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                     -1]
  ProximalPromoters$transcript_end[ProximalPromoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                   -1] + proximal_length
  colnames(ProximalPromoters)[c(2, 3)] <- c("proximal_promoter_start", "proximal_promoter_end")
  Downstreams <- Transcripts
  Downstreams$transcript_start[Downstreams$strand == 1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                        1]
  Downstreams$transcript_end[Downstreams$strand == 1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                      1] + downstream_length
  Downstreams$transcript_start[Downstreams$strand == -1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                           -1] - downstream_length
  Downstreams$transcript_end[Downstreams$strand == -1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                         -1]
  colnames(Downstreams)[c(2, 3)] <- c("downstream_start",
                                      "downstream_end")
  intronsGenerator <- function(exons_transcript) {
    IntronsTranscript <- exons_transcript[-1, ]
    IntronsTranscript$exon_chrom_start <- exons_transcript$exon_chrom_end[-dim(exons_transcript)[1]] +
      1
    IntronsTranscript$exon_chrom_end <- exons_transcript$exon_chrom_start[-1] -
      1
    return(IntronsTranscript)
  }
  ExonsSplit <- split(Exons, f = as.factor(Exons$ensembl_transcript_id))
  IntronsSplit <- base::lapply(ExonsSplit, intronsGenerator)
  Introns <- do.call(rbind, IntronsSplit)
  Structures <- list(ProximalPromoter = ProximalPromoters, Promoter = Promoters, X5UTR = X5UTRs,
                     CDS = CDSs, Intron = Introns, X3UTR = X3UTRs, Downstream = Downstreams,
                     TTS = TTS, TSS = TSS)
  StructureTracks <- list()
  StructuresNames <- c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                       "X3UTR", "Downstream", "TTS", "TSS")
  for (i in seq_along(Structures)) {
    structure <- Structures[[i]]
    structure <- structure[!(is.na(structure[, 2])), ]
    Widths <- structure[, 3] - structure[, 2] + 1
    structure <- structure[which(Widths > 0), ]
    chromosome_name <- structure$chromosome_name
    chromosome_name <- Wimtrap:::getRiddChr(chromosome_name)
    Track <- data.frame(seqnames = chromosome_name,
                        start = structure[,2],
                        end = structure[,3],
                        width = 1,
                        strand = structure$strand,
                        name = structure$ensembl_transcript_id
                        )
    Track$width <- Track$end - Track$start
    Track$strand[Track$strand==1] <- "+"
    Track$strand[Track$strand==-1] <- "-"
    Track <- GenomicRanges::makeGRangesFromDataFrame(Track, keep.extra.columns = TRUE)
    StructureTracks[[i]] <- Track
    names(StructureTracks)[[i]] <- StructuresNames[i]
  }
  return(StructureTracks)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Find the regions of 400bp that are potential regions bound by a considered TF
### based on a p-value matching threshold and that are located around a TSS.
### Annotate the matches with the different features considering different windows
### length
#========================================================================================#
getCandidatesRegions <- function(directInput_matches,
                                 Pwm,
                                 StructureTracks,
                                 FeatureRanges_list,
                                 genome,
                                 TFNames,
                                 pval_threshold,
                                 ChIP_regions,
                                 short_window = 20,
                                 medium_window = 400,
                                 long_window = 1000) {
  paths <- c()
  #@Define the regions that are around a TSS,
  #@from -promoter_length to +downstream_length
  structure_formatted <- lapply(seq_along(StructureTracks[seq(1, length(StructureTracks)-2)]),
                                function(i){
                                  x <- StructureTracks[[i]]
                                  x <- GenomicRanges::reduce(x)
                                  return(x)
                                })
  StructureTracks[seq(1, length(StructureTracks)-2)] <- structure_formatted
  rm(structure_formatted)
  TSSTrack <- as.data.frame(StructureTracks$TSS)
  TTSTrack <- as.data.frame(StructureTracks$TTS)
  colnames(TSSTrack)[6] <- "name"
  colnames(TTSTrack)[6] <- "name"
  TSSTrack <- TSSTrack[TSSTrack$name %in% TTSTrack$name,]
  TTSTrack <- TTSTrack[TTSTrack$name %in% TSSTrack$name,]
  TSSTrack <- TTSTrack[!(duplicated(TSSTrack$name)),]
  rownames(TSSTrack) <- TSSTrack$name
  TSSTrack <- TSSTrack[order(TSSTrack$name),]
  TTSTrack <- TTSTrack[order(TTSTrack$name),]
  TSSTrack$end <- TTSTrack$start
  DirectStranded <- TSSTrack[TSSTrack$strand == "+",]
  DirectStranded <- data.frame(seqnames = DirectStranded$seqnames,
                               start = DirectStranded$start - 5000,
                               end = DirectStranded$end + 5000,
                               strand = "+")
  ReverseStranded <- TSSTrack[TSSTrack$strand == "-",]
  ReverseStranded <- data.frame(seqnames = ReverseStranded$seqnames,
                                start = ReverseStranded$end - 5000,
                                end = ReverseStranded$start + 5000,
                                strand = "-")
  DRStranded <- rbind(DirectStranded,
                      ReverseStranded)
  rm(DirectStranded)
  rm(ReverseStranded)

  if (length(directInput_matches) == 0) {
    DRStranded <- DRStranded[DRStranded$seqnames %in% names(genome),]
    DRStranded <- split(DRStranded, f = DRStranded$seqnames)
    DRStranded <- lapply(DRStranded,
                         function(chrdata){
                           chrdata <- chrdata[which((chrdata$end < Biostrings::width(genome)[which(names(genome) == chrdata[1,1])]) & (chrdata$start > 0)),]
                           return(chrdata)
                         })
    DRStranded <- do.call(rbind, DRStranded)
    DRStranded <- GenomicRanges::makeGRangesFromDataFrame(DRStranded,
                                                          keep.extra.columns = TRUE)
    DRStranded <- GenomicRanges::reduce(DRStranded, ignore.strand = TRUE)
    message("Pattern matching...")
    #@ Scanning of the genome with the PWM
    Pwm <- Pwm[which(names(Pwm) %in% TFNames)]
    Matches <- list()
    for (pwm in names(Pwm)){
      Matches[[pwm]] <- matrixScan(pwm = Pwm[[pwm]], genome = genome, pval_threshold= pval_threshold)
      names(GenomicRanges::mcols(Matches[[pwm]])) <- c("matchScore", "matchLogPval")
    }
  } else {
    if (is.list(directInput_matches)) {
      Matches <- directInput_matches
      Pwm <- names(Matches)
      names(Pwm) <- names(Matches)
    } else {
      Matches <- list(directInput_matches)
      Pwm <- names(Matches)
      names(Pwm) <- names(Matches)
    }
  }
  Pwm <- Pwm[!duplicated(names(Pwm))]
  Matches <- lapply(Matches, function(x) {y <- x[rtracklayer::mcols(x)[,2] <= log10(pval_threshold)]; 
  if (length(y) == 0) {y <- x[which.min(rtracklayer::mcols(x)[,2])]};
  if (length(y) == 0) {x <- as.data.frame(x); x$matchLogPval <- log10(pval_threshold);
                       x <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE); y <- x};
  return(y)})
  for (i in seq_along(Matches)) {
    excluded <- c()
      if (length(Matches[[i]]) > 0){
      GenomeInfoDb::seqlevels(Matches[[i]]) <-
      getRiddChr(GenomeInfoDb::seqlevels(Matches[[i]]))
      message("Annotation of the matches...")

    #@ Annotate the matches associated to each TF
      print(paste("=>", names(Pwm)[i]))
      #Select only regions around TSS
      tmpa <- as.data.frame(Matches[[i]])
      tmpa$strand <- "*"
      tmpa <- GenomicRanges::makeGRangesFromDataFrame(tmpa)
          AllTracks <-
            c(FeatureRanges_list, Matches = tmpa)

      #@ Annotations with the structural features
      AnnotatedMatches <- Matches[[i]]
      for (structure in names(StructureTracks)[-c(length(StructureTracks)-1, length(StructureTracks))]){
        GenomicRanges::mcols(AnnotatedMatches)[structure] <- 0
       GenomicRanges::mcols(AnnotatedMatches)[GenomicRanges::findOverlaps(Matches[[i]], StructureTracks[[structure]])@from,structure] <- 1
      }
      AnnotatedMatches <- getClosest(AnnotatedMatches, StructureTracks, modeling = FALSE, TSS = FALSE)
      AnnotatedMatches <- getClosest(AnnotatedMatches, StructureTracks, modeling = FALSE, TSS = TRUE)
      if (length(StructureTracks) == 9){
        GenomicRanges::mcols(AnnotatedMatches)[,"Promoter"] <- 0
        promoter_length <- mean(GenomicRanges::width(StructureTracks$Promoter))
        downstream_length <- mean(GenomicRanges::width(StructureTracks$Downstream))
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS <=0) ,
                                              which(names(GenomicRanges::mcols(AnnotatedMatches)) %in% c("X5UTR", "CDS", "Intron", "X3UTR", "Downstream"))] <- 0
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS <=0 &
                                                         AnnotatedMatches$DistToClosestTSS > -promoter_length), "Promoter"] <- 1
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS > 0 &
                                                         AnnotatedMatches$DistToClosestTTS < downstream_length), "Downstream"] <- 1
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS > 0 &
                                                         AnnotatedMatches$DistToClosestTTS < 0),
                                                 "Downstream"] <- 0
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$Intron == 1 & AnnotatedMatches$CDS == 1), "Intron"] <- 0
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$X5UTR == 1 & AnnotatedMatches$Intron== 1), "Intron"] <- 0
        GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$X5UTR == 1 & AnnotatedMatches$CDS== 1), "CDS"] <- 0
      }

      #@ Annotations with the categorical features
      categorical <- c()
      for (feature in names(AllTracks)){
        if (length(AllTracks[[feature]]@elementMetadata@listData) > 0){
          if(is.factor(AllTracks[[feature]]@elementMetadata@listData[[1]])){
            categorical <- c(categorical, feature)
          }
        }
      }

      for (categoricalfeature in categorical){
       GenomicRanges::mcols(AnnotatedMatches)[categoricalfeature] <- "NA"
       GenomicRanges::mcols(AnnotatedMatches)[GenomicRanges::findOverlaps(Matches[[i]], AllTracks[[categoricalfeature]])@from,
                                               categoricalfeature] <- as.character(GenomicRanges::mcols(AllTracks[[categoricalfeature]])[GenomicRanges::findOverlaps(Matches[[i]], AllTracks[[categoricalfeature]])@to,1])
      }

      #@ Annotations with the overlapping and numeric features
      #@ considering different windows
      for (feature in names(AllTracks)){
        if (length(AllTracks[[feature]]@elementMetadata@listData) == 0){
         GenomicRanges::mcols(AllTracks[[feature]])[feature] <- 1
        }
      }

      Matches[[i]] <- Matches[[i]][as.character(GenomeInfoDb::seqnames(Matches[[i]])) %in% as.character(GenomeInfoDb::seqnames(StructureTracks$TSS)),]
      for (windowlength in c(short_window, medium_window, long_window)){
        queries <- GenomicRanges::resize(Matches[[i]], windowlength, fix = "center")
        for (feature in names(AllTracks)[!(names(AllTracks) %in% categorical)]){
          overlaps <- GenomicRanges::findOverlaps(queries, AllTracks[[feature]])
          peakcuts <- data.frame(start.query = GenomicRanges::start(queries[overlaps@from]),
                                 end.query = GenomicRanges::end(queries[overlaps@from]),
                                 start.subject = GenomicRanges::start(AllTracks[[feature]][overlaps@to]),
                                 end.subject = GenomicRanges::end(AllTracks[[feature]][overlaps@to]))
          peakcuts$start.subject[peakcuts$start.subject < peakcuts$start.query] <- peakcuts$start.query[peakcuts$start.subject < peakcuts$start.query]
          peakcuts$end.subject[peakcuts$end.subject > peakcuts$end.query] <- peakcuts$end.query[peakcuts$end.subject > peakcuts$end.query]
          peakcuts$contribution <- (GenomicRanges::mcols(AllTracks[[feature]])[overlaps@to,1]*(peakcuts$end.subject-peakcuts$start.subject+1))/windowlength
         GenomicRanges::mcols(AnnotatedMatches)[paste0(c(feature, "_", windowlength, "bp"), collapse = "")] <- 0
         GenomicRanges::mcols(AnnotatedMatches)[overlaps@from[!duplicated(overlaps@from)], paste0(c(feature, "_", windowlength, "bp"), collapse = "")] <- tapply(peakcuts$contribution, overlaps@from, sum)
        }
      }

      # Calculate the number of PWM matches in the 100bp and 400bp windows centered
      # on the considered PWM matches
      if (length(ncol(Pwm[[i]]))>0){
       GenomicRanges::mcols(AnnotatedMatches)["Matches_20bp"] <- NULL
       GenomicRanges::mcols(AnnotatedMatches)["Matches_1000bp"] <- (GenomicRanges::mcols(AnnotatedMatches)["Matches_1000bp"][,1]*1000)/ncol(Pwm[[i]])
       GenomicRanges::mcols(AnnotatedMatches)["Matches_400bp"] <- (GenomicRanges::mcols(AnnotatedMatches)["Matches_400bp"][,1]*400)/ncol(Pwm[[i]])
      }

      #Scale the predictive features
      GenomicRanges::mcols(AnnotatedMatches)[,-(which(colnames(GenomicRanges::mcols(AnnotatedMatches)) %in% c("ClosestTSS", "ClosestTTS")))] <-
        apply(GenomicRanges::mcols(AnnotatedMatches)[,-(which(colnames(GenomicRanges::mcols(AnnotatedMatches)) %in% c("ClosestTSS", "ClosestTTS")))],
                                                  2,
                                                  function(x){x <- as.numeric(x);
                                                              if (length(x[x>0])>0){
                                                                ## Outliers detection
                                                                if (stats::IQR(x[x>0], na.rm = TRUE) > 0){
                                                                  outlier_threshold_low <- stats::median(x[x>0], na.rm = TRUE)-1.5*stats::IQR(x[x>0], na.rm = TRUE)
                                                                  outlier_threshold_high <- stats::median(x[x>0], na.rm = TRUE)+1.5*stats::IQR(x[x>0], na.rm = TRUE)
                                                                  outliers_low <- which(x < outlier_threshold_low)
                                                                  x[outliers_low] <- outlier_threshold_low
                                                                  outliers_high <- which(x > outlier_threshold_high)
                                                                  x[outliers_high] <- outlier_threshold_high
                                                                  x <- (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
                                                                  } else {
                                                                    x[x>1] <- 1                                                            }
                                                                };
                                                               return(x)})

      # Export the annotated matches
      print(summary(as.data.frame(AnnotatedMatches)))
      path_considered <-paste0(getwd(), "/", names(Pwm)[i], "_", paste0(unlist(strsplit(as.character(Sys.time()), ":"))[c(2,3)], collapse = "_"), "_annotations.tsv")
      readr::write_tsv(as.data.frame(AnnotatedMatches), path = path_considered)
      paths <- c(paths, path_considered)
    } else {
      message(paste0("Warning: No match can cross the p-value threshold for: "), names(Matches)[i])
      excluded <- c(excluded, names(Matches)[i])
    }
  }
  if(length(paths) == 0){

  } else {
    names(paths) <- names(Pwm)[!(names(Pwm) %in% excluded)]
  }
  return(paths)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Implements a fast method of pattern-matching
#========================================================================================#
matrixScan <- function(pwm, genome, pval_threshold = 0.001){
  NbChromosomes <- length(genome)
  Start <- runif(200, 1, Biostrings::width(genome)[1] - 5000)
  Start <- Start[!(is.na(Start))]
  SequencesSample <- Biostrings::Views(genome[[1]], start = Start,
                                       end = Start + 5000)
  SequencesSample <- list(SequencesSample)
  if (NbChromosomes > 1) {
    for (i in seq(2, NbChromosomes)) {
      Start <- runif(200, 1, Biostrings::width(genome)[i] -
                       5000)
      Start <- Start[!(is.na(Start))]
      SequencesSample[[i]] <- Biostrings::Views(genome[[i]],
                                                start = Start, end = Start + 5000)
    }
  }
  RandomFwScan <- lapply(SequencesSample, Biostrings::matchPWM,
                         pwm = pwm, min.score = "0%", with.score = TRUE)
  RandomRwScan <- lapply(SequencesSample, Biostrings::matchPWM,
                         pwm = Biostrings::reverseComplement(pwm), min.score = "0%",
                         with.score = TRUE)
  RandomScores <- rtracklayer::mcols(RandomFwScan[[1]])
  if (length(genome) == 1) {
    init = 1
  } else {
    init = 2
  }
  for (i in seq(init, length(genome))) {
    RandomScores <- rbind(RandomScores, rtracklayer::mcols(RandomFwScan[[i]]))
  }
  for (i in seq_along(genome)) {
    RandomScores <- rbind(RandomScores, rtracklayer::mcols(RandomRwScan[[i]]))
  }
  RandomScores <- RandomScores@listData$score
  DistributionScores <- graphics::hist(RandomScores, breaks = 1e+06,
                                       plot = FALSE)
  PvaluesScores <- cumsum(DistributionScores$counts[seq(length(DistributionScores$counts),
                                                        1, -1)])/sum(DistributionScores$counts)
  PvaluesCentered <- abs(PvaluesScores - pval_threshold)
  ScoreThreshold <- DistributionScores$breaks[length(DistributionScores$breaks) -
                                                which.min(PvaluesCentered) + 1]
  ScanFw <- lapply(genome, Biostrings::matchPWM, pwm = pwm,
                   min.score = ScoreThreshold, with.score = TRUE)
  ScanRw <- lapply(genome, Biostrings::matchPWM, pwm = Biostrings::reverseComplement(pwm),
                   min.score = ScoreThreshold, with.score = TRUE)
  if ((sum(unlist(lapply(ScanFw, length))) + sum(unlist(lapply(ScanRw, length)))) < 10000){
    ScanFw <- lapply(genome, Biostrings::matchPWM, pwm = pwm, with.score = TRUE, min.score = 0.6*Biostrings::maxScore(pwm))
    ScanRw <- lapply(genome, Biostrings::matchPWM, pwm = Biostrings::reverseComplement(pwm), with.score = TRUE, min.score = 0.6*Biostrings::maxScore(pwm))
  }
  NbMatchesFw <- unlist(lapply(ScanFw, length))
  NbMatchesRw <- unlist(lapply(ScanRw, length))
  seqnamesMatches <- S4Vectors::Rle(rep(names(ScanFw), 2),
                                    c(NbMatchesFw, NbMatchesRw))
  rangesMatches <- IRanges::IRanges(ScanFw[[1]]@ranges)
  if (length(genome) > 1) {
    for (i in seq(2, length(ScanFw))) {
      rangesMatches <- c(rangesMatches, IRanges::IRanges(ScanFw[[i]]@ranges))
    }
  }
  for (i in seq_along(ScanRw)) {
    rangesMatches <- c(rangesMatches, IRanges::IRanges(ScanRw[[i]]@ranges))
  }
  strandMatches <- S4Vectors::Rle(c("+", "-"), c(sum(NbMatchesFw),
                                                 sum(NbMatchesRw)))
  scoreMatches <- rtracklayer::mcols(ScanFw[[1]])
  if (length(genome) > 1) {
    for (i in seq(2, length(ScanFw))) {
      scoreMatches <- rbind(scoreMatches, rtracklayer::mcols(ScanFw[[i]]))
    }
  }
  for (i in seq_along(ScanRw)) {
    scoreMatches <- rbind(scoreMatches, rtracklayer::mcols(ScanRw[[i]]))
  }
  names(scoreMatches) <- "matchScore"
  if ((nrow(rtracklayer::mcols(ScanRw[[i]]))) == 0){
    scoreMatches <- data.table::data.table(macthScore = rep(1,
                                                            (sum(unlist(lapply(ScanFw, length))) + sum(unlist(lapply(ScanRw, length))))))
    kept <- seq(1, length(scoreMatches))
    matchLogPval <- rep(log10(pval_threshold), length(kept))
    scores <- scoreMatches[,1]
  } else {
    x <- c(round(min(scoreMatches[, 1]), 5), round(RandomScores,
                                                   5), round(max(scoreMatches[, 1]), 5))
    if (length(which(is.infinite(x)))) {
      x <- x[-which(is.infinite(x))]
    }
    DistributionScores <- graphics::hist(x, breaks = c(0, seq(min(x,
                                                                  na.rm = TRUE) - pval_threshold/20, max(x, na.rm = TRUE) + pval_threshold/20,
                                                              pval_threshold/100)), ylim = c(0, length(x)), plot = FALSE)
    PvaluesScores <- cumsum(DistributionScores$counts[seq(length(DistributionScores$counts),
                                                          1, -1)])/sum(DistributionScores$counts)
    PvaluesTable <- data.frame(score = DistributionScores$breaks[2:(length(DistributionScores$breaks) -
                                                                      1)], p.value = PvaluesScores[seq((length(PvaluesScores) -
                                                                                                          1), 1, -1)])
    PvaluesTable <- rbind(PvaluesTable, c(round(max(scoreMatches[,
                                                                 1]), 4), 0))
    PvaluesTable <- rbind(PvaluesTable, c(round(min(scoreMatches[,
                                                                 1]), 4), max(PvaluesTable$p.value)))
    scores <- as.factor(round(scoreMatches[, 1], 4))
    PvaluesTable <- PvaluesTable[PvaluesTable$score %in% scores,
                                 ]
    PvaluesTable <- PvaluesTable[!duplicated(PvaluesTable$score),]
    kept <- which(scores %in% PvaluesTable$score)
    scores <- scores[scores %in% PvaluesTable$score]
    scores <- droplevels(scores)
    levels(scores) <- PvaluesTable$p.value
    matchLogPval <- log10(as.numeric(as.character(scores)))
    matchLogPval[is.infinite(matchLogPval)] <- min(matchLogPval[!(is.infinite(matchLogPval))])
  }
  if ((min(matchLogPval) >= log10(pval_threshold)) | length(kept) == 0){
    ScanFw <- lapply(genome, Biostrings::matchPWM, pwm = pwm, with.score = TRUE, min.score = 0.8*Biostrings::maxScore(pwm))
    ScanRw <- lapply(genome, Biostrings::matchPWM, pwm = Biostrings::reverseComplement(pwm), with.score = TRUE, min.score = 0.8*Biostrings::maxScore(pwm))
    
    NbMatchesFw <- unlist(lapply(ScanFw, length))
    NbMatchesRw <- unlist(lapply(ScanRw, length))
    seqnamesMatches <- S4Vectors::Rle(rep(names(ScanFw), 2),
                                      c(NbMatchesFw, NbMatchesRw))
    rangesMatches <- IRanges::IRanges(ScanFw[[1]]@ranges)
    if (length(genome) > 1) {
      for (i in seq(2, length(ScanFw))) {
        rangesMatches <- c(rangesMatches, IRanges::IRanges(ScanFw[[i]]@ranges))
      }
    }
    for (i in seq_along(ScanRw)) {
      rangesMatches <- c(rangesMatches, IRanges::IRanges(ScanRw[[i]]@ranges))
    }
    strandMatches <- S4Vectors::Rle(c("+", "-"), c(sum(NbMatchesFw),
                                                   sum(NbMatchesRw)))
    scoreMatches <- rtracklayer::mcols(ScanFw[[1]])
    if (length(genome) > 1) {
      for (i in seq(2, length(ScanFw))) {
        scoreMatches <- rbind(scoreMatches, rtracklayer::mcols(ScanFw[[i]]))
      }
    }
    for (i in seq_along(ScanRw)) {
      scoreMatches <- rbind(scoreMatches, rtracklayer::mcols(ScanRw[[i]]))
    }
    names(scoreMatches) <- "matchScore"
    matchLogPval <- log10(pval_threshold)
    ScanTrack <- GenomicRanges::GRanges(seqnames = seqnamesMatches,
                                        ranges = rangesMatches, strand = strandMatches, scoreMatches,
                                        matchLogPval = matchLogPval)
    overlaps <- GenomicRanges::findOverlaps(ScanTrack[GenomicRanges::strand(ScanTrack) == "+"], ScanTrack[GenomicRanges::strand(ScanTrack) == "-"],
                                            ignore.strand = TRUE)
    ScanTrack <- ScanTrack[-(length(ScanTrack[GenomicRanges::strand(ScanTrack) == "+"])+overlaps@to)]
  } else {

    ScanTrack <- GenomicRanges::GRanges(seqnames = seqnamesMatches[kept],
                                        ranges = rangesMatches[kept,], strand = strandMatches[kept], matchScore = scores,
                                        matchLogPval = matchLogPval
                                        )
    overlaps <- GenomicRanges::findOverlaps(ScanTrack[GenomicRanges::strand(ScanTrack) == "+"], ScanTrack[GenomicRanges::strand(ScanTrack) == "-"],
                                            ignore.strand = TRUE)
    if (length(overlaps@to) > 0){
    ScanTrack <- ScanTrack[-(length(ScanTrack[GenomicRanges::strand(ScanTrack) == "+"])+overlaps@to)]
    }
  }
  return(ScanTrack)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Allows to determine the transcript whose the TSS (TTS) is the closest to PWM matches
### and to calculte the distance to its TSS (TTS)
#========================================================================================#
getClosest <- function(GRanges_data, StructureTracks, modeling = TRUE,
                       strand_as_feature = FALSE,
                       TSS = TRUE){
  #Distance to TSS or TTS?
  if (TSS == TRUE){
    TSS_Track <- StructureTracks[[length(StructureTracks)]]
  } else {
    TSS_Track <- StructureTracks[[length(StructureTracks) - 1]]
  }
  GRanges_data <- GRanges_data[as.character(GenomicRanges::seqnames(GRanges_data)) %in%
                                 as.character(GenomeInfoDb::seqlevels(TSS_Track))]
  ClosestTSS_Indices <- GenomicRanges::nearest(GRanges_data,
                                               TSS_Track, ignore.strand = TRUE)
  ClosestTSS_Data <- TSS_Track[ClosestTSS_Indices, ]
  ClosestTSS <- rtracklayer::mcols(TSS_Track)[ClosestTSS_Indices,
                                              1]
  DistToClosestTSS <- GenomicRanges::distanceToNearest(GRanges_data,
                                                       TSS_Track, ignore.strand = TRUE)
  DistToClosestTSS <- rtracklayer::mcols(DistToClosestTSS)[,
                                                           1]
  DistToClosestTSS[as.logical((GenomicRanges::start(ClosestTSS_Data) >
                                 GenomicRanges::start(GRanges_data)) & (GenomicRanges::strand(ClosestTSS_Data) ==
                                                                          "+"))] <- -DistToClosestTSS[as.logical((GenomicRanges::start(ClosestTSS_Data) >
                                                                                                                    GenomicRanges::start(GRanges_data)) & (GenomicRanges::strand(ClosestTSS_Data) ==
                                                                                                                                                             "+"))]
  DistToClosestTSS[as.logical((GenomicRanges::start(ClosestTSS_Data) <
                                 GenomicRanges::start(GRanges_data)) & (GenomicRanges::strand(ClosestTSS_Data) ==
                                                                          "-"))] <- -DistToClosestTSS[as.logical((GenomicRanges::start(ClosestTSS_Data) <
                                                                                                                    GenomicRanges::start(GRanges_data)) & (GenomicRanges::strand(ClosestTSS_Data) ==
                                                                                                                                                             "-"))]
  if (TSS ==  TRUE){
    if (strand_as_feature) {
      if (modeling) {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          DistToClosestTSS = DistToClosestTSS, TranscriptOrientation = as.character(GenomicRanges::strand(ClosestTSS_Data)))
      }
      else {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          ClosestTSS = ClosestTSS, DistToClosestTSS = DistToClosestTSS,
                          TranscriptOrientation = as.character(GenomicRanges::strand(ClosestTSS_Data)))
      }
    } else {
      if (modeling) {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          DistToClosestTSS = DistToClosestTSS)
      }
      else {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          ClosestTSS = ClosestTSS, DistToClosestTSS = DistToClosestTSS)
      }
    }
  }
  else {
    if (strand_as_feature) {
      if (modeling) {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          DistToClosestTTS = DistToClosestTSS, TranscriptOrientation = as.character(GenomicRanges::strand(ClosestTSS_Data)))
      }
      else {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          ClosestTTS = ClosestTSS, DistToClosestTTS = DistToClosestTSS,
                          TranscriptOrientation = as.character(GenomicRanges::strand(ClosestTSS_Data)))
      }
    } else {
      if (modeling) {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          DistToClosestTTS = DistToClosestTSS)
      }
      else {
        MetaData <- cbind(as.data.frame(rtracklayer::mcols(GRanges_data)),
                          ClosestTTS = ClosestTSS, DistToClosestTTS = DistToClosestTSS)
      }
    }
  }
  GRanges_data <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(GRanges_data),
                                         ranges = IRanges::IRanges(start = GenomicRanges::start(GRanges_data),
                                                                   end = GenomicRanges::end(GRanges_data)), strand = GenomicRanges::strand(GRanges_data),
                                         MetaData)
  return(GRanges_data)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Allows to quickly download plant genome sequences. Function modified from biomaRt
#========================================================================================#
getsequence <- function (organism, release = NULL, type = "dna", id.type = "toplevel"){
  if (!is.element(type, c("dna", "cds", "pep", "ncrna"))) 
    stop("Please a 'type' argument supported by this function: \n                 'dna', 'cds', 'pep', 'ncrna'.")
  name <- NULL
  if (!suppressMessages(is.available.genome(organism = organism, 
                                            db = "ensemblgenomes", details = FALSE))) {
    warning("Unfortunately organism '", organism, "' is not available at ENSEMBLGENOMES. ", 
            "Please check whether or not the organism name is typed correctly or try db = 'ensembl'.", 
            " Thus, download of this species has been omitted. ", 
            call. = FALSE)
    return(FALSE)
  }
  else {
    taxon_id <- assembly <- accession <- NULL
    new.organism <- stringr::str_to_lower(stringr::str_replace_all(organism, 
                                                                   " ", "_"))
    ensembl_summary <- suppressMessages(as.data.frame(unlist(is.available.genome(organism = organism, 
                                                                                 db = "ensemblgenomes", details = TRUE))))
    if (nrow(ensembl_summary) == 0) {
      message("Unfortunately, organism '", organism, "' does not exist in this database. Could it be that the organism name is misspelled? Thus, download has been omitted.")
      return(FALSE)
    }
    new.organism <- paste0(stringr::str_to_upper(stringr::str_sub(ensembl_summary["name",], 
                                                                  1, 1)), stringr::str_sub(ensembl_summary["name",], 
                                                                                           2, nchar(ensembl_summary["name",])))
  }
  get.org.info <- ensembl_summary
  rest_url <- paste0("http://rest.ensembl.org/info/assembly/", 
                     ensembl_summary["name",], "?content-type=application/json")
  rest_api_status <- test_url_status(url = rest_url, organism = organism)
  if (is.logical(rest_api_status)) {
    return(FALSE)
  }
  else {
    release_api <- jsonlite::fromJSON("http://rest.ensembl.org/info/eg_version?content-type=application/json")
    if (!is.null(release)) {
      if (!is.element(release, seq_len(as.integer(release_api)))) 
        stop("Please provide a release number that is supported by ENSEMBLGENOMES.", 
             call. = FALSE)
    }
    if (is.null(release)) 
      core_path <- "http://ftp.ensemblgenomes.org/pub/current/"
    if (!is.null(release)) 
      core_path <- paste0("http://ftp.ensemblgenomes.org/pub/current", 
                          release, "/")
    ensembl.qry <- paste0(core_path, stringr::str_to_lower(stringr::str_replace(get.org.info$division[1], 
                                                                                "Ensembl", "")), "plants/fasta/",
                          stringr::str_to_lower(ensembl_summary["name",]), 
                          "/", type, "/", paste0(ensembl_summary["url_name",], ".", rest_api_status$default_coord_system_version, 
                                                 ".", "dna_rm", ifelse(id.type == "none", "", "."), 
                                                 ifelse(id.type == "none", "", id.type), ".fa.gz"))
    
    if (file.exists(file.path(path, paste0(new.organism, 
                                           ".", rest_api_status$default_coord_system_version, 
                                           ".", type, ifelse(id.type == "none", "", "."), ifelse(id.type == 
                                                                                                 "none", "", id.type), ".fa.gz")))) {
      path <- getwd()
      message("File ", file.path(path, paste0(new.organism, 
                                              ".", rest_api_status$default_coord_system_version, 
                                              ".", type, ifelse(id.type == "none", "", "."), 
                                              ifelse(id.type == "none", "", id.type), ".fa.gz")), 
              " exists already. Thus, download has been skipped.")
    }
    else {
      path <- getwd()
      file.fasta <- file.path(path, 
                              paste0(ensembl_summary["name",], ".", rest_api_status$default_coord_system_version, 
                                     ".", type, ifelse(id.type == "none", "", "."), 
                                     ifelse(id.type == "none", "", id.type), ".fa.gz"))
      utils::download.file(url = ensembl.qry, destfile = file.fasta)
    }
    return(file.fasta)
  }
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Checks whether the genome sequence of the considered organism can be quickly
### downloaded with getsequence. Function modified from biomaRt
#========================================================================================#
is.available.genome <- function (db = "refseq", organism, details = FALSE) 
{
  if (is.element(db, c("refseq", "genbank"))) {
    if (file.exists(file.path(tempdir(), paste0("AssemblyFilesAllKingdoms_", 
                                                db, ".txt")))) {
      suppressWarnings(AssemblyFilesAllKingdoms <- readr::read_tsv(file.path(tempdir(), 
                                                                             paste0("AssemblyFilesAllKingdoms_", db, ".txt")), 
                                                                   col_names = c("assembly_accession", "bioproject", 
                                                                                 "biosample", "wgs_master", "refseq_category", 
                                                                                 "taxid", "species_taxid", "organism_name", 
                                                                                 "infraspecific_name", "isolate", "version_status", 
                                                                                 "assembly_level", "release_type", "genome_rep", 
                                                                                 "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm", 
                                                                                 "paired_asm_comp", "ftp_path", "excluded_from_refseq"), 
                                                                   comment = "#", col_types = readr::cols(assembly_accession = readr::col_character(), 
                                                                                                          bioproject = readr::col_character(), biosample = readr::col_character(), 
                                                                                                          wgs_master = readr::col_character(), refseq_category = readr::col_character(), 
                                                                                                          taxid = readr::col_integer(), species_taxid = readr::col_integer(), 
                                                                                                          organism_name = readr::col_character(), infraspecific_name = readr::col_character(), 
                                                                                                          isolate = readr::col_character(), version_status = readr::col_character(), 
                                                                                                          assembly_level = readr::col_character(), release_type = readr::col_character(), 
                                                                                                          genome_rep = readr::col_character(), seq_rel_date = readr::col_date(), 
                                                                                                          asm_name = readr::col_character(), submitter = readr::col_character(), 
                                                                                                          gbrs_paired_asm = readr::col_character(), 
                                                                                                          paired_asm_comp = readr::col_character(), 
                                                                                                          ftp_path = readr::col_character(), excluded_from_refseq = readr::col_character())))
    }
    else {
      kgdoms <- biomartr::getKingdoms(db = db)(db = db)
      storeAssemblyFiles <- vector("list", length(kgdoms))
      for (i in seq_along(kgdoms)) {
        storeAssemblyFiles[i] <- list(getSummaryFile(db = db, 
                                                     kingdom = kgdoms[i]))
      }
      AssemblyFilesAllKingdoms <- dplyr::bind_rows(storeAssemblyFiles)
      readr::write_tsv(AssemblyFilesAllKingdoms, file.path(tempdir(), 
                                                           paste0("AssemblyFilesAllKingdoms_", db, ".txt")))
    }
    organism_name <- assembly_accession <- taxid <- NULL
    orgs <- stringr::str_replace_all(AssemblyFilesAllKingdoms$organism_name, 
                                     "\\(", "")
    orgs <- stringr::str_replace_all(orgs, "\\)", "")
    AssemblyFilesAllKingdoms <- dplyr::mutate(AssemblyFilesAllKingdoms, 
                                              organism_name = orgs)
    organism <- stringr::str_replace_all(organism, "\\(", 
                                         "")
    organism <- stringr::str_replace_all(organism, "\\)", 
                                         "")
    if (!is.taxid(organism)) {
      FoundOrganism <- dplyr::filter(AssemblyFilesAllKingdoms, 
                                     stringr::str_detect(organism_name, organism) | 
                                       assembly_accession == organism)
    }
    else {
      FoundOrganism <- dplyr::filter(AssemblyFilesAllKingdoms, 
                                     taxid == as.integer(organism))
    }
    if (nrow(FoundOrganism) == 0) {
      message("Unfortunatey, no entry for '", organism, 
              "' was found in the '", db, "' database. ", 
              "Please consider specifying ", paste0("'db = ", 
                                                    dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                     "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
              "' to check whether '", organism, "' is available in these databases.")
      return(FALSE)
    }
    if (nrow(FoundOrganism) > 0) {
      if (!details) {
        if (all(FoundOrganism$refseq_category == "na")) {
          message("Only a non-reference genome assembly is available for '", 
                  organism, "'.", " Please make sure to specify the argument 'reference = FALSE' when running any get*() function.")
        }
        else {
          message("A reference or representative genome assembly is available for '", 
                  organism, "'.")
        }
        if (nrow(FoundOrganism) > 1) {
          message("More than one entry was found for '", 
                  organism, "'.", " Please consider to run the function 'is.available.genome()' and specify 'is.available.genome(organism = ", 
                  organism, ", db = ", db, ", details = TRUE)'.", 
                  " This will allow you to select the 'assembly_accession' identifier that can then be ", 
                  "specified in all get*() functions.")
        }
        return(TRUE)
      }
      if (details) {
        if (all(FoundOrganism$refseq_category == "na")) {
          message("Only a non-reference genome assembly is available for '", 
                  organism, "'.", " Please make sure to specify the argument 'reference = FALSE' when running any get*() function.")
        }
        return(FoundOrganism)
      }
    }
  }
  if (db == "ensembl") {
    name <- accession <- accession <- assembly <- taxon_id <- NULL
    new.organism <- stringr::str_replace_all(organism, " ", 
                                             "_")
    ensembl.available.organisms <- get.ensembl.info()
    ensembl.available.organisms <- dplyr::filter(ensembl.available.organisms, 
                                                 !is.na(assembly))
    if (!is.taxid(organism)) {
      selected.organism <- dplyr::filter(ensembl.available.organisms, 
                                         stringr::str_detect(name, stringr::str_to_lower(new.organism)) | 
                                           accession == organism, !is.na(assembly))
    }
    else {
      selected.organism <- dplyr::filter(ensembl.available.organisms, 
                                         taxon_id == as.integer(organism), !is.na(assembly))
    }
    if (!details) {
      if (nrow(selected.organism) == 0) {
        message("Unfortunatey, no entry for '", organism, 
                "' was found in the '", db, "' database. ", 
                "Please consider specifying ", paste0("'db = ", 
                                                      dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                       "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
                "' to check whether '", organism, "' is available in these databases.")
        return(FALSE)
      }
      if (nrow(selected.organism) > 0) {
        message("A reference or representative genome assembly is available for '", 
                organism, "'.")
        if (nrow(selected.organism) > 1) {
          message("More than one entry was found for '", 
                  organism, "'.", " Please consider to run the function 'is.available.genome()' and specify 'is.available.genome(organism = ", 
                  organism, ", db = ", db, ", details = TRUE)'.", 
                  " This will allow you to select the 'assembly_accession' identifier that can then be ", 
                  "specified in all get*() functions.")
        }
        return(TRUE)
      }
    }
    if (details) 
      return(selected.organism)
  }
  if (db == "ensemblgenomes") {
    new.organism <- stringr::str_replace_all(organism, " ", 
                                             "_")
    name <- accession <- assembly <- taxon_id <- NULL
    selected.organism <- getensemblgenome.info(new.organism)
    if (!details) {
      if (length(selected.organism) == 0) {
        message("Unfortunatey, no entry for '", organism, 
                "' was found in the '", db, "' database. ", 
                "Please consider specifying ", paste0("'db = ", 
                                                      dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                       "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
                "' to check whether '", organism, "' is available in these databases.")
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
    if (details) 
      return(selected.organism)
  }
  if (db == "uniprot") {
    if (is.taxid(organism)) {
      unipreot_rest_url <- paste0("https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&taxid=", 
                                  as.integer(organism))
      rest_status_test <- curl::curl_fetch_memory(unipreot_rest_url)
      if (rest_status_test$status_code != 200) {
        message("Something went wrong when trying to access the API 'https://www.ebi.ac.uk/proteins/api/proteomes'.", 
                " Sometimes the internet connection isn't stable and re-running the function might help. Otherwise, could there be an issue with the firewall? ", 
                "Is it possbile to access the homepage 'https://www.ebi.ac.uk/' through your browser?")
      }
      uniprot_species_info <- tibble::as_tibble(jsonlite::fromJSON(unipreot_rest_url))
    }
    else {
      organism_new <- stringr::str_replace_all(organism, 
                                               " ", "%20")
      unipreot_rest_url_name <- paste0("https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&name=", 
                                       organism_new)
      rest_status_test_name <- curl::curl_fetch_memory(unipreot_rest_url_name)
      unipreot_rest_url_upid <- paste0("https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&upid=", 
                                       organism)
      rest_status_test_upid <- curl::curl_fetch_memory(unipreot_rest_url_upid)
      if ((rest_status_test_upid$status_code != 200) & 
          (rest_status_test_name$status_code != 200)) {
        message("Something went wrong when trying to access the API 'https://www.ebi.ac.uk/proteins/api/proteomes'.", 
                " Sometimes the internet connection isn't stable and re-running the function might help. Otherwise, could there be an issue with the firewall? ", 
                "Is it possbile to access the homepage 'https://www.ebi.ac.uk/' through your browser?")
      }
      if (rest_status_test_name$status_code == 200) {
        uniprot_species_info <- tibble::as_tibble(jsonlite::fromJSON(unipreot_rest_url_name))
      }
      if (rest_status_test_upid$status_code == 200) {
        uniprot_species_info <- tibble::as_tibble(jsonlite::fromJSON(unipreot_rest_url_upid))
      }
    }
    if (!details) {
      if (nrow(uniprot_species_info) == 0) {
        message("Unfortunatey, no entry for '", organism, 
                "' was found in the '", db, "' database. ", 
                "Please consider specifying ", paste0("'db = ", 
                                                      dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                       "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
                "' to check whether '", organism, "' is available in these databases.")
        return(FALSE)
      }
      if (nrow(uniprot_species_info) > 0) {
        message("A reference or representative genome assembly is available for '", 
                organism, "'.")
        if (nrow(uniprot_species_info) > 1) {
          message("More than one entry was found for '", 
                  organism, "'.", " Please consider to run the function 'is.available.genome()' and specify 'is.available.genome(organism = ", 
                  organism, ", db = ", db, ", details = TRUE)'.", 
                  " This will allow you to select the 'assembly_accession' identifier that can then be ", 
                  "specified in all get*() functions.")
        }
        return(TRUE)
      }
    }
  }
  if (details) 
    return(uniprot_species_info)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Outputs information about the storage of the genome sequence of the plant species
### of interest on Ensembl genomes.
#========================================================================================#
getensemblgenome.info <- function(organism){
  server <- "https://rest.ensembl.org"
  ext <- paste0("/info/genomes/", organism, "?")
  r <- httr::GET(paste(server, ext, sep = ""), httr::content_type("application/json"))
  httr::stop_for_status(r)

  dataorganism <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))
  return(dataorganism)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Hidden functions from biomartr
#========================================================================================#
get.ensembl.info <- function (update = FALSE) 
{
  if (file.exists(file.path(tempdir(), "ensembl_info.tsv")) && 
      !update) {
    suppressWarnings(ensembl.info <- readr::read_tsv(file.path(tempdir(), 
                                                               "ensembl_info.tsv"), col_names = TRUE, col_types = readr::cols(division = readr::col_character(), 
                                                                                                                              taxon_id = readr::col_integer(), name = readr::col_character(), 
                                                                                                                              release = readr::col_integer(), display_name = readr::col_character(), 
                                                                                                                              accession = readr::col_character(), common_name = readr::col_character(), 
                                                                                                                              assembly = readr::col_character())))
  }
  else {
    rest_url <- "http://rest.ensembl.org/info/species?content-type=application/json"
    rest_api_status <- curl::curl_fetch_memory(rest_url)
    if (rest_api_status$status_code != 200) {
      message("The API 'http://rest.ensembl.org' does not seem to\n                respond or work properly. Is the homepage 'http://rest.ensembl.org' currently available?", 
              " Could it be that there is a firewall issue on your side? Please re-run the function and check if it works now.")
    }
    ensembl.info <- tibble::as_tibble(jsonlite::fromJSON(rest_url)$species)
    aliases <- groups <- NULL
    ensembl.info <- dplyr::select(ensembl.info, -aliases, 
                                  -groups)
    readr::write_tsv(ensembl.info, file.path(tempdir(), 
                                             "ensembl_info.tsv"))
  }
  return(ensembl.info)
}

is.taxid <- function (x) 
{
  return(stringr::str_count(x, "[:digit:]") == nchar(x))
}

test_url_status <- function (url, organism) 
{
  test_status <- curl::curl_fetch_memory(url)
  if (test_status$status_code == 200) {
    json.qry.info <- jsonlite::fromJSON(url)
    return(json.qry.info)
  }
  else {
    message("Something went wrong when trying to access the url: ", 
            url, ". It seems like the organism '", organism, 
            "' does not exist in this database. Could it be that the organism name is misspelled? Thus, download has been omitted.")
    return(FALSE)
  }
}
