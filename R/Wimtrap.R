#' Imports genomic data into the R session
#' @export
#' @description
#' Imports genomic data that will allow to define contextual features around matches with the primary motifs of
#' transcrription factors. Data about gene structures might be optionally automatically downloaded from Biomart.
#' @param organism Binomial name of the organism. Can be set to NULL if you provide
#' the location of the transcription start sites (TSS), transcription termination site (TTS) and structures
#' of the protein-coding genes of the organism (see the arguments  *tss*, *tts* and *genomic_data*)
#' @param genomic_data A named character vector defining the local paths to BED files describing genomic features.
#' The vector has to be named according to the features described by the files indicated. All the data related to
#' the chromatin state have to be specific of the samed condition. The proporties of the BED files are the following,
#' depending on the type of feature: 'numeric' (the score field of the file is fulfilled or empty - if empty, the
#' score will automatically be set to '1') or 'categorical' (the score field of the file is empty while its name field
#' is fulfilled with the names of the categories) features.
#' @param biomart Logical. Should be automatically downloaded through biomart the location of the transcription
#' start sites (TSS), transcription termination site (TTS) and structures of the protein-coding genes of the organism?
#' Default is TRUE.
#' @param tss NULL (by default) or local path to a BED file defining the transcription stat site (TSS), name and
#' orientation of each protein-coding transcript of the organism. The default value allows to download automatically
#' these informations from biomart.
#' @param tts NULL (by default) or local path to a BED file defining the transcription termination site (TTS), name
#' and orientation of each protein-coding transcript of the organism. The default value allows to download automatically
#' these informations from biomart.
#' @param promoterLength Length of the promoters. By default, the promoter is defined as the region spanning the 2000bp
#' upstream of the transcription start site (TSS).
#' @param downstreamLength Length of the downstream regions. By default, the downstream region is defined as the region
#' spanning the 1000bp downstream of the transcription termnination site (TTS).
#' @param proximalLength Length of the proximal promoters. By default, the proximal promoter is defined as the region
#' spanning the 500bp upstream of the transcription start site (TSS).
#' @return A named list of GRanges objects. Each component of the list is named according to the nature of the genomic
#' data. The last two components describe respectively the position of the transcritpion termination site (TTS) and the
#' transcription start site (TSS) of each protein-coding transcripts of the organism considered.
#' @examples
#'## Without automatic download from Biomart of data related to gene structure
#' genomic_data.ex <- c(CE = system.file("conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex2,
#'                                               tss = system.file("tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("tts.bed", package = "Wimtrap"))
#'                                               genomic_data = genomic_data.ex)
#'
#' ##With automatic download from Biomart of data related to gene structure
#' genomic_data.ex <- c(CE = system.file("conserved_elements_example.bed", package = "Wimtrap"),
#'                      DGF = system.file("DGF_example.bed", package = "Wimtrap"),
#'                      DHS = system.file("DHS_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(organism = "Arabidopsis thaliana",
#'                                               biomart = TRUE,
#'                                               genomic_data = genomic_data.ex)
importGenomicData <- function(organism = NULL,
                              genomic_data,
                              biomart = TRUE,
                              tss = NULL,
                              tts = NULL,
                              promoterLength = 2000,
                              downstreamLength = 1000,
                              proximalLength = 500)
{
  #@ Defining all the parameters required by the function
  #@ _inputGenomicData()_
  if (biomart == TRUE) {
    proximalpromoter <- "biomart"
    promoter <- "biomart"
    x5utr <- "biomart"
    cds <- "biomart"
    intron <- "biomart"
    x3utr <- "biomart"
    downstream <- "biomart"
  } else {
    proximalpromoter <- NULL
    promoter <- NULL
    x5utr <- NULL
    cds <- NULL
    intron <- NULL
    x3utr <- NULL
    downstream <- NULL
  }
  if (length(tss) == 0){
    tss <- "biomart"
  }
  if (length(tts) == 0){
    tts <- "biomart"
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
      feature_paths <- genomic_data
      directInput_featureGRanges <- NULL
    } else {
      feature_paths <- NULL
      directInput_featureGRanges <- genomic_data
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

  if (UseOfBiomart) {
    message("Downloading structural genomic ranges...")
    #@ Load the biomart database
    Ensembl <-
      loadEnsembl(organism, NULL, NULL, NULL)
    #@ Query the database
    StructureBiomart <-
      getStructure(Ensembl, promoterLength, downstreamLength, proximalLength)
  } else {

  }

  #@ List the location of each genomic structures in "StructuralFeatures"
  #@ from the specified source (GRanges, source file, biomart)
  for (structure in seq_along(structuralFeatures)) {
    considered <- structuralFeatures[[structure]][1]
    if (length(considered) > 1) {
      #@ The input is a GRanges objectS
    } else {
      if (length(considered) == 1 & is.character(considered)) {
        if (considered == "biomart") {
          structuralFeatures[[structure]] <- StructureBiomart[[structure]]
        } else {
          if (structure < (length(structuralFeatures)-1)) {
            structuralFeatures[[structure]] <-
              importFeatures(considered)[[1]]
          } else {
            structuralFeatures$TTS <-
              rtracklayer::import(structuralFeatures[[(length(structuralFeatures)-1)]][1])
            structuralFeatures$TSS <-
              rtracklayer::import(structuralFeatures[[length(structuralFeatures)]][1])
          }
        }
      } else {

      }
    }
  }

  #@ Change the name of "structuralFeatures" to "StructureTracks"
  if (biomart){
    StructureTracks <- StructureBiomart
  } else {
    StructureTracks <- structuralFeatures
    GenomeInfoDb::seqlevels(StructureTracks$TSS) <- getRiddChr(GenomeInfoDb::seqlevels(StructureTracks$TSS))

    names(StructureTracks) <- c("ProximalPromoter",
                                "Promoter",
                                "X5UTR",
                                "CDS",
                                "Intron",
                                "X3UTR",
                                "Downstream",
                                "TTS",
                                "TSS")
    if (is.null(StructureTracks$Promoter)){
      PartA <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="+"])
      GenomicRanges::start(PartA) <- GenomicRanges::start(PartA) - promoterLength
      PartB <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="-"])
      GenomicRanges::end(PartB) <- GenomicRanges::start(PartB) + promoterLength
      StructureTracks$Promoter <- c(PartA, PartB)
    }
    if (is.null(StructureTracks$Downstream)){
      PartA <- GenomicRanges::granges(StructureTracks$TTS[GenomicRanges::strand(StructureTracks$TTS)=="+"])
      GenomicRanges::end(PartA) <- GenomicRanges::end(PartA) + downstreamLength
      PartB <- GenomicRanges::granges(StructureTracks$TTS[GenomicRanges::strand(StructureTracks$TTS)=="-"])
      GenomicRanges::start(PartB) <- GenomicRanges::start(PartB) - downstreamLength
      StructureTracks$Downstream <- c(PartA, PartB)
    }
    if (is.null(StructureTracks$ProximalPromoter)){
      PartA <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="+"])
      GenomicRanges::start(PartA) <- GenomicRanges::start(PartA) - proximalLength
      PartB <- GenomicRanges::granges(StructureTracks$TSS[GenomicRanges::strand(StructureTracks$TSS)=="-"])
      GenomicRanges::end(PartB) <- GenomicRanges::start(PartB) + proximalLength
      StructureTracks$ProximalPromoter <- c(PartA, PartB)
    }
    fulfilled <- lapply(StructureTracks, is.null)
    StructureTracks <- StructureTracks[unlist(fulfilled)==FALSE]
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


#' Builds a dataset that allows to train a model of transcription factor binding sites (TFBS)
#' @export
#' @importFrom rtracklayer mcols import
#' @importFrom GenomicRanges strand makeGRangesFromDataFrame GRanges reduce as.data.frame start end resize findOverlaps nearest distanceToNearest
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings readDNAStringSet width Views matchPWM reverseComplement letterFrequency
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom utils read.csv
#' @importFrom biomartr getENSEMBLGENOMESInfo getENSEMBLInfo getGenome
#' @importFrom biomaRt listMarts useMart listDatasets getBM
#' @importFrom methods as
#' @importFrom IRanges IRanges resize
#' @importFrom readr write_tsv
#' @importFrom data.table fread data.table
#' @importFrom graphics hist
#' @importFrom stats runif
#' @description
#' `getTFBSdata` outputs a balanced dataset characterizing the genomic context
#' around occurences of primary motifs that have been proved to be either TF-bound
#' or TF-not bound sites according to ChIP-seq data obtained in given experimental conditions.
#' @param  ChIPpeaks A named character vector defining the local paths to BED files encoding the
#' location of ChIP-peaks. The vector is named according to the transcription factors that are described
#' by the files indicated.
#' @param  pfm Path to a file in the raw Jaspar format. This file describes the Pseudo Frequency Matrices (PFMs)
#' of the primary motifs of the transcription factors whose ChIP-seq data are considered.
#' IMPORTANT: In the file indicated, the PFMs have to be named by their cognate transcription
#' factors. These names have to be consistent with those provided through the *ChIP-peaks* argument.
#' *pfm* can be set to NULL if you prefer to use pattern-matching results from elsewhere (see the
#' argument *matches*).
#' @param organism Binomial name of the organism. Can be set to NULL if you provide the genome sequence
#' (see the argument *genome_sequence*)
#' @param genome_sequence "getGenome" (by default) or local path to a FASTA file encoding the genomic sequence
#'  of the organism. The default value allows the automatic download of the genomic sequence (when *organism* is defined)
#'  from ENSEMBL or ENSEMBL GENOMES.
#' @param imported_genomic_data An object output by [importGenomicData()] and that includes data related to the chromatin
#' state that are specific to the condition in which the ChIP-seq data considered have been obtained.
#' @param matches NULL (by default) or a named list of GRanges objects. Each GRanges object defines the location along
#' the genome of the matches with the primary motif of a given transcription factor. It contains also a metadata column
#' named 'matchLogPval'that gives the p-value of the matches. The list input through *matches* has to be named according
#' to the names of the transcription factors considered. These names have to be consistent with those provided through
#' the *ChIP-peaks* argument. The default value allows to perform the pattern-matching analysis with the function encoded
#' by the Wimtrap package.
#' @param strandAsFeature Logical. Should be considered as feature the orientation of the matches in relation to the
#' direction of transcription of the closest transcript? Default is FALSE.
#' @param pvalThreshold P-value threshold to identify the matches with the primary motifs of the transcription factors.
#'  Default is set to 0.001.
#' @param short_window Integer (20 by default). Sets the length of the short-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param medium_window Integer (400 by default). Sets the length of the medium-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param long_window Integer (1000 by default). Sets the length of the long-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @seealso [importGenomicData()] for importing genomic data and [buildTFBSmodel()] to train a predictive model of
#' transcription factor binding sites.
#' @return A data.table describing the location of the potential binding sites in the 5 first columns, the
#' results of pattern-matching (raw score and/or p-value), the genomic features extracted on short-, medium-
#' and long-ranges windows centered on the potential binding sites, the label ('1' = "positive" = "ChIP-validated in
#' a given condition" or '0' = "negative") and the name of the cognate transcription factor.
#' @example
#' genomic_data.ex <- c(CE = system.file("conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex2,
#'                                               tss = system.file("tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("tts.bed", package = "Wimtrap"))
#'                                               genomic_data = genomic_data.ex)
#' TFBSdata.ex <- getTFBSdata(ChIPpeaks = c(PIF3 = system.file("PIF3_example.bed", package = "Wimtrap"),
#'                                          TOC1 = system.file("TOC1_example.bed", package = "Wimtrap")),
#'                            pfm = system.file("pfm_example.jaspar", package = "Wimtrap"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
getTFBSdata <- function(ChIPpeaks,
                        pfm = NULL,
                        organism = NULL,
                        genome_sequence = "getGenome",
                        imported_genomic_data,
                        matches = NULL,
                        strandAsFeature = FALSE,
                        pvalThreshold = 0.001)
{
  #@ Import the source data into the R session
  sourceData <- inputSourceData(ChIPpeaks,
                                pfm,
                                organism,
                                genome_sequence,
                                matches,
                                strandAsFeature)
  #@ Build the dataset
  DataSet <- getCandidatesRegions(sourceData$matches,
                                  sourceData$pwm,
                                  imported_genomic_data[names(imported_genomic_data) %in% c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                                                                                           "X3UTR", "Downstream", "TTS", "TSS")],

                                  imported_genomic_data[!(names(imported_genomic_data) %in% c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                                                                                            "X3UTR", "Downstream", "TTS", "TSS"))],
                                  sourceData$genome,
                                  TFNames = names(sourceData$ChIP),
                                  pvalThreshold,
                                  modeling = TRUE,
                                  sourceData$ChIP)
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
#' is fed by a dataset of potential binding sites of training TFs either labelled as 'positive' or 'negative'
#' (i.e. ChIP-seq 'validated' or 'not-validated' in a given condition).
#' @param TFBSdata A data.table object as output by the [getTFBSdata()] function.
#' @param TF4Validation NULL (by default) or character vector of names of training TFs. If NULL, the model performances will be
#' assessed from a validation dataset obtained by sampling randomly 20% of the potential binding sites of all the training TFs.
#' Otherwise the model will be validated with the data related to the training TFs named in the *TF4Validation* parameter and
#' trained using the remaining training TFs. Please note, that in the last case, the *TFBSdata* table must contain a
#' column named 'TF'.
#' @param model_assessment Logical. If TRUE (by default), 1) the imortance of features is represented and, based on the
#' validation dataset, 2) the confusion matrix is printed and 3) the ROC and AUC are plotted.
#' @return An object of xgb.Booster class.
#' @seealso [getTFBSdata()] for obtaining a training dataset and [predictTFBS()] to predict transcription factor
#' binding site location
#' @example
#' genomic_data.ex <- c(CE = system.file("conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex2,
#'                                               tss = system.file("tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("tts.bed", package = "Wimtrap"))
#'                                               genomic_data = genomic_data.ex)
#' TFBSdata.ex <- getTFBSdata(ChIPpeaks = c(PIF3 = system.file("PIF3_example.bed", package = "Wimtrap"),
#'                                          TOC1 = system.file("TOC1_example.bed", package = "Wimtrap")),
#'                            pfm = system.file("pfm_example.jaspar", package = "Wimtrap"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
#' TFBSmodel.ex <- buildTFBSmodel(TFBSdata = TFBSdata.ex, TF4Validation = "PIF3")

buildTFBSmodel <- function(TFBSdata, TF4Validation = NULL, model_assessment = TRUE){
  TFBSdata <- TFBSdata[,-seq(1,5), with = FALSE]
  #Split the dataset into a training and a validation datasets
  if(is.null(TF4Validation)){
    trainind <- sample(seq(1,nrow(filteredDescr)), nrow(filteredDescr)*0.8)
    testind <- seq(1,nrow(filteredDescr))[seq(1,nrow(filteredDescr)) != trainind]
  } else {
    trainind <- which(TFBSdata$TF %in% TF4Validation)
    testind <- which(!(TFBSdata$TF) %in% TF4Validation)
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

  dtrain <- xgboost::xgb.DMatrix(data = new_tr,label = labels)
  dtest <- xgboost::xgb.DMatrix(data = new_ts,label=ts_label)

  params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)
  xgbcv <- xgboost::xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, print.every.n = 10, early.stop.round = 20, maximize = F)

  xgb1 <- xgboost::xgb.train(params = params, data = dtrain, nrounds = xgbcv$best_iteration, watchlist = list(val=dtest,train=dtrain), print.every.n = 10, early.stop.round = 10, maximize = F , eval_metric = "error")
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
  return(xgb1)
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
#' @importFrom biomaRt listMarts useMart listDatasets getBM
#' @importFrom methods as
#' @importFrom IRanges IRanges resize
#' @importFrom readr write_tsv
#' @importFrom data.table fread data.table
#' @importFrom graphics hist
#' @importFrom stats runif
#' @description
#' Predicts the binding sites of a target transcription factor in a target condition
#' by applying a predictive model on potential binding sites
#' located by a prior pattern-matching analysis and annotated with genomic features extracted at their location.
#' @param TFBSmodel A xgb.Booster object as output by the function [buildTFBSmodel()] or any predictive model
#' that can be used by the function [stats::predict()].
#' @param pfm Path to a file in the raw Jaspar format. This file describes the Pseudo Frequency Matrices (PFMs),
#' among which that of the target Transcription factor.
#' @param targetTF  A character strinf defining the name of the target Transcription factor (the one used in the file
#' indicated by the *pfm* argument).
#' @param organism Binomial name of the organism. Can be set to NULL if you provide the genome sequence (see the argument
#' *genome_sequence*)
#' @param genome_sequence "getGenome" (by default) or local path to a FASTA file encoding the genomic sequence
#'  of the organism. The default value allows the automatic download of the genomic sequence (when *organism* is defined)
#'  from ENSEMBL or ENSEMBL GENOMES.
#' @param imported_genomic_data An object output by [importGenomicData()] and that includes chromatin state data
#' that are specific to the target condition>.
#' @param matches NULL (by default) or a named list of GRanges objects. Each GRanges object defines the location along
#' the genome of the matches with the primary motif of a given transcription factor. It contains also a metadata column
#' named 'matchLogPval'that gives the p-value of the matches. The default value allows to perform the pattern-matching analysis with the function encoded
#' by the Wimtrap package.
#' @param strandAsFeature Logical. Should be considered as feature the orientation of the matches in relation to the
#' direction of transcription of the closest transcript? Default is FALSE.
#' @param pvalThreshold P-value threshold to identify the matches with the primary motifs of the transcription factors.
#'  Default is set to 0.001.
#' @param short_window Integer (20 by default). Sets the length of the short-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param medium_window Integer (400 by default). Sets the length of the medium-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param long_window Integer (1000 by default). Sets the length of the long-ranges window centered on the potential binding
#' sites and on which the genomic features are extracted.
#' @param show_annotation Logical (FALSE by default). Should the annotations of the potential binding sites with genomic
#' features be output?
#' @details  Remark: Set the use of the strand or not as a predictive feature, the pattern-matching p-value threshold
#' and the lengths of the windows on which to extract the genomic feature in agreement as you did to
#' build the model.
#' @seealso [importGenomicData()] for importing genomic data and [buildTFBSmodel()] to train a predictive model of
#' transcription factor binding sites.
#' @return A data.table describing for the potential binding sites the location, the closest gene (relatively to the
#' transcript start site) and the prediction score. Higher the prediction score, higher the chance of the potential
#' binding site to be a binding site of the target condition in the target condition. Consider as minimal threshold,
#' the score of 0.5 and level it up for increasing the specificity. Optionally, the data.table might include the
#' genomic features.
#' genomic_data.ex <- c(CE = system.file("conserved_elements_example.bed", package = "Wimtrap"),
#'                       DGF = system.file("DGF_example.bed", package = "Wimtrap"),
#'                       DHS = system.file("DHS_example.bed", package = "Wimtrap"),
#'                       X5UTR = system.file("x5utr_example.bed", package = "Wimtrap"),
#'                       CDS = system.file("cds_example.bed", package = "Wimtrap"),
#'                       Intron = system.file("intron_example.bed", package = "Wimtrap"),
#'                       X3UTR = system.file("x3utr_example.bed", package = "Wimtrap")
#'                      )
#' imported_genomic_data.ex <- importGenomicData(biomart = FALSE,
#'                                               genomic_data = genomic_data.ex2,
#'                                               tss = system.file("tss_example.bed", package = "Wimtrap"),
#'                                               tts = system.file("tts.bed", package = "Wimtrap"))
#'                                               genomic_data = genomic_data.ex)
#' TFBSdata.ex <- getTFBSdata(ChIPpeaks = c(PIF3 = system.file("PIF3_example.bed", package = "Wimtrap"),
#'                                          TOC1 = system.file("TOC1_example.bed", package = "Wimtrap")),
#'                            pfm = system.file("pfm_example.jaspar", package = "Wimtrap"),
#'                            organism = NULL,
#'                            genome_sequence = system.file("genome_example.fa", package = "Wimtrap"),
#'                            imported_genomic_data = imported_genomic_data.ex)
#' TFBSmodel.ex <- buildTFBSmodel(TFBSdata = TFBSdata.ex, TF4Validation = "PIF3")
#' PIF3BS.predictions <- predictTFBS(TFBSmodel.ex,
#'                                   system.file("pfm_example.jaspar", package = "Wimtrap"),
#'                                   targetTF = "PIF3",
#'                                   genome_sequence = system.file("genome_example.fa", package = "Wimtrap"),
#'                                   imported_genomic_data = imported_genomic_data.ex)
#' ##To get the transcripts whose expression is potentially regulated by PIF3 using as prediction score
#' ##threshold of 0.5 to predict the binding sites in the target condition, do as follows:
#' PIF3_regulated.predictions <- as.character(PIF3BS.predictions$gene[!duplicated(PIF3BS.predictions$prediction.score >= 0.5)])
#' ###If you want to consider only the gene model,
#' ###then do as followq:
#' PIF3_regulated.predictions <- unlist(strsplit(PIF3_regulated.predictions, "[.]"))[seq(1, 2*length(PIF3_regulated.predictions),2)]
#' PIF3_regulated.predictions <- PIF3_regulated.predictions[!duplicated(PIF3_regulated.predictions)]

predictTFBS <- function(TFBSmodel,
                        pfm = NULL,
                        targetTF = NULL,
                        organism = NULL,
                        genome_sequence = "getGenome",
                        imported_genomic_data,
                        matches = NULL,
                        strandAsFeature = FALSE,
                        pvalThreshold = 0.001,
                        short_window = 20,
                        medium_window = 400,
                        long_window = 1000,
                        show_annotations = FALSE){
  #@ Import the source data into the R session
  sourceData <- inputSourceData(NULL,
                                pfm,
                                organism,
                                genome_sequence,
                                matches,
                                strandAsFeature,
                                targetTF)
  #@ Build the dataset
  DataSet <- getCandidatesRegions(sourceData$matches,
                                  sourceData$pwm,
                                  imported_genomic_data[names(imported_genomic_data) %in% c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                                                                                            "X3UTR", "Downstream", "TTS", "TSS")],

                                  imported_genomic_data[!(names(imported_genomic_data) %in% c("ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
                                                                                              "X3UTR", "Downstream", "TTS", "TSS"))],
                                  sourceData$genome,
                                  TFNames = names(sourceData$pwm),
                                  pvalThreshold,
                                  modeling = FALSE,
                                  NULL,
                                  short_window,
                                  medium_window,
                                  long_window)

  TFBSdata <- DataSet[,-seq(1,5), with = FALSE]
  #Pre-processing of the  dataset
  ##Remove columns that do not have to be taken into account
  if (length(grep(pattern = "matchScore", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "matchScore", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
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
  TFBSdata$matchLogPval[which(is.infinite(TFBSdata$matchLogPval))] <- max(TFBSdata$matchLogPval)
  ##Create dummy variables
  TFBSdata <- stats::model.matrix(~.+0, data = TFBSdata)
  NAs <- is.na(as.data.frame(TFBSdata))
  NAs <- apply(NAs, 1, function(x) {if (length(which(x==TRUE)) > 0 ) {return(TRUE)} else {return(FALSE)} })
  TFBSdata <- TFBSdata[!NAs,]
  TFBSdata <- TFBSdata[,colnames(TFBSdata)[colnames(TFBSdata) %in% TFBSmodel$feature_names]]
  TFBSdata <- data.table::as.data.table(TFBSdata)
  TFBSdata <- stats::model.matrix(~.+0,data = TFBSdata)

  TFBSdata <- xgboost::xgb.DMatrix(data = TFBSdata)
  xgbpred <- stats::predict(TFBSmodel,TFBSda)
  if (show_annotations){
    results <- cbind(DataSet, prediction.score = xgbpred)
    colnames(results)[which(colnames(results) == "ClosestTSS")] <- "gene"
  } else {
    results <- cbind(DataSet[,c(seq(1,5), which(colnames(DataSet)=="ClosestTSS")),with=FALSE], prediction.score = xgbpred)
    colnames(results)[which(colnames(results) == "ClosestTSS")] <- "gene"
    }
  return(results)
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
                            strandAsFeature = FALSE,
                            targetTF = NULL)
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
        if (strandAsFeature == FALSE) {
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

  message("Importing ChIP-peaks data and genomic-wide features...")

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
    if(!is.null(targetTF)){
      TFs <- targetTF
    } else {
      TFs <- NULL
    }
  }
  Pwm <- getPwm(pfm_path, directInput_pwmMatrices, organism, TFs)
  Pwm <- Pwm[names(Pwm) %in% TFs]

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
getPwm <- function(pfm_path = NULL, directInput_pwmMatrices = NULL, organism, TFs){
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
  return(list_pfm)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Download the genome sequence from ENSEMBL or ENSEMBL GENOMES
#========================================================================================#
getChromosomes <- function(organism)
{
  biomartr::getENSEMBLGENOMESInfo()
  biomartr::getENSEMBLInfo()
  wd <- getwd()
  LinneName <- organism
  filepath <- tryCatch(biomartr::getGenome(organism = LinneName,
                                           db = "ensembl", path = file.path(wd, "genome")), error = function(e) {
                                             warning("Genome not on the first db. We try with Ensembl Genomes db")
                                           })
  if (filepath == FALSE) {
    filepath <- tryCatch(biomartr::getGenome(organism = LinneName,
                                             db = "ensemblgenomes", path = file.path(wd, "genome")),
                         error = function(e) {
                           return("test")
                         })
  }
  else {
  }
  Genome <- Biostrings::readDNAStringSet(filepath)
  ChrNames <- GenomeInfoDb::seqlevels(Genome)
  SplitChrNames <- unlist(lapply(as.list(ChrNames), base::strsplit,
                                 split = " "))
  FieldNumber <- length(SplitChrNames)/length(ChrNames)
  IndexChrSelected <- 1 + (which(SplitChrNames == "dna:chromosome") -
                             2)/FieldNumber
  ChrNames <- SplitChrNames[which(SplitChrNames == "dna:chromosome") -
                              1]
  ChrNames <- getRiddChr(ChrNames)
  Genome <- Genome[IndexChrSelected]
  names(Genome) <- ChrNames
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
loadEnsembl <- function (organism, host = NULL, biomart = NULL, dataset = NULL)
{
  OrganismsList <- data.frame(organism = c(), host = c(),
                              dataset = c(), mart = c())
  if (missing(host) || is.null(host)) {
    Hosts <- c("ensembl.org", "bacteria.ensembl.org", "plants.ensembl.org",
               "fungi.ensembl.org", "protists.ensembl.org", "ensemblgenomes.org",
               "metazoa.ensembl.org")
  }
  else {
    Hosts = host
  }
  for (Host in Hosts) {
    ListMarts <- tryCatch(biomaRt::listMarts(host = Host),
                          error = function(e) {
                            return("No reachable")
                          })
    if (length(ListMarts) == 1 && ListMarts == "No reachable") {
    }
    else {
      if (missing(biomart) | is.null(biomart)) {
        Biomart <- ListMarts[1, 1]
      }
      else {
        Biomart <- biomart
      }
      ensembl <- biomaRt::useMart(biomart = Biomart, host = Host)
      ListDataSets <- biomaRt::listDatasets(ensembl)
      OrganismsInHost <- data.frame(organism = ListDataSets$description,
                                    host = rep(Host, length(ListDataSets$description)),
                                    dataset = ListDataSets$dataset, mart = rep(Biomart,
                                                                               length(ListDataSets$description)))
      OrganismsList <- rbind(OrganismsList, OrganismsInHost)
      if (Host == "ensembl.org") {
        Biomart <- ListMarts[2, 1]
        ensembl <- biomaRt::useMart(biomart = Biomart,
                                    host = Host)
        ListDataSets <- biomaRt::listDatasets(ensembl)
        OrganismsInHost <- data.frame(organism = ListDataSets$description,
                                      host = rep(Host, length(ListDataSets$description)),
                                      dataset = ListDataSets$dataset, mart = rep(Biomart,
                                                                                 length(ListDataSets$description)))
        OrganismsList <- rbind(OrganismsList, OrganismsInHost)
      }
      else {
      }
    }
  }
  OrganismIndex <- base::grep(organism, OrganismsList$organism,
                              ignore.case = TRUE)
  if (length(OrganismIndex) == 0) {
    warning(base::paste0(c("Biomart and/or host parameters are not appropriate.\nNo dataset related to",
                           organism, "has been found.\nThe current parameters are:\nbiomart:",
                           base::paste(biomart, ","), "\nhost:", base::paste(Hosts,
                                                                             ","))))
    stop()
  }
  biomart <- ifelse(missing(biomart) || is.null(biomart),
                    as.character(OrganismsList$mart[OrganismIndex]), biomart)
  dataset <- ifelse(missing(dataset) || is.null(dataset),
                    OrganismsList$dataset[OrganismIndex], dataset)
  host <- ifelse(missing(host) || is.null(host), as.character(OrganismsList$host[OrganismIndex]),
                 host)
  Ensembl <- biomaRt::useMart(biomart, dataset, host)
  return(Ensembl)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Retrieve the structure of protein-coding genes from biomart database
#========================================================================================#
getStructure <- function (ensembl, promoterLength = 2000, downstreamLength = 1000, proximalLength = 500)
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
                                                                                      1] - promoterLength
  Promoters$transcript_end[Promoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                    1]
  Promoters$transcript_start[Promoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                     -1]
  Promoters$transcript_end[Promoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                   -1] + promoterLength
  colnames(Promoters)[c(2, 3)] <- c("promoter_start", "promoter_end")
  ProximalPromoters <- Transcripts
  ProximalPromoters$transcript_start[ProximalPromoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                      1] - proximalLength
  ProximalPromoters$transcript_end[ProximalPromoters$strand == 1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                    1]
  ProximalPromoters$transcript_start[ProximalPromoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                     -1]
  ProximalPromoters$transcript_end[ProximalPromoters$strand == -1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                   -1] + proximalLength
  colnames(ProximalPromoters)[c(2, 3)] <- c("proximal_promoter_start", "proximal_promoter_end")
  Downstreams <- Transcripts
  Downstreams$transcript_start[Downstreams$strand == 1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                        1]
  Downstreams$transcript_end[Downstreams$strand == 1] <- Transcripts$transcript_end[Transcripts$strand ==
                                                                                      1] + downstreamLength
  Downstreams$transcript_start[Downstreams$strand == -1] <- Transcripts$transcript_start[Transcripts$strand ==
                                                                                           -1] - downstreamLength
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
    chromosome_name <- getRiddChr(chromosome_name)
    Track <- GenomicRanges::GRanges(seqnames = methods::as(chromosome_name,
                                                           "Rle"), ranges = IRanges::IRanges(start = structure[,
                                                                                                               2], end = structure[, 3]), strand = structure[,
                                                                                                                                                             5])
    if (i < (length(Structures)-1)) {
      Track <- GenomicRanges::reduce(Track)
    }
    else {
      if (i == (length(Structures)-1)){
        Track <- GenomicRanges::GRanges(seqnames = methods::as(chromosome_name,
                                                               "Rle"), ranges = IRanges::IRanges(start = structure[,
                                                                                                                   3], end = structure[, 3]), strand = structure[,
                                                                                                                                                                 5], TTS = methods::as(TTS$ensembl_transcript_id,
                                                                                                                                                                                       "Rle"))

      } else {
        Track <- GenomicRanges::GRanges(seqnames = methods::as(chromosome_name,
                                                               "Rle"), ranges = IRanges::IRanges(start = structure[,
                                                                                                                   3], end = structure[, 3]), strand = structure[,
                                                                                                                                                                 5], TSS = methods::as(TSS$ensembl_transcript_id,                                                                                                                                                                                "Rle"))
      }

    }
    Track <- GenomicRanges::as.data.frame(Track)
    Track <- GenomicRanges::makeGRangesFromDataFrame(Track,
                                                     keep.extra.columns = TRUE)
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
                                 pvalThreshold,
                                 modeling,
                                 ChIP_regions,
                                 short_window = 20,
                                 medium_window = 400,
                                 long_window = 1000) {
  paths <- c()
  #@Define the regions that are around a TSS,
  #@from -promoterLength to +downstreamLength
  TSSTrack <- as.data.frame(StructureTracks$TSS)
  TTSTrack <- as.data.frame(StructureTracks$TTS)
  colnames(TSSTrack)[6] <- "name"
  colnames(TTSTrack)[6] <- "name"
  TSSTrack <- TSSTrack[TSSTrack$name %in% TTSTrack$name,]
  TTSTrack <- TTSTrack[TTSTrack$name %in% TSSTrack$name,]
  TSSTrack <- TTSTrack[!(duplicated(TSSTrack$name)),]
  rownames(TSSTrack) <- TSSTrack$name
  TSSTrack <- TSSTrack[TTSTrack$name,]
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

  if (length(directInput_matches) == 0) {
    message("Pattern matching...")
    #@ Scanning of the genome with the PWM
    Pwm <- Pwm[which(names(Pwm) %in% TFNames)]
    Matches <- list()
    for (pwm in names(Pwm)){
      Matches[[pwm]] <- matrixScan(pwm = Pwm[[pwm]], genome = genome, pvalThreshold= pvalThreshold)
    }
  } else {
    if (is.list(directInput_matches)) {
      Matches <- directInput_matches
    } else {
      Matches <- list(directInput_matches)
      Pwm <- names(Matches)
      names(Pwm) <- names(Matches)
    }
  }
  Matches <- lapply(Matches, function(x) {x <- x[rtracklayer::mcols(x)[,2] <= log10(pvalThreshold)]; return(x)})
  for (i in seq_along(Matches)) {
    GenomeInfoDb::seqlevels(Matches[[i]]) <-
      getRiddChr(GenomeInfoDb::seqlevels(Matches[[i]]))
  }



  message("Annotation of the matches...")

  #@ Annotate the matches associated to each TF

  AnnotatedScanResults <- list()
  for (i in seq_along(Matches)) {
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
    GenomicRanges::mcols(AnnotatedMatches)[,"Promoter"] <- 0
    promoterLength <- mean(GenomicRanges::width(StructureTracks$Promoter))
    DownstreamLength <- mean(GenomicRanges::width(StructureTracks$Downstream))
     GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS <=0) ,
                                             5:(which(names(GenomicRanges::mcols(AnnotatedMatches))=="DistToClosestTTS")-2)] <- 0
     GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS <=0 &
                                                     AnnotatedMatches$DistToClosestTSS > -promoterLength), 4] <- 1
     GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS > 0 &
                                                     AnnotatedMatches$DistToClosestTTS < downstreamLength),
                                             (which(names(GenomicRanges::mcols(AnnotatedMatches))=="DistToClosestTTS")-2)] <- 1
     GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS > 0 &
                                                     AnnotatedMatches$DistToClosestTTS < 0),
                                             (which(names(GenomicRanges::mcols(AnnotatedMatches))=="DistToClosestTTS")-2)] <- 0
     GenomicRanges::mcols(AnnotatedMatches)[which(AnnotatedMatches$DistToClosestTSS > 0 &
                                                     apply(GenomicRanges::mcols(AnnotatedMatches)[5:(which(names(GenomicRanges::mcols(AnnotatedMatches))=="DistToClosestTTS")-3)], 1, sum) > 0),
                                             (which(names(GenomicRanges::mcols(AnnotatedMatches))=="DistToClosestTTS")-2)] <- 0


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

    # Identify the true and the false binding sites
    if(modeling){
     GenomicRanges::mcols(AnnotatedMatches)["ChIP.peak"] <- 0
     GenomicRanges::mcols(AnnotatedMatches)[GenomicRanges::findOverlaps(Matches[[i]], ChIP_regions[[names(Pwm)[i]]])@from,"ChIP.peak"] <- 1
    }

    # Calculate the number of PWM matches in the 100bp and 400bp windows centered
    # on the considered PWM matches
   GenomicRanges::mcols(AnnotatedMatches)["Matches_1000bp"] <- (GenomicRanges::mcols(AnnotatedMatches)["Matches_1000bp"][,1]*1000)/ncol(Pwm[[i]])
   GenomicRanges::mcols(AnnotatedMatches)["Matches_400bp"] <- (GenomicRanges::mcols(AnnotatedMatches)["Matches_400bp"][,1]*400)/ncol(Pwm[[i]])

    # Export the annotated matches
    readr::write_tsv(as.data.frame(AnnotatedMatches), path = paste0(names(Pwm)[i], "_annotations.tsv"))
    paths <- c(paths, paste0(names(Pwm)[i], "_annotations.tsv"))
  }
  annotations_files <- paths
  DataSet <- data.frame()
  for (file in annotations_files){
    TF <- unlist(strsplit(file, "_annotations"))[1]
    considered <- data.table::fread(file,
                                    stringsAsFactors = TRUE)
    considered$TF <- TF
    if(modeling){
    NbTrueBs <- nrow(considered[considered$ChIP.peak == 1,])
    DataSet <- rbind(DataSet,
                     considered[considered$ChIP.peak == 1,],
                     considered[sample(which(considered$ChIP.peak == 0), NbTrueBs),])
    } else {
      DataSet <- rbind(DataSet, considered)
    }
  }
  DataSet <- DataSet[, -c("Matches_20bp"), with = FALSE]
  return(DataSet)
}

#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Implements a fast method of pattern-matching
#========================================================================================#
matrixScan <- function (pwm, genome, pvalThreshold = 0.001){
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
  PvaluesCentered <- abs(PvaluesScores - pvalThreshold)
  ScoreThreshold <- DistributionScores$breaks[length(DistributionScores$breaks) -
                                                which.min(PvaluesCentered) + 1]
  ScanFw <- lapply(genome, Biostrings::matchPWM, pwm = pwm,
                   min.score = ScoreThreshold, with.score = TRUE)
  ScanRw <- lapply(genome, Biostrings::matchPWM, pwm = Biostrings::reverseComplement(pwm),
                   min.score = ScoreThreshold, with.score = TRUE)
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
  x <- c(round(min(scoreMatches[, 1]), 5), round(RandomScores,
                                                 5), round(max(scoreMatches[, 1]), 5))
  if (length(which(is.infinite(x)))) {
    x <- x[-which(is.infinite(x))]
  }
  DistributionScores <- graphics::hist(x, breaks = c(0, seq(min(x,
                                                                na.rm = TRUE), max(x, na.rm = TRUE) + pvalThreshold/10,
                                                            pvalThreshold/100)), ylim = c(0, length(x)), plot = FALSE)
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
  ScanTrack <- GenomicRanges::GRanges(seqnames = seqnamesMatches[kept],
                                      ranges = rangesMatches[kept,], strand = strandMatches[kept], scoreMatches[kept,],
                                      matchLogPval)
  overlaps <- GenomicRanges::findOverlaps(ScanTrack[GenomicRanges::strand(ScanTrack) == "+"], ScanTrack[GenomicRanges::strand(ScanTrack) == "-"],
                                          ignore.strand = TRUE)
  ScanTrack <- ScanTrack[-(length(ScanTrack[GenomicRanges::strand(ScanTrack) == "+"])+overlaps@to)]
  return(ScanTrack)
}


#________________________________________________________________________________________#
#HIDDEN FUNCTION
#========================================================================================#
### Allows to determine the transcript whose the TSS (TTS) is the closest to PWM matches
### and to calculte the distance to its TSS (TTS)
#========================================================================================#
getClosest <- function(GRanges_data, StructureTracks, modeling = TRUE,
                       strandAsFeature = FALSE,
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
    if (strandAsFeature) {
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
    if (strandAsFeature) {
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
