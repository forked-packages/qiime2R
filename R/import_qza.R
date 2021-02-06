#' Read QIIME2 Artifacts (.qza)
#'
#' Reads QIIME2 artifacts (.qza) file and extracts embedded data and object
#' metadata into an R object.
#'
#' @param file (string) Path to the input file,
#'        e.g.: `file = "~/data/moving_pictures/table.qza"`.
#' @param tmp (string) A temporary directory that the object will be
#'        decompressed to (`default = "tempdir()"`).
#' @param rm (`TRUE`|`FALSE`)
#'        Should the decompressed object be removed at completion of function
#'        (default = `TRUE`).
#'
#' @return An R object with data and metadata in QIIME artifact file (.qza) as
#'        `list`, `data.frame`, or `matrix` (depending on the type of artifact).
#'        The object has given 2 additional classes:
#'  - Class `qiime2_qza`;
#'  - Other class which coincides with QIIME2 semantic type, e.g.
#'   `FeatureData[Sequence]`.
#'
#' The object contains attribute `metadata`, which is a named list of the
#' following fields:
#'
#' - `data` -- the raw data, e.g., OTU table as matrix or tree in phylo format;
#' - `uuid` -- the unique identifer of the artifact;
#' - `type` -- the semantic type of the object (e.g., `FeatureData[Sequence]`);
#' - `format` -- the format of the qiime artifact;
#' - `provenance` -- information tracking how the object was created;
#' - `contents` -- a table of all the files contained within the artifact and
#'    their file size;
#' - `version` -- the reported version for the artifact, a warning error may be
#'    thrown if a new version is seen;
#'
#' @md
#' @examples
#' \dontrun{
#' asv <- import_qza("data/table.qza")
#'
#' # Get metadata
#' attr(asv, "metadata")
#' }
#'
#' @export


import_qza <- function(file, tmp = tempdir(), rm = TRUE) {

  if (missing(file)) {
    stop("Path to artifact (.qza) not provided.")
  }

  if (!file.exists(file)) {
    stop("Input artifact (", file, ") not found. ",
      "Please check path and/or use list.files() ",
      "to see files in current working directory.")
  }

  if (!grepl("qza$", file)) {
    stop("Provided file is not qiime2 artifact (.qza).")
  }

  unzip(file, exdir = tmp) # Unzfips to tmp dir
  unpacked <- unzip(file, exdir = tmp, list = TRUE) # Lists contents but does not unzip

  # Start by loading in the metadata not assuming it will be first file listed
  metadata <- yaml::read_yaml(
    paste0(tmp, "/", paste0(gsub("/..+", "", unpacked$Name[1]), "/metadata.yaml"))
  )

  # Construnct ath to unzipped metadata directory
  qza_dir <- function(...) {
    paste0(tmp, "/", metadata$uuid, ...)
  }


  metadata$contents <- data.frame(files = unpacked)
  metadata$contents$size <- sapply(paste0(tmp, "/", metadata$contents$files), file.size)
  metadata$version <- read.table(qza_dir("/VERSION"))


  # Get data dependent on format
  switch(
    metadata$format,

    "NewickDirectoryFormat" = {
      qza_data <- read.tree(qza_dir("/data/tree.nwk"))
    },

    "DistanceMatrixDirectoryFormat" = {
      qza_data <- as.dist(read.table
        (qza_dir("/data/distance-matrix.tsv"),
          header = TRUE, row.names = 1, fill = TRUE)
      )
    },

    "TSVTaxonomyDirectoryFormat" = {
      qza_data <- read.table(
        qza_dir("/data/taxonomy.tsv"),
        sep = "\t", header = TRUE, quote = "", comment.char = ""
      )
    },

    "OrdinationDirectoryFormat" = {
      qza_data <- suppressWarnings(
        readLines(qza_dir("/data/ordination.txt"))
      )
      metadata <- parse_ordination(metadata, tmp)
    },

    "DNASequencesDirectoryFormat" = {
      qza_data <- Biostrings::readDNAStringSet(
        qza_dir("/data/dna-sequences.fasta")
      )
    },

    "AlignedDNASequencesDirectoryFormat" = {
      qza_data <-
        Biostrings::readDNAMultipleAlignment(
          qza_dir("/data/aligned-dna-sequences.fasta"),
          format = "fasta"
        )
    },

    "EMPPairedEndDirFmt" = ,
    "EMPSingleEndDirFmt" = ,
    "FastqGzFormat" = ,
    "MultiplexedPairedEndBarcodeInSequenceDirFmt" = ,
    "MultiplexedSingleEndBarcodeInSequenceDirFmt" = ,
    "PairedDNASequencesDirectoryFormat" = ,
    "SingleLanePerSamplePairedEndFastqDirFmt" = ,
    "SingleLanePerSampleSingleEndFastqDirFmt" = {

      qza_data <- data.frame(files = list.files(qza_dir("/data")))

      qza_data$size <- format(sapply(qza_data$files, function(x) {
        file.size(qza_dir("/data/", x))
      }, simplify = TRUE))

    },

    "AlphaDiversityDirectoryFormat" = {
      qza_data <- read.table(qza_dir("/data/alpha-diversity.tsv"))
    },

    "DifferentialDirectoryFormat" = {
      defline <- suppressWarnings(
        readLines(qza_dir("/data/differentials.tsv"))[2]
      )
      defline <- strsplit(defline, split = "\t")[[1]]
      defline[grep("numeric", defline)] <- "double"
      defline[grep("categorical|q2:types", defline)] <- "factor"

      coltitles <- strsplit(
        suppressWarnings(
          readLines(qza_dir("/data/differentials.tsv"))[1]
        ),
        split = "\t")[[1]]

      qza_data <-
        read.table(qza_dir("/data/differentials.tsv"),
          header = F, col.names = coltitles, skip = 2, sep = "\t",
          colClasses = defline, check.names = FALSE)

      colnames(qza_data)[1] <- "Feature.ID"

    },
    # else
    {
      if (grepl("BIOMV", metadata$format)) {
        qza_data <- read_q2biom(qza_dir("/data/feature-table.biom"))

      } else if (grepl("StatsDirFmt", metadata$format)) {
        # Can be tsv or csv

        if (paste0(metadata$uuid, "/data/stats.csv") %in% metadata$contents$files.Name) {
          qza_data <- read.csv(
            qza_dir("/data/stats.csv"), header = TRUE, row.names = 1)
        }

        if (paste0(metadata$uuid, "/data/stats.tsv") %in% metadata$contents$files.Name) {
          qza_data <- read.table(
            qza_dir("/data/stats.tsv"), header = TRUE, row.names = 1, sep = "\t"
          )
        }

      } else {
        warning(
          "Format '", metadata$format, "' not supported. ",
          "Only a list of internal files and provenance are being imported."
        )

        qza_data <- list.files(qza_dir("/data"))
      }
    }
  )

  # Add Provenance
  pfiles <-
    paste0(tmp, "/", grep("..+provenance/..+action.yaml", unpacked$Name, value = TRUE))

  metadata$provenance <- lapply(pfiles, yaml::read_yaml)

  names(metadata$provenance) <-
    grep("..+provenance/..+action.yaml", unpacked$Name, value = TRUE)

  if (rm == TRUE) {
    unlink(qza_dir(), recursive = TRUE)
  }

  structure(
    qza_data,
    class = c(metadata$type, "qiime2_qza", class(qza_data)),
    metadata = structure(metadata, class = c("qza_metadata", "list"))
  )
}

#' @rdname import_qza
#' @export
print.qiime2_qza <- function(x, ...) {
  cat("QIIME 2 artifact (.qza) object \n")
  # cat("     uuid: ", attr(x, "metadata")$uuid, " \n")
  # cat("   Format: ", attr(x, "metadata")$format, " \n")
  # cat("     Type: ", attr(x, "metadata")$type, " \n")
  cat("R classes: ", paste(class(x), collapse = ", "), sep = "")

  if (!is.data.frame(x) && !is.matrix(x)) {
    # print(ls.str(x))
    cat("\n")
    cat(str(x, max.level = 1, give.attr = FALSE))

  } else {
    cat("\nnrow:", nrow(x), "\nncol:", ncol(x))
  }
}

#' @rdname import_qza
#' @export
print.qza_metadata <- function(x, ...) {
  cat("Metadata of QIIME 2 artifact. \n")
  cat(str(x, max.level = 1, give.attr = FALSE))
}
