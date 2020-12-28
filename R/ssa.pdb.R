#########################################################################
## Project: Secondary Structure Analysis (SSA)
## Script purpose: Calculate Secondary Structures of Structures
## Date: Dec 28, 2020
## Author: Umut Gerlevik
#########################################################################

#' @title ggstride
#'
#' @description ssa.pdb() returns a data frame that contains
#' residue identifiers (resname_chain_resid; e.g., MET_A_1)
#' as column names and STRIDE assignments of each residue as
#' the row
#'
#' @param pdb: the pdb file including same atoms with the trajectory. bio3d::read.pdb() should be used to read.
#'
#' @return dataframe
#'
#' @examples
#' pdb <- read.pdb("example/wt.pdb")
#' ssa_df <- ssa.pdb(pdb)
#'
#' package::ssa(pdb)
#'
#' @export
#'
#' @import "stringi", "readr", "bio3d", "stats"

# Requirements:
# R packages: "stringi", "readr", "bio3d"
# In system: "stride" (http://webclu.bio.wzw.tum.de/stride/)
## Stride should be executable in your terminal as "stride".
## Simply, you can add the environmental variables.

# Run R or RStudio as administrator (or root with sudo) to run stride
# from the system terminal from R.

# ssa.pdb() will return a ssa data frame that contains
# residue identifiers as colnames and STRIDE
# secondary structure results as the row

# pdb: the pdb file including same atoms with the trajectory,
# bio3d::read.pdb() should be used to read

ssa.pdb <- function(pdb) {
  require(stringi)
  require(readr)
  require(bio3d)

  tmp <- stri_rand_strings(1, 5)
  dir.create(paste0("stride_", tmp))
  out <- matrix(ncol = length(pdb[pdb$calpha]))
  colnames(out) <- paste(pdb$atom$resid[pdb$calpha],
                         pdb$atom$chain[pdb$calpha],
                         pdb$atom$resno[pdb$calpha], sep = "_")
  out <- as.data.frame(na.omit(out))

  operating_system <- Sys.info()
  operating_system <- operating_system['sysname']

  write.pdb(pdb, file = paste0("stride_", tmp, "/", tmp, ".pdb"))

  if(operating_system == "Windows") {
    exe <- system.file("bin/stride_WIN.exe", package = "ggstride")
    res <- suppressWarnings(system(paste0(exe, " stride_", tmp, "/", tmp, ".pdb"), intern = TRUE))
  } else if(operating_system == "Linux") {
    exe <- system.file("bin/stride_LINUXAMD64", package = "ggstride")
    res <- suppressWarnings(system(paste0(exe, " stride_", tmp, "/", tmp, ".pdb"), intern = TRUE))
  } else {
    tryCatch({
      res <- suppressWarnings(system(paste0("stride", " stride_", tmp, "/", tmp, ".pdb"), intern = TRUE))
    },
    error = function(e) {
      print("Error: Define 'stride' as an environmental variable!")
    })
  }

    res <- read_table(as.character(res[grep("ASG", res)]), col_names = FALSE)[, c(2, 3, 4, 6, 7)]
    out <- rbind(out, res$X6)

  colnames(out) <- paste(pdb$atom$resid[pdb$calpha],
                         pdb$atom$chain[pdb$calpha],
                         pdb$atom$resno[pdb$calpha], sep = "_")

  unlink(paste0("stride_", tmp), recursive = TRUE)

  return(as.data.frame(out))
}
