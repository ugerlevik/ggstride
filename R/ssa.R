#########################################################################
## Project: Secondary Structure Analysis (SSA)
## Script purpose: Calculate Secondary Structures in Trajectory
## Date: Dec 19, 2020
## Author: Umut Gerlevik
#########################################################################

#' @title ggstride
#'
#' @description ssa() returns a data frame that contains
#' residue identifiers (resname_chain_resid; e.g., MET_A_1)
#' as column names and STRIDE assignments of each residue as
#' each row for each frame of the trajectory.
#'
#' @param pdb,traj,start,last
#'
#' @return dataframe
#'
#' @examples
#' pdb_WT <- read.pdb("example/wt.pdb")
#' dcd_WT <- read.dcd(trjfile = "example/wt.dcd")
#' ssa_WT <- ssa(pdb_WT, dcd_WT)
#'
#' package::ssa(pdb, traj, start = 1, last = nrow(traj))
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

# ssa() will return a ssa data frame that contains
# residue identifiers as colnames and STRIDE
# secondary structure results at each row for
# each frame of the trajectory

# pdb: the pdb file including same atoms with the trajectory,
	# bio3d::read.pdb() should be used to read
# traj: the trajectory file, xyz format of bio3d should be used,
	# could be any type of trajectory supported by bio3d
# start (default = 1): the number of first frame you will analyze
# last (default = number of frames in your trajectory):
	# the number of last frame you will analyze

ssa <- function(pdb, traj, start = 1, last = nrow(traj)) {
  require(stringi)
  require(readr)
  require(bio3d)

  tmp <- stri_rand_strings(1, 5)
  dir.create(paste0("_outputs/stride_", tmp))
  out <- matrix(ncol = length(pdb[pdb$calpha]))
  colnames(out) <- paste(pdb$atom$resid[pdb$calpha],
                         pdb$atom$chain[pdb$calpha],
                         pdb$atom$resno[pdb$calpha], sep = "_")
  out <- as.data.frame(na.omit(out))

  for(i in 1:nrow(traj)) {
    write.pdb(pdb, file = paste0("_outputs/stride_", tmp, "/", tmp, "_", i, ".pdb"), xyz = traj[i,])
    res <- suppressWarnings(system(paste0("stride ", "_outputs/stride_", tmp, "/", tmp, "_", i, ".pdb"), intern = TRUE))
    res <- read_table(as.character(res[grep("ASG", res)]), col_names = FALSE)[, c(2, 3, 4, 6, 7)]
    out <- rbind(out, res$X6)
  }
  colnames(out) <- paste(pdb$atom$resid[pdb$calpha],
                         pdb$atom$chain[pdb$calpha],
                         pdb$atom$resno[pdb$calpha], sep = "_")

  unlink(paste0("_outputs/stride_", tmp), recursive = TRUE)

  rownames(out) <- start:last

  return(as.data.frame(out))
}
