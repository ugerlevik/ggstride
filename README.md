# <img src="https://github.com/umutgerlevik/ggstride/blob/master/vignettes/logo.png?raw=true" align="left" height=65/> ggstride: An R Package To Analyze and Visualize 2D Structures of Proteins

Three functions are included in this package:  
__ssa()__ returns a data frame that contains residue identifiers (resname_chain_resid; e.g., MET_A_1) as column names and STRIDE assignments of each residue as each row for each frame of the trajectory.  

__ssa_plot()__ returns a plot that contains percentage of secondary structures observed along the trajectories per residue by comparing the two ssa data frames obtained by ssa(). 

__ssa.pdb()__ returns a data frame that contains residue identifiers (resname_chain_resid; e.g., MET_A_1) as column names and STRIDE assignments of each residue of the pdb as a row.

## Installation
You can install ggstride from GitHub with:

```{r}
# install.packages("devtools")
devtools::install_github("ugerlevik/ggstride")
```

## Usage

Load required libraries to your R session:
```{r}
library(ggstride)
library(bio3d)
```

Read your trajectory files and the related pdb files:
```{r}
pdb_WT <- read.pdb(system.file("extdata/wt.pdb", package = "ggstride"))
dcd_WT <- read.dcd(trjfile = system.file("extdata/wt.dcd", package = "ggstride"))

pdb_mutant <- read.pdb(system.file("extdata/mutant.pdb", package = "ggstride"))
dcd_mutant <- read.dcd(trjfile = system.file("extdata/mutant.dcd", package = "ggstride"))
```

Use ssa() to calculate secondary structures:
```{r}
# Note: If your trajectory has two or more parts, you can use rbind: rbind(dcd_WT_part1, dcd_WT_part2)
ssa_WT <- ssa(pdb_WT, dcd_WT)
ssa_mutant <- ssa(pdb_mutant, dcd_mutant) 
```

Use ssa_plot() to visualize your results:
```{r}
# Plot all:
ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_mutant,
         name1 = "Wild-type",  name2 = "Mutant",
         color_number1 = 1, color_number2 = 2)
```

![SSA Plot](https://github.com/umutgerlevik/ggstride/blob/master/vignettes/ssa_plot_all.png?raw=true "SSA Plot")

You can assign the plot to a variable:
```{r}
plot_ssa_all <- ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_mutant,
                         name1 = "Wild-type", color_number1 = 1,
                         name2 = "Mutant", color_number2 = 2)
```

Plot with a focus between residues 520 and 530:
```{r}
plot_ssa_520_530 <- ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_mutant,
                             name1 = "Wild-type",  name2 = "Mutant",
                             resid1 = 520, resid2 = 530,
                             color_number1 = 1, color_number2 = 3)
plot_ssa_520_530
```

![SSA Plot Focused](https://github.com/umutgerlevik/ggstride/blob/master/vignettes/ssa_plot_520_530.png?raw=true "SSA Plot Focused")


Get the pdf outputs:
```{r}
pdf("ssa_plot_all.pdf", width = 20, height = 20)
plot_ssa_all
dev.off()

pdf("ssa_plot_520_530.pdf", width = 20, height = 20)
plot_ssa_520_530
dev.off()
```

### References

Please cite this paper:
Gerlevik U, Ergoren MC, Sezerman OU, Temel SG. 2022. Structural analysis of M1AP variants associated with severely impaired spermatogenesis causing male infertility. PeerJ 10:e12947 https://doi.org/10.7717/peerj.12947


You might also cite the studies below:

STRIDE method (http://webclu.bio.wzw.tum.de/stride/):
Frishman D, Argos P. Knowledge-based protein secondary structure assignment. Proteins. 1995 Dec;23(4):566-79. doi: 10.1002/prot.340230412. PMID: 8749853.

STRIDE executables (https://www.ks.uiuc.edu/Research/vmd/):
Humphrey W, Dalke A, Schulten K. VMD: visual molecular dynamics. J Mol Graph. 1996 Feb;14(1):33-8, 27-8. doi: 10.1016/0263-7855(96)00018-5. PMID: 8744570.

bio3d R package (http://thegrantlab.org/bio3d/):
Grant BJ, Rodrigues AP, ElSawy KM, McCammon JA, Caves LS. Bio3d: an R package for the comparative analysis of protein structures. Bioinformatics. 2006 Nov 1;22(21):2695-6. doi: 10.1093/bioinformatics/btl461. Epub 2006 Aug 29. PMID: 16940322.

ggplot2 R package:
Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org.

ggpubr R package:
Alboukadel Kassambara (2020). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.3.0.
https://CRAN.R-project.org/package=ggpubr

