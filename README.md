# ggstride
Protein secondary structure analysis (with STRIDE) and visualization (with ggplot2) in R.

Two functions are included in this package:  
__ssa()__ returns a data frame that contains residue identifiers (resname_chain_resid; e.g., MET_A_1) as column names
and STRIDE assignments of each residue as each row for each frame of the trajectory.  

__ssa_plot()__ returns a plot that contains percentage of secondary structures observed along the trajectories per residue 
by comparing the two ssa data frames obtained by ssa()   

## Requirements
__R packages:__ "stringi", "readr", "bio3d", "ggplot2", "ggpubr", "ggsci", "reshape2"  
__In system:__ STRIDE (http://webclu.bio.wzw.tum.de/stride/). Stride should be executable in your terminal as "stride".  


## Installation
You can install ggstride from GitHub with:

```{r}
# install.packages("devtools")
devtools::install_github("umutgerlevik/ggstride")
```

## Usage

Load required libraries to your R session:
```{r}
library(ggstride)
library(bio3d)
```

Read your trajectory files and the related pdb files:
```{r}
pdb_WT <- read.pdb(system.file("example/wt.pdb", package = "ggstride"))
dcd_WT <- read.dcd(trjfile = system.file("example/wt.dcd", package = "ggstride"))

pdb_mutant <- read.pdb(system.file("example/mutant.pdb", package = "ggstride"))
dcd_mutant <- read.dcd(trjfile = system.file("example/mutant.dcd", package = "ggstride"))
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
pdf("example/ssa_plot_all.pdf", width = 20, height = 20)
plot_ssa_all
dev.off()

pdf("example/ssa_plot_520_530.pdf", width = 20, height = 20)
plot_ssa_520_530
dev.off()
```


