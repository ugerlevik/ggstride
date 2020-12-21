# ggstride
Protein secondary structure analysis (with STRIDE) and visualization (with ggplot2) in R.

Two functions are included in this package:  
__ssa()__ returns a data frame that contains residue identifiers (resname_chain_resid; e.g., MET_A_1) as column names
and STRIDE assignments of each residue as each row for each frame of the trajectory.  

__ssa_plot()__ returns a plot that contains percentage of secondary structures observed along the trajectories per residue 
by comparing the two ssa data frames obtained by ssa()   

## Requirements:
__R packages:__ "stringi", "readr", "bio3d", "ggplot2", "ggpubr", "ggsci", "reshape2"  
__In system:__ STRIDE (http://webclu.bio.wzw.tum.de/stride/). Stride should be executable in your terminal as "stride".  

## Usage:
Open your R or RStudio as administrator (or root/sudo) to be able to run "stride" from your system terminal.

Load required libraries to your R session:
```{r}
library(stringi)
library(readr)
library(bio3d)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)
```

Source the functions to define ssa() and ssa_plot():
```{r}
source("functions/ssa.R")
source("functions/ssa_plot.R")
```

Read your trajectory files and the related pdb files:
```{r}
pdb_WT <- read.pdb("example/wt.pdb")
dcd_WT <- read.dcd(trjfile = "example/wt.dcd")

pdb_S50P <- read.pdb("example/S50P_mutant.pdb")
dcd_S50P <- read.dcd(trjfile = "example/S50P_mutant.dcd")
```

Use ssa() to calculate secondary structures:
```{r}
ssa_WT <- ssa(pdb_WT, dcd_WT) # Note: If your trajectory has two parts, you can use rbind: rbind(dcd_WT_r1, dcd_WT_r2).
ssa_S50P <- ssa(pdb_S50P, dcd_S50P) 
```

Use ssa_plot() to visualize your results:
```{r}
ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_S50P,
         name1 = "Wild-type", color_number1 = 1,
         name2 = "Ser50Pro", color_number2 = 2)
         
# You can assign the plot to a variable:
plot_ssa_WT_S50P <- ssa_plot(ssa1 = ssa_WT, ssa2 = ssa_S50P,
                             name1 = "Wild-type", color_number1 = 1,
                             name2 = "Ser50Pro", color_number2 = 2)
```
