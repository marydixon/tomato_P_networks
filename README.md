# Tomato-P-Networks

## Overview

This repository contains the R code and data associated with the following manuscript:

Dixon M M, Afkairin A, Manter D K, Vivanco J M. Rhizosphere microbiome co-occurrence network analysis across a tomato domestication gradient

When plant available phosphorus (P) is lost from soil solution, it often accumulates in the soil as a pool of unavailable legacy P. We grew wild tomato, traditional tomato, and modern tomato throughout their vegetative developmental stage. Co-occurrence network analysis revealed that as tomato progressed along stages of domestication, the rhizosphere microbiome complexity changed, wherein there was a decline in complexity with traditional tomato and subsequent increase in network complexity with modern tomato. Further, traditional tomato showcased a unique rhizosphere by harboring keystone taxa (Peribacillus muralis) capable of performing many beneficial soil functions. By illustrating these changing patterns of network complexity in the tomato rhizosphere microbiome, we can further understand how plant domestication and breeding has shaped plant-microbe interactions.

## File descriptions

1.  "Networks.Rmd" includes code to generate network files to upload into Gephi.
    -   Data inputs:
        -   "P.count.RDS": Phyloseq object that includes count data.\
    -   Data outputs:
        -   "networkModernSpecies.gexf": Network file for modern tomato.
        -   "networkTraditionalSpeciesgexf": Network file for traditional tomato.
        -   "networkWildSpecies.gexf": Network file for wild tomato.
        -   "WholeNetworkProject.gephi": Network file for all tested tomato.
        -   Results corresponding to figure 1
        -   Results corresponding to figure 4
2.  "keystonetaxa.Rmd" includes code to determine differences in the presence of identified keystone taxa.
    -   Data inputs:
        -   "EMU_database.GIBBs.KO.PICRUST2.xlsx": Shows information from PICRUST database that shows functional abundance of mapped bacteria. "P.rel.RDS": Phyloseq object showing relative abundance of bacteria.
    -   Outputs:
        -   "Figure2.png" Predicted functions of keystone taxa
        -   "Figure3.png" Square-root-transformed relative abundance of the sum of keystone bacteria in the rhizosphere of wild, traditional, and modern tomato
        -   Results corresponding to figure 2
        -   Results corresponding to figure 3
3.  "Microbiome_Read_In_Data.R" includes code to convert ".log" demultiplex files to excel worksheets. It also includes codes to convert EMU files for indivial sequencing runs and compiles them into one phyloseq object.

-   Function inputs:
    -   "Microbiome_Functions.R": Includes custom functions that allow for conversion from EMU to phyloseq ojects.
    -   Data inputs:
        -   "demultiplex_plateA.log": Includes sequence reads for plate A.\
        -   "demultiplex_plateB.log": Includes sequence reads for plate B.\
        -   "demultiplex_plateC.log": Includes sequence reads for plate C.\
        -   "EMU_Relative_Abundance_plateA.csv": Taxonomic relative abundance information for plate A.
        -   "EMU_Relative_Abundance_plateB.csv": Taxonomic relative abundance information for plate B.
        -   "EMU_Relative_Abundance_plateC.csv": Taxonomic relative abundance information for plate C.
    -   Outputs:
        -   "demultiplex_plateA.xlsx" Sequence reads for plate A.
        -   "demultiplex_plateB.xlsx" Sequence reads for plate B.
        -   "demultiplex_plateC.xlsx" Sequence reads for plate C.
        -   "Phyloseq.Plate.A.RDS" Phyloseq object that includes sequence reads and relative abundance information for plate A. - "Phyloseq.Plate.B.RDS" Phyloseq object that includes sequence reads and relative abundance information for plate B.
        -   "Phyloseq.Plate.C.RDS" Phyloseq object that includes sequence reads and relative abundance information for plate C.
        -   "Phyloseq.RelativeAbundance.RDS" merged phyloseq object that combines "Phyloseq.Plate.A.RDS", "Phyloseq.Plate.B.RDS", "Phyloseq.Plate.C.RDS"
        -   "Phyloseq.Count.Data.RDS" A phyloseq object that includes count data for microbial abundance, combined for all plates.
