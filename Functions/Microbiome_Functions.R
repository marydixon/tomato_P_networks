## ---------------------------
##
## Script name: Custom R functions
##
## Author: Dan Manter
##
## Date Created: March 6 2022
##
## ---------------------------
##
## Notes: List of functions
##   tax.clean
##   emu_to_phyloseq 
##
## ---------------------------



# this function changes blanks in a taxonomy table to the lowest common ancestor
tax.clean <- function (tax) {
  for (i in 1:nrow(tax)) {
    if (tax[i,7] == "") {
      if (tax[i,6] == "") {
        if (tax[i,5] == "") {
          if (tax[i,4] == "") {
            if (tax[i,3] == "") {
              if (tax[i,2] == "") {
                if (tax[i,1] == "") {
                  tax[i,7] <- "<NA>"
                } else {
                  tax[i,7] <- paste0("s_of_", tax[i,1])
                }
              } else {
                tax[i,7] <- paste0("s_of_", tax[i,2])
              }
            } else {
              tax[i,7] <- paste0("s_of_", tax[i,3])
            }
          } else {
            tax[i,7] <- paste0("s_of_", tax[i,4])
          }
        } else{
          tax[i,7] <- paste0("s_of_", tax[i,5])
        }
      } else {
        tax[i,7] <- paste0("s_of_", tax[i,6])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,6] == "") {
      if (tax[i,5] == "") {
        if (tax[i,4] == "") {
          if (tax[i,3] == "") {
            if (tax[i,2] == "") {
              if (tax[i,1] == "") {
                tax[i,6] <- "<NA>"
              } else {
                tax[i,6] <- paste0("g_of_", tax[i,1])
              }
            } else {
              tax[i,6] <- paste0("g_of_", tax[i,2])
            }
          } else {
            tax[i,6] <- paste0("g_of_", tax[i,3])
          }
        } else {
          tax[i,6] <- paste0("g_of_", tax[i,4])
        }
      } else {
        tax[i,6] <- paste0("g_of_", tax[i,5])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,5] == "") {
      if (tax[i,4] == "") {
        if (tax[i,3] == "") {
          if (tax[i,2] == "") {
            if (tax[i,1] == "") {
              tax[i,5] <- "<NA>"
            } else {
              tax[i,5] <- paste0("f_of_", tax[i,1])
            }
          } else {
            tax[i,5] <- paste0("f_of_", tax[i,2])
          }
        } else {
          tax[i,5] <- paste0("f_of_", tax[i,3])
        }
      } else {
        tax[i,5] <- paste0("f_of_", tax[i,4])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,4] == "") {
      if (tax[i,3] == "") {
        if (tax[i,2] == "") {
          if (tax[i,1] == "") {
            tax[i,4] <- "<NA>"
          } else {
            tax[i,4] <- paste0("o_of_", tax[i,1])
          }
        } else {
          tax[i,4] <- paste0("o_of_", tax[i,2])
        }
      } else {
        tax[i,4] <- paste0("o_of_", tax[i,3])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,3] == "") {
      if (tax[i,2] == "") {
        if (tax[i,1] == "") {
          tax[i,3] <- "<NA>"
        } else {
          tax[i,3] <- paste0("c_of_", tax[i,1])
        }
      } else {
        tax[i,3] <- paste0("c_of_", tax[i,2])
      }
    }
  }
  
  for (i in 1:nrow(tax)) {
    if (tax[i,2] == "") {
      if (tax[i,1] == "") {
        tax[i,2] <- "<NA>"
      } else {
        tax[i,2] <- paste0("p_of_", tax[i,1])
      }
    }
  }
  return(tax)  
}


# custom function to import emu output files into phyloseq
emu_to_phyloseq <- function (RA_file=NULL, meta_file=NULL, sheet=NULL, 
                             range=NULL, sample_names=NULL, run_name=NULL) {
  df <- read.csv(RA_file, header=T,  
                 row.names='tax_id')
  tax <- df[,2:8]
  tax <- tax.clean(tax)
  TAX <- tax_table(as.matrix(tax))
  
  otu <- df[,9:ncol(df)]
  otu[is.na(otu)] <- 0

  meta <- data.frame(read_excel(meta_file, sheet=sheet, range=range))
  common <- Reduce(intersect, list(meta$sample, names(otu)))
  meta <- meta[meta$sample %in% common,]
  SAMP <- sample_data(meta)
  sample_names(SAMP) <- meta[,1]
  
  otu <- otu[,common]
  OTU <- otu_table(otu, taxa_are_rows=T)
  
  phy_obj <- phyloseq(SAMP, OTU, TAX )
  sample_names(phy_obj) <- meta[,sample_names]
  
  return(phy_obj)
}

# custom function to import emu output files into phyloseq
emu_to_phyloseq_edited <- function (RA_file=NULL, meta_file=NULL, sheet=NULL, 
                             range=NULL, run_name=NULL) {
  df <- read.csv(RA_file, header=T,  
                 row.names='tax_id')
  tax <- df[,2:8]
  tax <- tax.clean(tax)
  TAX <- tax_table(as.matrix(tax))
  
  otu <- df[,9:ncol(df)]
  otu[is.na(otu)] <- 0
  
  meta <- data.frame(read_excel(meta_file, sheet=sheet, range=range))
  common <- Reduce(intersect, list(meta$sample, names(otu)))
  meta <- meta[meta$sample %in% common,]
  SAMP <- sample_data(meta)
  sample_names(SAMP) <- meta[,1]
  
  otu <- otu[,common]
  OTU <- otu_table(otu, taxa_are_rows=T)
  
  phy_obj <- phyloseq(SAMP, OTU, TAX )
  sample_names(phy_obj) <- meta[,sample_names]
  
  return(phy_obj)
}


