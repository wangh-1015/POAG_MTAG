### LDSC  SNP A1 A2 N Z
# more detailed instruction could be found in website:"https://github.com/mglev1n/ldscr"
# GWAS summary data used in our analysis: GBMI_POAG, MVP_POAG, Finn_POAG, 
# prev_MTAG(MTAG(Han et al.,2023))

# install.packages("devtools")
devtools::install_github("mglev1n/ldscr")

library(dplyr)
library(data.table)
library(ldscr)

### genetic correlation
ldsc <- ldsc_rg(
  munged_sumstats = list(
  "GBMI"=GBMI_POAG,
  "MVP" = MVP_POAG, 
  "Finn" = Finn_POAG,
  "prev_MTAG"=prev_MTAG,
  "pres_MTAG"=our_MTAG
                      ),
  ancestry = "EUR"
  )

h2 <- ldsc$h2
rg <- ldsc$rg
h2
rg


