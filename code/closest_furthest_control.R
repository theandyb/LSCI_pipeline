library(tidyverse)
library(data.table)

population <- "SAS"

for(subtype in c("AT_CG","AT_GC", "AT_TA",
                 "GC_AT", "GC_TA", "GC_CG",
                 "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")){
  DT <- data.table(read_csv(paste0("output/controls/", population,"/", subtype, ".csv"),
                            col_names = c("chr", "pos", "s_ref", "c_ref", "window", "dist", "spos", "motif")))
  min_df <- DT[ , .SD[which.min(dist)], by = pos]
  max_df <- DT[ , .SD[which.max(dist)], by = pos]

  out_min <- paste0("output/controls/", population,"/", subtype, ".csv.min")
  out_max <- paste0("output/controls/", population,"/", subtype, ".csv.max")

  fwrite(min_df, out_min)
  fwrite(max_df, out_max)
}
