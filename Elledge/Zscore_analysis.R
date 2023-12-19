library(virScanR)
library(mmR)
library(dplyr)
library(data.table)

counts <- mmR::mm.fastread("path_to_counts_file")

propTrim = .04

convert_params1 <-
  virScanR::vs.set_params_convert(stat  = "Z_score",  #two other pertinent: NLP_nBinom", "NLP_pois"
                                  makeGroups = TRUE,
                                  makeGroupSize = 300,
                                  makeGroup_by_col = "input",
                                  idCol = "id",
                                  groupTogetherIfGreaterThanGroupSize = TRUE,
                                  splitHighVarMeanGrps = TRUE,
                                  cols_to_remove = NULL,
                                  cols_to_not_evaluate = c("id","input"), 
                                  propExtremesToRemove = propTrim,
                                  removeTail = "both",
                                  coresAcrossGroups = parallel::detectCores()-2,
                                  returnAs = "data.table", 
                                  returnParams = TRUE)

datZ = 
  virScanR::vs.convert(
    data = counts, 
    paramsList = convert_params1
  )


mmR::mm.fastwrite(datZ$out, path = "path_to_output_file")
