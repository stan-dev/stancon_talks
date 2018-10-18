# test hypothesis of equal accurary for two sequences of predictions (rps can be seen as error)
# use two sided test: ie eiher method can be worse

# 
AddDMtest <- function(result, rps_data, ref_model_nr = 12){
  require(forecast)
  require(data.table)
  if(!is.data.table(result)) stop('not a data.table')
  if(!is.data.table(rps_data)) stop('not a data.table')
  result <- result[, DMstat := 0]
  result <- result[, DMpval := 0]
  # make sure matches are sorted properly
  setkey(rps_data, row_id)
  e1 <- rps_data[model_nr == ref_model_nr & !is.na(rps_vec), ]$rps_vec 
  # run dm.test
  for(i in result$model_nr){
    e2 <- rps_data[model_nr == i & !is.na(rps_vec), ]$rps_vec 
    if(i != ref_model_nr) {dm_res <- dm.test(e1, e2, h = 1, power = 1, alternative = "two.sided")
      result <- result[model_nr == i, DMstat := round(dm_res$statistic, 3)]
      result <- result[model_nr == i, DMpval := round(dm_res$p.value, 3)]
      }
  }
  result <- result[model_nr == ref_model_nr, DMpval := NA]
  result <- result[model_nr == ref_model_nr, DMstat := NA]
  setnames(result, "DMpval", paste("DMpval", ref_model_nr, sep = '_'))
  setnames(result, "DMstat", paste("DMstat", ref_model_nr, sep = '_'))
  result
}
