devtools::document()
# 1. Prepare input
data("experimental_data") # experimental data with two column (feature and rt)
data("experimental_car") # car info in experimental data, it comes from experimental data
# 2. Prepare Standard car
data("standard_car") # standard car or QC standard car (when user use a workflow different from Tanglab)
# 3. Predict RT for experimental_car
ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
experimental_car <- ResList$lc_tibble
ResList$p
# 4. Predict RT for standard_car
ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
standard_car <- ResList$lc_tibble
ResList$p
# 5. Add dead time
experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
# 6. Calculate retention index
experimental_data$ri <- sapply(experimental_data$rt, function(x) {
  rt2ri(x, experimental_lc = experimental_car)
})
# 7. Retention time correction
# Make the experimental lc_tibble and standard lc_tibble have the same number of lc.
ResList <- compare_lc_tibble(experimental_lc = experimental_car, standard_lc = standard_car)
experimental_car <- ResList$experimental_lc
standard_car <- ResList$standard_lc
# Calculate new retention time
experimental_data <- orignRT2newRT(experimental_data = experimental_data,
                                   experimental_lc = experimental_car,
                                   standard_lc = standard_car, thread = 1)
