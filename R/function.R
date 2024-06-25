#' @title lc2rt
#' @description
#' Use lc_tibble to predict rt.
#'
#' @param lc_tibble A lc_tibble.
#' @param targetRt The maximum retention time user want to overwrite.
#'
#' @return A list.
#' @export
#'
#' @examples
#' data("experimental_data")
#' data("experimental_car")
#' data("standard_car")
#' ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
#' experimental_car <- ResList$lc_tibble
#' ResList$p
lc2rt <- function(lc_tibble, targetRt){
  #maxC0 <- max(lc_tibble$lc)
  lc_tibble$lnC <- log(lc_tibble$lc)
  predicted_lc <- lc_tibble %>%
    dplyr::filter(predicted) %>%
    dplyr::filter(lc > 0)
  model_lc <- lc_tibble %>%
    dplyr::filter(!predicted) %>%
    dplyr::filter(lc > 0)
  lmRes <- lm(formula = rt ~ lnC, data = model_lc)
  intercept <- as.numeric(coefficients(lmRes)[1])
  slope <- as.numeric(coefficients(lmRes)[2])
  r_squared <- summary(lmRes)$r.squared
  sign <- ifelse(intercept < 0, "-", "+")
  text <- paste0("y = ", round(slope, 4), "lnC ", sign, " ", abs(round(intercept, 4)), "\n",
                 "R2 = ", round(r_squared, 4))
  predicted_lc$rt <- slope * predicted_lc$lnC + intercept
  lc_tibble_tmp <- lc_tibble %>% dplyr::filter(lc <= 0)
  lc_tibble <- rbind(lc_tibble_tmp, model_lc, predicted_lc) %>%
    dplyr::arrange(ri)
  maxRt <- max(lc_tibble$rt)
  if(targetRt < maxRt) targetRt <- maxRt
  maxC <- max(lc_tibble$lc)
  while(maxRt < targetRt){
    maxC <- maxC + 1
    maxRt <- slope * log(maxC) + intercept
    lc_tibble_tmp <- dplyr::tibble(id = NA, name = NA, rt = maxRt, lc = maxC, ri = maxC * 100,
                                    exactmass = NA, formula = NA, smiles = NA, inchi = NA, inchikey = NA, predicted = TRUE,
                                   lnC = log(maxC))
    lc_tibble <- rbind(lc_tibble, lc_tibble_tmp)
  }
  df <- lc_tibble %>%
    dplyr::filter(lc > 0) %>%
    dplyr::arrange(lnC)
  # df <- lc_tibble %>%
  #   dplyr::filter(lc > 0) %>%
  #   dplyr::arrange(lnC) %>%
  #   dplyr::filter(lc < maxC0 | lc == max(lc))
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = lnC, y = rt, color = predicted), shape = 19, size = 3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
    ggplot2::theme_bw() +
    ggplot2::annotate("text", x = min(df$lnC), y = max(df$rt),
                      label = text, color = "red", hjust = 0, vjust = 1, size = 6) +
    ggplot2::geom_abline(data = df,
                         ggplot2::aes(intercept = intercept, slope = slope),
                         color = "red", linewidth = 1, linetype = 2)
  lc_tibble <- lc_tibble %>%
    dplyr::select(id, name, rt, lc, ri, exactmass, formula, smiles, inchi, inchikey, predicted)
  # lc_tibble <- lc_tibble %>%
  #   dplyr::select(id, name, rt, lc, ri, exactmass, formula, smiles, inchi, inchikey, predicted) %>%
  #   dplyr::filter(lc < maxC0 | lc == max(lc))
  return(list(slope = slope, intercept = intercept, r_squared = r_squared, p = p, lc_tibble = lc_tibble))
}

#' @title deadTime_add
#' @description
#' Give car_tibble a dead time.
#'
#' @param car_tibble A car_tibble.
#' @param deadTime dead time.
#'
#' @return A car_tibble.
#' @export
#'
#' @examples
#' data("experimental_data")
#' data("experimental_car")
#' data("standard_car")
#' ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
#' experimental_car <- ResList$lc_tibble
#' ResList$p
#' ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
#' standard_car <- ResList$lc_tibble
#' ResList$p
#' experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
#' standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
deadTime_add <- function(lc_tibble, deadTime = 0){
  if(deadTime >= min(lc_tibble$rt)) stop("deadTime wrong!")
  car_tibble_tmp <- dplyr::tibble(id = NA, name = NA, rt = deadTime, lc = -1,
                                  ri = -100, exactmass = NA, formula = NA,
                                  smiles = NA, inchi = NA, inchikey = NA, predicted = FALSE)
  lc_tibble <- rbind(car_tibble_tmp, lc_tibble)
  return(lc_tibble)
}

#' @title compare_lc_tibble
#' @description
#' The minimum and maximum lc of two lc_tibble should be kept consistent when comparing them.
#'
#' @param experimental_lc experimental_lc.
#' @param standard_lc standard_lc.
#'
#' @return A list.
#' @export
#'
#' @examples
#' data("experimental_data")
#' data("experimental_car")
#' data("standard_car")
#' ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
#' experimental_car <- ResList$lc_tibble
#' ResList$p
#' ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
#' standard_car <- ResList$lc_tibble
#' ResList$p
#' experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
#' standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
#' ResList <- compare_lc_tibble(experimental_lc = experimental_car, standard_lc = standard_car)
#' experimental_car <- ResList$experimental_lc
#' standard_car <- ResList$standard_lc
compare_lc_tibble <- function(experimental_lc, standard_lc){
  experimental_maxlc <- max(experimental_lc$lc)
  standard_maxlc <- max(standard_lc$lc)
  maxlc <- ifelse(experimental_maxlc > standard_maxlc, experimental_maxlc, standard_maxlc)
  if(maxlc != experimental_maxlc){
    ResList <- lc2rt(lc_tibble = experimental_lc, targetRt = 0)
    slope <- ResList$slope;intercept <- ResList$intercept
    C_vec <- (experimental_maxlc + 1):maxlc
    rt_vec <- sapply(C_vec, function(x) {
      slope * log(x) + intercept
    })
    ri_vec <- C_vec * 100
    experimental_lc_tmp <- dplyr::tibble(id = NA, name = NA, rt = rt_vec, lc = C_vec, ri = ri_vec,
                                         exactmass = NA, formula = NA, smiles = NA, inchi = NA, inchikey = NA, predicted = TRUE)
    experimental_lc <- rbind(experimental_lc, experimental_lc_tmp) %>%
      dplyr::arrange(lc)
  }
  if(maxlc != standard_maxlc){
    ResList <- lc2rt(lc_tibble = standard_lc, targetRt = 0)
    slope <- ResList$slope;intercept <- ResList$intercept
    C_vec <- (standard_maxlc + 1):maxlc
    rt_vec <- sapply(C_vec, function(x) {
      slope * log(x) + intercept
    })
    ri_vec <- C_vec * 100
    standard_lc_tmp <- dplyr::tibble(id = NA, name = NA, rt = rt_vec, lc = C_vec, ri = ri_vec,
                                         exactmass = NA, formula = NA, smiles = NA, inchi = NA, inchikey = NA, predicted = TRUE)
    standard_lc <- rbind(standard_lc, standard_lc_tmp) %>%
      dplyr::arrange(lc)
  }
  return(list(experimental_lc = experimental_lc, standard_lc = standard_lc))
}
#' Title
#'
#' @param rt
#' A vector of experimental data rt.
#' @param experimental_lc
#' This is the lc_tibble of the correction metabolite corresponding to the experimental data.
#'
#' @return A vector of ri.
#' @export
#'
#' @examples
#' data("experimental_data")
#' data("experimental_car")
#' data("standard_car")
#' ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
#' experimental_car <- ResList$lc_tibble
#' ResList$p
#' ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
#' standard_car <- ResList$lc_tibble
#' ResList$p
#' experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
#' standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
#' ResList <- compare_lc_tibble(experimental_lc = experimental_car, standard_lc = standard_car)
#' experimental_car <- ResList$experimental_lc
#' standard_car <- ResList$standard_lc
#' experimental_data$ri <- sapply(experimental_data$rt, function(x) {
#'   rt2ri(x, experimental_lc = experimental_car)
#' })
rt2ri <- function(rt, experimental_lc){
  idx <- which(experimental_lc$rt == rt)
  if(length(idx) == 1){
    ri <- experimental_lc[idx, ]$ri
    return(ri)
  }
  experimental_lc_m <- experimental_lc[experimental_lc$rt < rt, ]
  experimental_lc_n <- experimental_lc[experimental_lc$rt > rt, ]
  RIm <- experimental_lc_m[nrow(experimental_lc_m), ]$ri
  RIn <- experimental_lc_n[1, ]$ri
  RTm <- experimental_lc_m[nrow(experimental_lc_m), ]$rt
  RTn <- experimental_lc_n[1, ]$rt
  ri <- RIm + (RIn - RIm) * ((rt - RTm) / (RTn - RTm))
  return(ri)
}
#' @title orignRT2newRT
#' @description
#' Convert experimental retention time to standard retention time.
#'
#' @param experimental_data experimental_data tibble.
#' @param experimental_lc
#' This is the lc_tibble of the correction metabolite corresponding to the experimental data.
#' @param standard_lc
#' This is the lc_tibble of the correction metabolite corresponding to the standard chromatographic condition.
#' @param thread The number of parallel thread.
#'
#' @return A new experimental data.
#' @export
#'
#' @examples
#' data("experimental_data")
#' data("experimental_car")
#' data("standard_car")
#' ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
#' experimental_car <- ResList$lc_tibble
#' ResList$p
#' ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
#' standard_car <- ResList$lc_tibble
#' ResList$p
#' experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
#' standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
#' ResList <- compare_lc_tibble(experimental_lc = experimental_car, standard_lc = standard_car)
#' experimental_car <- ResList$experimental_lc
#' standard_car <- ResList$standard_lc
#' experimental_data$ri <- sapply(experimental_data$rt, function(x) {
#'   rt2ri(x, experimental_lc = experimental_car)
#' })
#' experimental_data <- orignRT2newRT(experimental_data = experimental_data,
#'                                    experimental_lc = experimental_car,
#'                                    standard_lc = standard_car, thread = 1)
orignRT2newRT <- function(experimental_data, experimental_lc, standard_lc, thread = 1){
  feature_number <- nrow(experimental_data)
  message(paste0("Number of features: ", feature_number))

  loop <- function(i){
    tmp_tibble <- experimental_data[i, ]
    orign_rt <- tmp_tibble$rt
    ri <- tmp_tibble$ri
    experimental_lc_tibble_m <- experimental_lc[experimental_lc$rt <= orign_rt, ]
    experimental_lc_tibble_n <- experimental_lc[experimental_lc$rt > orign_rt, ]
    RIm <- experimental_lc_tibble_m[nrow(experimental_lc_tibble_m), ]$ri
    lc_m <- experimental_lc_tibble_m[nrow(experimental_lc_tibble_m), ]$lc
    RIn <- experimental_lc_tibble_n[1, ]$ri
    lc_n <- experimental_lc_tibble_n[1, ]$lc
    standard_lc_tibble_m <- standard_lc[standard_lc$lc == lc_m, ]
    standard_lc_tibble_n <- standard_lc[standard_lc$lc == lc_n, ]
    standard_RTm <- standard_lc_tibble_m[1, ]$rt
    standard_RTn <- standard_lc_tibble_n[1, ]$rt
    new_rt <- standard_RTm + (standard_RTn - standard_RTm) * ((ri - RIm) / (RIn - RIm))
    return(new_rt)
  }

  pb <- utils::txtProgressBar(max = feature_number, style = 3)
  if(thread == 1){
    new_rt_list <- lapply(1:feature_number, function(i) {
      utils::setTxtProgressBar(pb, i)
      loop(i)
    })
  }else if(thread > 1){
    cl <- snow::makeCluster(thread)
    doSNOW::registerDoSNOW(cl)
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb,
                                                                 n))
    new_rt_list <- foreach::`%dopar%`(foreach::foreach(i = 1:feature_number,
                                                       .options.snow = opts,
                                                       .packages = c("dplyr")))
    snow::stopCluster(cl)
    gc()
  }else stop("thread is wrong!")

  new_rt_vec <- purrr::list_c(new_rt_list)
  colnames(experimental_data)[2] <- "orign_rt"
  experimental_data$new_rt <- new_rt_vec
  return(experimental_data)
}
