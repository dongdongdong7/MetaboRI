#' @title predictCar
#' @description
#' Construct acylcarnitine prediction model and predict.
#'
#' @param car_tibble A car_tibble.
#' @param targetRt Target retention time, the maximum retention time you want to predict.
#'
#' @return A list.
#' @export
#'
#' @examples
#' data("experimental_car", package = "MetaboRI")
#' predictCar_res <- predictCar(car_tibble = experimental_car)
#' predictCar_res$slope
#' predictCar_res$intercept
#' predictCar_res$r_squared
#' predictCar_res$car_tibble
#' predictCar_res$p
predictCar <- function(car_tibble, targetRt){
  car_tibble$lnC <- log(as.integer(stringr::str_replace(car_tibble$carnitine, "C", "")))
  predicted_car <- car_tibble %>%
    dplyr::filter(predicted) %>%
    dplyr::filter(carnitine != "C0")
  model_car <- car_tibble %>%
    dplyr::filter(!predicted) %>%
    dplyr::filter(carnitine != "C0")
  lmRes <- lm(formula = rt ~ lnC, data = model_car)
  intercept <- as.numeric(coefficients(lmRes)[1])
  slope <- as.numeric(coefficients(lmRes)[2])
  r_squared <- summary(lmRes)$r.squared
  sign <- ifelse(intercept < 0, "-", "+")
  text <- paste0("y = ", round(slope, 4), "lnC ", sign, " ", abs(round(intercept, 4)), "\n",
                 "R2 = ", round(r_squared, 4))
  predicted_car$rt <- slope * predicted_car$lnC + intercept
  df <- rbind(model_car, predicted_car) %>%
    dplyr::arrange(lnC)
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = lnC, y = rt, color = predicted), shape = 19, size = 3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
    ggplot2::theme_bw() +
    ggplot2::annotate("text", x = min(df$lnC), y = max(df$rt),
                      label = text, color = "red", hjust = 0, vjust = 1, size = 6) +
    ggplot2::geom_abline(data = df,
                         ggplot2::aes(intercept = intercept, slope = slope),
                         color = "red", linewidth = 1, linetype = 2)
  car_tibble <- car_tibble %>% dplyr::filter(carnitine == "C0")
  car_tibble <- rbind(car_tibble, model_car, predicted_car) %>%
    dplyr::select(id, name, rt, carnitine, ri, exactmass, formula, smiles, inchi, inchikey, predicted) %>%
    dplyr::arrange(ri)
  maxRt <- max(car_tibble$rt)
  C <- 26
  while(maxRt < targetRt){
    C <- C + 1
    maxRt <- slope * log(C) + intercept
    car_tibble_tmp <- dplyr::tibble(id = NA, name = NA, rt = maxRt, carnitine = paste0("C", C), ri = C * 100,
                                    exactmass = NA, formula = NA, smiles = NA, inchi = NA, inchikey = NA, predicted = TRUE)
    car_tibble <- rbind(car_tibble, car_tibble_tmp)
  }
  df <- car_tibble %>%
    dplyr::filter(carnitine != "C0") %>%
    dplyr::mutate(lnC = log(as.integer(stringr::str_replace(carnitine, "C", "")))) %>%
    dplyr::arrange(lnC)
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_point(ggplot2::aes(x = lnC, y = rt, color = predicted), shape = 19, size = 3) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
    ggplot2::theme_bw() +
    ggplot2::annotate("text", x = min(df$lnC), y = max(df$rt),
                      label = text, color = "red", hjust = 0, vjust = 1, size = 6) +
    ggplot2::geom_abline(data = df,
                         ggplot2::aes(intercept = intercept, slope = slope),
                         color = "red", linewidth = 1, linetype = 2)
  return(list(slope = slope, intercept = intercept, r_squared = r_squared, p = p, car_tibble = car_tibble))
}
rt2ri_old <- function(rt, car_tibble){
  idx <- which(car_tibble$rt == rt)
  if(length(idx) == 1){
    ri <- car_tibble[idx, ]$ri
    return(ri)
  }
  car_tibble_m <- car_tibble[car_tibble$rt < rt, ]
  car_tibble_n <- car_tibble[car_tibble$rt > rt, ]
  RIm <- car_tibble_m[nrow(car_tibble_m), ]$ri
  RIn <- car_tibble_n[1, ]$ri
  RTm <- car_tibble_m[nrow(car_tibble_m), ]$rt
  RTn <- car_tibble_n[1, ]$rt
  ri <- RIm + (RIn - RIm) * ((rt - RTm) / (RTn - RTm))
  return(ri)
}
#' Title
#'
#' @param experimental_data test.
#' @param experimental_car test.
#' @param standard_car test.
#' @param thread test.
#'
#' @return test.
#' @export
#'
#' @examples "A test"
orignRT2newRT_old <- function(experimental_data, experimental_car, standard_car, thread = 1){
  feature_number <- nrow(experimental_data)
  message(paste0("Number of features: ", feature_number))

  loop <- function(i){
    tmp_tibble <- experimental_data[i, ]
    orign_rt <- tmp_tibble$rt
    ri <- tmp_tibble$ri
    experimental_car_tibble_m <- experimental_car[experimental_car$rt <= orign_rt, ]
    experimental_car_tibble_n <- experimental_car[experimental_car$rt > orign_rt, ]
    RIm <- experimental_car_tibble_m[nrow(experimental_car_tibble_m), ]$ri
    carnitine_m <- experimental_car_tibble_m[nrow(experimental_car_tibble_m), ]$carnitine
    RIn <- experimental_car_tibble_n[1, ]$ri
    carnitine_n <- experimental_car_tibble_n[1, ]$carnitine
    standard_car_tibble_m <- standard_car[standard_car$carnitine == carnitine_m, ]
    standard_car_tibble_n <- standard_car[standard_car$carnitine == carnitine_n, ]
    standard_RTm <- standard_car_tibble_m[1, ]$rt
    standard_RTn <- standard_car_tibble_n[1, ]$rt
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
