#' @title predictCar
#' @description
#' Construct acylcarnitine prediction model and predict.
#'
#' @param car_tibble A car_tibble.
#'
#' @return A list.
#' @export
#'
#' @examples
#' predictCar_res <- predictCar(car_tibble = standard_car)
predictCar <- function(car_tibble){
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
  return(list(slope = slope, intercept = intercept, r_squared = r_squared, car_tibble = car_tibble))
}
