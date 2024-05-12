library(magrittr)
standard_car <- dplyr::as_tibble(openxlsx::read.xlsx("D:/fudan/Projects/2023/MetaboRI/Progress/database/demo/tidy/standard_car.xlsx", sheet = 1))
standard_car$lnC <- log(as.integer(stringr::str_replace(standard_car$carnitine, "C", "")))
predicted_car <- standard_car %>%
  dplyr::filter(predicted) %>%
  dplyr::filter(carnitine != "C0")
model_car <- standard_car %>%
  dplyr::filter(!predicted) %>%
  dplyr::filter(carnitine != "C0")
ggplot2::ggplot(model_car) +
  ggplot2::geom_point(ggplot2::aes(x = lnC, y = rt))
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
ggplot2::ggplot(df) +
  ggplot2::geom_point(ggplot2::aes(x = lnC, y = rt, color = predicted), shape = 19, size = 3) +
  ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
  ggplot2::theme_bw() +
  ggplot2::annotate("text", x = min(df$lnC), y = max(df$rt),
                    label = text, color = "red", hjust = 0, vjust = 1, size = 6) +
  ggplot2::geom_abline(data = df,
                       ggplot2::aes(intercept = intercept, slope = slope),
                       color = "red", linewidth = 1, linetype = 2)

experimental_car <- dplyr::as_tibble(openxlsx::read.xlsx("D:/fudan/Projects/2023/MetaboRI/Progress/database/demo/tidy/experimental_car.xlsx", sheet = 1))
