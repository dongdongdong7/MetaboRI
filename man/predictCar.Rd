% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/endoRI.R
\name{predictCar}
\alias{predictCar}
\title{predictCar}
\usage{
predictCar(car_tibble, targetRt)
}
\arguments{
\item{car_tibble}{A car_tibble.}

\item{targetRt}{Target retention time, the maximum retention time you want to predict.}
}
\value{
A list.
}
\description{
Construct acylcarnitine prediction model and predict.
}
\examples{
data("experimental_car", package = "MetaboRI")
predictCar_res <- predictCar(car_tibble = experimental_car)
predictCar_res$slope
predictCar_res$intercept
predictCar_res$r_squared
predictCar_res$car_tibble
predictCar_res$p
}
