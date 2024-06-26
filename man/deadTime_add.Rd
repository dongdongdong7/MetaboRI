% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function.R
\name{deadTime_add}
\alias{deadTime_add}
\title{deadTime_add}
\usage{
deadTime_add(lc_tibble, deadTime = 0)
}
\arguments{
\item{deadTime}{dead time.}

\item{car_tibble}{A car_tibble.}
}
\value{
A car_tibble.
}
\description{
Give car_tibble a dead time.
}
\examples{
data("experimental_data")
data("experimental_car")
data("standard_car")
ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
experimental_car <- ResList$lc_tibble
ResList$p
ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
standard_car <- ResList$lc_tibble
ResList$p
experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
}
