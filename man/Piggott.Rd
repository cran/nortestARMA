\name{Piggott}
\alias{Piggott}
\alias{resid.temp}
\docType{data}
\title{
Daily gas demand study by Piggott
}
\description{
Two series consisting of daily temperature values and the (square root
transformed) wind speed series obtained from Piggott (1980). This data
set was first studied by Shea (1987) and then by Ducharme et al. (2004).
}
\usage{data(Piggott)}
\format{
  A data frame with 366 observations on the following 2 variables.
  \describe{
    \item{\code{temp}}{numeric vector of daily temperatures.}
    \item{\code{wind.trans}}{numeric vector of the (square root
transformed) wind speed series.}
  }
}
\references{
Ducharme G., Lafaye de Micheaux P., (2004). Goodness-of-fit tests of normality for the
innovations in ARMA models. \emph{Journal of Time Series Analysis} 25 (3), 373--395.

Duchesne P., Lafaye de Micheaux P., Tagne J. (2016). Neyman Smooth Test of Normality for ARMA Time Series
Models with Unknown Mean. Submitted.

Piggott J. L., (1980). The use o f Box-Jenkins modelling for the forecasting of daily gas demand.
Paper presented to the \emph{Royal Statistical Society}.

Shea B. L., (1987). Estimation of Multivariate Time
Series. \emph{Journal of Time Series Analysis} 8 (1), 95--109.
}
\author{
Duchesne P., Lafaye de Micheaux P.
}
\examples{
data(Piggott)
}
\keyword{datasets}
