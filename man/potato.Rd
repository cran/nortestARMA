\name{potato}
\alias{potato}
\docType{data}
\title{
Prince Edward Island; Average Yield, Potatoes (Hundredweight per harvested acres (1957-2014))
}
\description{
We retrieved the data from the Canadian Socio-Economic Information Management System
(CANSIM), which is a computerized database of Statistics Canada. The time series
consists of annual data on the productivity of potatoes (per acre) in Prince Edward
Island for the time period 1957-2014 (CANSIM Series V47152).
}
\usage{data(potato)}
\format{
  A data frame with 58 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{numeric vector of the years of harvest}
    \item{\code{potatoeyield}}{numeric vector of the productivity of
      potatoes (Hundredweight per harvested acres)}
  }
}
\references{
Cheng H., (2005). Competitive relationship among potato production areas in northeastern
america. \emph{Journal of Food Distribution Research Proceedings} 36 (1), 27--32.

Ducharme G., Lafaye de Micheaux P., (2004). Goodness-of-fit tests of normality for the
innovations in ARMA models. \emph{Journal of Time Series Analysis} 25 (3), 373--395.

Duchesne P., Lafaye de Micheaux P., Tagne J. (2016). Neyman Smooth Test of Normality for ARMA Time Series
Models with Unknown Mean. Submitted.
}
\author{
Duchesne P., Lafaye de Micheaux P.
}
\examples{
data(potato)
}
\keyword{datasets}
