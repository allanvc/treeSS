\name{fl_deaths}
\alias{fl_deaths}
\docType{data}
\title{Florida General Mortality Data (2016)}
\description{
Death counts by county and ICD-10 cause of death for the state of Florida,
USA, 2016. This is \strong{raw data} -- the user builds the ICD-10 tree
and cases matrix in the analysis script, making it a pedagogical example
of the full workflow.
}
\format{
A data.frame with 3,066 rows and 6 columns:
\describe{
  \item{county_fips}{5-digit FIPS code (character), e.g. \code{"12086"}.}
  \item{county_name}{County name, e.g. \code{"Miami-Dade County"}.}
  \item{icd10_code}{ICD-10 cause of death code, e.g. \code{"I25.1"}.}
  \item{icd10_desc}{Description of the ICD-10 code, e.g.
    \code{"Atherosclerotic heart disease"}.}
  \item{deaths}{Number of deaths (integer).}
  \item{population}{County population estimate (integer).}
}
}
\details{
Covers 65 of Florida's 67 counties, 253 ICD-10 codes, and 157,000 total
deaths. All ages, all causes.

County polygons and centroids can be obtained via
\code{tigris::counties(state = "FL", cb = TRUE)}.

The ICD-10 descriptions in \code{icd10_desc} can be used to build a
lookup table: \code{unique(fl_deaths[, c("icd10_code", "icd10_desc")])}.
}
\source{
CDC WONDER Compressed Mortality File 1999--2016
(\url{https://wonder.cdc.gov/cmf-icd10.html}).
}
\examples{
data(fl_deaths)
head(fl_deaths)
cat("Counties:", length(unique(fl_deaths$county_fips)), "\n")
cat("ICD-10 codes:", length(unique(fl_deaths$icd10_code)), "\n")
cat("Total deaths:", sum(fl_deaths$deaths), "\n")
}
\keyword{datasets}
