% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_libraries.R
\name{read_lib}
\alias{read_lib}
\title{Read msp/mgf mass spectral libraries}
\usage{
read_lib(
  file,
  format = "msp",
  type = "MS2",
  remove_ri = TRUE,
  remove_rt = TRUE
)
}
\arguments{
\item{file}{Mass spectral library in \code{msp} or \code{mgf} format.}

\item{format}{The format of the library, either \code{msp} or \code{mgf}.}

\item{type}{The type of the library, either \code{EI} or \code{MS2}.}

\item{remove_ri}{A logical scalar only used in case of EI mass spectral
library. Should retention index (RI) be removed? \code{TRUE} or \code{FALSE}}

\item{remove_rt}{A logical scalar only used in case of MS2 mass spectral
library. Should retention time (RT) be removed? \code{TRUE} or \code{FALSE}}
}
\value{
A \code{list} with each spectral entry as a list element for further
processing.
}
\description{
\code{read_lib} offers a way to read mass spectral libraries into R
for further processing.
}
\details{
This is a generic function to read either EI or MS2 mass spectral libraries.
The library can be either in \code{msp} or \code{mgf} form. For this reason,
it is required to set the format and the type of the input library. The
default is \code{MS2} in \code{msp} format. In the case of EI mass spectral
library, an additional Boolean parameter \code{remove_ri} can be set to
remove or keep the retention index (RI). In the case of MS2 mass spectral
library, an additional Boolean parameter \code{remove_rt} can be set to
remove or keep the retention time (RT). This function supports parallel
computing making use of the \pkg{future.apply}. Please see the vignette
for more details.
}
\examples{
# The first 2 lines only indicate the location where the example files are
# stored. You might not need them.
EI_file <- system.file("EI.msp", package = "mspcompiler")
MS2_mgf_file <- system.file("MS2.mgf", package = "mspcompiler")

EI <- read_lib(file = EI_file, format = "msp", type = "EI", remove_ri = FALSE)
MS2_mgf <- read_lib(file = MS2_mgf_file, format = "mgf", type = "MS2")
}
