% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organize_libraries_MS2.R
\name{separate_polarity}
\alias{separate_polarity}
\title{Separate positive and negative modes in MS2 library}
\usage{
separate_polarity(lib, polarity)
}
\arguments{
\item{lib}{A MS2 library mixed with positive and negative modes.}

\item{polarity}{The polarity, can be either "pos" or "neg"}
}
\value{
A \code{list}, being positive or negative
}
\description{
\code{separate_polarity} offers a way to separate a MS2 library based on
polarity.
}
\details{
Some libraries, e.g., NIST and GNPS have both positive and negative MS2
records mixed in a singled file. However, in practice, a MS2 library should
be either positive or negative based on the polarity used in the experiment.
Therefore, this function provides a way to separate positive and negative
modes in MS2 library.
}
