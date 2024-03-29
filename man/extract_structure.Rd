% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_structures.R
\name{extract_structure}
\alias{extract_structure}
\title{Extract SMILES from the sdf file generated by \code{\link{combine_mol2sdf}}}
\usage{
extract_structure(input, output)
}
\arguments{
\item{input}{The sdf file generated by \code{\link{combine_mol2sdf}}, e.g.,
"/home/exported.sdf".}

\item{output}{The location where the structure information will be stored
and its name, e.g., "/home/exported.txt".}
}
\value{
A data.frame and creates a *.txt file.
}
\description{
\code{extract_structure} offers a way to retrieve SMILES from the sdf file.
}
\details{
The function is a wrapper of the \code{convertFormatFile} function from the
\pkg{ChemmineOB} package. As InChI and InChIKey are not supported in Windows-
based systems, this function will automatically determine which type of
operating system you are working with. Only \strong{name} and \strong{SMILES}
will be retrieved if you work with Windows, while \strong{InChI} and
\strong{InChIKey} will be exported as well in Linux-based or Mac OS systems.
}
