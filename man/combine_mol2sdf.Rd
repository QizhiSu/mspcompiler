% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_structures.R
\name{combine_mol2sdf}
\alias{combine_mol2sdf}
\title{Combine multiple mol files into a single sdf file and remove duplicates}
\usage{
combine_mol2sdf(input, output, use_filename = FALSE)
}
\arguments{
\item{input}{The location of the exported *.MOL folder from Lib2NIST,
e.g., "/home/nist.MOL".}

\item{output}{The location where the sdf file will be stored and its name,
e.g., "/home/exported.sdf".}

\item{use_filename}{In case you want to use the file name as the Molecule_Name
in the sdf file, please use \code{use_filename = TRUE}. This is useful when
you draw your own chemicals which might not have Molecule_Name in the .MOL files.
With this option, you can use the name of the .MOL files.}
}
\value{
It will return no value but only creates a sdf file.
}
\description{
\code{combine_mol2sdf} offers a way to combine multiple mol files in to
a single sdf file removing duplicates.
}
\details{
The msp file exported from NIST format by the \strong{Lib2NIST} software
has no SMILES which is used for viewing the structure in MS-DIAL and is
crucial as well if you want to predict retention index (RI) for compounds
with no experimental RI (in the case of EI library). Chemical structure of
each compound in the NIST library can be exported as a mol file. Every
library entry will have one corresponding mol file. All these mol files will
be stored in a folder with ".MOL" suffix. This function aims to combine all
these mol files into a single sdf file which can be then used for retrieve
SMILES for each entry. This function supports parallel computing.
}
