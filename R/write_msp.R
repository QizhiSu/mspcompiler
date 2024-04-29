#' Write EI library into a msp file
#'
#' \code{write_EI_msp} offers a way to write the organized  EI library into a
#' msp file.
#'
#' @param lib The organized EI library
#' @param filename The location and filename of the msp file to be exported,
#'   e.g., "/home/exported.msp".
#'
#' @return No return but create a msp file
#' @export
write_EI_msp <- function(lib, filename) { # nolint
  zz <- file(description = filename, open = "w")
  for (i in 1:length(lib)) {
    cat(paste("Name: ", lib[[i]]$Name, "\n", sep = ""), file = zz)
    cat(paste("InChIKey: ", lib[[i]]$InChIKey, "\n", sep = ""), file = zz)
    cat(paste("SMILES: ", lib[[i]]$Smiles, "\n", sep = ""), file = zz)
    cat(paste("Formula: ", lib[[i]]$Formula, "\n", sep = ""), file = zz)
    cat(paste("MW: ", lib[[i]]$"Molecular weight", "\n", sep = ""), file = zz)
    cat(paste("RI: ", lib[[i]]$RI, "\n", sep = ""), file = zz)
    cat(paste("Comment: ", lib[[i]]$Comment, "\n", sep = ""), file = zz)
    cat(paste("Num peaks: ", lib[[i]]$"Number of peaks", "\n", sep = ""), file = zz)
    cat(paste(paste(lib[[i]]$Spectra$mz, lib[[i]]$Spectra$ins, sep = " "), "", sep = ""),
      file = zz, sep = "\n"
    )
    cat("\r\n", file = zz)
  }
  close(zz)
}


#' Write MS2 library into a msp file
#'
#' \code{write_EI_msp} offers a way to write the organized  MS2 library into a
#' msp file.
#'
#' @param lib The organized MS2 library
#' @param filename The location and filename of the msp file to be exported,
#'   e.g., "/home/exported.msp".
#'
#' @return No return but create a msp file
#' @export
write_MS2_msp <- function(lib, filename) {
  zz <- file(description = filename, open = "w")
  for (i in 1:length(lib)) {
    cat(paste("Name: ", lib[[i]]$Name, "\n", sep = ""), file = zz)
    cat(paste("PrecursorMZ: ", lib[[i]]$PrecursorMZ, "\n", sep = ""), file = zz)
    cat(paste("PrecursorType: ", lib[[i]]$PrecusorType, "\n", sep = ""), file = zz)
    cat(paste("IonMode: ", lib[[i]]$IonMode, "\n", sep = ""), file = zz)
    cat(paste("Formula: ", lib[[i]]$Formula, "\n", sep = ""), file = zz)
    cat(paste("SMILES: ", lib[[i]]$Smiles, "\n", sep = ""), file = zz)
    cat(paste("InChIKey: ", lib[[i]]$InChIKey, "\n", sep = ""), file = zz)
    cat(paste("RetentionTime: ", lib[[i]]$RetentionTime, "\n", sep = ""), file = zz)
    cat(paste("CCS: ", lib[[i]]$CCS, "\n", sep = ""), file = zz)
    cat(paste("CollisionEnergy: ", lib[[i]]$CollisionEnergy, "\n", sep = ""), file = zz)
    cat(paste("InstrumentType: ", lib[[i]]$InstrumentType, "\n", sep = ""), file = zz)
    cat(paste("Comment: ", lib[[i]]$Comment, "\n", sep = ""), file = zz)
    cat(paste("Num peaks: ", lib[[i]]$"Number of peaks", "\n", sep = ""), file = zz)
    cat(paste(paste(lib[[i]]$Spectra$mz, lib[[i]]$Spectra$ins, sep = " "), "", sep = ""),
      file = zz, sep = "\n"
    )
    cat("\r", file = zz)
  }
  close(zz)
}
