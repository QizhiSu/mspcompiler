#' Read msp/mgf mass spectral libraries
#'
#' \code{read_lib} offers a way to read mass spectral libraries into R
#' for further processing.
#'
#' This is a generic function to read either EI or MS2 mass spectral libraries.
#' The library can be either in \code{msp} or \code{mgf} form. For this reason,
#' it is required to set the format and the type of the input library. The
#' default is \code{MS2} in \code{msp} format. In the case of EI mass spectral
#' library, an additional Boolean parameter \code{remove_ri} can be set to
#' remove or keep the retention index (RI). In the case of MS2 mass spectral
#' library, an additional Boolean parameter \code{remove_rt} can be set to
#' remove or keep the retention time (RT). This function supports parallel
#' computing making use of the \pkg{future.apply}. Please see the vignette
#' for more details.
#'
#' @param file Mass spectral library in \code{msp} or \code{mgf} format.
#' @param format The format of the library, either \code{msp} or \code{mgf}.
#' @param type The type of the library, either \code{EI} or \code{MS2}.
#' @param remove_ri A logical scalar only used in case of EI mass spectral
#'   library. Should retention index (RI) be removed? \code{TRUE} or \code{FALSE}
#' @param remove_rt A logical scalar only used in case of MS2 mass spectral
#'   library. Should retention time (RT) be removed? \code{TRUE} or \code{FALSE}
#'
#' @return A \code{list} with each spectral entry as a list element for further
#'   processing.
#'
#' @import stringr
#' @import future.apply
#'
#' @export
#'
#' @examples
#' # The first 2 lines only indicate the location where the example files are
#' # stored. You might not need them.
#' EI_file <- system.file("EI.msp", package = "mspcompiler")
#' MS2_mgf_file <- system.file("MS2.mgf", package = "mspcompiler")
#'
#' EI <- read_lib(file = EI_file, format = "msp", type = "EI", remove_ri = FALSE)
#' MS2_mgf <- read_lib(file = MS2_mgf_file, format = "mgf", type = "MS2")
read_lib <-
  function(file, format = "msp", type = "MS2",
           remove_ri = TRUE, remove_rt = TRUE) {
    tmp <- readLines(file)
    tmp <- gsub("ï»¿", "", tmp) # fix files that have different encoding
    # Individual compounds are recognized differently
    # depending on the format
    if (format == "msp") {
      start_line <- grep('^name:', tmp, ignore.case = TRUE)
    } else {
      start_line <- grep('^begin ions', tmp, ignore.case = TRUE)
    }
    num_line <- diff(c(start_line, length(tmp) + 1))
    split_factor <- rep(1:length(start_line), num_line)
    cmp_list <- split(tmp, split_factor)

    # For EI spectral libraries
    if (type == "EI") {
      get_msp <- function(cmp) {
        name <- cmp[grep('^name:', cmp, ignore.case = TRUE)]
        name <- gsub('^name: ', '', name, ignore.case = TRUE)
        smiles <- cmp[grep('^smiles:', cmp, ignore.case = TRUE)]
        smiles <- gsub('^smiles: ', '', smiles, ignore.case = TRUE)
        inchikey <- cmp[grep('^inchikey:', cmp, ignore.case = TRUE)]
        inchikey <- gsub('inchikey: ', '', inchikey, ignore.case = TRUE)
        formula <- cmp[grep('^formula:', cmp, ignore.case = TRUE)]
        formula <- gsub('^formula: ', '', formula, ignore.case = TRUE)
        mw <- cmp[grep('^mw:',  cmp, ignore.case = TRUE)]
        mw <- gsub('^mw: ', '', mw, ignore.case = TRUE)
        comment <- cmp[grep('^comments?:', cmp, ignore.case = TRUE)]
        comment <- gsub('^comments?: ', '', comment, ignore.case = TRUE)
        if (remove_ri) {
          RI <- NA
        } else {
          RI <- cmp[grep('^retentionindex:', cmp, ignore.case = TRUE)]
          RI <- gsub('^retentionindex: ', '', RI, ignore.case = TRUE)
          RI <- round(as.numeric(RI))
        }
        # Dealing with the spectrum
        peak_numbers <- cmp[grep('^num peaks:', cmp, ignore.case = TRUE)]
        peak_number <- gsub('^num peaks: ?', '', peak_numbers, ignore.case = TRUE)
        # matrix of masses and intensities
        if (as.numeric(peak_number) > 0) {
          # Determine position of mass intensity pairs
          mass_inten_posi <- which(grepl('^[0-9]', cmp) & !grepl(': ', cmp))
          # Turn mass intensity pairs to a numeric vector
          mass_inten <- str_remove_all(cmp[mass_inten_posi], '\n')
          mass_inten <- str_remove(mass_inten, '".*"$')
          mass_inten <- unlist(strsplit(mass_inten, '\t| '))
          mass_inten <-
            as.numeric(mass_inten[grep('^[0-9].*[0-9]$|^[0-9]$', mass_inten)])
          # Extract mz and intensity
          mz <- mass_inten[seq(1, length(mass_inten), 2)]
          intensity <-  mass_inten[seq(2, length(mass_inten), 2)]
          spectra <- cbind.data.frame(mz = mz, ins = intensity)
          return(
            list(
              Name = name,
              InChIKey = inchikey,
              Smiles = smiles,
              Formula = formula,
              'Molecular weight' = mw,
              RI = RI,
              Comment = comment,
              'Number of peaks'  =  peak_number,
              Spectra = spectra
            )
          )
        } else {
          return(
            list(
              Name = name,
              InChIKey = inchikey,
              Smiles = smiles,
              'Molecular weight' = mw,
              RI = RI,
              Comment = comment,
              Formula = formula,
              'Number of peaks' = peak_number
            )
          )
        }
      }
      # For tandem spectral libraries, treatment varies depending on the format
      # For msp format
    } else {
      get_msp <- function(cmp) {
        if (format == "msp") {
          name <- cmp[grep('^name:', cmp, ignore.case = TRUE)]
          name <- gsub('^name: ', '', name, ignore.case = TRUE)
          smiles <- cmp[grep('^smiles:', cmp, ignore.case = TRUE)]
          smiles <- gsub('^smiles: ', '', smiles, ignore.case = TRUE)
          inchikey <- cmp[grep('^inchikey:',  cmp, ignore.case = TRUE)]
          inchikey <- gsub('inchikey: ', '', inchikey, ignore.case = TRUE)
          formula <- cmp[grep('^formula:', cmp, ignore.case = TRUE)]
          formula <- gsub('^formula: ', '', formula, ignore.case = TRUE)
          precursor_ion <- cmp[grep('^precursormz:',  cmp, ignore.case = TRUE)]
          precursor_ion <- gsub('^precursormz: ', '',
                                precursor_ion, ignore.case = TRUE)
          precursor_type <-cmp[grep('^precursor_?type:',
                                    cmp, ignore.case = TRUE)]
          precursor_type <-gsub('^precursor_?type: ', '',
                                precursor_type, ignore.case = TRUE)
          ion_mode <- cmp[grep('ion_?mode:', cmp, ignore.case = TRUE)]
          ion_mode <- gsub('ion_?mode: ', '', ion_mode, ignore.case = TRUE)
          if (remove_rt) {
            retention_time <- NA
          } else {
            retention_time <- cmp[grep('^retention_?time:',
                                       cmp, ignore.case = TRUE)]
            retention_time <- gsub('^retention_?time: ', '',
                                   retention_time, ignore.case = TRUE)
          }
          ccs <- cmp[grep('^ccs:', cmp, ignore.case = TRUE)]
          ccs <- gsub('^ccs: ', '', ccs, ignore.case = TRUE)
          collision_energy <- cmp[grep('^collision_?energy:',
                                       cmp, ignore.case = TRUE)]
          collision_energy <- gsub('^collision_?energy: ', '',
                                   collision_energy, ignore.case = TRUE)
          instrument_type <- cmp[grep('^instrument_?type:',
                                      cmp, ignore.case = TRUE)]
          instrument_type <- gsub('^instrument_?type: ', '',
                                  instrument_type, ignore.case = TRUE)
          comment <- cmp[grep('^comments?:', cmp, ignore.case = TRUE)]
          comment <- gsub('^comments?: ', '', comment, ignore.case = TRUE)
          peak_numbers <- cmp[grep('^num peaks:', cmp, ignore.case = TRUE)]
          peak_number <- gsub('^num peaks: ', '',
                              peak_numbers, ignore.case = TRUE)

        # For mgf format
        } else {
          name <- cmp[grep('^name=', cmp, ignore.case = TRUE)]
          name <- gsub('^name=', '', name, ignore.case = TRUE)
          name <- gsub(' \\[?[0-9]?M(-|\\+).*$', '', name, ignore.case = TRUE)
          smiles <- cmp[grep('^smiles=', cmp, ignore.case = TRUE)]
          smiles <- trimws(gsub('^smiles=', '', smiles, ignore.case = TRUE))
          inchikey <- cmp[grep('^inchiaux=', cmp, ignore.case = TRUE)]
          inchikey <- gsub('^inchiaux=', '', inchikey, ignore.case = TRUE)
          formula <- NA
          precursor_ion <- cmp[grep('^pepmass=',  cmp, ignore.case = TRUE)]
          precursor_ion <- gsub('^pepmass=', '',
                                precursor_ion, ignore.case = TRUE)
          ion_mode <- cmp[grep('^ionmode=', cmp, ignore.case = TRUE)]
          ion_mode <- gsub('^ionmode=', '', ion_mode, ignore.case = TRUE)
          retention_time <- NA
          ccs <- NA
          collision_energy <- NA
          instrument_type <- NA
          precursor_type <- str_extract(name, "\\[?[0-9]?M(-|\\+).*$")
          if (grepl('\\](-\\+)+$', precursor_type)) {
            precursor_type <- precursor_type
          } else{
            if (grepl('positive', ion_mode, ignore.case = TRUE)) {
              if (grepl('\\[', precursor_type)) {
                precursor_type <- paste0(precursor_type, "+")
              } else{
                precursor_type <- paste0("[", precursor_type, "]+")
              }
            } else{
              if (grepl('\\[', precursor_type)) {
                precursor_type <- paste0(precursor_type, "-")
              } else{
                precursor_type <- paste0("[", precursor_type, "]-")
              }
            }
          }
          PI <- cmp[grep('PI=', cmp, ignore.case = TRUE)]
          collector <- cmp[grep('^datacollector=', cmp, ignore.case = TRUE)]
          submit_user <- cmp[grep('^submituser=', cmp, ignore.case = TRUE)]
          spectrum_type <- cmp[grep('^mslevel', cmp, ignore.case = TRUE)]
          library_quality <- cmp[grep('^libraryquality=',
                                      cmp, ignore.case = TRUE)]
          comment <- paste(PI,
                           collector,
                           submit_user,
                           spectrum_type,
                           library_quality,
                           sep = "; ")
          peak_number <- length(which(grepl('^[0-9]', cmp)))
        }

        # Manipulation of peak matrix are the same for both msp and mgf
        if (as.numeric(peak_number) > 0) {
          # Determine position of mass intensity pairs
          mass_inten_posi <- which(grepl('^[0-9]', cmp) & !grepl(': ', cmp))
          # Turn mass intensity pairs to a numeric vector
          mass_inten <- str_remove_all(cmp[mass_inten_posi], '\n')
          mass_inten <- str_remove(mass_inten, '".*"$')
          mass_inten <- unlist(strsplit(mass_inten, '\t| '))
          mass_inten <-
            as.numeric(mass_inten[grep('^[0-9].*[0-9]$|^[0-9]$', mass_inten)])
          # Extract mz and intensity
          mz <- mass_inten[seq(1, length(mass_inten), 2)]
          intensity <-  mass_inten[seq(2, length(mass_inten), 2)]
          spectra <- cbind.data.frame(mz = mz, ins = intensity)
          return(
            list(
              Name = name,
              PrecursorMZ = precursor_ion,
              PrecusorType = precursor_type,
              IonMode = ion_mode,
              Formula = formula,
              Smiles = smiles,
              InChIKey = inchikey,
              RetentionTime =  retention_time,
              CCS = ccs,
              CollisionEnergy = collision_energy,
              InstrumentType = instrument_type,
              Comment = comment,
              'Number of peaks' = peak_number,
              Spectra = spectra
            )
          )
        } else {
          return(
            list(
              Name = name,
              PrecursorMZ = precursor_ion,
              PrecusorType = precursor_type,
              IonMode = ion_mode,
              Formula = formula,
              Smiles = smiles,
              InChIKey = inchikey,
              RetentionTime =  retention_time,
              CCS = ccs,
              CollisionEnergy = collision_energy,
              InstrumentType = instrument_type,
              Comment = comment,
              'Number of peaks' = peak_number
            )
          )
        }
      }
    }

    # Apply get_msp to the cmp_list using multiple sessions
    cmp_list <- future.apply::future_lapply(cmp_list, get_msp)

    # Return the organized list
    return(cmp_list)
  }


#' A wrapper to read multiple msp files at a time
#'
#' \code{read_multilibs} offers a way to read multiple msp files at a time and
#' combine them into a single file.
#'
#' When you are building your in-house libraries, you may probably have multiple
#' msp files at hand (e.g., one msp for one group of compounds). To avoid
#' empolying \code{read_lib} several times, this function provides a way to read
#' all these files at once.
#'
#' @param folder The folder that contains multiple msp files.
#'
#' @return A single \code{list} combining all msp files
#' @export
#'
read_multilibs <- function(folder) {
  all_files <- list.files(path = folder, pattern = "*.msp", full.names = TRUE)
  do.call(c, lapply(all_files, read_lib))
}

