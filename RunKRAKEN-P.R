#' @title Run KRAKEN 2
#' @description Executes taxonomy classification with Kraken 2.
#' @param patients_dir path that indicates the directory where fastq files from a patient are saved. It can be the path of a directory where other folders are stored (each one corresponding to one sample) or the folder of one patient.
#' @param de_host can be "Bowtie", "BWA", "RSubread", "" for no dehosting.
#' @param path_kraken2 path to the executable kraken2
#' @export

RunKRAKEN <- function(patients_dir, de_host) {

  libPath <- dirname(system.file(package = "PipelineBiota"))
  pipe_sof <- sprintf("%s/PipelineBiota-Softwares", libPath)
  softwares <- readLines(sprintf("%s/path_to_soft.txt", pipe_sof))
  linea_software <- grep("(?i)KRAKEN2", softwares, ignore.case = TRUE, value = TRUE)

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "sin"
  }

  out_kraken2 <- downloadKRAKEN2()
  path_kraken2 <- out_kraken2[[1]]
  db_kraken2 <- out_kraken2[[2]]
  
  list_dirs <- list.dirs(patients_dir, full.names = TRUE, recursive = FALSE)
  list_files <- list.files(patients_dir, full.names = TRUE, recursive = FALSE)

  #Process of ONE SINGLE SAMPLE  ---------------------------------------
  if( any(grepl("fastq.gz", list_files)) ) {
    message("The process of one single sample will be executed.")
    patient <- basename(patients_dir)
    print(patient)
    patient_dir <- paste(patients_dir, "/trimmed", sep="")
    list_files <- list.files(patient_dir, full.names = TRUE)

    fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq.gz", patient, de_host_file), list_files)]
    fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq.gz", patient, de_host_file), list_files)]

    if(de_host == "") {
      fileR1 <- list_files[grepl("T_S04_L001_R1_001.fastq.gz", list_files)]
      fileR2 <- list_files[grepl("T_S04_L001_R2_001.fastq.gz", list_files)]
    }

    dir.create(sprintf("%s/Resultados_KRAKEN", patient_dir))

    if(!(file.exists(sprintf("%s/Resultados_KRAKEN/report_%s.sequences", patient_dir, de_host_file)))) {
      
      args <- c(
        sprintf("--db %s", db_kraken2),
        "--threads 15",
        "--use-names",
        sprintf("--output %s/Resultados_KRAKEN/output_%s.txt", patient_dir, de_host_file),
        sprintf("--report %s/Resultados_KRAKEN/report_%s.sequences", patient_dir, de_host_file),
        sprintf("--paired %s %s", fileR1, fileR2)
      )

      start_time_kraken <- Sys.time()
      system2(path_kraken2, args = args)
      end_time_kraken <- Sys.time()
      time_kraken <- end_time_kraken - start_time_kraken
      print(paste("KRAKEN Execution Time:", time_kraken))


    } else {
      return(message(sprintf("KRAKEN report has already been generated for this patient: %s", patient)))
    }

    #Process of MULTIPLE SAMPLES  -------------------------------
  } else {
    message("The process of multiple samples will be executed.")
    for (patient in list_dirs) {
      id <- basename(patient)
      patient <- paste(patient, "/trimmed", sep= "")
      print(patient)
      list_files <- list.files(patient, full.names = TRUE)

      fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq", id, de_host_file), list_files)]
      fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq", id, de_host_file), list_files)]

      if(de_host == "") {
        fileR1 <- list_files[grepl("T_S04_L001_R1_001.fastq", list_files)]
        fileR2 <- list_files[grepl("T_S04_L001_R2_001.fastq", list_files)]
      }

      gz <- ifelse(grepl(".gz", fileR1), "--gzip-compressed", "")
      dir.create(sprintf("%s/Resultados_KRAKEN", patient))

      if(!(file.exists(sprintf("%s/Resultados_KRAKEN/report_%s.sequences", patient, de_host_file)))) {
        args <- c(
          sprintf("--db %s", db_kraken2),
          "--threads 15",
          "--use-names",
          sprintf("--output %s/Resultados_KRAKEN/output_%s.txt", patient, de_host_file),
          sprintf("--report %s/Resultados_KRAKEN/report_%s.sequences", patient, de_host_file),
          sprintf("%s --paired %s %s", gz, fileR1, fileR2)
        )

        system2(path_kraken2, args = args)
      } else {
        print("This patient's report has already been generated.")
      }

    }
  }

}


