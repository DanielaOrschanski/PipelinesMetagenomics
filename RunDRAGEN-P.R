# This includes: 
#  - downloadBS
#  - Upload to BS
#  - Executes DRAGEN Metagenomics Pipeline
#  - Download DRAGEN's Results

#' @title Download BS (Base Space)
#' @description downloads and installs the executable of bs for DRAGEN metagenomics processes by CLI Base Space aplication
#' @param dir_bs path where to store the executable Base Space
#' @export
#' @return path of the bs executable
downloadBS <- function(dir_bs) {
  
  if(!file.exists(sprintf("%s/bs", dir_bs))) {
    setwd(dir_bs)
    system("wget 'https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs'")
    system(sprintf("chmod u+x %s/bs", dir_bs))
    system(sprintf("export PATH=$PATH:%s/bs", dir_bs))
    system(sprintf("%s/bs auth", dir_bs))
  }
  
  bs_path  <- sprintf("%s/bs", dir_bs)
  return(bs_path)
}



#' @title upload to BS (Base Space)
#' @description uploads fastq files to BS
#' @param patient_dir path that indicates the directory where fastq files from a patient are saved.
#' @param bs_path path to base space executable
#' @param de_host can be "Bowtie", "BWA", "RSubread" or "" for none previous dehosting.
#' @export
uploadtoBS <- function(patient_dir, bs_path, de_host = "") {

  id <- basename(patient_dir)

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
    project_id <- "438230795"
    carpeta_bs <- "dehostBowtie"
    nombre_p <- sprintf("%sDHBo", id)

  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
    project_id <- "436873438"
    carpeta_bs <- "dehostBWA"
    nombre_p <- sprintf("%sDHbwa", id)
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
    project_id = "436780349"
    carpeta_bs ="dehostRsubread"
    nombre_p <- sprintf("%sDHRs", id)
  } else if(de_host == "") {
    de_host_file <- "sin"
    project_id <- "436799364"
    carpeta_bs <- "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
  }

  patient <- basename(patient_dir)
  #All the results are stored in "trimmed" folder within each patients folder.
  patient_dir <- paste(patient_dir, "/trimmed", sep="")

  nuevo_path_report <- sprintf("%s/Resultados_DRAGEN/%s%s.DRAGEN-report.tsv", patient_dir, patient, de_host_file)

  if(file.exists(nuevo_path_report)) {
    return("This patient has already been processed with DRAGEN")
  }

  #Select the correspondant R1 y R2:
  print(patient)
  list_files <- list.files(patient_dir, full.names = TRUE)

  fileR1 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R1_001.fastq.gz", patient, de_host_file), list_files)]
  fileR2 <- list_files[grepl(sprintf("%sDH%s_S04_L001_R2_001.fastq.gz", patient, de_host_file), list_files)]

  if(de_host == "") {
    fileR1 <- list_files[grepl("T_S04_L001_R1_001.fastq", list_files)]
    fileR2 <- list_files[grepl("T_S04_L001_R2_001.fastq", list_files)]
  }


  #Upload to base space ----------------------------------------
  command <- sprintf("%s list biosample", bs_path)
  output <- capture.output(system(command, intern = TRUE))
  output <- output[!grepl("^\\+|^\\s*$", output)]
  output <- gsub("^\\|\\s*|\\s*\\|$", "", output)
  data_list <- strsplit(output, "\\s*\\|\\s*")
  df_output <- do.call(rbind, data_list)

  colnames(df_output) <- trimws(df_output[2,])
  df_output <- df_output[-c(1:3),]

  biosamples <- df_output[,2]
  command <- sprintf("%s upload dataset --name %s -p %s %s %s", bs_path, carpeta_bs, project_id, fileR1, fileR2)
  system(command)

}


#' @title Run DRAGEN
#' @description Executes DRAGEN Metagenomics Pipeline using CLI Base Space
#' @param patient_dir path that indicates the directory where fastq files from a patient are saved.
#' @param reference "hg38" as default to indicate the genome reference.
#' @param bs_path path to executable Base Space
#' @param de_host can be "Bowtie", "BWA", "RSubread", "" for no previous dehosting but enables DRAGEN's dehosting, or "sin_PDH" for neither previous de-hosting, nor DRAGEN's de-hosting.
#' @export
#' @import tidyr
#' @import readr
#' @import stringr

RunDRAGEN <- function(patient_dir, bs_path, reference = "hg38", de_host) {

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
    project_id = "438230795"
    carpeta_bs = "dehostBowtie"
    nombre_p <- sprintf("%sDHBo", id)
    dsin_dehost = TRUE

  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
    project_id = "436873438"
    carpeta_bs = "dehostBWA"
    nombre_p <- sprintf("%sDHbwa", id)
    dsin_dehost = TRUE

  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
    project_id = "436780349"
    carpeta_bs ="dehostRsubread"
    nombre_p <- sprintf("%sDHRs", id)
    dsin_dehost = TRUE

  } else if(de_host == "") {
    de_host_file <- "sin"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = FALSE #executes dehosting with dragen for trimmed files

  } else if(de_host == "sinDH_PD") { # neither previous de-hosting, nor DRAGEN de-hosting
    de_host_file <- "sinDH_PD"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = TRUE
  }

  patient <- basename(patient_dir)
  patient_dir <- paste(patient_dir, "/trimmed", sep="")
  id <- patient

  if (!(file.exists(sprintf("%s/Resultados_DRAGEN/%s%s_DRAGEN-report.tsv", patient_dir, patient, de_host_file)))) {

 
    #Avoid reprocessing same sample ----
    command <- sprintf("%s list appsession --project-id %s", bs_path, project_id)

    output <- capture.output(system(command, intern = TRUE))
    df_output <- as.data.frame(do.call(rbind, strsplit(output, " +", perl = TRUE)))
    df_output <- df_output[which(df_output$V3 == sprintf("DRAGENHG38_%s", de_host_file) ),]
    df_output <- df_output[, c(3,4,5,6)]
  
    if( (patient %in% df_output[,3]) | (patient %in% df_output[,4])| (patient %in% df_output[,2]) ) {
      return(message(sprintf("The DRAGEN aplication has already been launched for this patient: %s", patient)))
    }

    #---------------------

    #Identification of the id provided by Base Space to execute the analysis --------
    
    command <- sprintf("%s list biosamples --project-id %s", bs_path, project_id)
    output <- capture.output(system(command, intern = TRUE))
    df_output <- as.data.frame(do.call(rbind, strsplit(gsub("[|]", "", output), " +")))
    df_output <- df_output[,-c(1,2)]
    colnames(df_output) <- df_output[2,]
    df_output <- df_output[-c(1,2,3, nrow(df_output)),]
    df_output <- df_output[, 1:3]
    colnames(df_output)[1] <- "BioSampleName1"

    id_bs <- df_output$Id[which(df_output$BioSampleName == nombre_p)]
    if(length(nchar(id_bs)) == 0 ){
      id_bs <- df_output$BioSampleName[which(df_output$BioSampleName1 == nombre_p)]
    }

    # The table can appear with different formats:
    if(length(nchar(id_bs))==0) {
      id_bs <- df_output$`|`[which(df_output$`"|` == patient)]
    }

    if( length(nchar(id_bs)) == 0) {
      output <- system(command, intern = TRUE)
      output_text <- paste(output, collapse = "\n")
      output_lines <- unlist(strsplit(output_text, "\n"))
      data_lines <- output_lines[-which(grepl("-", output_lines))]
      data_lines <- data.frame(data_lines)
      data_lines <- as.data.frame(lapply(data_lines, function(x) gsub("^\\||\\|$", "", x)))
      data_lines <- data_lines %>%
        separate(data_lines, into = c("BioSampleName", "Id", "ContainerName", "ContainerPosition", "Status"), sep = "\\|")
      id_bs <- data_lines$Id[which(grepl(patient, data_lines$BioSampleName))]
      id_bs <- gsub(" ", "", id_bs)
    }

    #If there are no coincidences:
    if(length(nchar(id_bs)) == 0) {
      print(system(command, intern = TRUE))
      print(patient)
      id_bs <- readline(prompt = "Enter the id of your sample: ")
    }

    #---------------------------------------------------------------
    
    ref <- ifelse(reference == "hg19", "hg19-altaware-cnv-anchor.v8", "hg38-altaware-cnv-anchor.v8")

    if(dsin_dehost == TRUE) { #DRAGEN WITHOUT DEHOSTING
      command <- sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 -o app-session-name:"DRAGENHG38_%s %s" -l "DRAGENHG38_%s %s" -o project-id:%s -o sample-id:%s -o ht-ref:%s -o dehost-checkbox:false -o basespace-labs-disclaimer:Accepted -o db:minikraken20200312',
                         bs_path,
                         de_host_file,
                         patient,
                         de_host_file,
                         patient,
                         project_id,
                         id_bs,
                         ref)


      system(command)
    } else { # DRAGEN WITH DEHOSTING ENABLED
      command <- sprintf('%s launch application -n "DRAGEN Metagenomics Pipeline" --app-version 3.5.12 -o app-session-name:"DRAGENHG38_%s %s" -l "DRAGENHG38_%s %s" -o project-id:%s -o sample-id:%s -o ht-ref:%s -o basespace-labs-disclaimer:Accepted -o db:minikraken20200312',
                         bs_path,
                         de_host_file,
                         patient,
                         de_host_file,
                         patient,
                         project_id,
                         id_bs,
                         ref)
      system(command)
    }

  } else {
    return(message(sprintf("The DRAGEN report has already been generated for this patient: %s", patient)))
  }

}

#' @title Download DRAGEN Report
#' @description downloads DRAGEN Metagenomics Pipeline outputs
#' @param patient_dir path that indicates the directory where fastq files from a patient are saved.
#' @param bs_path path to executable Base Space
#' @param de_host can be "Bowtie", "BWA", "RSubread", "" for no previous dehosting but enables DRAGEN's dehosting, or "sin_PDH" for neither previous de-hosting, nor DRAGEN's de-hosting.
#' @examples path_report <- download_DRAGENReport(patient_dir ="~/Biota/Nuevas13/Muestras/33" )
#' @export
#' @import tidyr
download_DRAGENReport <- function(patient_dir, bs_path, de_host) {

  patient <- basename(patient_dir)
  print(patient)
  id <- patient
  patient_dir <- paste(patient_dir, "/trimmed", sep="")

  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
    project_id = "438230795"
    carpeta_bs = "dehostBowtie"
    nombre_p <- sprintf("%sDHBo", id)
    dsin_dehost = TRUE
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
    project_id = "436873438"
    carpeta_bs = "dehostBWA"
    nombre_p <- sprintf("%sDHbwa", id)
    dsin_dehost = TRUE
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
    project_id = "436780349"
    carpeta_bs ="dehostRsubread"
    nombre_p <- sprintf("%sDHRs", id)
    dsin_dehost = TRUE
  } else if(de_host == "") {
    de_host_file <- "sin"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = FALSE
  } else if(de_host == "sinDH_PD") { 
    de_host_file <- "sinDH_PD"
    project_id = "436799364"
    carpeta_bs = "MuestrasTrimmed"
    nombre_p <- sprintf("%sT", id)
    dsin_dehost = TRUE
  }

  nuevo_path_report <- sprintf("%s/Resultados_DRAGEN/%s%s.DRAGEN-report.tsv", patient_dir, patient, de_host_file)

  if(file.exists(nuevo_path_report)) {
    return(message(sprintf("The DRAGEN report was already downloaded for patient %s", patient)))
  }
  
  #Find appsesion id:
  
  command <- sprintf("%s list appsession --project-id %s", bs_path, project_id)

  output <- capture.output(system(command, intern = TRUE))
  df_output <- as.data.frame(do.call(rbind, strsplit(output, " +", perl = TRUE)))

  df_output <- df_output[which(df_output$V3 == sprintf("DRAGENHG38_%s", de_host_file) | df_output$V4 == sprintf("DRAGENHG38_%s", de_host_file)),]
  df_output <- df_output[, c(3,4,5,6,7)]

  appsession_id <- unique(df_output[which(df_output[,2] == patient), 4])
  if(length(nchar(appsession_id)) == 0) {
    appsession_id <- unique(df_output[which(df_output[, 3] == patient), 5])
  }

  if(length(nchar(appsession_id)) == 0) {
    df_output <- df_output[, c(1,2,4)]
    colnames(df_output) <- c(sprintf("DRAGENHG38_%s", de_host_file), "SAMPLE", "ID")
    appsession_id <- df_output$ID[which(df_output$SAMPLE == patient)]
  }

  if(length(nchar(appsession_id)) == 0) {
    output <- system(command, intern = TRUE)
    output_text <- paste(output, collapse = "\n")
    output_lines <- unlist(strsplit(output_text, "\n"))
    data_lines <- output_lines[-which(grepl("-", output_lines))]
    data_lines <- data.frame(data_lines)
    data_lines <- as.data.frame(lapply(data_lines, function(x) gsub("^\\||\\|$", "", x)))
    data_lines <- data_lines %>%
      separate(data_lines, into = c("BioSampleName", "Id", "ContainerName", "ContainerPosition", "Status"), sep = "\\|")
    appsession_id <- data_lines$Id[which(grepl(patient, data_lines$BioSampleName))]
    appsession_id <- gsub(" ", "", appsession_id)

  }

  if(length(nchar(appsession_id)) == 0) {
    print(system(command, intern = TRUE))
    print(patient)
    appsession_id <- readline(prompt = "Enter the id of your appsession: ")
  }

  if(length(appsession_id)>1) {
    appsession_id <- appsession_id[1]
  }

  #------------------------------------
  
  dir.create(sprintf("%s/DRAGEN_Reports_%s", patient_dir, de_host_file))
  command <- sprintf("%s download appsession -i %s -o %s",
                     bs_path,
                     appsession_id,
                     sprintf("%s/DRAGEN_Reports_%s", patient_dir, de_host_file))
  system(command)

  folder_report <- list.dirs( sprintf("%s/DRAGEN_Reports_%s", patient_dir, de_host_file), full.names = TRUE, recursive = FALSE)
  file_report <- sprintf("%s/%s.microbe-classification-report.tsv", folder_report, nombre_p)
  nuevo_path_report <- sprintf("%s/Resultados_DRAGEN/%s%s.DRAGEN-report.tsv", patient_dir, patient, de_host_file)
  file.rename(from = file_report, to = nuevo_path_report)

  message(sprintf("The DRAGEN report was downloaded for patient %s", patient))

  return(nuevo_path_report)
}

