#' @title Run HUMAnN
#' @description Executes concat of R1 and R2 fastq files, executes HUMAnN 3.0 and generates CPM.
#' @param patient_dir path to the folder of the patient that will be analyzed.
#' @param de_host indicates which alignment-based de-hosting software will be used. It can be "Bowtie", "BWA", "RSubread", "" for no dehosting.
#' @export

RunHuman <- function(patient_dir, de_host) {
  
  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  } else if(de_host == "") {
    de_host_file <- "T"
  } else {
    stop("de_host must be Bowtie, BWA, RSubread or empty string")
  }

  id <- basename(patient_dir)
  print(id)
  
  #1. Concatenate R1 and R2:
  concat_file <- sprintf("%s/Vias/%s/%s%s_concatR1R2.fastq.gz", patient_dir, de_host_file, id, de_host_file)

  if(!file.exists(concat_file) | file.info(concat_file)$size == 0) {
    message("Concatenando ...")
    R1 <- sprintf("%s/trimmed/%sDH%s_S04_L001_R1_001.fastq.gz", patient_dir, id, de_host_file)
    R2 <- sprintf("%s/trimmed/%sDH%s_S04_L001_R2_001.fastq.gz", patient_dir, id, de_host_file)

    if(de_host == "") {
      R1 <- sprintf("%s/trimmed/%s%s_S04_L001_R1_001.fastq.gz", patient_dir, id, de_host_file)
      R2 <- sprintf("%s/trimmed/%s%s_S04_L001_R2_001.fastq.gz", patient_dir, id, de_host_file)
    }

    dir.create(sprintf("%s/Vias", patient_dir))
    dir.create(sprintf("%s/Vias/%s", patient_dir, de_host_file))

    start_time_cat <- Sys.time()
    system(sprintf("cat %s %s > %s", R1, R2, concat_file))
    end_time_cat <- Sys.time()
    cat_duration <- end_time_cat - start_time_cat

  }

  #2. Execute HUMAnN:
  vias_dir <-  sprintf("%s/Vias/%s", patient_dir, de_host_file)

  path_human <- "/home/daniela/miniconda3/envs/biobakery3/bin/humann"

  command <- paste(
    "bash -c 'source activate biobakery3 &&",  
    path_human,                                
    sprintf("--input %s", concat_file),
    sprintf("--output %s", vias_dir),
    "--threads 15'",
    sep = " "
  )

  tsv_file <- sprintf("%s/%s%s_concatR1R2_pathabundance.tsv", vias_dir, id, de_host_file)
  if(!file.exists(tsv_file)) {
    message("Executing HUMAnN ...")
    start_time_human <- Sys.time()
    system(command = command, intern = FALSE)
    end_time_human <- Sys.time()
    human_duration <- end_time_human - start_time_human

  } else {
    message("HUMAnN was already executed before!")
  }

  #Add CPM: ---------------------------
  rpk_file <- sprintf("%s/Vias/%s/%s%s_concatR1R2_pathabundance.tsv", patient_dir, de_host_file, id, de_host_file)
  dir.create(sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s", de_host_file))
  cpm_file <- sprintf("/home/daniela/Daniela/Biota/Muestras/SubsetPathways/%s/%s%s_CPM_pathabundance.tsv", de_host_file, id, de_host_file)


  command <- paste(
    "bash -c 'source activate biobakery3 &&",  
    "humann_renorm_table",                                
    sprintf("--input %s", rpk_file),
    sprintf("--output %s", cpm_file),
    "--units cpm",
    "--update-snames'",
    sep = " "
  )
  if(!file.exists(cpm_file)) {
    start_time_cpm <- Sys.time()
    system(command = command, intern = FALSE)
    end_time_cpm <- Sys.time()
    cpm_duration <- end_time_cpm - start_time_cpm
  } else {
    message("CPM file has already been generated for this patient!")
  }

  return(cpm_file)

}
