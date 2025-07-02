#' @title De Hosting
#' @import reshape2
#' @import openxlsx
#' @import Rsubread
#' @description Removes the reads that map to the genome reference aplying a specific approach.
#' @param patient_dir path that indicates the directory where fastq files from a patient are saved.
#' @param de_host. It can be "Bowtie", "BWA" or "RSubread".
#' @export
deHosting <- function(patient_dir, de_host, compatible_DRAGEN = FALSE) {
  
  if(de_host  == "Bowtie") {
    de_host_file <- "Bo"
  } else if( de_host == "BWA") {
    de_host_file <- "bwa"
  } else if(de_host == "RSubread") {
    de_host_file <- "Rs"
  }
  
  id <- basename(patient_dir)
  patient_dir_trim <- paste0(patient_dir, "/trimmed", sep="")
  file_list_trimmed <- list.files(patient_dir_trim, full.names = TRUE, recursive = FALSE)
  
  log_file <- path.expand(sprintf("%s/%s_%s_Log.txt", patient_dir_trim, id, de_host_file))
  error_file <- path.expand(sprintf("%s/%s_%s_Error.txt", patient_dir_trim, id, de_host_file))
  
  
  if(length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("DH%s_S04_L001_R1_001.fastq.gz", de_host_file))])) == 0 | length(nchar(file_list_trimmed[endsWith(file_list_trimmed, sprintf("DH%s_S04_L001_R2_001.fastq.gz", de_host_file))])) == 0 ) {
    message(sprintf("Removing host DNA with %s...", de_host))
    
    #Uses trimmed fastq files as input:
    r1_trim <- sprintf("%s/%sT_S04_L001_R1_001.fastq.gz", patient_dir_trim, id)
    r2_trim <- sprintf("%s/%sT_S04_L001_R2_001.fastq.gz", patient_dir_trim, id)
    out_r1_hrmv <- sprintf("%s/%sDH%s_S04_L001_R1_001.fastq", patient_dir_trim, id, de_host_file)
    out_r2_hrmv <- sprintf("%s/%sDH%s_S04_L001_R2_001.fastq", patient_dir_trim, id, de_host_file)
    
    out <- downloadHG38()
    hg38Fasta <- out[[2]]
    #hg38Fasta  = Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
    
    # DEHOSTING WITH BOWTIE ----
    if(de_host == "Bowtie") {
    
      B <- downloadBowtie2()
      Bowtie2 <- B[[1]]
      hg38_indexBowtie2 <- B[[2]]
      
      bowtie_output_bam <- sprintf("%s/%sDH%s.bam", patient_dir_trim, id, de_host_file)
      start_time_bowtie <- Sys.time()
      
      #Step 1: execute bowtie2
      system2(
        Bowtie2,
        args = c(
          "-x", hg38_indexBowtie2,
          "-1", r1_trim,
          "-2", r2_trim,
          "-q", "--phred33", "--sensitive", "--end-to-end", "-p", "10",
          "-S", bowtie_output_bam 
        ),
        stdout = log_file,
        stderr = error_file
      )
      end_time_bowtie <- Sys.time()
      bowtie_duration <- end_time_bowtie - start_time_bowtie
      
   
      # Step 2: Process Bowtie2 output with Samtools y Bedtools
      start_time_samtools <- Sys.time()
      Samtools <- downloadSamtools()
      system(sprintf("%s --version", Samtools))
      system(sprintf(
        "%s view -bh -h %s | %s view -bh -h -f 12 -F 256 - | %s sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
        Samtools, bowtie_output_bam, Samtools, Samtools, out_r1_hrmv, out_r2_hrmv
      ))
      end_time_samtools <- Sys.time()
      samtools_duration <- end_time_samtools - start_time_samtools
      time_bo <- samtools_duration + bowtie_duration
      print(paste("Execution time:", time_bo))
      
      # DEHOSTING WITH BWA ---------
    } else if (de_host == "BWA") {
      
      BWA_out <- downloadBWA()
      BWA <- BWA_out[[1]]
      indexBWA <- BWA_out[[2]]
      
      if(!file.exists(indexBWA) | length(list.files(indexBWA)) == 0) {
        # Generate BWA index
        system2(
          BWA,
          args = c(
            "index",
            "-p", file.path(indexBWA, "Homo_sapiens.GRCh38.dna_sm.primary_assembly"), 
            hg38Fasta))
      }
      
      start_time_bwa <- Sys.time()
      system2(
        BWA,
        args = c(
          "mem",                                   
          "-t", "10",                              
          file.path(indexBWA, "Homo_sapiens.GRCh38.dna_sm.primary_assembly"), 
          r1_trim,                                
          r2_trim,                                
          sprintf(
            "| samtools view -bh -h - | samtools view -bh -h -f 12 -F 256 - | samtools sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
            out_r1_hrmv, out_r2_hrmv
          )  ))
      end_time_bwa <- Sys.time()
      time_bwa <- end_time_bwa - start_time_bwa
      print(paste("BWA Execution Time:", time_bwa))
      
      
      #DEHOSTING WITH RSUBREAD ---------------
    } else if (de_host == "RSubread") {
      
      library(Rsubread)
      
      if(length(list.files("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread")) == 0) {
        # Generate Rsubread index
        dir.create("/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread")
        buildindex(
          basename = "/home/daniela/R/x86_64-pc-linux-gnu-library/4.1/PipelineBiota-Softwares/HG38/indexRsubread/Homo_sapiens.GRCh38.dna_sm.primary_assembly",
          reference = hg38Fasta)
      }
      
      # Alignment with Rsubread
      start_time_rs <- Sys.time()
      out_bam <- paste0(patient_dir_trim, "/", id ,"_R1R2_Rsubread.bam")
      align(
        index = indexRsubread,
        readfile1 = r1_trim,
        readfile2 = r2_trim,
        input_format = "FASTQ",
        output_file = out_bam,
        phredOffset = 33,
        nthreads = 10,
        PE_orientation = "fr", # "forward-reverse"
        type = "dna", 
        useAnnotation = FALSE
      )
      
      # Filter and conversion to FASTQ
      system(sprintf("samtools view -bh -h %s | samtools view -bh -h -f 12 -F 256 - | samtools sort -n - | bedtools bamtofastq -i - -fq %s -fq2 %s",
                     out_bam, out_r1_hrmv, out_r2_hrmv))
      end_time_rs <- Sys.time()
      time_rs <- end_time_rs - start_time_rs
      print(paste("Rsubread Execution Time:", time_rs))
      
      file.remove(out_bam)
      file.remove(paste0(patient_dir_trim, "/", id ,"_R1R2_Rsubread.bam.indel.vcf"))
    }
    
    
    system2("gzip", args = c(out_r1_hrmv, out_r2_hrmv))
    
    out_r1_hrmv <- sprintf("%s/%sDH%s_S04_L001_R1_001.fastq.gz", patient_dir_trim, id, de_host_file)
    out_r2_hrmv <- sprintf("%s/%sDH%s_S04_L001_R2_001.fastq.gz", patient_dir_trim, id, de_host_file)
    
    # Modify headers for compatibility formats to be upleaded to BS Illumina: ------
    if(compatible_DRAGEN == "TRUE") {
      cmd <- paste(
        "zcat", out_r1_hrmv,
        "| sed 's/\\/1$/ 1:N:0:CGAGGCTG+CTCCTTAC/'",
        "| gzip >",
        paste0(out_r1_hrmv, "_fixed")
      )
      system(cmd)
      
      cmd <- paste(
        "zcat", out_r2_hrmv,
        "| sed 's/\\/2$/ 2:N:0:CGAGGCTG+CTCCTTAC/'",
        "| gzip >",
        paste0(out_r2_hrmv, "_fixed")
      )
      system(cmd)
      
      file.remove(out_r1_hrmv)
      file.remove(out_r2_hrmv)
      
      old_r1_hrmv <- sprintf("%s/%sDH%s_S04_L001_R1_001.fastq.gz_fixed", patient_dir_trim, id, de_host_file)
      old_r2_hrmv <- sprintf("%s/%sDH%s_S04_L001_R2_001.fastq.gz_fixed", patient_dir_trim, id, de_host_file)
      
      
      file.rename(old_r1_hrmv, out_r1_hrmv)
      file.rename(old_r2_hrmv, out_r2_hrmv)
      
    }
    
  } else {
    message("De-hosting already done.")
  }
}