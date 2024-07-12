library(limma)
# library(Rbwa)
library(dplyr)
library(tools)
library(optparse)

fn_exists_fasta <- function(folder_fasta) {
  extension = unlist(lapply(list.files(folder_fasta, pattern = "fa"), function(x)
    file_ext(x)))
  extension_fa <- extension[grep("^fa$", extension)]
  extension_fasta <- extension[grep("^fasta$", extension)]
  ## Es fa o fasta, y existe ?
  if (length(extension_fa) == 0 && length(extension_fasta) == 0) {
    stop("No existe archivo de referencia")
    
  } else if (length(extension_fa) != 0 ||
             length(extension_fasta) != 0) {
    if (length(extension_fa) != 0) {
      extension <- extension_fa
    } else if (length(extension_fasta) != 0) {
      extension <- extension_fasta
    }
    
    
  }
  ## el archivo fasta ?
  fasta_file <-
    list.files(folder_fasta,
               pattern = paste0(".", extension, "$"),
               full.names = T)
  return(fasta_file)
  
  
}

index_fasta_samtools <- function(folder_fasta) {
  ## Vemos si existe el archivo con patron fa
  
  fasta_file <- fn_exists_fasta(folder_fasta)
  if (length(list.files(folder_fasta, pattern = "fai$")) == 0) {
    command <- paste("samtools faidx", fasta_file)
    print(command)
    system(command, intern = T)
  } else{
    message("Ya esta el index fai")
  }
  
  
  
  
}

index_bwa <- function(folder_fasta) {
  fasta_file <- list.files(folder_fasta, pattern = ".fasta$", full.names =
                             T)
  
  extension <- c("amb", "ann", "bwt", "pac", "sa")
  
  if (length(file.exists(list.files(folder_fasta, extension[1], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[2], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[3], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[4], full.names = T))) ==
      0 |
      length(file.exists(list.files(folder_fasta, extension[5], full.names = T))) ==
      0) {
    ## no existen bwa index
    print("Creando ficheros índices para bwa mem...")
    
    comando <- paste("bwa index", fasta_file)
    system(command = comando, intern = T)
    
  } else if (length(file.exists(list.files(folder_fasta, extension[1], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[2], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[3], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[4], full.names = T))) !=
             0 &
             length(file.exists(list.files(folder_fasta, extension[5], full.names = T))) !=
             0) {
    message("Ya se han creado los ficheros para el alineamiento")
    
  }
  
  
  
}



bwamem <- function(fastq_dir ,
                   folder_fasta,
                   threads,
                   folder_gatk,
                   output_dir) {
  ### conseguimos el archivo fasta
  
  fasta_file <- list.files(folder_fasta, pattern = ".fasta$", full.names =
                             T)
  
  ####fai index exist ?######
  
  if (!length(file.exists(list.files(dirname(fasta_file), "fai")))) {
    print("### generating fai index...")
    
    index_fasta_samtools(folder_fasta)
    
  }
  
  ## buscamos los archivos fastq
  fastq_files <- list.files(fastq_dir, full.names = T)
  print(fastq_files)
  ## Archivos fastq concatenados
  fastq_full_path_files <-
    c(fastq_1 = fastq_files[1], fastq_2 = fastq_files[2])
  ## Creamos directorio de mapeo
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  
  output_file_name <-
    file_path_sans_ext(output_file_name[length(output_file_name)])
  ### Creando el directorio de mapeo
  
  mapping_output_dir <- file.path(output_dir, "mapping_output", output_file_name)
  if (!dir.exists(mapping_output_dir)) {
    dir.create(mapping_output_dir, recursive = T)
  }
  ### out file name sam file
  output_file_sam <-
    file.path(mapping_output_dir, paste0(output_file_name, ".sam"))
  ### out file name bam file
  
  
  output_file_bam <-
    file.path(mapping_output_dir, paste0(output_file_name, ".bam"))
  ### out file name bam sort file
  
  output_file_sorted_bam <-
    file.path(mapping_output_dir,
              paste0(output_file_name, ".sorted.bam"))
  ### Si no esta mapepeado, mapear
  print(output_file_bam)
  if (!file.exists(output_file_bam)) {
    print("#### MAPPING...#####")
    
    comando <-
      paste(
        "bwa mem -HMpP -v 3 -t",
        threads,
        fasta_file,
        fastq_full_path_files[1],
        fastq_full_path_files[2],
        ">",
        output_file_sam
      )
    
    
    system(command = comando, intern = T)
    ### Si ya se ha mapeado pero no esta el bam, crearlo
    
    print("#### SAM TO BAM")
    
    command_sam_to_bam <-
      paste("samtools view -S -b -h -@",
            threads,
            output_file_sam,
            "-o",
            output_file_bam)
    
    print(command_sam_to_bam)
    system(command_sam_to_bam, intern = T)
    ### Si ya esta el bam, sortearlo
    
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        folder_gatk,
        "SortSam -CREATE_INDEX true -INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
    print(command_out_bam_sorted)
    system(command_out_bam_sorted, intern = T)
    ### AHORA SE PROCEDERIA A MERGE LOS BAMS EN CASO DE TRIO
    
  } else if (file.exists(output_file_bam) &&
             !file.exists(output_file_sam) &&
             !file.exists(output_file_sorted_bam)) {
    print("#### SAM TO BAM")
    
    command_sam_to_bam <-
      paste("samtools view -S -b -@",
            threads,
            output_file_sam,
            "-o",
            output_file_bam)
    
    print(command_sam_to_bam)
    system(command_sam_to_bam, intern = T)
    ### Si ya esta el bam, sortearlo
    
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        folder_gatk,
        "SortSam -CREATE_INDEX true -R",
        fasta_file,
        "-INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
    print(command_out_bam_sorted)
    system(command_out_bam_sorted, intern = T)
  }
  else if (file.exists(output_file_sam) &&
           file.exists(output_file_bam) &&
           !file.exists(output_file_sorted_bam)) {
    print("### BAM to sorted BAM")
    command_out_bam_sorted <-
      paste(
        folder_gatk,
        "SortSam -CREATE_INDEX true -R",
        fasta_file,
        "-INPUT",
        output_file_bam,
        "-OUTPUT",
        output_file_sorted_bam,
        "-SORT_ORDER coordinate -VALIDATION_STRINGENCY STRICT"
      )
    
  } else if (file.exists(output_file_sam) &&
             file.exists(output_file_bam) &&
             file.exists(output_file_sorted_bam)) {
    print("YA SE HA MAPEADO")
    
  } else{
    stop(message("No se ha mapeado bien", call = T))
  }
  
  
}




markdups <- function(output_dir, fastq_dir , gatk_file) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  
  ## Creamos directorio de mapeo, renombramos al archivo
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  ## nombramos al arhivo bam sorteado
  bam_file <-
    file.path(mapping_output_dir,
              output_file_name,
              paste0(output_file_name, ".sorted.bam"))
  ## nombramos al archivo bam marcado con duplciados
  mark_file <-
    file.path(
      mapping_output_dir,
      output_file_name,
      paste0(output_file_name, ".sorted.mark_dup.bam")
    )
  ## nombramos al archivo con las meteicas
  metrics_file <-
    file.path(
      mapping_output_dir,
      output_file_name,
      
      paste0(output_file_name, ".sorted.mark_dup.txt")
    )
  ## Si no existe el archivo crearlo
  if (!file.exists(mark_file)) {
    command <- paste(
      gatk_file,
      " MarkDuplicates -CREATE_INDEX true -INPUT",
      bam_file,
      "-VALIDATION_STRINGENCY STRICT -OUTPUT",
      mark_file,
      "-M",
      metrics_file
    )
    print("AQUIIIIIIIIIIIII******************************")
    print(command)
    system(command = command, intern = T)
    
  } else{
    message("Ya se han marcado los duplicados")
  }
  ## llamamos al comando
  
  
}


create_dict <- function(folder_fasta, gatk_file) {
  ## llamamos al archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  fasta_dict <- paste0(file_path_sans_ext(fasta_file), ".dict")
  if (!file.exists(fasta_dict)) {
    ## creamos el diccionario del fasta
    command <-
      paste(gatk_file, " CreateSequenceDictionary -R", fasta_file)
    system(command, intern = T)
    
  } else{
    message("Ya se ha creado el diccionario fasta")
  }
  
  
}

creacion_readgroup <-
  function(output_dir , fastq_dir , gatk_file) {
    ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
    mapping_output_dir <- file.path(output_dir, "mapping_output")
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    codigo <-
      unlist(strsplit(gsub("R[12]", "sample", fastq_files[1]), "/"))
    
    ## quitamos la extension del archivo
    output_file_name <- file_path_sans_ext(output_file_name)
    ## obtenemos el archivo marcado con duplicados
    mark_file <-
      file.path(
        mapping_output_dir,
        output_file_name,
        paste0(output_file_name, ".sorted.mark_dup.bam")
      )
    ## quitamos la extension y renombramos al archivo de salida
    out_file <- paste0(file_path_sans_ext(mark_file), "_RG.bam")
    
    ## si ya se ha creado el archivo con grupo
    
    if (!file.exists(out_file)) {
      command <-
        paste(
          gatk_file,
          "AddOrReplaceReadGroups I=",
          mark_file,
          "O=",
          out_file,
          paste0(
            "RGID=",
            codigo,
            " RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=1"
          )
        )
      system(command = command, intern = T)
      
    } else{
      message("Ya estan los grupos")
    }
    
    
    
    
  }

base_recalibrator <-
  function(folder_fasta,
           output_dir,
           folder_data_gatk,
           fastq_dir,
           gatk_file) {
    ## llamamos al knwon sites file vcf
    known_sites_file <-
      list.files(folder_data_gatk,
                 pattern = ".vcf$",
                 full.names = T)
    ## llamamos al archivo fasta
    fasta_file <- fn_exists_fasta(folder_fasta)
    ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
    mapping_output_dir <- file.path(output_dir, "mapping_output")
    fastq_files <- list.files(fastq_dir, full.names = F)
    output_file_name <-
      unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
    ## quitamos la extension del archivo
    output_file_name <- file_path_sans_ext(output_file_name)
    ## obtenemos el archivo marcado con duplicados
    RG_file <-
      file.path(
        mapping_output_dir,
        output_file_name,
        
        paste0(output_file_name, ".sorted.mark_dup_RG.bam")
      )
    ## quitamos la extension y renombramos al archivo de salida
    out_file <- file.path(dirname(RG_file), "recal_data.table")
    ## si no existe la tabla de recalibracion, la calculamos
    
    if (!file.exists(out_file)) {
      command <-
        paste(
          gatk_file,
          "BaseRecalibrator -I",
          RG_file,
          " -R",
          fasta_file,
          " --known-sites",
          known_sites_file,
          " -O"  ,
          out_file
        )
      print(command)
      system(command = command, intern = T)
    } else{
      message("ya existe la tabla de base recalculator")
    }
    
    
  }


applybqsr <- function(folder_fasta,
                      output_dir,
                      fastq_dir,
                      gatk_file) {
  ## llamamos al archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  ## quitamos la extension del archivo
  output_file_name <- file_path_sans_ext(output_file_name)
  ## obtenemos el archivo marcado con duplicados
  RG_file <-
    file.path(
      mapping_output_dir,
      output_file_name,
      
      paste0(output_file_name, ".sorted.mark_dup_RG.bam")
    )
  
  recal_data.table <-
    file.path(dirname(RG_file), "recal_data.table")
  
  out_file <-
    file.path(
      dirname(recal_data.table),
      paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam")
    )
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        gatk_file,
        "ApplyBQSR -I",
        RG_file ,
        "-R",
        fasta_file,
        " --bqsr-recal-file",
        recal_data.table,
        " -O",
        out_file
      )
    print(command)
    system(command = command, intern = T)
  } else{
    message("Ya se ha aplicado el bsqr")
  }
  
}

haplotype_caller <- function(output_dir,
                             folder_fasta,
                             fastq_dir,
                             gatk_file) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## recuperamos el ultimo archivo bam
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  
  output_file_name <- file_path_sans_ext(output_file_name)
  out_dir <- file.path(output_dir, "variantCalling", output_file_name)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  
  ## bam file
  bam_file <-
    file.path(
      mapping_output_dir,
      output_file_name,
      
      paste0(output_file_name, ".sorted.mark_dup_RG_bqsr.bam")
    )
  
  out_file <-
    file.path(out_dir, paste0(basename(output_file_name), ".g.vcf.gz"))
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        gatk_file,
        "HaplotypeCaller -I",
        bam_file,
        "-R",
        fasta_file,
        "-ERC GVCF -O",
        out_file
      )
    system(command, intern = T)
  } else{
    message("Ya se han llamado a las variantes")
  }
  
  
}

compute_depth <- function(folder_fasta, output_dir) {
  ## el directorio del mapeo esta creado pero lo tenemos que guardar de nuevo
  ## recuperamos los archivos bam
  mapping_output_dir <- file.path(output_dir, "mapping_output")
  dirs <- list.dirs(mapping_output_dir, recursive = F)
  dirs.l <- length(dirs)
  dirs <- dirs[-dirs.l]
  ## recuperamos el que tenemos
  full_files <- sapply(dirs[-5], function(x)
    list.files(x, full.names = T, pattern = "_bqsr.bam$"))
  
  ## creamos el direcotrio de cobertura
  
  dir_coverage <- file.path(output_dir, "coverage_and_stats")
  if (!dir.exists(dir_coverage)) {
    dir.create(dir_coverage)
  }
  
  samples_names <- strsplit2(dirs, "/")
  samples_names <- samples_names[, ncol(samples_names)]
  ## buscamos los nombres de los archibvos y para cada archivo
  ## se computa la cobertura
  for (s in 1:dirs.l) {
    sample_name <- samples_names[s]
    coverage_sample_dir <- file.path(dir_coverage, sample_name)
    bam_file <- full_filess[s]
    if (!dir.exists(coverage_sample_dir)) {
      dir.create(coverage_sample_dir)
    }
    outfile_coverage <- file.path(coverage_sample_dir, "coverage.txt")
    print(paste("computando cobertura muestra: ", sample_name, " ...."))
    
    comando <- paste("./compute_depth.sh", bam_file, ">", outfile_coverage)
    print(comando)
    system(comando, intern = T)
    
  } else{
    print("la cobertura ya esta computada")
    next
  }
  
  
}
fun_reheader <- function(output_base, output_dir, fastq_dir) {
  fastq_files <- list.files(fastq_dir, full.names = F)
  codigo <-
    file_path_sans_ext(unlist(strsplit(
      gsub("R[12]", "sample", fastq_files[1]), "/"
    )))
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  
  output_file_name <- file_path_sans_ext(output_file_name)
  out_dir <- file.path(output_dir, "variantCalling", output_file_name)
  out_file <-
    file.path(out_dir, paste0(basename(output_file_name), ".g.vcf.gz"))
  
  out_reheader <- file.path(out_dir, paste0(output_file_name, "_reheader.g.vcf.gz"))
  if (!file.exists(out_reheader)) {
    command <- paste("echo", codigo, ">", file.path(out_dir, "./sample.txt"))
    
    system(command = command, intern = T)
    print("Reheadeando ...")
    command <- paste(
      "bcftools reheader -s",
      file.path(out_dir, "sample.txt"),
      "-o",
      out_reheader,
      out_file
    )
    print(command)
    system(command = command, intern = T)
    print("tabixeando...")
    command <- paste("tabix -p vcf", out_reheader)
    print(command)
    system(command = command, intern = T)
    
  } else{
    print("ya esta reheado")
  }
  
  
}







fun_merge <- function(output_dir,
                      fastq_dir,
                      file_cohort,
                      file_gatk,
                      folder_fasta,
                      samples,
                      output_cohort) {
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  fastq_files <- list.files(fastq_dir, full.names = F)
  output_file_name <-
    unlist(strsplit(gsub("R[12]", "map", fastq_files[1]), "/"))
  
  output_file_name <- file_path_sans_ext(output_file_name)
  out_dir <- file.path(output_dir, "variantCalling", output_file_name)
  
  
  if (!is.null(file_cohort)) {
    file_reheader <- list.files(
      out_dir,
      pattern = "reheader.g.vcf.gz$",
      full.names = T,
      recursive = T
    )
    
    writeLines(file_reheader, "./samples_to_merge.list")
    
    
    command <- paste(
      file_gatk,
      "CombineGVCFs -R --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' ",
      fasta_file,
      "--variant",
      file_cohort,
      
      "--variant ./samples_to_merge.list --read-validation-stringency LENIENT -O",
      output_cohort
    )
    
    
    
  } else{
    files_to_merge <- vector(mode = "character", length = length(samples))
    
    for (s in 1:length(samples)) {
      sample <- samples[s]
      files_to_merge[s] <- file.path(
        "output",
        "variantCalling",
        paste0(sample, "_map"),
        paste0(sample, "_map_reheader.g.vcf.gz")
      )
      
    }
    print("aqui")
    print(output_cohort)
    print(files_to_merge)
    print(paste(files_to_merge, sep = "\n"))
    cat(paste(files_to_merge, collapse = "\n"), file = "samples_to_merge.list")
    command <- paste(
      file_gatk,
      "CombineGVCFs --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -R ",
      fasta_file,
      
      "--variant samples_to_merge.list --read-validation-stringency LENIENT  --output",
      output_cohort
    )
    
  }
  print(command)
  if (!file.exists(output_cohort)) {
    system(command = command, intern = TRUE)
    
  } else{
    print("ya se han mezclado los archivos")
  }
  
  
  
}

genotypeGVCF <- function(folder_fasta,
                         output_dir,
                         cohort_file,
                         file_gatk) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  
  output_file_name <- file_path_sans_ext(basename(cohort_file))
  ## indicamos el archivo de entrada
  dir_out <- file.path(output_dir, "variantCallingPostMerge")
  if (!dir.exists(dir_out)) {
    dir.create(dir_out)
  }
  file_out <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz"))
  
  if (!file.exists(file_out)) {
    command <-
      paste(file_gatk,
            " GenotypeGVCFs -R",
            fasta_file,
            "-V",
            cohort_file,
            "-O",
            file_out)
    system(command = command, intern = T)
  } else{
    message("Ya se ha calculado la probabilidad posterior de alelo no referente")
  }
  
  
}

variantRecallibrator <-
  function(folder_fasta,
           folder_data_gatk,
           output_dir,
           file_gatk,
           cohort_file) {
    fasta_file <- fn_exists_fasta(folder_fasta)
    dir_out <- file.path(output_dir, "variantCallingPostMerge")
    output_file_name <- file_path_sans_ext(basename(cohort_file))
    file_out <-
      file.path(dir_out,
                paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz"))
    
    
    
    in_file <-
      file.path(dir_out,
                paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz"))
    
    snps_recal_file <-
      file.path(dir_out,
                paste0(output_file_name, "_apply_genotypeGVCF.recal"))
    
    tranches_file <-
      file.path(dir_out,
                paste0(output_file_name, "_apply_genotypeGVCF.g.tranches"))
    
    if (!file.exists(snps_recal_file) |
        !file.exists(tranches_file)) {
      command <-
        paste(
          file_gatk,
          " VariantRecalibrator -V",
          in_file,
          " --resource:hapmap,known=false,training=true,truth=true,prior=15",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
          ),
          "--resource:omni,known=false,training=true,truth=true,prior=12",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
          ),
          "--resource:1000G,known=false,training=true,truth=false,prior=10",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
          ),
          "--resource:dbsnp,known=true,training=false,truth=false,prior=7",
          file.path(
            folder_data_gatk,
            "resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
          ),
          "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP",
          "-O",
          snps_recal_file,
          "--tranches-file",
          tranches_file
        )
      system(command = command, intern = T)
      
    } else{
      message("Ya se ha hecho el pre recalibrado de variantes")
    }
    
    
  }
applyVQSR <- function(folder_fasta, output_dir, cohort_file) {
  ## recuperamos el archivo fasta
  fasta_file <- fn_exists_fasta(folder_fasta)
  ## obtenemos el nombre del archivo
  dir_out <- file.path(output_dir, "variantCallingPostMerge")
  output_file_name <- file_path_sans_ext(basename(cohort_file))
  in_file <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.g.vcf.gz"))
  
  snps_recal_file <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.recal"))
  
  tranches_file <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.g.tranches"))
  
  out_file <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf"))
  
  if (!file.exists(out_file)) {
    command <-
      paste(
        "~/tools/gatk-4.3.0.0/gatk ApplyVQSR -R",
        fasta_file,
        "-V",
        in_file,
        "-O",
        out_file,
        "--truth-sensitivity-filter-level 99.0 --tranches-file",
        tranches_file,
        "--recal-file",
        snps_recal_file,
        "-mode SNP"
      )
    system(command = command, intern = T)
    
  } else{
    message("Ya se ha aplicado VQSR")
  }
  
}
analysisReady <- function(output_dir, cohort_file) {
  ## obtenemos el nombre del archivo
  dir_out <- file.path(output_dir, "variantCallingPostMerge")
  output_file_name <- file_path_sans_ext(basename(cohort_file))
  
  
  in_file <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.vqsr.vcf"))
  out_file <-
    file.path(dir_out,
              paste0(output_file_name, "_apply_genotypeGVCF.vqsr.pass.vcf"))
  
  if (!file.exists(out_file)) {
    command <-
      paste("bcftools view -f 'PASS,.' -O v -o", out_file, in_file)
    system(command = command, intern = T)
    
  } else{
    "ya se ha hecho el PASS filter"
  }
  
  
}

anotation <-
  function(path_snpeff,
           output_dir,
           clinvar_db,
           dbsnp_db,
           gwas_db,
           dbnsfp_db,
           cohort_file,
           path_gene_sets) {
    ## obtenemos el nombre del archivo
    output_file_name <- file_path_sans_ext(basename(cohort_file))
    
    anotacion_dir <- file.path(output_dir, "anotation")
    if (!dir.exists(anotacion_dir)) {
      dir.create(anotacion_dir)
    }
    
    in_file <-
      file.path(
        output_dir,
        "variantCallingPostMerge",
        paste0(output_file_name, "_apply_genotypeGVCF.vqsr.pass.vcf")
      )
    output_file_anno1 <-
      file.path(anotacion_dir, paste0(output_file_name, "anno_snpeff.vcf"))
    output_file_anno1_1 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2.vcf"))
    output_file_anno2 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2_clinvar.vcf"))
    
    output_file_anno3 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2_clinvar_freqs.vcf"))
    
    output_file_anno4 <-
      file.path(anotacion_dir,
                paste0(output_file_name, "anno_snpeff2_clinvar_freqs_gwas.vcf"))
    output_file_anno5 <-
      file.path(
        anotacion_dir,
        paste0(output_file_name, "anno_snpeff2_clinvar_freqs_gwas_2.vcf")
      )
    output_file_anno6 <- file.path(anotacion_dir, paste0(output_file_name, "_anno_final.vcf"))
    
    print(file.path(path_snpeff, "snpEff.jar"))
    print("arribaaaa")
    
    comando1 <-
      paste(
        "java ",
        paste0("-Xmx", threads, "g"),
        "-jar",
        file.path(path_snpeff, "snpEff.jar"),
        "hg38 -v ",
        in_file,
        " > ",
        output_file_anno1
      )
    comando2 <-
      paste(
        "java ",
        paste0("-Xmx", threads, "g"),
        "-jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " varType -v ",
        output_file_anno1,
        " > ",
        output_file_anno1_1
      )
    comando3 <-
      paste(
        "java ",
        paste0("-Xmx", threads, "g"),
        "-jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " annotate -v ",
        clinvar_db,
        output_file_anno1_1,
        " > ",
        output_file_anno2
      )
    comando4 <-
      paste(
        "java ",
        paste0("-Xmx", threads, "g"),
        "-jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " annotate -v ",
        dbsnp_db,
        output_file_anno2,
        ">",
        output_file_anno3
      )
    comando5 <-
      paste(
        "java ",
        paste0("-Xmx", threads, "g"),
        "-jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " gwasCat -db",
        gwas_db,
        output_file_anno3,
        ">",
        output_file_anno4
      )
    campos <-
      "aaref,aaalt,rs_dbSNP151,HGVSc_snpEff,HGVSp_snpEff,APPRIS,M-CAP_pred,CADD_phred,GTEx_V8_gene,GTEx_V8_gene,GTEx_V8_tissue,Geuvadis_eQTL_target_gene,Reliability_index"
    comando6 <-
      paste(
        "java",
        paste0("-Xmx", threads, "g"),
        "-jar",
        file.path(path_snpeff, "SnpSift.jar"),
        " dbnsfp  -v -db",
        dbnsfp_db,
        "-f",
        campos,
        output_file_anno4,
        ">",
        output_file_anno5
      )
    comando7 <- paste(
      "java",
      paste0("-Xmx", threads, "g"),
      "-jar",
      file.path(path_snpeff, "SnpSift.jar"),
      "geneSets -v",
      path_gene_sets,
      output_file_anno5,
      ">",
      output_file_anno6
    )
    
    if (!file.exists(output_file_anno1)) {
      print(comando1)
      system(command = comando1, intern = T)
      
    }
    if (!file.exists(output_file_anno1_1)) {
      print(comando2)
      system(command = comando2, intern = T)
      
      
    }
    
    
    if (!file.exists(output_file_anno2)) {
      print(comando3)
      system(command = comando3, intern = T)
      
      
    }
    if (!file.exists(output_file_anno3)) {
      print(comando4)
      system(command = comando4, intern = T)
    }
    
    if (!file.exists(output_file_anno4)) {
      print(comando5)
      system(command = comando5, intern = T)
    }
    
    if (!file.exists(output_file_anno5)) {
      print(comando6)
      system(command = comando6, intern = T)
    }
    if (!file.exists(output_file_anno6)) {
      print(comando7)
      system(command = comando7, intern = T)
    }
    
    
    ##
    print("ya esta anotado, falta limpiar el VCF")
    
    
    
  }



library(limma)
# library(Rbwa)
library(dplyr)
library(tools)
library(optparse)



fun_mergeSample <- function(df_clean, codigo) {
  columnas.soi <- df_clean[, grep(codigo, colnames(df_clean))]
  columna_gt <- columnas.soi[, grep("GT", colnames(columnas.soi))]
  columnas_merged <- data.frame(DP_GQ_PL_BLOQUE_DE_FASE = apply(columnas.soi[, -grep("GT", colnames(columnas.soi))], 1, function(x)
    paste0(x, collapse = "_")))
  
  columnas_return <- cbind(columna_gt, columnas_merged)
  values_to_grep <- (columnas.soi[, grep("GT", colnames(columnas.soi))])
  
  aux <- gsub("/", "|", columnas_return[, 1])
  aux.matrix <- strsplit2(aux, "|")
  aux.vector <- ifelse(aux.matrix[, 1] == aux.matrix[, 3], "HOM", "HET")
  aux2 <- ifelse(grepl("0", columnas_return[, 1]), "REF", "ALT")
  aux3.matrix <- cbind(aux.vector, aux2)
  aux3.matrix[, 2] <- ifelse(is.na(aux3.matrix[, 1]), yes = NA, no = aux3.matrix[, 2])
  
  aux3.paste <- apply(aux3.matrix, 1, function(x)
    paste0(x, collapse = "_"))
  columnas_return$CIGOSIDAD <- aux3.paste
  columnas_return$CIGOSIDAD <- ifelse(grepl("\\.|\\.", aux), NA, columnas_return$CIGOSIDAD)
  columnas_return <- apply(columnas_return, 1, function(x)
    paste(x, collapse = "_"))
  columnas_return.split <- strsplit2(columnas_return, "_")[, 1]
  columnas_return <- ifelse(grepl("\\.|\\.", columnas_return), NA, columnas_return)
  columnas_return <- data.frame(GT_DP_GQ_PL_BLOQUE_DE_FASE = columnas_return)
  return(columnas_return)
}

buscar_herencia <- function(df) {
  vector_hpo <- df$hpo  # Vector de la columna "hpo"
  vector_valores <-
    strsplit(vector_hpo, ";")  # Separar los valores por el delimitador ";"
  vector_resultado1 <-
    sapply(vector_valores, function(x)
      paste(unique(x[grep("inheritance", ignore.case = T, x = x)]), sep = ",", collapse = ","))
  X <- bind_cols(herencia = vector_resultado1, df)
  return(X)
}


fun_post_process <- function(hpo_file, sois, output_dir) {
  ### DEBUG
  dir_out <- file.path(output_dir, "postProcess")
  if (!dir.exists(dir_out)) {
    dir.create(dir_out)
  }
  
  if (!file.exists(file.path(output_dir, "anotation", "anno_final_validated.vcf"))) {
    vcf_file <- list.files(file.path(output_dir, "anotation"),
                           full.names = T,
                           pattern = "_anno_final.vcf")
    print("vcf Scan, needs to be cleaned the vcf file..")
    vcf <- VariantAnnotation::readVcf(vcf_file, genome = "hg38")
    print("write the cleaned vcf file ...")
    VariantAnnotation::writeVcf(vcf,
                                file.path(output_dir, "anotation", "anno_final_validated.vcf"))
    print("convert to tsv file for processing")
    command <- paste(
      "vk vcf2tsv wide --print-header --ANN",
      file.path(output_dir, "anotation", "anno_final_validated.vcf"),
      ">",
      file.path(dir_out, "pre_tsv.tsv")
    )
    system(command = command, intern = T)
    print("process tsv file to ready analysis file")
  } else{
    print("ya esta limpio y listo para post procesdado")
  }
  print("reading tsv ...")
  files_to_exists <- list.files(file.path(output_dir, "postProcess"),
                                recursive = T,
                                pattern = "*.csv")
  if (length(files_to_exists) == 0) {
    tsv <- read.delim(file.path(dir_out, "pre_tsv.tsv"), na.strings = ".")
    print("reading done .")
    df_clean <- tsv[, colSums(is.na(tsv)) != nrow(tsv)]
    rm(list = "tsv")
    vcf_file <- file.path(output_dir, "anotation", "anno_final_validated.vcf")
    print(
      "scan vcf looking for LOF, errors, distance to feature, protein position and nt position"
    )
    vcf_lines <- readLines(vcf_file)
    
    # Keep only the header lines (those starting with #)
    vcf_header <- vcf_lines[grep("^#", vcf_lines)]
    vcf_header.columnas <- tail(vcf_header, 1)
    vcf_grep_LOF <- vcf_lines[grep("LOF", vcf_lines)]
    vcf_read <- read.delim(vcf_file,
                           header = T,
                           skip = length(vcf_header) -
                             1)
    formate <- strsplit2(vcf_read$INFO, ";")
    looklof <- vcf_read[grep("LOF", vcf_read$INFO), c("POS", "INFO"), ]
    looklof <- cbind(looklof$POS,
                     looklof$INFO,
                     limma::strsplit2(looklof[, "INFO"], "LOF="))
    looklof <- looklof[, c(1, 4)]
    looklof <-
      cbind(looklof[, 1], limma::strsplit2(looklof[, 2], ";"))[, 1:2]
    looklof[, 2] <- gsub("=", "", looklof[, 2])
    colnames(looklof) <- c("POS", "LOF")
    ## errores
    distance_to_feature <- df_clean$error
    protein_position <- df_clean$distance_to_feature
    df_clean$cDNA_position.cDNA_len <- df_clean$protein_position
    errors <- df_clean$LOF
    df_clean$error <- errors
    df_clean$LOF <- NA
    df_clean$distance_to_feature <- distance_to_feature
    df_clean$protein_position <- protein_position
    
    looklof <- as.data.frame(looklof)
    looklof$POS <- as.numeric(looklof$POS)
    ##joint
    df_clean <- right_join(looklof, df_clean, by = c("POS", "LOF"))
    
    df_clean <- df_clean[, -grep("feature_id", colnames(df_clean))]
    df_clean <- df_clean[!duplicated(df_clean), ]
    df_clean <- df_clean[, -grep("^SAMPLE.", colnames(df_clean))]
    df_clean <- df_clean[, -grep("PGT$", colnames(df_clean))]
    df_clean <- df_clean[, -grep("_PID$", colnames(df_clean))]
    df_clean <- df_clean[which(!is.na(df_clean$FILTER)), ]
    df_clean <- df_clean[which(df_clean$error == "" |
                                 grepl("^INFO", df_clean$error)), ]
    ## codigo de muestra analizar
    
    codigos <- colnames(df_clean)[grep("^DX", colnames(df_clean))]
    codigos.vector <- unique(strsplit2(codigos, "_")[, 1])
    codigos.lista <- vector(mode = "list", length = length(codigos.vector))
    for (cod in 1:length(codigos.vector)) {
      codigo_i <- codigos.vector[cod]
      codigos.lista[[cod]] <- fun_mergeSample(df_clean = df_clean, codigo = codigo_i)
      
    }
    
    names(codigos.lista) <- codigos.vector
    
    codigos.df <- Reduce("cbind", codigos.lista)
    colnames(codigos.df) <- paste(codigos.vector, colnames(codigos.df), sep =
                                    "_")
    df_quasi <- as.data.frame(bind_cols(codigos.df, df_clean[, -grep("^DX", colnames(df_clean))]))
    
    grep_all <- df_quasi[, grep("^DX", colnames(df_quasi))]
    N <- rowSums(!is.na(grep_all))
    
    df_final <- bind_cols(N = N, df_quasi)
    
    column_names_non_na <- apply(codigos.df, 1, function(row) {
      colnames(codigos.df)[!is.na(row)]
    })
    
    samples <- unlist(lapply(lapply(column_names_non_na, function(x)
      strsplit2(x, "_")[, 1]), function(x)
        paste(x, collapse = ";")))
    
    X <- as.data.frame(bind_cols(samples = samples, df_final))
    
    
    hpo <- read.delim(hpo_file, skip = 1, header = F)[, c(2, 4)]
    
    hpo_ <- aggregate(V4 ~ V2, hpo, FUN = paste, collapse = ";")
    colnames(hpo_) <- c("gene_name", "hpo")
    X <- left_join(X, hpo_, by = "gene_name")[, -1]
    X <- buscar_herencia(X)
    X <- bind_cols(samples = samples, X)
    ### frecuencias
    ## computamos frecuencias
    frecuencias <-
      as.data.frame(lapply(X[, c(grep("^AF_", colnames(X)), grep("_AF$", colnames(X)))], as.numeric))
    freqs <- rowMeans(frecuencias, na.rm = T)
    X <- X[, -which(colnames(X) %in% colnames(frecuencias))]
    X <- bind_cols(freq = freqs, X)
    
    
    sois <- unlist(lapply(sois, function(x)
      gsub("-", ".", x)))
    print("for each sample to analyze:")
    print(sois)
    for (s in 1:length(sois)) {
      soi <- sois[s]
      print(soi)
      soi.df <- X[complete.cases(X[, grep(soi, colnames(X))]), ]
      soi.df <- as.data.frame(bind_cols(SAMPLE_code = soi, soi.df))
      canonica.chrom <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
      hap <- soi.df[!soi.df$CHROM %in% canonica.chrom, ]
      hap_any <- hap[which(!is.na(hap$ALLELEID)), ]
      canonical <- soi.df[soi.df$CHROM %in% canonica.chrom, ]
      if (dim(hap_any)[1] != 0) {
        final <- as.data.frame(bind_rows(canonical, hap_any))
      } else{
        final <- as.data.frame(canonical)
      }
      
      final <- final[which(!grepl("HOM_REF", final[, grep(soi, colnames(final))])), ]
      print(
        "write final vcf of the sample of interest \n
		without reference variants on this sample"
      )
      
      
      grep_samples <- final[, grep("^DX", colnames(final))]
      
      aux_unique_fun <- function(X) {
        terminos <- c("HOM_REF")
        idx <- unlist(lapply(terminos, function(x)
          grepl(x, X, perl = T)))
        X[idx] <- "CUBIERTO PERO NO HAY VARIANTES"
        
        return(X)
      }
      
      grep_samples_aux <- as.data.frame(apply(grep_samples, 2, function(columna)
        aux_unique_fun(columna)))
      
      final[, grep("^DX", colnames(final))] <- grep_samples_aux
      
      n_hom_alt <- colSums(apply(grep_samples_aux, 1, function(x)
        grepl("HOM_ALT", x)), na.rm = T)
      
      n_het_ref <- colSums(apply(grep_samples_aux, 1, function(x)
        grepl("HET_REF", x)), na.rm = T)
      
      n_het_alt <- colSums(apply(grep_samples_aux, 1, function(x)
        grepl("HET_ALT", x)), na.rm = T)
      
      final <- as.data.frame(
        bind_cols(
          N_HOM_ALT = n_hom_alt,
          N_HET_REF = n_het_ref,
          N_HET_ALT = n_het_alt,
          final
        )
      )
      
      dir.sample <- file.path(dir_out, soi)
      if (!dir.exists(dir.sample)) {
        dir.create(dir.sample)
      }
      
      write.csv(
        final,
        file = file.path(dir_out, soi, "post_process_canonical_filtered.csv"),
        row.names = F
      )
      
      ## unicas
      
      print("write unique variants of the sample of interest")
      
      final.unicas <- final[final$N == 1, ]
      
      write.csv(
        final.unicas,
        file = file.path(
          dir_out,
          soi,
          "post_process_canonical_filtered_unique.csv"
        ),
        row.names = F
      )
      
      rm(list = "final")
    }
    
  } else{
    print("ya estan procesados los archivos")
  }
  
  ##end for
  
}




wrapper_fun <- function(folder_fasta_,
                        folders_fastq_,
                        folder_gatk_,
                        output_dir_,
                        threads_,
                        folder_data_gatk_,
                        sample_or_folder,
                        cohort_file,
                        cohort_output_,
                        dbnsfp_db_,
                        gwas_db_,
                        dbsnp_db_,
                        clinvar_db_,
                        path_snpeff_,
                        path_gene_sets_,
                        hpo_file_,
                        regions_) {
  index_bwa(folder_fasta_)
  folders_fastq_ <- list.dirs(folders_fastq_)[-1]
  number_of_samples <- length(folders_fastq_)
  ### loop for samples
  for (s in 1:number_of_samples) {
    folder_sample <- folders_fastq_[s]
    
    salida <- bwamem(
      fastq_dir = folder_sample,
      folder_fasta = folder_fasta_,
      threads = threads_,
      folder_gatk = folder_gatk_,
      output_dir = output_dir_
    )
    
    salida <- markdups(output_dir = output_dir_,
                       fastq_dir  = folder_sample,
                       gatk_file = folder_gatk_)
    salida <- create_dict(folder_fasta = folder_fasta_, gatk_file = folder_gatk_)
    
    salida <- creacion_readgroup(output_dir = output_dir_,
                                 fastq_dir = folder_sample,
                                 gatk_file = folder_gatk_)
    
    salida <- base_recalibrator(
      folder_fasta = folder_fasta_,
      output_dir = output_dir_,
      folder_data_gatk = folder_data_gatk_,
      fastq_dir = folder_sample,
      gatk_file = folder_gatk_
    )
    
    salida <- applybqsr(
      folder_fasta = folder_fasta_,
      output_dir = output_dir_,
      fastq_dir = folder_sample,
      gatk_file = folder_gatk_
    )
    
    salida <- haplotype_caller(
      output_dir = output_dir_,
      folder_fasta = folder_fasta_,
      fastq_dir = folder_sample,
      gatk_file = folder_gatk_
    )
    
    salida <- fun_reheader(output_dir = output_dir_, fastq_dir = folder_sample)
    
    
    
  }
  
  
  
  salida <- fun_merge(
    output_dir = output_dir_,
    fastq_dir = folder_sample,
    file_cohort = cohort_file,
    file_gatk = folder_gatk_,
    folder_fasta = folder_fasta_,
    samples = sample_or_folder,
    output_cohort = cohort_output_
  )
  
  salida <- genotypeGVCF(
    folder_fasta = folder_fasta_,
    output_dir = output_dir_,
    cohort_file = cohort_output_,
    file_gatk = folder_gatk_
  )
  
  salida <- variantRecallibrator(
    folder_fasta = folder_fasta_,
    folder_data_gatk = folder_data_gatk_,
    output_dir = output_dir_,
    file_gatk = folder_gatk_,
    cohort_file = cohort_output_
  )
  
  
  salida <- applyVQSR(
    folder_fasta = folder_fasta_,
    output_dir = output_dir_,
    cohort_file = cohort_output_
  )
  
  salida <- analysisReady(output_dir = output_dir_, cohort_file = cohort_output_)
  
  salida <- anotation(
    path_snpeff = path_snpeff_,
    output_dir = output_dir_,
    clinvar_db = clinvar_db_,
    dbsnp_db = dbsnp_db_,
    gwas_db = gwas_db_,
    dbnsfp_db = dbnsfp_db_,
    cohort_file = cohort_output_,
    path_gene_sets = path_gene_sets_
  )
  
  salida <-  fun_post_process(hpo_file = hpo_file_,
                              soi = sample_or_folder,
                              output_dir = output_dir_)
  
  compute_depth(folder_fasta = folder_fasta_,
                regions = regions_,
                output_dir = output_dir_)
  
  
}


# Define options
# Define the option list
option_list <- list(
  make_option(c("-i", "--folder_input"), type = "character", help = "Comma-separated list of folder names for input data"),
  make_option(c("-f", "--folder_fasta"), type = "character", help = "Path to store fasta files"),
  make_option(c("-g", "--folder_gatk"), type = "character", help = "Folder GATK"),
  make_option(c("-t", "--threads"), type = "numeric", help = "Threads for run"),
  make_option(c("-o", "--output"), type = "character", help = "Path to store alignment results"),
  make_option(c("-d", "--data"), type = "character", help = "Path to gatk data"),
  make_option(c("-c", "--cohort"), type = "character", help = "Cohort file"),
  make_option(c("-O", "--output_cohort"), type = "character", help = "Cohort output file"),
  make_option(c("-p", "--path_snpeff"), type = "character", help = "Path to snpEff"),
  make_option(c("-b", "--dbsnp_db"), type = "character", help = "Path to dbsnp"),
  make_option(c("-n", "--clinvar_db"), type = "character", help = "Path to clinvar"),
  make_option(c("-w", "--gwas_db"), type = "character", help = "Path to gwas"),
  make_option(c("-q", "--dbnsfp_db"), type = "character", help = "Path to dbnsfp"),
  make_option(c("-G", "--path_gene_sets"), type = "character", help = "Path to gene sets"),
  make_option(c("-H", "--hpo_file"), type = "character", help = "hpo file genes to phenotype"),
  make_option(c("-R", "--regions"), type = "character", help = "regions covered")
  
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, positional_arguments = TRUE)

# Print parsed options and positional arguments
print(opt$options) # Named list of options
print(opt$args)    # Positional arguments

# Use the parsed arguments
folder_input <- opt$options$folder_input
folder_fasta <- opt$options$folder_fasta
folder_gatk <- opt$options$folder_gatk
threads <- opt$options$threads
output <- opt$options$output
data <- opt$options$data
cohort <- opt$options$cohort
output_cohort <- opt$options$output_cohort
samples <- opt$args # Positional arguments (sample codes)
path_snpeff <- opt$options$path_snpeff
clinvar_db = opt$options$clinvar_db
gwas_db = opt$options$gwas_db
dbnsfp_db = opt$options$dbnsfp_db
dbsnp_db = opt$options$dbsnp_db
path_gene_sets = opt$options$path_gene_sets
hpo_file = opt$options$hpo_file
regions = opt$options$regions
# Print the values for verification
cat("Folder input:", folder_input, "\n")
cat("Folder fasta:", folder_fasta, "\n")
cat("Folder GATK:", folder_gatk, "\n")
cat("Threads:", threads, "\n")
cat("Output path:", output, "\n")
cat("Data path:", data, "\n")
cat("Cohort file:", cohort, "\n")
cat("Output cohort file:", output_cohort, "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")
cat("Path to snpeff:", paste(path_snpeff, collapse = ", "), "\n")
cat("path to clinvar:", paste(clinvar_db, collapse = ", "), "\n")
cat("path to gwas:", paste(gwas_db, collapse = ", "), "\n")
cat("path to dbnsfp:", paste(dbnsfp_db, collapse = ", "), "\n")
cat("path to dbsnp:", paste(dbsnp_db, collapse = ", "), "\n")
cat("path to genesets:", paste(path_gene_sets, collapse = ", "), "\n")
cat("path to hpo_file:", paste(hpo_file, collapse = ", "), "\n")

print(args)
# Check for missing arguments
# Split folder_input (assuming comma-separated)

# Access arguments

wrapper_fun(
  folder_fasta_ = folder_fasta,
  folders_fastq_ = folder_input,
  folder_gatk_ = folder_gatk,
  output_dir_ = output,
  threads_ = threads,
  folder_data_gatk_ = data,
  sample_or_folder = samples,
  cohort_file = cohort,
  cohort_output_ = output_cohort,
  dbnsfp_db_ = dbnsfp_db,
  gwas_db_ = gwas_db,
  dbsnp_db_ = dbsnp_db,
  clinvar_db_ = clinvar_db,
  path_snpeff_ = path_snpeff,
  path_gene_sets_ = path_gene_sets,
  hpo_file_ = hpo_file,
  regions_ = regions
  
)
