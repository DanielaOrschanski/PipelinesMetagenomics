GitHub repository containing the code and processed data necessary to reproduce the analyses and figures presented in manuscript: 

**“Dermatological Implications of De-hosting and Bioinformatics Pipelines on Shotgun Microbiome Analysis”** 

- Daniela Orschanski¹, Leonardo Néstor Rubén Dandeu², Martín Nicolás Rivero³, Vivian Labovsky², Elmer Andrés Fernández¹,⁴*.


¹ Fundación para el Progreso de la Medicina, Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET), Córdoba, Argentina  
² Laboratorio de Inmunohematología, Instituto de Biología y Medicina Experimental (IBYME), Fundación IBYME, Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET), Ciudad Autónoma de Buenos Aires, Buenos Aires, Argentina
³ BiotaLife SKIN, Ciudad Autónoma de Buenos Aires, Buenos Aires, Argentina
⁴ Universidad Argentina de la Empresa (UADE), Buenos Aires, Argentina


In this repository you will find:

- Execution scripts:
  1. DeHosting.R: Contains all the coding necessary for applying any studied de-hosting method (Bowtie2, Rsubread and BWA)
  2. RunDRAGEN-P.R
  3. RunKRAKEN-P.R
  4. RunHuman-P.R
 
- Scripts for reproducibility of figures:
1. Figure2.R: makes use of Fig2a-MappedReads.xlsx, Fig2c-GC.xlsx and Fig2d-HomoSapiens.xlsx.
2. Figure3.R: makes use of Fig3a-MicroorgPerLevel.xlsx, Fig3-All_Taxa_Kraken.xlsx and Fig3-data_list.RData.
3. Figure4.R: makes use of Fig4-RA6genera.xlsx.
4. Figure5.R: makes use of Fig5-DFProteobacteria.xlsx, Fig5-Genera_DifSign_Age.xlsx, Fig5-Genera_DifSign_Sex.xlsx 
5. Figure6.R: makes use of Fig6a-NumberPathways.xlsx, Fig6b-DataList-Pathways.RData, Fig6c-Pathways-Sex.xlsx, Fig6d-Pathways-Age.xlsx, TableS2, TableS3, TableS4.
