install.packages("usethis")
library("usethis")
install.packages("rlang")
library("rlang")
install.packages("devtools")
library("devtools")
install_github("compgenomics/MeTDiff")
library(MeTDiff)
gtf <-"/path of bam file/Annotation.gtf"

C1_m6A <-"/path of bam file/ctrl1.m6A.bam"
C2_m6A <- "/path of bam file/ctrl2.m6A.bam"
C3_m6A <- "/path of bam file/ctrl3.m6A.bam"
C4_m6A <-"/path of bam file/ctrl4.m6A.bam"
C5_m6A <- "/path of bam file/ctrl5.m6A.bam"
C6_m6A <- "/path of bam file/ctrl6.m6A.bam"
C7_m6A <- "/path of bam file/ctrl7.m6A.bam"
C1_input <-"/path of bam file/ctrl1.input.bam"
C2_input <- "/path of bam file/ctrl2.input.bam"
C3_input <- "/path of bam file/ctrl3.input.bam"
C4_input <-"/path of bam file/ctrl4.input.bam"
C5_input <- "/path of bam file/ctrl5.input.bam"
C6_input <- "/path of bam file/ctrl6.input.bam"
C7_input <- "/path of bam file/ctrl7.input.bam"
T1_m6A <-"/path of bam file/case1.m6A.bam"
T2_m6A <- "/path of bam file/case2.m6A.bam"
T3_m6A <- "/path of bam file/case3.m6A.bam"
T4_m6A <-"/path of bam file/case4.m6A.bam"
T5_m6A <- "/path of bam file/case5.m6A.bam"
T6_m6A <- "/path of bam file/case6.m6A.bam"

T1_input <-"/path of bam file/case1.input.bam"
T2_input <- "/path of bam file/case2.input.bam"
T3_input <- "/path of bam file/case3.input.bam"
T4_input <-"/path of bam file/case4.input.bam"
T5_input <- "/path of bam file/case5.input.bam"
T6_input <- "/path of bam file/case6.input.bam"

treated_ip <- "/path of bam file/treated_IP1.bam"
treated_input <- "/path of bam file/treated_Input1.bam"

IP_BAM <- c(C1_m6A,C2_m6A,C3_m6A,C4_m6A,C5_m6A,C6_m6A,C7_m6A)
INPUT_BAM <- c(C1_input,C2_input,C3_input,C4_input,C5_input,C6_input,C7_input)
TREATED_IP_BAM <- c(T1_m6A,T2_m6A,T3_m6A,T4_m6A,T5_m6A,T6_m6A)
TREATED_INPUT_BAM <- c(T1_input,T2_input,T3_input,T4_input,T5_input,T6_input)

metdiff(GENE_ANNO_GTF=gtf,IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM,
        TREATED_IP_BAM = TREATED_IP_BAM,TREATED_INPUT_BAM=TREATED_INPUT_BAM,
        EXPERIMENT_NAME="MeTDiff_diff_site")