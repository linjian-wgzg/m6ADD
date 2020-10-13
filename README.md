"# m6ADD" 
"# m6ADD-code" 
m6ADD is a database containing manually collected experimentally confirmed m6A-disease data and data obtained from high-throughput disease m6A modification profiles, aimed at exploring the association between m6A modified gene disorders and diseases. The m6ADD database contains 222 experimentally confirmed m6A-disease relationship pairs (human 185, mouse 37). We screened out differential m6A data from 30 sets of sequencing data of 16 diseases by two calculation methods, and provided a statistical evaluation result. The m6A-disease data includes m6A genomic location, disease name, m6A protein, regulatory mode, tissue / cell line, experimental method, and data source. In addition, we have developed a ppi network tool for obtaining differential m6A genes to show the function of these genes and a tool to predict m6A regulatory proteins associated with 24 types of cancers.
First:
High-throughput sequencing data of m6ADD
The test data and the reference genome annotation file have been uploaded to Baidu cloud disk（Link：https://pan.baidu.com/s/1h6aFgn8pBgGOnwp57kutzg Extraction code：8hat）. 
After downloading the data, 
1. Use fastq software to convert the .sra file to a .fastq file; 
2. Use HISAT2 to compare the sequencing data to the reference genome; 
3. Use samtools to convert .sam files into .bam files and complete the sorting;
4. Use RADAR and MeTDiff to identify differential methylation sites respectively; （MeTDiff_diff.R、RADAR_diff.R）
5. Calculate the score value of the intersection to generate the final result file.（Form_result.R）
Second:
Expanded data of m6ADD
1. We obtain general gene expression data from TCGA, and identify differentially expressed genes and m6A modified regulatory protein genes;
2. Obtain protein interaction data from STRING and HPRD databases;
3. The first step is to obtain differentially expressed genes, use Pearson correlation coefficient to construct a correlation matrix, and combine the obtained protein interaction data to construct a one-step neighbor network of differential genes. Integrate the above results to reconstruct a neighbor network, and finally perform network module mining to predict the function of differentially expressed genes.
