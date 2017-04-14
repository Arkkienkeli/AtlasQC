# CueSea

### INTRODUCTION
This tool is designed for quality control of Illumina  microarrays genotyping data.


### PREREQUISITES
To run this tool on your machine, you need to have Python 2.7 or 3 installed.
Besides that, you need several packages, which can be installed using pip.
-scikit-learn
-plotly


### HOW TO RUN
1. Most common scenario:
python3 cuesea.py Genome_Studio_Report.txt --snps=snps.txt
This command will run QC using your Genome_Studio_Report.txt as data file and
snps.txt as list of snps to analyze. X and Y chromosome snps treated as autosomal.
2. Gender
python3 cuesea.py Genome_Studio_Report.txt --snps=snps.txt --gender gender.txt
This command will run QC using your Genome_Studio_Report.txt as data file and
snps.txt as list of snps to analyze. This will include proper analyze of X and Y
chromosome snps as described in HOW IT WORKS paragraph.
3. Angle
python3 cuesea.py Genome_Studio_Report.txt --angle 8.0
The --angle key specifies max angle between points of each cluster. The default
value is 8.0 degrees, but you can slightly increase or dicrease it.

### INPUT FILE FORMATS
1. Genome Studio Report
This tool uses default Genome Studio report, and we use columns:
- SNP Name (which is Illumina Name)
- Sample ID
- Allele 1 - Top
- Allele 2 - Top
- Chr
- X
- Y
2. Snps list
List of snps, one per line.
3. Gender
File with two tab-separated columns (Sample ID, gender) no header.
F - for women, M - for men. Other symbols treated as undefined gender.
For example:
0000-1111	F

### OUTPUT FILES
Tool outputs 3 types of reports:
1. Sample Report - list of snps, which have not pass the QC.
This file contains columns:
SNP - Sample ID - error code.
Error code 0 - snp fails  intensity check 
Error code 1 - snp failed  clusterization step
Error code 2 - snp failed gender check
Error code 3 - snp’s genotype distribution doesn’t fit to HWE
2. Overall report - list of snps with statistics.
This file contains columns:
SNP - Number of samples - Excluded by QC - Excluded by R(intensity) - Excluded by clusterization
Or if it was gender file on input:
SNP - № of samples - Excluded by QC - Excluded by R(intensity) - Excluded by clusterization - Excluded from X\Y
3. Genotype Report - information about SNPs with changed genotype class
This file contains columns:
SNP - Sample ID - old GT class - new GT class
4. Plots
Samples marked green - passes QC
Samples marked blue - not passed QC
Samples marked red - changed GT class

### AUTHORS
Daria Iakovishina, PhD

Nikita Moshkov

### LICENSE
GNU GPL v.3
