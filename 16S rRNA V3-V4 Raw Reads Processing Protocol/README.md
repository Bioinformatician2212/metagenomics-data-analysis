

This pipeline processes raw paired-end 16S rRNA V3-V4 sequencing data. It includes quality control, trimming, merging, filtering, chimera removal, taxonomic classification, and BIOM file generation.

---

## Prerequisites and Installation

### 1. **FastQC**
FastQC is a quality control tool for high-throughput sequence data.

**Installation**:  
```bash
sudo apt-get install fastqc
```

For detailed instructions, visit the [FastQC website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

---

### 2. **MultiQC**
MultiQC aggregates results from FastQC and other quality tools.

**Installation**:  
```bash
pip install multiqc
```

For more details, visit the [MultiQC GitHub page](https://github.com/ewels/MultiQC).

---

### 3. **TrimGalore**
TrimGalore is used for adapter trimming.

**Installation**:  
```bash
wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.10.tar.gz
tar -xzvf 0.6.10.tar.gz
cd TrimGalore-0.6.10/
chmod +x trim_galore
```

For more details, visit the [TrimGalore GitHub page](https://github.com/FelixKrueger/TrimGalore).

---

### 4. **VSEARCH**
VSEARCH performs merging, filtering, chimera removal, and deduplication.

**Installation**:  
```bash
sudo apt-get install vsearch
```

Alternatively, download the latest binaries from the [VSEARCH GitHub page](https://github.com/torognes/vsearch).

---

### 5. **Kraken2**
Kraken2 is used for taxonomic classification.

**Installation**:  
```bash
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
```

To download the reference database:  
```bash
kraken2-build --standard --db NCBI_refseq
```

For more details, visit the [Kraken2 GitHub page](https://github.com/DerrickWood/kraken2).

---

### 6. **kraken-biom**
kraken-biom converts Kraken2 reports to BIOM files.

**Installation**:  
```bash
pip install kraken-biom
```

For more details, visit the [kraken-biom GitHub page](https://github.com/smdabdoub/kraken-biom).

---

### 7. **R (biomformat and tidyverse libraries)**
R is used to process the BIOM file into an OTU table.

**Installation**:
- Install R:  
  ```bash
  sudo apt-get install r-base
  ```
- Install required R libraries:
  ```r
  install.packages("biomformat")
  install.packages("tidyverse")
  ```

For more details, visit the [R project website](https://www.r-project.org/).

---

## Pipeline Steps

### Step 1: Quality Check of Raw Reads

Run FastQC on the raw reads and summarize the reports with MultiQC.

```bash
#!/bin/bash
mkdir RAWDATA_FastQC_results
fastqc RAWDATA/*.fastq.gz --outdir RAWDATA_FastQC_results
multiqc RAWDATA_FastQC_results/*_fastqc*
```

---

### Step 2: Trimming Adapters and Barcodes

Trim adapters and barcodes using TrimGalore.

```bash
#!/bin/bash
mkdir trimmed_reads
for R1 in *1.fastq.gz; do
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    ./trim_galore --paired $R1 $R2 --output_dir trimmed_reads
done
echo "All samples are processed."
```

**Note**: Adjust the `R1` and `R2` file suffixes based on your file naming convention.

---

### Step 3: Quality Check of Trimmed Reads

Run FastQC and MultiQC on the trimmed reads.

```bash
#!/bin/bash
mkdir trimmed_reads_FastQC_results
fastqc trimmed_reads/*.fq.gz --outdir trimmed_reads_FastQC_results
multiqc trimmed_reads_FastQC_results/*_fastqc*
```

---

### Step 4: Merging Paired-End Reads

Merge paired-end reads using VSEARCH.

```bash
#!/bin/bash
mkdir merged_files
for R1 in trimmed_reads/*_1_val_1.fq.gz; do
    R2="${R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"
    Base_Name=$(basename "$R1" _1_val_1.fq.gz)
    vsearch --fastq_mergepairs "$R1" --reverse "$R2" --fastq_qmax 46 --fastqout "merged_files/${Base_Name}_merged.fq.gz"
done
```

**Note**: Adjust the `R1` and `R2` file suffixes based on your file naming convention.

---

### Step 5: Quality Filtering of Merged Reads

Filter reads by quality, reducing error rates.

```bash
#!/bin/bash
mkdir q_filtered
for R1 in merged_files/*_merged.fq.gz; do
    base_name=$(basename "$R1" "_merged.fq.gz")
    vsearch --fastq_filter "$R1" --fastq_maxee 2.0 --fastq_qmax 46 --fastqout "q_filtered/${base_name}__merged_qual.fastq.gz"
done
echo "All files are filtered."
```

---

### Step 6: Removing Duplicate Reads

Remove duplicate reads using VSEARCH.

```bash
#!/bin/bash
mkdir unique_reads
for R1 in q_filtered/*_merged_qual.fastq.gz; do
    base_name=$(basename "$R1" "_merged_qual.fastq.gz")
    vsearch --fastx_uniques "$R1" --fastqout "unique_reads/${base_name}__merged_unique.fastq.gz" -sizeout -uc uc_out
done
echo "All duplicate reads are removed."
```

---

### Step 7: Removing Chimeric Reads

Detect and remove chimeric sequences.

```bash
#!/bin/bash
mkdir -p chimeras_removed
echo "Aborted files:" > aborted_files.log
for R1 in unique_reads/*___merged_unique.fastq.gz; do
    base_name=$(basename "$R1" "___merged_unique.fastq.gz")
    vsearch --chimeras_denovo "$R1" --nonchimeras "chimeras_removed/${base_name}_merged_unique_non_chimeric.fastq.gz"
    if [ $? -ne 0 ]; then
        echo "$R1" >> aborted_files.log
    fi
done
echo "All files processed. Aborted files listed in aborted_files.log."
```

---

### Step 8: Taxonomic Classification

Classify sequences against a 16S reference database using Kraken2.

```bash
#!/bin/bash
mkdir kraken_reports
krakendb=NCBI_refseq
for R1 in chimeras_removed/*_merged_unique_non_chimeric.fastq.gz; do
    base_name=$(basename "$R1" "_merged_unique_non_chimeric.fastq.gz")
    kraken2 --use-names --db $krakendb --threads 16 "$R1" --output kraken_reports/${base_name}.kraken --report kraken_reports/${base_name}.report
done
echo "All files classified."
```

---

### Step 9: Generating BIOM File

Generate a BIOM file from Kraken2 reports.

```bash
kraken-biom kraken_reports/*.report --fmt json -o filename.biom
```

---

### Step 10: Convert BIOM File to OTU Table in R

```r
# Load necessary libraries
library(biomformat)
library(tidyverse)

# Read the BIOM file
biom_data <- read_biom("filename.biom")

# Extract OTU table and taxonomy
otu_table <- as.data.frame(as.matrix(biom_data(biom_data)))
taxonomy <- as.data.frame(observation_metadata(biom_data))

# Combine taxonomy and OTU table
otu_with_taxonomy <- taxonomy %>%
  mutate(genus_species = paste(gsub("g__", "", taxonomy6), gsub("s__", "", taxonomy7), sep = "_")) %>%
  inner_join(otu_table, by = c("row.names" = "row.names"))

# Save the result
write.csv(otu_with_taxonomy, "otu_table_with_taxonomy.csv", row.names = FALSE)
```




