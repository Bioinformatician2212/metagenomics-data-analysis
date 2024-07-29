# Metagenomics Data Analysis

This project provides a comprehensive analysis of metagenomics data with a focus on alpha and beta diversity. The analysis is performed using various statistical tests and visualizations to understand the microbial community structure. Below is an outline of the workflow and the tools used in this analysis.

# Workflow

 1. Data Preparation
-OTU Table: The Operational Taxonomic Unit (OTU) table contains the abundance of different OTUs across samples. The rows represent OTUs, and the columns represent samples.
- Metadata Table: The metadata table includes sample-specific information such as environmental conditions, sample type, or any other relevant variables.

 2. Alpha Diversity Analysis
Alpha diversity measures the diversity within a single sample. We use the following indices:
- Shannon Index: Measures the entropy or uncertainty in the species distribution.
- Simpson Index: Measures the probability that two randomly selected individuals belong to the same species.

*Statistical Tests for Alpha Diversity:*
- Kruskal-Wallis Test: Non-parametric test used to compare differences in alpha diversity indices between multiple groups.
- Generalized Linear Model (GLM): Used to model the relationship between alpha diversity and explanatory variables.

 3. Beta Diversity Analysis
Beta diversity assesses the diversity between different samples or groups. We utilize the following methods:
- Bray-Curtis Dissimilarity: Measures the dissimilarity between pairs of samples based on species abundance.
- Principal Coordinates Analysis (PCoA): Visualizes the differences between samples in a reduced dimensional space.

*Statistical Test for Beta Diversity:*
- Adonis Test: A permutation-based test to assess the significance of the observed differences in beta diversity.

# Libraries Used
- ggplot2: For creating advanced visualizations.
- vegan: For performing diversity analysis and statistical tests.
- phyloseq: For handling and analyzing microbiome data.
- rstatix: For additional statistical tests and data manipulation.

# Example Analysis
# Load necessary libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(rstatix)

# Define file paths
otu_file <- "C:\\Users\\Jayesh Punde\\Desktop\\analysis\\new assignment\\otu_table.csv"

meta_file <- "C:\\Users\\Jayesh Punde\\Desktop\\analysis\\new assignment\\metadata.csv"

# Read in the OTU table and metadata
otu_data <- read.csv(otu_file, row.names = 1, check.names = FALSE)
meta_data <- read.csv(meta_file, row.names = 1)

# Create the phyloseq object
otu_table <- otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)
sample_data <- sample_data(meta_data)
physeq <- phyloseq(otu_table, sample_data)

# Define alpha diversity metric and grouping variable
alpha_metric <- "Simpson"  # Change to "Chao1" or "Simpson" as needed
group_var <- "ANEMIA_STATUS"  # Replace with the column name in your metadata

# Calculate alpha diversity
alpha_div <- estimate_richness(physeq, measures = c("Chao1", "Shannon", "Simpson"))
alpha_div$group <- factor(sample_data(physeq)[[group_var]])

# Add covariates to alpha_div data frame
covariates <- c("AGE", "SEX")  # Replace with the column names of your covariates
alpha_div <- cbind(alpha_div, sample_data(physeq)[, covariates])

# Create boxplot for alpha diversity
ggplot(alpha_div, aes(x = group, y = get(alpha_metric))) +
  geom_boxplot(aes(color = group)) +
  geom_jitter(width = 0.2, aes(color = group)) +
  labs(y = alpha_metric, x = group_var) +
  theme_minimal()

# Perform Kruskal-Wallis test
kruskal_test <- kruskal_test(alpha_div, as.formula(paste0(alpha_metric, " ~ group")))
print(kruskal_test)

# Create formula for GLM
covariates_formula <- paste(covariates, collapse = " + ")
formula <- as.formula(paste0(alpha_metric, " ~ group + ", covariates_formula))

# Perform GLM
glm_result <- summary(glm(formula, data = alpha_div, family = Gamma(link = "log")))
print(glm_result$coefficients)
 
# Define beta diversity metric
beta_metric <- "bray"  # Change to "jaccard" if needed

# Calculate distance matrix
dist_matrix <- vegdist(t(otu_table(physeq)), method = beta_metric)

# Perform PCoA
pcoa <- cmdscale(dist_matrix, k = 2)
pcoa_df <- as.data.frame(pcoa)
# Add grouping information from metadata
pcoa_df$group <- factor(sample_data(physeq)[[group_var]])
colnames(pcoa_df) <- c("PCoA1", "PCoA2", "group")

# Create PCoA plot
ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = group)) +
  geom_point(size = 3) +
  labs(x = "PCoA 1", y = "PCoA 2", color = group_var) +
  theme_minimal() +
  ggtitle("PCoA of Bray-Curtis Distance Matrix")

# Perform PERMANOVA (adonis2) test
adonis_result <- adonis2(dist_matrix ~ ANEMIA_STATUS, data = as(sample_data(physeq), "data.frame"), permutations = 99)
print(adonis_result)

# Data Tables

# OTU Table
- Format: CSV
- Description: Contains OTU counts across samples.

# Metadata Table
- Format: CSV
- Description: Contains sample-specific metadata.

Please refer to the `data/` directory for the raw OTU and metadata files used in this analysis.



Feel free to adjust any details or add additional sections as needed!
















                         
