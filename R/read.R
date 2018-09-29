# see README in data/hitchip
# library(devtools)
# # Install the HITChipDB from Github
# install_github("microbiome/HITChipDB")
library(HITChipDB)
library(here)
# Reading the L2 (genus) level OTU table from the folder
otu <- read.table(here("data/hitchip/l2-rpa.tab"))

# Reading L3 (species) level data as phyloseq object
species <- import_hitchip(data.dir = here("data/hitchip/"), method = "rpa", detection.threshold = 1e-6, verbose = F)
# Summarize species to genus level
genus <- aggregate_taxa(species, level = "L2")

# Compare the phyloseq version and the version that was in the folder.
# Before this we need to make sure that the genus names are compatible..
otu2 <- abundances(genus)
rownames(otu2) <- gsub(" ", "_", rownames(otu2))
rownames(otu2) <- gsub("\\.", "", rownames(otu2))
rownames(otu2)[which(rownames(otu2) == "Clostridium_\\(sensu_stricto\\)")] <- "Clostridium_sensu_stricto"
otu2 <- otu2[rownames(otu),]

# Correlation is good, data is essentially the same, except some difference in scale but this is not relevant
# if you use relative abundances and perhaps not relevant otherwise either (we are not interested in absolute HITChip abundances)
# cors <- c(); for (k in 1:130) {cors[[k]] <- cor.test(unlist(otu[k,]), unlist(otu2[k,]))$estimate}
# hist(cors)


# The phyloseq object you can then use as usual. Note that HITChip data is not NGS data and
# there is no guarantee that all distributional properties are the same. This may or may not limit the use of some methods.
# library(microbiome)
# core(transform(genus, "compositional"), detection = 0.1/100, prevalence = 50/100)
# plot(density(log10(abundances(genus)["Dialister",])))



