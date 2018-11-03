library(vegan)
library(HITChipDB)
library(microbiome)
library(tidyverse)
library(here)

# load data and helper functions
source(here("R/bayesian_helper.R"))
source(here("R/mb_helper.R"))
load(here("data/data_transfer.RData"))
source(here("R/read.R"))

# take over the meta variables I created in other docs
meta_new <- data_transfer[, 1:9] 
# for ordination, I create bf versus no bf
qplot(data = meta_new, bf_ratio)
# based on the plot maybe it is sensible to make 3 categories

# just for handson purposes I create categories from the bf_ratio.
meta_new <- meta_new %>%
  mutate(
      groups = ifelse(time == "pre" & cc == "no", "noCCpre", ifelse(
          time == "pre" & cc == "yes", "CCpre", ifelse(
              time == "post" & cc == "no", "noCCpost", "CCpost"))),
      bf = ifelse(bf_ratio <= 0.25, "lowBF", ifelse(
          bf_ratio <0.75, "mediumBF", "highBF"))) %>% 
  mutate(groups = as.factor(groups), bf = as.factor(bf))

# create new pseq object (from read.R results a pseq object Leo created)
otu <- otu_to_df(genus)
otu <- otu %>% select(species, meta_new$sample_id) %>% df_to_otu()
pseq <- phyloseq(otu, df_to_sd(meta_new), tax_table(genus))

pseq
genus

# add diversity indeces to sample data
diversities <- 
    global(pseq, index = "all") %>% 
    select(contains("diversities")) %>% 
    rownames_to_column("sample_id")
colnames(diversities) <- gsub("diversities_", "", colnames(diversities))

sample_data(pseq) <- 
    sd_to_df(pseq) %>% 
    left_join(diversities, by = "sample_id") %>%
    df_to_sd()

# use relative abundance data
pseq_c <- microbiome::transform(pseq, "compositional") 

# PCoA bray
cc_pcoa <- ordinate(pseq_c, method = "PCoA", distance = "bray")
e_values <- cc_pcoa$values$Eigenvalues
plot_ordination(pseq_c, cc_pcoa, color = "groups") +
    geom_point(size = 3) +
    coord_fixed(sqrt(e_values[2] / e_values[1]))

# # PCoA weighted unifrac (phy_tree slot is empty)
# cc_pcoa <- ordinate(pseq_c, method = "PCoA", distance = "wunifrac")
# e_values <- cc_pcoa$values$Eigenvalues
# plot_ordination(pseq_c, cc_pcoa, color = "groups") +
#     geom_point(size = 3) +
#     coord_fixed(sqrt(e_values[2] / e_values[1]))


# MDS
cc_mds <- ordinate(pseq_c, "MDS", "bray")
e_values <- cc_mds$values$Eigenvalues
plot_ordination(pseq_c, cc_mds, color = "groups") +
    geom_point(size = 3) +
    coord_fixed(sqrt(e_values[2] / e_values[1]))

pseq.cca <- ordinate(pseq_c, "CCA", formula = pseq ~ groups)
plot_ordination(pseq_c, pseq.cca, type = "samples", color = "bf") +
  geom_point(size = 3)

# prepare some objects as in the tutorial for plotting
tax <- 
    tax_table(pseq)@.Data %>% 
        data.frame(stringsAsFactors = FALSE)
ps_scores <- vegan::scores(pseq.cca)
sites <- data.frame(ps_scores$sites)
sites$sample_id <- rownames(sites)
sites <- sites %>%
    left_join(sd_to_df(pseq_c), by = "sample_id")
species <- data.frame(ps_scores$species) %>% 
    rownames_to_column("L2") %>%
    left_join(tax, by = "L2")
species$otu_id <- 1:dim(species)[2]

evals_prop <- 100 * pseq.cca$CCA$eig[1:2] / sum(pseq.cca$CCA$eig)

library(ggrepel)
# 
ggplot() +
    geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
    geom_point(data = species, aes(x = CCA1, y = CCA2), size = 0.5) +
    geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
                    size = 1.5, segment.size = 0.1) +
    facet_grid(. ~ groups) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
         y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
    scale_color_brewer(palette = "Set2") +
    coord_fixed(sqrt(pseq.cca$CCA$eig[2] / pseq.cca$CCA$eig[1])*0.33) +
    theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

library(caret)

dataMatrix <- data.frame(groups = sample_data(pseq_c)$groups, t(otu_table(pseq_c)))
pls.fit <- train(groups ~., data = na.omit(dataMatrix), method = "pls", preProc = "center")

# create biplot for interpretation
pls_biplot <- list("loadings" = loadings(pls.fit$finalModel),
                  "scores" = scores(pls.fit$finalModel))
class(pls_biplot$scores) <- "matrix"
pls_biplot$scores <- data.frame(sample_data(pseq_c), pls_biplot$scores)
tax <- tax_table(pseq)@.Data %>%
    data.frame(stringsAsFactors = FALSE)
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)


ggplot() +
    geom_point(data = pls_biplot$scores, 
               aes(x = Comp.1, y = Comp.2), shape = 2) +
    geom_point(data = pls_biplot$loadings,
              aes(x = 25 * Comp.1, y = 25 * Comp.2),
               size = 0.3, alpha = 0.6) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "Axis1", y = "Axis2", col = "bf") +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    facet_grid(~groups) +
    theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

rf.fit <- train(groups ~., data = na.omit(dataMatrix), method = "rf", preProc = "center", proximity = T)
rf.fit

# cc*time vs div
ggplot(sd_to_df(pseq), aes(x = cc, shannon)) +
    geom_violin() +
    geom_jitter(width = 0.05) +
    stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) +
    facet_wrap(.~time) +
    geom_hline(aes(yintercept = mean(shannon)), linetype = "dashed") +
    coord_flip()

# age vs div
ggplot(sd_to_df(pseq), aes(x = age_d_s, shannon, col = cc)) +
    geom_point() +
    geom_smooth(method = "lm", se = F)

# bf vs div
ggplot(sd_to_df(pseq), aes(x = bf_count_s, shannon, col = cc)) +
    geom_point() +
    geom_smooth(method = "lm", se = F)

library(afex)

df <- sd_to_df(pseq)
# check contrasts for interpretation
contrasts(df$cc)
contrasts(df$time)

# mixed function in afex uses contr.sum
mfit <- mixed(shannon ~ cc*time + age_d_s  + (1|subject_id),
            data = df, method = "KR")
# output
mfit$anova_table
summary(mfit)

# alternatively we can use lme4
lme4.fit <- lme4::lmer(shannon ~ cc*time + age_d_s  + (1|subject_id),
            data = df)
summary(lme4.fit)
car::Anova(lme4.fit, test = "F", type = 2)

newdata <- select(df, subject_id, cc, time, age_d_s)
newdata$pred <- predict(lme4.fit, newdata = newdata)

library(DESeq2)
# how to avoid NA formation here??
# ok it seems that the following is the problem:
# R is restricted to integers that are smaller than 2147483648. 
# therefore I apply disvision by 1000, then it works
pseq2 <- pseq
otu_table(pseq2) <- otu_table(pseq2)/1000
ds2 <- phyloseq_to_deseq2(pseq2, ~ time + cc + time:cc + (1|subject_id))

varianceStabilizingTransformation(ds2, blind = TRUE, fitType = "parametric")

ds2 <- estimateSizeFactors(ds2)
ds2 <- estimateDispersions(ds2)
abund <- getVarianceStabilizedData(ds2)

dds <- DESeq(ds2)
res <- results(dds) %>% 
    as.data.frame() %>%
    rownames_to_column("taxon") %>%
    arrange(padj, log2FoldChange)



res

# Running the DESeq2 analysis
ds2 <- phyloseq_to_deseq2(pseq, ~ nationality)
dds <- DESeq(ds2)
res <- results(dds)
df <- as.data.frame(res)
df$taxon <- rownames(df)
df <- df %>% arrange(log2FoldChange, padj)

library(Rtsne)

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"
ps <- microbiome::transform(pseq_c, trans)
dm <- vegdist(otu_table(ps), distance)
tsne_out <- Rtsne(dm, dims = 2)
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(ps))
plot_landscape(proj, legend = T, size = 1)



pseq %>% sample_data()
