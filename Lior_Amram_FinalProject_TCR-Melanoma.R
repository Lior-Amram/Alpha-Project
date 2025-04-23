#### 1 Loading the packages: ####
library(igraph)
library(dplyr)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(stringdist)
library(caret)
library(randomForest)
library(visNetwork)
library(tidyr)
library(immunarch)
library(ggpubr)
library(vegan)
library(pheatmap)
library(ComplexHeatmap)

#### 2 Loading the dataset: ####
TCR_melanoma <- repLoad('r_files/adaptive_melanoma')

#### 3 Analysis ####
#### 3.1 Show all samples with number of CDR3.aa sequences ####
# before subsampling
melanoma_clntps <- repExplore(TCR_melanoma$data,.method = 'volume') %>% left_join(.,TCR_melanoma$meta,by='Sample')

#set criteria for discarding samples:
#a remove sample with M6 (no other samples to compare with)
#b remove samples under N clonotypes
#c remove duplicated samples

TCR_melanoma$data <- TCR_melanoma$data[-79]
TCR_melanoma$meta <- TCR_melanoma$meta %>% filter(Sample != '37_6m')
TCR_melanoma$data <- TCR_melanoma$data[-which(names(TCR_melanoma$data) %in% melanoma_clntps[25,1])]
TCR_melanoma$meta <- TCR_melanoma$meta %>% filter(!Sample == melanoma_clntps[25,1])

### subsampling ###
TCR_melanoma_sub <- list(data = repSample(TCR_melanoma$data,'downsample'),
                         meta = TCR_melanoma$meta)
TCR_melanoma_sub$meta <- TCR_melanoma_sub$meta %>% mutate(response_new = ifelse(str_detect(Response ,'CR|PR'),'Responders','Non-responders'))
melanoma_clntps_sub <- repExplore(TCR_melanoma_sub$data,.method = 'volume') %>% left_join(.,TCR_melanoma_sub$meta,by='Sample')

#plot (fig1ab - in results)

num_of_sequences_raw <- ggplot(melanoma_clntps,aes(x=Sample,y=Volume)) +
  theme_bw() +
  geom_bar(stat='identity',) +
  scale_fill_brewer(palette="Set1")+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x='Melanoma blood Samples',y='Number of unique CDR3 amino acid sequences')
num_of_sequences_sub <- ggplot(melanoma_clntps,aes(x=Sample,y=Volume)) +
  theme_bw() +
  geom_bar(stat='identity',) +
  scale_fill_brewer(palette="Set1")+
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x='Melanoma Samples',y='Number of unique CDR3 amino acid sequences')

num_of_sequences_raw+num_of_sequences_sub

##### richness analysis + shannon diversity ####

repExplore(TCR_melanoma_sub$data) %>% 
  left_join(., TCR_melanoma_sub$meta, by = 'Sample') %>% 
  ggplot(aes(x = response_new, y = Volume)) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(~TP, scales = 'free_x') + 
  stat_compare_means(aes(group = response_new), 
                     method = "wilcox.test", # Wilcoxon test for two groups
                     label = "p.signif") + 
  theme_bw() + 
  labs(x='',y='Number of Clonotypes (Unique CDR3-beta amino acid sequences)') +
  ggtitle('Richness CDR3-beta AA comparison between Responders to Non-responders at different timepoints')

Shannon_div_DF <- data.frame(Shannon_div = unlist(lapply(TCR_melanoma_sub$data,function(x) diversity(x$Proportion)))) %>% 
  mutate(Sample = rownames(.)) %>% left_join(.,TCR_melanoma_sub$meta,by='Sample')

ggplot(Shannon_div_DF,aes(x = response_new, y = Shannon_div)) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(~TP, scales = 'free_x') + 
  stat_compare_means(aes(group = response_new), 
                     method = "wilcox.test", # Wilcoxon test for two groups
                     label = "p.signif") + 
  theme_bw() + 
  labs(x='',y='Shannon-diversity') +
  ggtitle('Shannon-diversity CDR3-beta AA comparison between Responders to Non-responders at different timepoints')


# repertoire overlap analysis - Jaccard
melanoma_jaccard_mtx <- repOverlap(TCR_melanoma_sub$data,'jaccard',.col = 'aa') %>% vis(.text.size = 0)
melanoma_jaccard_mtx

### NETWORK ANALYSIS ###
pub_melanoma_no_singletons <- publicRepertoire(TCR_melanoma_sub$data,'aa','prop') %>% as.data.frame() %>% filter(Samples > 1)
# create levenshtein distance matrix
pub_melanoma_no_singletons_sdist <- stringdistmatrix(pub_melanoma_no_singletons$CDR3.aa,pub_melanoma_no_singletons$CDR3.aa,'lv',useNames = T)
# load antigen database - filter to Human TRB sequences and Cancer
mcpas = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB") %>% filter(Category == 'Cancer') %>% filter(CDR3.beta.aa != 'NA')
# keep only sequences with lv-ED1
pub_melanoma_no_singletons_sdist[pub_melanoma_no_singletons_sdist>1]=NA # rmv dist > 1 
pub_melanoma_no_singletons_sdist[pub_melanoma_no_singletons_sdist==0]=NA # rmv dist = 0
pairs_pub_melanoma_no_singletons_sdist <- unique(melt(pub_melanoma_no_singletons_sdist,na.rm = TRUE,varnames = c('seq1','seq2'),as.is = T)) # melt mtx and get unique pairs
#creating lv CDR3.aa.beta network
melanoma_sub_G <- graph_from_data_frame(edge_list_TCR_melanoma_sub,directed = F)
V(melanoma_sub_G)$seq_len <- nchar(V(melanoma_sub_G)$name)


# Convert relevant data frames to data.table
pub_dt <- as.data.table(pub_melanoma_no_singletons)
mcpas_dt <- as.data.table(mcpas)
meta_dt <- as.data.table(TCR_melanoma_sub$meta)

# Set keys for fast filtering
setkey(pub_dt, CDR3.aa)
setkey(meta_dt, Sample)
setkey(mcpas_dt, CDR3.beta.aa)

# Extract all vertex names in advance
vertex_names <- V(melanoma_sub_G)$name
num_vertices <- length(vertex_names)

# Precompute `cnames`
cnames_list <- lapply(vertex_names, function(v) {
  tmp <- pub_dt[J(v), nomatch = 0L]
  if (nrow(tmp) == 0) return(character(0))
  valid_cols <- names(tmp)[colSums(!is.na(tmp)) > 0] # Keep non-empty columns
  valid_cols[3:length(valid_cols)] # Exclude first two columns
})

# Precompute `pids`, `tms`, and `response`
meta_filtered_list <- lapply(cnames_list, function(cnames) {
  if (length(cnames) == 0) return(NULL)
  meta_dt[Sample %in% cnames]
})

pids_list <- lapply(cnames_list, function(cnames) unique(str_extract(cnames, "[^_]+")))
tms_list <- lapply(cnames_list, function(cnames) unique(str_extract(cnames, "(?<=_)[0-9]+")))
response_list <- lapply(meta_filtered_list, function(meta) {
  if (is.null(meta) || nrow(meta) == 0) return(NA)
  responses <- unique(meta$response_new)
  if (length(responses) > 1) return('Mixed')
  return(responses)
})

# Precompute `ismcpas` and `epitope`
ismcpas_list <- vertex_names %in% mcpas_dt$CDR3.beta.aa
epitope_list <- lapply(vertex_names, function(v) {
  if (v %in% mcpas_dt$CDR3.beta.aa) {
    return(mcpas_dt[J(v), Pathology, nomatch = NA][[1]])  # Extract first match
  }
  return(NA)
})

# Assign computed values to vertices
V(melanoma_sub_G)$sids <- cnames_list
V(melanoma_sub_G)$pids <- pids_list
V(melanoma_sub_G)$tps <- tms_list
V(melanoma_sub_G)$res <- unlist(response_list)
V(melanoma_sub_G)$ismcpas <- ismcpas_list
V(melanoma_sub_G)$epitope <- epitope_list

# add vertex attributes - clustering by louvain Identify clusters in the global network
clusters <- cluster_louvain(melanoma_sub_G)

# Assign cluster membership to nodes
V(melanoma_sub_G)$cluster <- clusters$membership




# Function to check if a node is cancer-associated
is_cancer_associated <- function(seq, mcpas_sequences) {
  min_dist <- min(stringdist(seq, mcpas_sequences, method = "lv"))
  return(ifelse(min_dist <= 1, 1, 0))
}

# Assign cancer association attribute
V(melanoma_sub_G)$cancer_assoc <- sapply(V(melanoma_sub_G)$name, is_cancer_associated, mcpas$CDR3.beta.aa)

# function to Extract node attributes for classification
extract_features <- function(network) {
  data.frame(
    name = names(V(network)),  
    degree = degree(network),
    betweenness = betweenness(network),
    closeness = closeness(network),
    #pagerank = page_rank(network)$vector,
    cluster = membership(cluster_louvain(network)),
    cancer_assoc = V(network)$cancer_assoc,  
    sample_count = sapply(V(network)$sids, length),  # Count samples per node
    patient_count = sapply(V(network)$pids, length),  # Count patients per node
    timepoint_count = sapply(V(network)$tps, length),  # Count timepoints per node
    epitope = ifelse(is.na(V(network)$epitope), 0, 1),
    response = as.factor(V(network)$res)  
  )
}


# Extract features
data_ml <- extract_features(melanoma_sub_G)
ext_features <- data_ml


#ML - RF - task1 response classification
# Split into train and test sets
set.seed(42)
trainIndex <- createDataPartition(data_ml$response, p = 0.7, list = FALSE)
trainData <- data_ml[trainIndex, ]
testData <- data_ml[-trainIndex, ]

# Train a Random Forest model
rf_model <- randomForest(response ~ ., data = trainData, ntree = 500)

# Predict on test data
predictions <- predict(rf_model, testData)

# Evaluate performance
confusionM_T1 <- confusionMatrix(predictions, testData$response)
confusionM_T1 <- confusionM_T1$table

bal_acc_T1 <- mean(diag(confusionM_T1) / rowSums(confusionM_T1))


### ML - RF - task2 - cancer associated classification

# Convert categorical variables to factors
ext_features$response <- as.factor(ext_features$response)
ext_features$cancer_assoc <- factor(ext_features$cancer_assoc, levels = c(0, 1), labels = c("Not Cancer Associated", "Cancer Associated"))

# One-hot encoding response for Task 2

dummy_vars <- model.matrix(~ response - 1, data = ext_features)
ext_features <- cbind(ext_features, dummy_vars)
ext_features$response <- NULL


# Split Data
set.seed(42)
colnames(ext_features)[12] <- 'response_NonResponders'
trainIndex <- createDataPartition(ext_features$cancer_assoc, p = 0.7, list = FALSE)
trainData <- ext_features[trainIndex, ]
testData <- ext_features[-trainIndex, ]


rf_model2 <- randomForest(cancer_assoc ~ ., data = trainData, ntree = 500)

# Predict on test data
predictions <- predict(rf_model2, testData)

# Evaluate performance
confusionM_T2 <- confusionMatrix(predictions, testData$cancer_assoc)
confusionM_T2 <- confusionM_T2$table

bal_acc_T2 <- mean(diag(confusionM_T2) / rowSums(confusionM_T2))



#### plot network ####
plot(clusters, melanoma_sub_G, vertex.label=NA, vertex.size=5)
