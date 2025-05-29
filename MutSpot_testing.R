getwd()
setwd("")
rm(list = ls())

##install.packages("devtools")
library(devtools)
# install_github("skandlab/MutSpot", subdir="MutSpot_Rpackage")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


## install.packages("nloptr", type = "source", verbose = TRUE)
 # install.packages("ragg", type = "source")
## install.packages("textshaping", type = "source")
 # install.packages("Matrix")
 # library(nloptr)
 # library(quantreg)

library(data.table) ## faster than numpy and pandas usually
library(BSgenome.Hsapiens.UCSC.hg38) ## full human genome reference
library(Biostrings) ## allow DNA string manipulation and sequence extraction
library(GenomicRanges) ## genomic intervals
library(rtracklayer)
library(MutSpot)

# mutspot_dir <- "/Users/jj75/Desktop/Singapore Coop/MutSpot_repo/MutSpot/MutSpot_Rpackage/R"
# r_files <- list.files(mutspot_dir, pattern = "\\.R$", full.names = TRUE)
# sapply(r_files, source)


## detect and remove out-of-bound mutations
# mut <- fread("mutation_file/gastric_RF_nonMSI_prefiltered.MAF/gastric_RF_nonMSI_prefiltered_new.MAF", header = FALSE)
# setnames(mut, c("chr", "start", "end", "ref", "alt", "sample_id", "cohort"))
chr_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38) ## retrieve the actual length of each chromosome
# chr_lengths
# valid_chr <- intersect(mut$chr, names(chr_lengths))
mut_filtered <- fread("cleaned.MAF")
setnames(mut_filtered, c("chr", "start", "end", "ref", "alt", "sample_id", "cohort"))
# head(mut_filtered)
## fwrite(mut_filtered, "cleaned.MAF", sep = "\t", col.names = FALSE)

# 5-mer centered at each mutation
gr5 <- GRanges(seqnames = mut_filtered$chr, ranges = IRanges(mut_filtered$start - 2, mut_filtered$start + 2))

# Character vector of 5bp sequences
seqs5 <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr5))
head(seqs5)

# head(seqs5)
## get 5-mer centered at the mutation site from the reference genome
## seqs 5 is a vector of DNA strings like: TAAGG

## build sequencing context feature matrix: return a 6-column matrix of binary indicators
seq_mat <- data.frame(
  seq_CG  = as.integer(substr(seqs5, 2, 3) == "CG"),
  seq_AAG = as.integer(substr(seqs5, 2, 4) == "AAG"),
  seq_AA  = as.integer(substr(seqs5, 2, 3) == "AA"),
  seq_GA  = as.integer(substr(seqs5, 3, 4) == "GA"),
  seq_CA  = as.integer(substr(seqs5, 2, 3) == "CA"),
  seq_A   = as.integer(substr(seqs5, 3, 3) == "A")
)

# colnames(feature_mat) <- paste0("SEQ_", colnames(feature_mat))

# colnames(feature_mat) <- make.names(colnames(feature_mat), unique = TRUE)



## Prepare GRanges for mutated sites
gr_mut <- GRanges(seqnames = mut_filtered$chr, ranges = IRanges(mut_filtered$start, mut_filtered$end))

## Load Bed and compute epigenetic overlaps
epi_meta <- fread("fixed_epigenetic.txt")
epi_matrix <- matrix(0, nrow = length(gr_mut), ncol = nrow(epi_meta))
## a matrix of shape [number of mutations, number of epigenetic features]
colnames(epi_matrix) <- epi_meta$feature_name

for (i in seq_len(nrow(epi_meta))) {
  bed_gr <- import(epi_meta$file_path[i]) # import() from rtracklayer will create a GRanges object
  seqlevelsStyle(bed_gr) <- seqlevelsStyle(gr_mut) # standardize naming style of chromosomes
  
  # remove chromosomes from bed_gr that do not exist in gr_mut
  bed_gr <- keepSeqlevels(bed_gr, intersect(seqlevels(bed_gr), seqlevels(gr_mut)), pruning.mode = "coarse")
  
  ## building a binary matrix indicating if a mutation lies within an epigenetic feature region
  ## findOverlaps will return two columns:
  ## first is named quaerHits, and is the index of mutation site that is contained within some epigenetic feature
  ## second is named subjectHits, and is the index/name of the epigenetic feature
  overlap_idx <- queryHits(findOverlaps(gr_mut, bed_gr))
  epi_matrix[overlap_idx, i] <- 1
}

# bed_gr <- import(epi_meta$file_path[1])
# seqlevelsStyle(bed_gr)
# seqlevelsStyle(gr_mut)

## merge sequencing context and epigenetic matrix

## Combines everything into a final feature matrix for mutated sites.
final_mut_matrix <- cbind(seq_mat, epi_matrix)



## non-mutated sites creation

## !!!! below method is commented out for now because it can reach vector memory limit of 16 GB and we will do it on the server

# sample_background_sites <- function(mut_gr, genome = BSgenome.Hsapiens.UCSC.hg38,
#                                     n_sample = 1e6, seed = 42) {
#   set.seed(seed)
#   
#   # 1. Get full genome
#   chr_lengths <- seqlengths(genome)
#   genome_gr <- GRanges(seqnames = names(chr_lengths),
#                        ranges = IRanges(start = 1, end = chr_lengths))
#   
#   # 4. Exclude real mutated positions
#   genome_gr <- setdiff(genome_gr, mut_gr)
#   
#   # 5. Unlist to per-base resolution
#   ## genome_gr is a GRanges object like: chr1 start = 1000, end = 1005
#   ## tile function splits it into 5 1bp sites
#   ## then this returns a GRangeList object
#   ## unlist will break that GRangeList and return a single GRange object with each row representing one site
#   genome_1bp <- unlist(tile(genome_gr, width = 1))
#   
#   message("Total positions: ", length(genome_1bp))
#   
#   # 6. Randomly sample background sites
#   if (length(genome_1bp) < n_sample) {
#     warning("Fewer mappable sites than requested samples. Returning all available.")
#     sampled_sites <- genome_1bp
#   } else {
#     sampled_sites <- genome_1bp[sample(seq_along(genome_1bp), n_sample)]
#   }
#   
#   return(sampled_sites)
# }
# 
# 
# bg_sites <- sample_background_sites(gr_mut,
#                                     n_sample = 1e6)

bg_sites <- GRanges(seqnames = mut_filtered$chr, ranges = IRanges(mut_filtered$start + 10000L, width = 1)) ## offset each mutation by 10000 bp to generate non-mutated sites
bg_sites <- bg_sites[!bg_sites %over% gr_mut] ## Remove any background site that overlaps a true mutation site

gr5_bg <- resize(bg_sites, width = 5, fix = "center") ## prepare into 5-mer window
in_bounds <- start(gr5_bg) >= 1 & end(gr5_bg) <= chr_lengths[as.character(seqnames(gr5_bg))]
gr5_bg <- gr5_bg[in_bounds]
bg_sites <- bg_sites[in_bounds]


seqs_bg <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr5_bg))
# head(seqs_bg)

feature_bg <- data.frame(
  CG_G  = as.integer(substr(seqs_bg, 2, 3) == "CG"),   # G is mutation site (C before it)
  AAG_middle = as.integer(substr(seqs_bg, 2, 4) == "AAG"),  # middle A is mutation site
  AA_right  = as.integer(substr(seqs_bg, 2, 3) == "AA"),   # second A is mutated
  GA_G  = as.integer(substr(seqs_bg, 3, 4) == "GA"),   # G is mutation site (A follows)
  CA_right  = as.integer(substr(seqs_bg, 2, 3) == "CA"),   # A is mutation site (C before it)
  A_single   = as.integer(substr(seqs_bg, 3, 3) == "A")     # 1-mer: mutation itself is A
)

# colnames(feature_bg) <- paste0("SEQ_", colnames(feature_bg))
colnames(feature_bg) <- make.names(colnames(feature_bg), unique = TRUE)

## !! how does the model know if for GA, mutation site is in G or A:
## it knows because we only indicate 1 if substring(3,4) is GA. which means, only indicate 1 if G is the mutation site

## epigenetic features for non-mutated
epi_bg_matrix <- matrix(0, nrow = length(bg_sites), ncol = nrow(epi_meta))
colnames(epi_bg_matrix) <- epi_meta$feature_name

for (i in seq_len(nrow(epi_meta))) {
  bed_gr <- import(epi_meta$file_path[i])
  seqlevelsStyle(bed_gr) <- seqlevelsStyle(bg_sites)
  bed_gr <- keepSeqlevels(bed_gr, intersect(seqlevels(bed_gr), seqlevels(bg_sites)), pruning.mode = "coarse")
  overlap_idx <- queryHits(findOverlaps(bg_sites, bed_gr))
  epi_bg_matrix[overlap_idx, i] <- 1
}

## merge non-mutated features and save
final_bg_matrix <- cbind(feature_bg, epi_bg_matrix)

# nrow(final_bg_matrix)
# nrow(final_mut_matrix)

## local mutation rate
bin_size <- 100000L ## 100kb

## 100kb bins across genome
bins_list <- lapply(names(chr_lengths), function(chr) {
  chr_len <- chr_lengths[[chr]]
  starts <- seq(1, chr_len, by = bin_size)
  ends <- pmin(starts + bin_size - 1, chr_len)
  gr <- GRanges(seqnames = chr, ranges = IRanges(starts, ends))
  # seqlevels(gr) <- chr  # make sure it’s only this chr
  # seqlengths(gr) <- chr_len
  gr
})



all_bins <- unlist(GRangesList(bins_list))


## mutation in each bin
# mut_gr <- GRanges(seqnames = mut_filtered$chr, ranges = IRanges(mut_filtered$start, mut_filtered$end))
bin_mut_counts <- countOverlaps(all_bins, gr_mut)
mutation_density <- bin_mut_counts / bin_size

## assign to each mutation and background site
## mut_bin_idx[i] returns the index of all_bins that gr_mut[i] belongs to
mut_bin_idx <- findOverlaps(gr_mut, all_bins, select = "first")
bg_bin_idx <- findOverlaps(bg_sites, all_bins, select = "first")

## comment below two lines to bin the local mutation rate for now to avoid vector limit issue
final_mut_matrix$local_mut_rate <- mutation_density[mut_bin_idx]
final_bg_matrix$local_mut_rate <- mutation_density[bg_bin_idx]

# # Assign mutation rate
# final_mut_matrix$local_mut_rate <- mutation_density[mut_bin_idx]
# final_bg_matrix$local_mut_rate  <- mutation_density[bg_bin_idx]
# 
# # Combine to compute bins on entire range
# combined_local_mut_rate <- c(final_mut_matrix$local_mut_rate, final_bg_matrix$local_mut_rate)
# 
# # Create 10 quantile-based bins
# nbins <- 10
# bin_breaks <- quantile(combined_local_mut_rate, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
# bin_breaks[1] <- bin_breaks[1] - 1e-8  
# 
# # Assign bin labels (1 to nbins)
# mut_bins <- cut(final_mut_matrix$local_mut_rate, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)
# bg_bins  <- cut(final_bg_matrix$local_mut_rate, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)
# 
# # Add dummy variables
# for (i in 1:nbins) {
#   final_mut_matrix[[paste0("local_mut_bin_", i)]] <- as.integer(mut_bins == i)
#   final_bg_matrix[[paste0("local_mut_bin_", i)]]  <- as.integer(bg_bins == i)
# }
# 
# # Drop original continuous column
# final_mut_matrix$local_mut_rate <- NULL
# final_bg_matrix$local_mut_rate  <- NULL

# Bin local mutation rate into 10 quantiles and assign bin means
combined_local_mut_rate <- c(final_mut_matrix$local_mut_rate, final_bg_matrix$local_mut_rate)
nbins <- 10
bin_breaks <- quantile(combined_local_mut_rate, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE) ## bin edges of the 10 bins
bin_breaks[1] <- bin_breaks[1] - 1e-8
mut_bins <- cut(final_mut_matrix$local_mut_rate, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE) ## assign each site to a bin index
bg_bins <- cut(final_bg_matrix$local_mut_rate, breaks = bin_breaks, labels = FALSE, include.lowest = TRUE)
bin_means <- tapply(combined_local_mut_rate, cut(combined_local_mut_rate, breaks = bin_breaks, include.lowest = TRUE), mean)
final_mut_matrix$local_mut_rate <- bin_means[mut_bins]
final_bg_matrix$local_mut_rate <- bin_means[bg_bins]

# Label and combine for response matrix
final_mut_matrix$mut.count <- 1
final_mut_matrix$nonmut.count <- 0
final_bg_matrix$mut.count <- 0
final_bg_matrix$nonmut.count <- 1
colnames(final_bg_matrix) <- colnames(final_mut_matrix)
combined_matrix <- rbind(final_mut_matrix, final_bg_matrix)

# Aggregate by unique covariate combinations
response_agg <- aggregate(cbind(mut.count, nonmut.count) ~ ., data = combined_matrix, FUN = sum)
response_agg <- response_agg[!(response_agg$mut.count == 0 & response_agg$nonmut.count == 0), ] ## exclude zero sites

# Extract matrices
response_matrix <- response_agg[, c("mut.count", "nonmut.count")]
covariate_matrix <- response_agg[, !(colnames(response_agg) %in% c("mut.count", "nonmut.count"))]

dir.create("results", showWarnings = FALSE)
saveRDS(covariate_matrix, "results/mutCovariate-sparse-p1.RDS")
saveRDS(response_matrix, "results/mutCovariate-sparse-p2.RDS")
saveRDS(colnames(seq_mat), "results/nucleotide_selected.RDS")
saveRDS(bg_sites, "results/sampled.sites.snv.RDS")

# rm(list = ls())
features <- readRDS("results/mutCovariate-sparse-p1.RDS")
responses <- readRDS("results/mutCovariate-sparse-p2.RDS")

features$local_mut_rate <- as.numeric(features$local_mut_rate)

mutfreq.aggregated <- cbind(responses, features)

LRmodel <- glm(cbind(mut.count, nonmut.count) ~ ., 
               family = binomial(link = "logit"), 
               data = mutfreq.aggregated)

LRmodel$coefficients
## try copy paste MutSpot's codes for plot of variable importance
t = summary(LRmodel)
df = coef(t)
df = as.data.frame(df)
df$feat = rownames(df)
df = df[-1, ]
colnames(df) = c("estimate", "std.error", "z.value", "p.value", "feat")
df$feat = gsub("`", "", gsub("`1", "", df$feat))
# df$feat = ifelse(!df$feat %in% feature.names, substr(df$feat, 1, nchar(df$feat) - 1), df$feat)

df = df[order(df$z.value, decreasing = TRUE), ]
df$feat = factor(df$feat, levels = df$feat)
colnames(df)[3] = "imp"
pdf(paste("results/", "SNV", "_features.pdf", sep = ""))
print(ggplot2::ggplot(df, ggplot2::aes(x = feat, y = imp)) + ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black"))+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
        ggplot2::xlab("Features") + 
        ggplot2::ylab("Variable importance"))
dev.off()
## MutSpot's codes produce the same result as mine below codes 


## try to calculate McFadden's pseudo-R squared

# Null model (intercept only)
model_null <- glm(cbind(mut.count, nonmut.count) ~ 1, 
                  family = binomial(link = "logit"), data = mutfreq.aggregated)

# McFadden's pseudo-R²
R2_McFadden <- 1 - (as.numeric(logLik(LRmodel)) / as.numeric(logLik(model_null)))
print(R2_McFadden)

# logLik(LRmodel)
# logLik(model_null)

# colnames(final_mut_matrix) <- make.names(colnames(final_mut_matrix), unique = TRUE)
# colnames(final_bg_matrix)  <- make.names(colnames(final_bg_matrix), unique = TRUE)
# 
# n = min(nrow(final_bg_matrix), nrow(final_mut_matrix), 50000) ## the matrix for muted and non-muted must be same
# 
# print(n)
# 
# final_bg_matrix <- final_bg_matrix[1:n, ]
# final_mut_matrix <- final_mut_matrix[1:n, ]
# bg_sites <- bg_sites[1:n, ]

# combined_names <- colnames(feature_mat)
# epi_names <- colnames(epi_matrix)
# duplicated_names <- intersect(combined_names, epi_names)
# print(duplicated_names)
# any(duplicated(colnames(final_mut_matrix)))
# any(duplicated(colnames(final_mut_matrix)))
# any(duplicated(colnames(final_bg_matrix)))

# x <- readRDS("results/mutCovariate-sparse-p1.RDS")
# any(duplicated(colnames(x)))


# combined_matrix <- rbind(final_mut_matrix, final_bg_matrix)
# 
# saveRDS(as.data.frame(combined_matrix), "results/mutCovariate-sparse-p1.RDS")
# 
# n <- nrow(final_mut_matrix)
# p2 <- data.frame(
#   mut.count = c(rep(1, n), rep(0, n)),
#   nonmut.count = c(rep(0, n), rep(1, n))
# )
# saveRDS(p2, "results/mutCovariate-sparse-p2.RDS")
# 
# 
# 
# ## save metadata
# saveRDS(colnames(feature_mat), "results/nucleotide_selected.RDS")


# saveRDS(bg_sites, "results/sampled.sites.snv.RDS")
# 
# 
# # nrow(final_bg_matrix) == length(bg_sites)
# 
# # setwd("/Users/jj75/Desktop/Singapore Coop/MutSpot_testing")
# getwd()

# y <- readRDS("results/mutCovariate-sparse-p2.RDS")
# head(colnames(y))

# 
# MutSpot(run.to = 6, genome.build = "Ch38", cores = 3)

## there are exceeding vector limit problem when run genome-wide
## now running 1000000 size instead
 
# features <- readRDS("results/mutCovariate-sparse-p1.RDS")
# responses <- readRDS("results/mutCovariate-sparse-p2.RDS")
# str(features$local_mut_rate)
# mm <- model.matrix(~ ., data = features)
# colnames(mm)

# mat <- matrix(rnorm(5000000 * 50), nrow = 5000000, ncol = 50)
# print(object.size(mat), units = "auto")
# 
# getwd()
# model_path <- "./results/snv-LRmodel"
# load(model_path)
# print(head(LRmodel$coef, 50))
# length(LRmodel$coef)
# max(LRmodel$coef)
# 
# sum(responses$mut.count)
# sum(responses$nonmut.count)        
