setwd("/home/ericcumunel/Documents/Aster/")

real_quantifications <- read.table(file = 'quantif/transcript_real_quantif')
transcript_names <- read.table(file = "quantif/transcript_names")
real_quantifications["transcript_names"] <- transcript_names$V1

filtered_countings <- read.table(file = 'bam/spike_in/BYK_GDB_ONT_1_PAD67469_Aflipflop.sorted.filtered_80QC.tab')
new_primary_countings <- read.table(file = "bam/spike_in/BYK_GDB_ONT_1_PAD67469_Aflipflop.filter_80QC_new_primary.sorted.tab")
filtered_countings <- read.table(file = 'bam/spike_in/countings.txt')


new_primary_df <- merge(new_primary_countings, real_quantifications, by.x ="V2", by.y ="transcript_names")
colnames(new_primary_df) <- c("transcript_name","nb_of_read", "real_quantif")
cor(new_primary_df$nb_of_read, new_primary_df$real_quantif, method = "spearman")

filtered_df <- merge(filtered_countings, real_quantifications, by.x ="V2", by.y ="transcript_names")
colnames(filtered_df) <- c("transcript_name","nb_of_read", "real_quantif")
cor(filtered_df$nb_of_read, filtered_df$real_quantif, method = "spearman")


test <- read.table(file = 'bam/spike_in/new_counting.tsv')
cor(test$V2, test$V3, method = "spearman")
