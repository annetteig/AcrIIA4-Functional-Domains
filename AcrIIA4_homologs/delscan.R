library(Biostrings)
library(msa)

#There are 3390 unique deletions (if we keep the methionine and stop codon). There are 3741 if you don't remove duplicates.
acr <- "MNINDLIREIKNKDYTVKLSGTDSNSITQLIIRVNNDGNEYVISESENESIVEKFISAFKNGWNQEYEDEEEFYNDMQTITLKSELN"
length(unique(unlist(sapply(2:87, function (x) sapply((x + 1):88, function (y) paste0(substr(acr, 1, x), substr(paste0(acr, "*"), y, 88)))))))
length(unlist(sapply(2:87, function (x) sapply((x + 1):88, function (y) paste0(substr(acr, 1, x), substr(paste0(acr, "*"), y, 88))))))

#Comparing to the chosen gRNA, the closest match in the BY genome that is followed by a PAM only matches at 14 bases.
load("~/path/to/file/BYGuideTable_smaller.RData")
grna <- strsplit("ACGTATCGGTACATGCGCAT", "")[[1]]
max(sapply(cut.sites$guide, function (x) length(which(strsplit(x, "")[[1]][1:20] == grna))))



acriia4s <- readAAStringSet("~/path/seqdump_acriia4s.txt")
acriia4s <- acriia4s[which(sapply(acriia4s, function (x) substr(as.character(x), 1, 1)) == "M")]
acriia4s <- acriia4s[which(nchar(acriia4s) > 70)] #I want to remove acriia4 genes that are short due to sequencing errors, pseudogenes, etc. Basically, if they are implausibly short (70 aa or less when the canonical AcrIIA4 is 87), I will remove them.

#Left with 141 AcrIIA4s. As one of these is the canonical AcrIIA4, I would say we have found 140 homologs.

identify.gaps <- function (query, target) {
  tempmsa <- msa(c(query, as.character(target)), type = "protein")
  which(strsplit(as.character(tempmsa@unmasked[2]), "")[[1]][which(strsplit(as.character(tempmsa@unmasked[1]), "")[[1]] %in% LETTERS)] == "-")
}


all.deletions <- sapply(acriia4s, function (x) identify.gaps(acriia4s[2], x))
writeXStringSet(acriia4s[which(lengths(all.deletions) > 0)], filepath = "~/path/to/file/acriia4s_with_dels.fasta")


table(sapply(all.deletions, function (x) paste0(x, collapse = "_"))) 

#Checking deletions present at the start or end of the gene for whether they were artefactually caused by issues with genome sequencing and assembly.
#The 1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16 deletion isn't real - it is caused by being present on a contig missing the start codon (see https://www.ncbi.nlm.nih.gov/nuccore/DAAIPA010000038.1) or annotators missing the start codon (see https://www.ncbi.nlm.nih.gov/nuccore/1092873060 or https://www.ncbi.nlm.nih.gov/nuccore/NZ_CYWP01000005.1?from=222494&to=223101&report=fasta)
#Same for the 1_2 deletion (see https://www.ncbi.nlm.nih.gov/nuccore/AAAUUU010000083.1 and https://www.ncbi.nlm.nih.gov/nuccore/AAATBZ010000094.1)
#The first deletion in the 1_2_3_4_5_22_23_24_25_26_27 allele (1-5) is similarly caused by the contig missing the start codon (see https://www.ncbi.nlm.nih.gov/nuccore/AAFFKY010000083.1?from=4499&to=4753&report=fasta). The latter deletion (22-27) is represented by a different allele on its own.

#The 84_85_86_87 deletion is from the end of a contig (https://www.ncbi.nlm.nih.gov/protein/EAC3591559.1/ and EAC2201986.1)
#Same for 78_79_80_81_82_83_84_85_86_87 (https://www.ncbi.nlm.nih.gov/protein/EAF0101949.1)
#And 82_83_84_85_86_87 and 18_19_20_21_22_23_82_83_84_85_86_87 (EAG8194425.1, EAD8129392.1, EAD3168424.1, EAD8602602.1)
#And 85_86_87 (EAC6112265.1 and EAH1632085.1)

#The "81" deletions are from intact contigs (WP_185605294.1, EAD7580492.1)
#At least one of the "84" deletions is from an intact contig, so I won't check the rest (EAE7080009.1)
#86_87 is from an intact contig (https://www.ncbi.nlm.nih.gov/nuccore/AAALLZ010000003.1?from=56498&to=56759&report=fasta)

deletions.intact.contigs <- all.deletions[-which(sapply(all.deletions, function (x) paste0(x, collapse = "_")) %in% c("1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16", "1_2", "1_2_3_4_5_22_23_24_25_26_27", "84_85_86_87", "78_79_80_81_82_83_84_85_86_87", "82_83_84_85_86_87", "18_19_20_21_22_23_82_83_84_85_86_87", "85_86_87"))]
acriia4s.intact.contigs <- acriia4s[-which(sapply(all.deletions, function (x) paste0(x, collapse = "_")) %in% c("1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16", "1_2", "1_2_3_4_5_22_23_24_25_26_27", "84_85_86_87", "78_79_80_81_82_83_84_85_86_87", "82_83_84_85_86_87", "18_19_20_21_22_23_82_83_84_85_86_87", "85_86_87"))]

writeXStringSet(acriia4s.intact.contigs, filepath = "~/path/to/file/acriia4.homologs.fasta")

table(sapply(deletions.intact.contigs, function (x) paste0(x, collapse = "_")))

#I took the deletions marked as 25_26_27_28_36 and 26_27_28_29_36 and aligned them to each other and to AcrIIA4 in clustalw.
#The alignment found that the deletion event is in common. I would call it 26_27_28_29_35.
which(sapply(all.deletions, function (x) paste0(x, collapse = "_")) %in% c("26_27_28_29_36", "25_26_27_28_36"))

#As above, the deletion patterns 16_17_30_31_47_48_49_50_51 and 16_17_30_31_50_51_52_53_54 are in fact shared when aligned.
#I would probably call it 19_20_30_31_47_48_49_50_51, though 19_20 could swap for 17_18 and 47_48_49_50_51 for 45_46_47_48_49.

#The deletions 18_19_20_21_22_23, 22_23_24_25_26_27, and 23_24_25_26_27_28 also look shared to me. I would call them 23_24_25_26_27_28.

#The alleles marked as having a single codon missing at 81 or 84 are also clearly related upon alignment.


#ClinVar deletions

clinvar <- fread("~/path/variant_summary.txt", select = c("Type", "Name", "ClinicalSignificance")) #Downloaded Jan 31, 2024 from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
sort(table(clinvar$ClinicalSignificance)) #Decided the top 8 classifications are the most relevant ones
#sapply(names(sort(table(clinvar$ClinicalSignificance), decreasing = T)[1:8]), function (x) sort(table(clinvar$Type[which(clinvar$ClinicalSignificance == x)])/length(clinvar$Type[which(clinvar$ClinicalSignificance == x)])))
#sapply(names(sort(table(clinvar$ClinicalSignificance), decreasing = T)[1:8]), function (x) sort(table(clinvar$Type[which(clinvar$ClinicalSignificance == x)])))
clinvar.table <- t(sapply(names(sort(table(clinvar$ClinicalSignificance), decreasing = T)[1:8]), function (x) sapply(c("single nucleotide variant", "Deletion", "Insertion", "Duplication", "Microsatellite", "copy number gain", "copy number loss", "Indel", "Inversion", "Translocation", "Complex"), function (y) length(which(clinvar$Type[which(clinvar$ClinicalSignificance == x)] == y)))))
clinvar.table <- cbind.data.frame(clinvar.table, Other = sapply(names(sort(table(clinvar$ClinicalSignificance), decreasing = T)[1:8]), function (x) length(which(clinvar$ClinicalSignificance == x))) - apply(clinvar.table, 1, sum))
write.csv(clinvar.table, file = "~/path/to/file/clinvar_table.csv")
