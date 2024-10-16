setwd("/data/SBGE/simone/9-23-24")

library(tidyverse)
library(stringr)

#read in output file of foldx stability (concatenated), format
df <- read.table("relaxed_rank_001/all_results.txt")
df_filtered <- df %>%
  select(V1, V2) %>% 
  rename(pdb = V1, free_energy = V2) %>%  
  mutate(a3m = str_remove(str_extract(pdb, "(?<=\\./).*(?=_relaxed)"), "_relaxed"))  

#read in fasta file that contains all deletions used to build model, format
fasta_lines <- readLines("all_deletions.fasta")

fasta_table <- data.frame(start = numeric(), end = numeric(), length = numeric(), fasta_seq = character(), stringsAsFactors = FALSE)

# data frame from the fasta file
fasta_table <- tibble(lines = fasta_lines) %>%
  mutate(group = cumsum(str_starts(lines, ">"))) %>%
  mutate(lines = ifelse(str_starts(lines, ">"), str_remove(lines, "> "), lines)) %>%
  group_by(group) %>%
  summarise(metadata = first(lines[str_starts(lines, "start")]),
            fasta_seq = paste(lines[!str_starts(lines, "start")], collapse = "")) %>%
  separate(metadata, into = c("start_label", "start", "end_label", "end", "length_label", "length"), sep = " ") %>%
  select(start, end, length, fasta_seq) %>%
  mutate(across(c(start, end, length), as.numeric))

#make row number a column (to assign a3m corresponding value)
fasta_table <- fasta_table %>%
  mutate(row_number = row_number()-1)
fasta_table <- fasta_table %>%
  mutate(row_number = as.character(row_number))

#merge fasta file (key) for deciphering rows of foldx output
merged_df <- df_filtered %>%
  inner_join(fasta_table, by = c("a3m" = "row_number"))

df_filtered %>%
  anti_join(fasta_table, by = c("a3m" = "row_number")) 

#create column containing deletion end (start + length)
merged_df$end_value <- merged_df$start + merged_df$length
merged_df$delta_delta <- merged_df$free_energy - 13.0565


#ggplot(merged_df, aes(x = start, y = end_value, color = delta_delta)) +
  geom_point(size = 3) +  # Adjust point size as needed
  scale_color_gradient(low = "blue", high = "red") +  
  labs(
       x = "deletion start",
       y = "deletion end",
       color = "delta delta G") +
  scale_x_continuous(breaks = seq(10, max(merged_df$start), by = 10)) +  
  scale_y_continuous(breaks = seq(10, max(merged_df$end_value), by = 10)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

#omit wt (./0 ...)
df_mutants <- merged_df[!grepl("./0_relaxed_rank_001_alphafold2_ptm_model_3_seed_000", merged_df$pdb), ]

#create heatmap using annette's colors
ggplot(df_mutants, aes(x = start, y = end_value, fill = delta_delta)) +
  geom_tile() +  
  scale_fill_gradientn(colors = c("#ffff00", "#a020f0", "#a020f0"),  
                       values = scales::rescale(c(-2, 19, 100), to = c(0, 1)),  
                       limits = c(-2, 100)) + 
  labs(title = "mutant ddG",
    x = "deletion start",
    y = "deletion end",
    color = "delta delta G") +
  scale_x_continuous(breaks = seq(10, max(df_mutants$start), by = 10)) +  # X-axis labels every 10 units
  scale_y_continuous(breaks = seq(10, max(df_mutants$end_value), by = 10)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# save df_mutants as a CSV file
write.csv(df_mutants, "del_insilico_stability.csv", row.names = FALSE)



