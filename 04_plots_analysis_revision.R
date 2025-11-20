# Rhythm Analysis of DoReCo languages
# Involved: Susanne Fuchs, Ludger Paschen, Lara S. Burchardt
# R Codes from: Lara S. Burchardt 

# Script  6 of 6

# Plots and Analyses for rhythm analysis
# Version after Revisions 1 from ANYAS

# Major changes: Sümi filtered out as outlier, code cleaned


###############################################################################

# 00: preparations -----

## 00a: load packages -----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse", "psych", "tidygeocoder", "countrycode", "devtools",
             "lme4", "maps", "effsize","praatpicture", "grid", "ggplotify",
             "magick", "scales", "cowplot")

## 00b: prepare themes, color palettes, etc. ----

# color blind friendly palette with 28 colors for language families

colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
            "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1f77b4", "#aec7e8", 
            "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", 
            "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", 
            "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d")

species_comp_color <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728")

my_custom_theme <- theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

# 01: load data ----- 

# for rhythm analysis results per sequence
doreco_rhythm_results_complete <- read_rds("rhythm_results_doreco_ioi_meta_complete.rds")

# for pause duration and speaking duration
ioi_data_alternative <- read_delim("unsplit_ioi_data_for_rhythm_analysis_including_meta_data_run_Oct24.csv", delim = ",")

# for iois 
ioi_data <- read_delim("rhythm_analysis_results/revision/iois_doreco_rerun_ir.csv", delim = ",")

# ir values calculated and saved in script 03_post_meta
all_results <- read_delim("doreco_interger_ratios_all_pairs_all_files.csv")

# meta data languages
meta_languages <- read_delim("doreco_languages_metadata_2_0.csv", delim = ",")

# 02: data wrangling ----

## 02a: filtering Sümi ----

# for Sümi we only have very few data points, during the revision process and due to reviewer feedback, 
# we decided to not include Sümi in the results, so we filter it out

doreco_rhythm_results_complete <- doreco_rhythm_results_complete %>% 
  filter(Language != "Sümi")

ioi_data_alternative <- ioi_data_alternative %>% 
  filter( glottocode != "sumi1235")

meta_from_complete <- doreco_rhythm_results_complete %>% 
  select(filename,Language, speaker, speaker_age, speaker_sex, synthesis, tone)

ioi_data_meta <- left_join(ioi_data, meta_from_complete,
                           by = "filename")

ioi_data_meta <- left_join(ioi_data_meta, meta_languages, by = "Language", multiple = "any")

ioi_data_meta <- ioi_data_meta %>% 
  filter(Language != "Sümi")

all_results <- all_results %>% 
  filter(Language != "Sümi")

## 02b: sort median ioi beat ----

# for visualization we also produce a dataset where languages and language family are ordered by median ioi-beat

# Calculate median ioi_beat per language
median_ioi_lan <- doreco_rhythm_results_complete %>%
  group_by(Language) %>%
  summarize(median_ioi_beat = median(ioi_beat, na.rm = TRUE))

# Calculate median ioi_beat per language family
median_ioi_fam <- doreco_rhythm_results_complete %>%
  group_by(Family) %>%
  summarize(median_ioi_beat = median(ioi_beat, na.rm = TRUE))

# copy dataframe
doreco_rhythm_results_complete_ordered <- doreco_rhythm_results_complete 

# Reorder the levels of Language based on median ioi_beat
doreco_rhythm_results_complete_ordered$Language <- factor(doreco_rhythm_results_complete_ordered$Language,
                                                          levels = median_ioi_lan$Language[order(median_ioi_lan$median_ioi_beat)])
doreco_rhythm_results_complete_ordered$Family <- factor(doreco_rhythm_results_complete_ordered$Family,
                                                        levels = median_ioi_fam$Family[order(median_ioi_fam$median_ioi_beat)])


## 02c: summarize by file ----

# after splitting files because of long silences, we summarize by file, to avoid pseudoreplication and inflation of n

### 02c_1: results ----
doreco_rhythm_results_complete_summarized_file <- doreco_rhythm_results_complete %>%
  filter(!is.na(ioi_beat)) %>%
  group_by(Family, file) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), 
            across(where(is.character), first), 
            across(where(is.factor), first))

### 02c_2: iois ----

# for plots: summarize ioi per file to get readable plots

ioi_data_meta_summarized_file <- ioi_data_meta %>% 
  group_by(filename) %>% 
  summarize(median_ioi = round(median(ioi, na.rm = TRUE), digits = 2))

ioi_data_meta_summarized_file <- left_join(ioi_data_meta_summarized_file, meta_from_complete, by = "filename")

# 03: statistics -----

## 03a: file & sequence numbers ----

# summarize to get number of files, sequences, speakers, IOIs per Language (Table 2 in Manuscript)

doreco_rhythm_results_complete %>%
  group_by(speaker) %>% 
  select('ioi_beat', 'unbiased_cv', 'npvi', 'n_elements') %>% 
  summary(na.rm = TRUE)

summary_by_language_median <- doreco_rhythm_results_complete %>% 
  group_by(Language) %>% 
  summarize(median_ioi = round(median(ioi_beat, na.rm = TRUE), digits = 2),
            median_elements = median(n_elements, na.rm = TRUE),
            speaker_nr = length(unique(speaker)),
            file_nr = length(unique(file)),
            sequence_nr = length(unique(filename)))

summary_by_speaker <- doreco_rhythm_results_complete %>% 
  group_by(speaker) %>% 
  summarize(median_ioi_speaker = round(median(ioi_beat, na.rm = TRUE), digits = 2),
            filename_nr = length(unique(filename)),
            file_nr = length(unique(file)),
            age = median(speaker_age))

summary_by_language_iois <- ioi_data_meta %>% 
  group_by(Glottocode) %>% 
  summarize(nr_io_units = length(Glottocode),
            median_ioi_dur = median(ioi, na.rm = TRUE))


## 03b: correlation & effect size ----
### 03b_1: ioi beat----


#### effect sizes ----
# (cohen's D, psych package)
# effect size gender

men <- doreco_rhythm_results_complete %>% 
  filter(speaker_sex == "m")
women <- doreco_rhythm_results_complete %>% 
  filter(speaker_sex == "f")

d_gender <- cohen.d(men$ioi_beat, women$ioi_beat)
# d = 0.11 --> negligible

t.test(men$ioi_beat, women$ioi_beat, var.equal = FALSE)
# p = 0.00007*** --> significant but negligible

# effect size tone language

toneyes <- doreco_rhythm_results_complete %>% 
  filter(tone == "yes")
toneno <- doreco_rhythm_results_complete %>% 
  filter(tone == "no")

d_tone <- cohen.d(toneyes$ioi_beat, toneno$ioi_beat)
# d = 0.08 --> negligible

t.test(toneyes$ioi_beat, toneno$ioi_beat, var.equal = FALSE)
# p = 0.015 --> significant but negligible

#### correlations ----

# correlation with median age

cor_age <- corr.test(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$speaker_age)
r_age <- cor_age$r
R_squared_age <- r_age^2

# ~ 0.003 % of variance are explained by age --> negligible

# correlation with morphological synthesis

cor_morph <- corr.test(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$synthesis)
r_morph <- cor_morph$r
R_squared_morph <- r_morph^2

# 1.14% of variance explained, so even though statistically significant, very small effect size


##03b_2:raw ioi -----

#### effect sizes ----
# (cohen's D, psych package)
# effect size gender

men_ioi <- ioi_data_meta %>% 
  filter(speaker_sex == "m") %>% 
  drop_na(ioi)
women_ioi <- ioi_data_meta %>% 
  filter(speaker_sex == "f") %>% 
  drop_na(ioi)

d_gender_ioi <- cohen.d(men_ioi$ioi, women_ioi$ioi)
# d = -0.10 --> negligible

t.test(men_ioi$ioi, women_ioi$ioi, var.equal = FALSE)
# p = 0.00000 --> significant but negligible

# effect size tone language

toneyes_ioi <- ioi_data_meta %>% 
  filter(tone == "yes") %>% 
  drop_na(ioi)
toneno_ioi <- ioi_data_meta %>% 
  filter(tone == "no") %>% 
  drop_na(ioi)

d_tone_ioi <- cohen.d(toneyes_ioi$ioi, toneno_ioi$ioi)
# d = -0.06 --> negligible

t.test(toneyes_ioi$ioi, toneno_ioi$ioi, var.equal = FALSE)
# p = 0.00000 --> significant but negligible

#### correlations ----

# correlation with median age

cor_age_ioi <- corr.test(ioi_data_meta$ioi, ioi_data_meta$speaker_age)
r_age_ioi <- cor_age_ioi$r
R_squared_age_ioi <- r_age_ioi^2

# ~ 0.022 % of variance are explained by age --> negligible

# correlation with morphological synthesis

cor_morph_ioi <- corr.test(ioi_data_meta$ioi, ioi_data_meta$synthesis)
r_morph_ioi <- cor_morph_ioi$r
R_squared_morph_ioi <- r_morph_ioi^2

# 0.3% of variance explained, so even though statistically significant, very small effect size


## 03c: beat precision --------

# beat precision values need to be transformed, to account for symmetry around 0.5 when summarizing over a sequence
# (see https://doi.org/10.31219/osf.io/69v5a_v1, equation 6 and 7)

doreco_rhythm_results_complete$ugof_ioi_sym <- pmin(doreco_rhythm_results_complete$ugof_ioi, 1 - doreco_rhythm_results_complete$ugof_ioi)

# overview of bp values

median(doreco_rhythm_results_complete$ugof_ioi_sym)
mean(doreco_rhythm_results_complete$ugof_ioi_sym)
max(doreco_rhythm_results_complete$ugof_ioi_sym)
min(doreco_rhythm_results_complete$ugof_ioi_sym)


# 04: plots ----

## 04a: map plot ----

# map plot, where are languages spoken?

doreco_rhythm_results_complete_summarized_file$Longitude<-sapply(doreco_rhythm_results_complete_summarized_file$Longitude,function(x) ifelse(x<(-25),x + 360,x))
world <- map_data('world', interior=F, wrap=c(-25,335), ylim=c(-54,79))

map_doreco <- ggplot() +
  geom_polygon(
    data=world,
    aes(x=long,y=lat,group=group),
    colour="grey",fill="#F2DDC1" ,linewidth=0.2, 
  ) + 
  geom_jitter(
    data = doreco_rhythm_results_complete_summarized_file, aes(x = Longitude, y = Latitude), color = "red", size = 2
  ) +
  #scale_x_continuous(limits = c(0, 360)) +  # Set x limits
  coord_fixed(ratio = 1) +  # Ensure aspect ratio is fixed
  #scale_fill_viridis_d(option="D") +
  # geom_label_repel(box.padding=0.5, point.padding=0.5,
  #                  data=languages, aes(Longitude, Latitude, label=Name), 
  #                  min.segment.length=unit(0, 'lines'),
  #                  size=5, max.overlaps=99) +
  scale_x_continuous(name=NULL, breaks=NULL) +
  scale_y_continuous(name=NULL, breaks=NULL) +
  theme_minimal() +
  theme(legend.position="none") 

ggsave("map_plot_fig1_part1.jpg", dpi = 300,
       width = 12,
       height = 6)

## 04b: ioi plots -----

### raw interval distribution----

hist_ioi_raw <- ioi_data_meta %>% 
  ggplot(aes(x= ioi))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 color = "white", fill = "darkblue", bins = 100)+
  theme_minimal()+
  coord_cartesian(xlim = c(0,10))+
  annotate("text", x = 5, y = 3.5, label = paste("n = ", nrow(ioi_data_meta)), size = 6)+
  ylab("Percentage [%]")+
  xlab("IOI [sec]")+
  #annotate("text", x = 0.6, y = 9, label = "n = 1535")+
  my_custom_theme
print(hist_ioi_raw)


### ioi beat plots ----
hist_ioi_beat <- doreco_rhythm_results_complete %>% 
  ggplot(aes(x= ioi_beat))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 color = "white", fill = "darkblue")+
  theme_minimal()+
  ylab("Percentage [%]")+
  xlab("IOI Beat [Hz]")+
  annotate("text", x = 1.5, y = 9, label = paste("n =", nrow(doreco_rhythm_results_complete)), size = 6)+
  my_custom_theme
print(hist_ioi_beat)

hist_cv <- doreco_rhythm_results_complete %>% 
  ggplot(aes(x= unbiased_cv))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 color = "white", fill = "darkblue", binwidth = 0.1)+
  theme_minimal()+
  ylab("Percentage [%]")+
  xlab("CV of IOI Durations per Sequence")+
  annotate("text", x = 1.0, y = 15, label = paste("n =", nrow(doreco_rhythm_results_complete)), size = 6)+
  my_custom_theme
print(hist_cv)

### ioi boxplot per language ----

median_raw_ioi_lan <- ioi_data_meta %>%
  group_by(Language) %>%
  summarize(median_ioi = median(ioi, na.rm = TRUE))

ioi_data_meta$Language <- factor(ioi_data_meta$Language,
                               levels = median_raw_ioi_lan$Language[order(median_raw_ioi_lan$median_ioi)])

# boxplot 

boxplot_languages_ioi <- ioi_data_meta %>% 
  dplyr::filter(is.na(ioi) == FALSE) %>% 
  group_by(Language) %>% 
  ggplot(aes(y = ioi, x = Language, fill = Family))+
  #geom_boxplot(outliers = FALSE)+
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.6)+
  #geom_jitter(alpha = 0.2, size = 0.5, shape = 01)+
  geom_hline(yintercept = 2.015, linetype="dotted", linewidth = 0.7, color = "black")+  
  theme_minimal()+
  scale_fill_manual(values = colors)+
  coord_cartesian(ylim = c(0,8))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, 'cm'))+
  xlab("Languages")+
  ylab("IOI [sec]")+
  my_custom_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))


### ioi density plots ----

# additional plot for revision density of raw iois 

density_ioi <- 
  ioi_data_meta %>% 
  ggplot(aes(x = ioi))+
  geom_density()+
  #geom_vline(xintercept = min_outlier, linetype="dotted", linewidth = 2 )+  # zscore > 2, corresponding ioi beat
  #geom_vline(xintercept = max_outlier, linetype="dotted", linewidth = 2)+  # zscore < -2, corresponding ioi beat
  geom_jitter(aes(y = -1, color = filename), alpha = 0.5, size = 0.5)+
  #annotate("text", x = max_outlier-0.5, y = -0.3, label = paste("n =", n_outliers_bottom), size = 5)+
  #annotate("text", x = max_outlier-0.5, y = 2, label = expression("Z-Score "<="-1.96"), size = 5)+
  #annotate("text", x = 0.5, y = -0.3, label = paste("n =",n_elements_95 ), size = 5)+
  #annotate("text", x = 1.5, y = -0.3, label = paste("n =", n_outliers_top), size = 5)+
  #annotate("text", x = 1.5, y = 2, label = expression("Z-Score ">="1.96"), size = 5)+  
  #coord_cartesian(xlim = c(-0.5, 3))+
  ylab("Density")+
  xlab("IOI [sec]")+
  my_custom_theme

## 04b: CV plots ----
### density plot ----

# addition for review density_cv
density_cv <- 
  doreco_rhythm_results_complete %>% 
  ggplot(aes(x = unbiased_cv))+
  geom_density()+
  #geom_vline(xintercept = min_outlier, linetype="dotted", linewidth = 2 )+  # zscore > 2, corresponding ioi beat
  #geom_vline(xintercept = max_outlier, linetype="dotted", linewidth = 2)+  # zscore < -2, corresponding ioi beat
  geom_jitter(aes(y = -1, color = Language), alpha = 0.4, size = 0.2)+
  #annotate("text", x = max_outlier-0.5, y = -0.3, label = paste("n =", n_outliers_bottom), size = 5)+
  #annotate("text", x = max_outlier-0.5, y = 2, label = expression("Z-Score "<="-1.96"), size = 5)+
  #annotate("text", x = 0.5, y = -0.3, label = paste("n =",n_elements_95 ), size = 5)+
  #annotate("text", x = 1.5, y = -0.3, label = paste("n =", n_outliers_top), size = 5)+
  #annotate("text", x = 1.5, y = 2, label = expression("Z-Score ">="1.96"), size = 5)+  
  #coord_cartesian(xlim = c(-0.5, 3))+
  ylab("Density")+
  xlab("Coefficient of Variation")+
  my_custom_theme

ggsave("plot_density_unbiasedCV.jpg", dpi = 300,
       width = 12,
       height = 4)

## 04c: ioi beat plots ----

### ioi beat density ----
density_ioibeat <- 
  doreco_rhythm_results_complete %>% 
  ggplot(aes(x = ioi_beat))+
  geom_density()+
  #geom_vline(xintercept = min_outlier, linetype="dotted", linewidth = 2 )+  # zscore > 2, corresponding ioi beat
  #geom_vline(xintercept = max_outlier, linetype="dotted", linewidth = 2)+  # zscore < -2, corresponding ioi beat
  geom_jitter(aes(y = -1, color = Language), alpha = 0.5, size = 0.5)+
  #annotate("text", x = max_outlier-0.5, y = -0.3, label = paste("n =", n_outliers_bottom), size = 5)+
  #annotate("text", x = max_outlier-0.5, y = 2, label = expression("Z-Score "<="-1.96"), size = 5)+
  #annotate("text", x = 0.5, y = -0.3, label = paste("n =",n_elements_95 ), size = 5)+
  #annotate("text", x = 1.5, y = -0.3, label = paste("n =", n_outliers_top), size = 5)+
  #annotate("text", x = 1.5, y = 2, label = expression("Z-Score ">="1.96"), size = 5)+  
  #coord_cartesian(xlim = c(-0.5, 3))+
  ylab("Density")+
  xlab("IOI Beat [Hz]")+
  my_custom_theme




## what are the outliers? is there a trend here? age? 
doreco_rhythm_results_complete %>% 
  #ggplot(aes(x= unbiased_cv, y = z_scores$z_scores))+
  ggplot(aes(x= ugof_ioi, y = npvi))+
  geom_point()#+
geom_smooth()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= speaker_age, y = avg_height))+
  geom_point()


## 04d: ioi beat language/ language family plots ----

#languages_ioi

boxplot_languages <- doreco_rhythm_results_complete_ordered %>% 
  dplyr::filter(is.na(ioi_beat) == FALSE) %>% 
  group_by(speaker, Language) %>% 
  ggplot(aes(y = ioi_beat, x = Language, fill = Family))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.2, size = 0.7)+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  coord_cartesian(ylim = c(0,1.75))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, 'cm'))+
  xlab("Languages")+
  ylab("IOI Beat [Hz]")+
  my_custom_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

# boxplot languages cv (added for revision)

boxplot_languages_cv <- doreco_rhythm_results_complete_ordered %>% 
  dplyr::filter(is.na(unbiased_cv) == FALSE) %>% 
  group_by(speaker, Language) %>% 
  ggplot(aes(y = unbiased_cv, x = Language, fill = Family))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.2, size = 0.7)+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  coord_cartesian(ylim = c(0,1.75))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, 'cm'))+
  xlab("Languages")+
  ylab("unbiased CV")+
  my_custom_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

# not in publication, sorted by Language Family, note "Isolates" are not included
doreco_rhythm_results_complete_ordered %>%
  dplyr::filter(is.na(ioi_beat) == FALSE) %>%
  dplyr::filter(Family != "Isolate") %>% 
  #group_by(speaker,Family) %>%
  group_by(Family) %>%
  ggplot(aes(y = ioi_beat, x = Family))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.2)+
  theme_minimal()+
  theme(#legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, 'cm'))+
  xlab("Language Family")+
  ylab("IOI Beat [Hz]")


## 04e: sex differences ----

# original submission - ioi beat
box_gender <- doreco_rhythm_results_complete_summarized_file %>%
  drop_na(speaker_sex) %>%
  group_by(speaker) %>% 
  ggplot(aes(y = ioi_beat, x = speaker_sex))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.3, size = 0.8)+
  theme_minimal()+
  xlab("Gender")+
  ylab("IOI Beat [Hz]")+
  my_custom_theme+
  scale_x_discrete(labels = c("f" = "Female", "m" = "Male"))+
  annotate("text", x = 1.5, y = 0.8, label = paste("Cohen's d = ", round(d_gender$estimate, 2)), size = 5)

doreco_rhythm_results_complete_summarized_file %>%
  drop_na(speaker_sex) %>%
  group_by(speaker) %>%  
  ggplot(aes(y = unbiased_cv, x = speaker_sex))+
  geom_boxplot(outliers = FALSE)+
  coord_cartesian(ylim = c(0.15,1.5))+
  geom_jitter(alpha = 0.3, size = 0.8)+
  theme_minimal()+
  xlab("Gender")+
  ylab("Coefficient of Variation")+
  my_custom_theme+
  scale_x_discrete(labels = c("f" = "Female", "m" = "Male"))

doreco_rhythm_results_complete_summarized_file %>%
  drop_na(speaker_sex) %>%
  group_by(speaker) %>% 
  ggplot(aes(y = npvi, x = speaker_sex))+
  geom_boxplot(outliers = FALSE)+
  #coord_cartesian(ylim = c(0,1))+
  geom_jitter(alpha = 0.3, size = 0.8)+
  theme_minimal()+
  xlab("Gender")+
  ylab("nPVI")+
  my_custom_theme+
  scale_x_discrete(labels = c("f" = "Female", "m" = "Male"))

#### revision - raw ioi ----

box_gender_ioi <- ioi_data_meta_summarized_file %>%
  drop_na(speaker_sex) %>%
  #group_by(speaker) %>% 
  ggplot(aes(y = median_ioi, x = speaker_sex))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.3, size = 0.6, shape = 01)+
  theme_minimal()+
  xlab("Gender")+
  ylab("IOI [sec]")+
  my_custom_theme+
  scale_x_discrete(labels = c("f" = "Female", "m" = "Male"))+
  annotate("text", x = 1.5, y = 5, label = paste("Cohen's d = ", round(d_gender_ioi$estimate, 2)), size = 5)


## 04f age differences -----

# medium ioi beat per language vs. medium age per language
median_data <- aggregate(cbind(speaker_age, ioi_beat) ~ Language, data = doreco_rhythm_results_complete_summarized_file, median)

# Create a scatter plot with median ioi_beat against median age
### LP: I would remove this plot, don't see the need to average per language
###LSB: we discussed this, to avoid seeing a bias of differences in median age spoken per language

ggplot(median_data, aes(x = speaker_age, y = ioi_beat)) +
  geom_point() +
  geom_smooth(method = "lm", color = "#0072B2")+
  labs(x = "Median Age per Language",
       y = "Median IOI Beat [Hz] per Language") +
  theme_minimal()+
  my_custom_theme


# age against ioi beat directly (potentially biased because of differences in age structure per language)
scatter_age <- doreco_rhythm_results_complete_summarized_file %>% 
  ggplot(aes(x= speaker_age, y = ioi_beat ))+
  geom_jitter(alpha = 0.8, size = 0.8)+
  geom_smooth(method = "lm", color = "#0072B2")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab('Speaker Age ')+
  ylab("IOI beat [Hz]")+
  my_custom_theme+
  annotate("text", x = 80, y = 0.8, label = paste("R² =", round(R_squared_age, 3)), size = 5)

#### revision - raw ioi ----

scatter_age_ioi <- ioi_data_meta_summarized_file %>% 
  ggplot(aes(x= speaker_age, y = median_ioi ))+
  geom_jitter(alpha = 0.3, size = 0.5)+
  geom_smooth(method = "lm", color = "#0072B2")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab('Speaker Age ')+
  ylab("IOI [sec]")+
  my_custom_theme+
  annotate("text", x = 80, y = 6, label = paste("R² =", round(R_squared_age_ioi, 3)), size = 5)

## 04g: tone and morphological complexity ----

box_tone <- doreco_rhythm_results_complete_summarized_file %>%
  filter(is.na(tone) == FALSE) %>%
  group_by(speaker) %>% 
  ggplot(aes(x= tone, y = ioi_beat ))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(size= 0.5, alpha = 0.5)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab(' Tone Language')+
  ylab("IOI beat [Hz]")+
  my_custom_theme+
  scale_x_discrete(labels = c("no" = "No", "yes" = "Yes"))+
  annotate("text", x = 1.5, y = 0.8, label = paste("Cohen's d = ", round(d_tone$estimate, 2)), size = 5)


###SF: Can you explain how morphological synsthesis is calculated. For each speaker here?
### This means the lower the morphological complexity the higher the IOI beat? 

cor_morph <- cor(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$synthesis, use = "complete.obs")

scatter_morph <- doreco_rhythm_results_complete_summarized_file %>%
  group_by(speaker) %>% 
  ggplot(aes(x= synthesis, y = ioi_beat ))+
  geom_point()+
  geom_smooth(method = "lm", color = "#0072B2")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab('Morphological Synthesis')+
  ylab("IOI beat [Hz]")+
  my_custom_theme+
  annotate("text", x = 3, y = 0.8, label = paste("R² =", round(R_squared_morph, 2)), size = 5)

#### revision- war ioi -----

box_tone_ioi <- ioi_data_meta_summarized_file %>%
  filter(is.na(tone) == FALSE) %>%
  group_by(speaker) %>% 
  ggplot(aes(x= tone, y = median_ioi ))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(size= 0.5, alpha = 0.3, shape = 01)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab(' Tone Language')+
  ylab("IOI [sec]")+
  my_custom_theme+
  scale_x_discrete(labels = c("no" = "No", "yes" = "Yes"))+
  annotate("text", x = 1.5, y = 5, label = paste("Cohen's d = ", round(d_tone_ioi$estimate, 2)), size = 5)

# Synthesis

scatter_morph_ioi <- ioi_data_meta_summarized_file %>%
   #group_by(speaker) %>%
   #summarize(
  #   synthesis = mean(synthesis, na.rm = TRUE),
  #   median_ioi = mean(median_ioi, na.rm = TRUE)
  # ) %>%
  ggplot(aes(x= synthesis, y = median_ioi))+
  geom_jitter(alpha = 0.3, size = 0.5)+
  geom_smooth(method = "lm", color = "#0072B2")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab('Morphological Synthesis')+
  ylab("IOI [sec]")+
  my_custom_theme+
  #annotate("text", x = 3, y = 3.5, label = paste("R² =", round(R_squared_morph_ioi_speaker, 2)), size = 5)
  annotate("text", x = 3, y = 6, label = paste("R² =", round(R_squared_morph_ioi, 2)), size = 5)


# ioi_per_speaker <- ioi_data_meta_summarized_file %>%
#   group_by(speaker) %>%
#   summarize(
#     synthesis = mean(synthesis, na.rm = TRUE),
#     median_ioi = mean(median_ioi, na.rm = TRUE)
#   )
# 
# cor_morph_ioi_per_speaker <- corr.test(ioi_per_speaker$median_ioi, ioi_per_speaker$synthesis)
# r_morph_ioi_speaker <- cor_morph_ioi_per_speaker$r
# R_squared_morph_ioi_speaker <- r_morph_ioi_speaker^2
## R² = 0.05 (without summarizing: R² = 0.003, see line 268)


## 04h: different continents ----

doreco_rhythm_results_complete_summarized_file %>% 
  ggplot(aes(x= Area, y = ioi_beat ))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(size= 0.5, alpha = 0.5)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab(' Area')+
  ylab("IOI beat [Hz]")+
  my_custom_theme

## 04i: height information -----

doreco_rhythm_results_complete_summarized_file %>% 
  ggplot(aes(x= avg_height, y = ioi_beat ))+
  geom_jitter(alpha = 0.5, aes(color = Area))+
  geom_smooth(method = "lm")+
  theme_minimal()+
  xlab('Average Height per Country [cm]')+
  ylab("IOI beat [Hz]")+
  my_custom_theme+
  scale_color_manual(values = colors)

## 04j: n element information ----

hist_n_element <- doreco_rhythm_results_complete %>% 
  ggplot(aes(x= n_elements))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 color = "white", fill = "darkblue", 
                 boundary = 4)+
  theme_minimal()+
  #coord_cartesian(xlim = c(4,270))+
  ylab("Percentage [%]")+
  xlab("n Elements per Sequence")+
  annotate("text", x = 100, y = 20, label = paste("n =", nrow(doreco_rhythm_results_complete)), size = 6)+
  my_custom_theme
print(hist_n_element)


# n elements vs. cv 

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= n_elements, y = unbiased_cv))+
  geom_point(alpha = 0.5, shape = 21)+
  #geom_smooth(method = "lm")+
  theme_minimal()+
  #coord_cartesian(xlim = c(4,270))+
  ylab("CV")+
  xlab("n Elements per Sequence")+
  #annotate("text", x = 100, y = 20, label = paste("n =", nrow(doreco_rhythm_results_complete)), size = 6)+
  my_custom_theme

cor(doreco_rhythm_results_complete$n_elements, doreco_rhythm_results_complete$unbiased_cv)
# r = 0.02 --> negligible, but higher variability for shorter sequences 

## 04k: Annotation Plot from Praat ------

png("annotation_plot_fig1_part2.png", width = 12, height = 6, units = "in",  res = 300)

annotation_plot <- praatpicture(
  'BeAM_199X_HumanInLandOfDeath_flk_fragment.wav',
  start = 0, 
  end = 3.3,
  frames = c('sound', 'spectrogram', 'TextGrid'),
  proportion = c(30, 40, 30),
  spec_freqRange = c(0, 16000),
  
  tg_tiers = c('text-dolgan', 'text-english', 'word-dolgan', 'word-english'),
  tg_color = c('blue', 'black', 'black', 'black'),
  tg_focusTier = 'all',
  tg_focusTierLineType = c('solid', 'solid', 'dashed', 'dashed'),
  tg_focusTierColor = c('black', 'black', 'black', 'black'),
  
  draw_arrow = c('spectrogram', 0.478, 8500, 2.321, 8500, code = 3,
                 length = 0.15, 
                 angle = 20, col = 'blue', lwd = 2, lty = 'solid'),
  annotate = list(
    c('spectrogram', 1.3, 7000, col = 'blue', font = 2, cex = 1.8, labels = 'Inter-Onset-Interval (IOI)')#,
    #c('spectrogram', 2.0, 5000, col = 'red', font = 2, cex = 1.8, labels = 'Second Annotation')
  )
)

dev.off()

##04l: pause & ipu histograms for supps (after revision)----

## sprachdauer/ipu

mean(ioi_data_alternative$sprach_dauer, na.rm = TRUE)
#1.43 sec

hist_data <- ggplot_build(
  ggplot(ioi_data_alternative, aes(x = sprach_dauer)) +
    geom_histogram(binwidth = 0.1, na.rm = TRUE)
)$data[[1]]

# Compute percentage
hist_data <- hist_data %>%
  mutate(percentage = round(count / sum(count) * 100, digits = 2))

# Plot with precomputed values
hist_sprachdauer <- ggplot(hist_data, aes(x = xmin, y = percentage)) +
  geom_bar(
    stat = "identity",
    aes(),#text = paste0("Percentage: ", round(percentage, 2), "%")),
    width = diff(hist_data$x)[1],  # Ensures correct bar width
    color = "white",
    fill = "darkblue"
  ) +
  coord_cartesian(ylim= c(0,10), xlim = c(0,7.5))+
  xlab("IPU [sec]") +
  ylab("Percentage [%]") +
  theme_minimal()

hist_sprachdauer

## pause duration

mean(ioi_data_alternative$pause_duration, na.rm = TRUE)
# 0.83 sec

hist_data <- ggplot_build(
  ggplot(ioi_data_alternative, aes(x = pause_duration)) +
    geom_histogram(binwidth = 0.1, na.rm = TRUE)
)$data[[1]]

# Compute percentage
hist_data <- hist_data %>%
  mutate(percentage = round(count / sum(count) * 100, digits = 2))

# Plot with precomputed values
hist_pauseduration <- ggplot(hist_data, aes(x = xmin, y = percentage)) +
  geom_bar(
    stat = "identity",
    aes(),#text = paste0("Percentage: ", round(percentage, 2), "%")),
    width = diff(hist_data$x)[1],  # Ensures correct bar width
    color = "white",
    fill = "darkblue"
  ) +
  geom_vline(xintercept = 3.36, linetype="dotted", linewidth = 1, color = "darkred" )+
  annotate("text", x = 5, y = 4, label = "Cutoff Criterion, 99th Quantile", size = 4, color = "darkred")+
  xlab("Pause Duration [sec]") +
  ylab("Percentage [%]") +
  coord_cartesian(ylim= c(0,10), xlim = c(0,7.5))+
  theme_minimal()

hist_pauseduration

## for reviewer: pause duration per language
## the threshold to split a sequence is 3.36 seconds of pause duration
## for all but one langugage such long pause durations are only found in the outliers, 
## therefore we deam this approach suitable to not split languages with very different 
## pause properties (i.e. durations)

boxplot_languages_pause <-ioi_data_alternative %>% 
  #dplyr::filter(is.na(ugof_ioi_sym) == FALSE) %>% 
  group_by(speaker, glottocode) %>% 
  ggplot(aes(y = pause_duration, x = glottocode))+ #, fill = Family))+
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.4, outlier.color = "darkred")+
  #geom_jitter(alpha = 0.2, size = 0.7)+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  geom_hline(yintercept = 3.36, linetype="dotted", linewidth = 1, color = "black" )+  
  #coord_cartesian(ylim = c(0,0.55))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, 'cm'))+
  xlab("Languages")+
  ylab("Pause Duration [sec]")+
  my_custom_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

## plot ioi beat vs cv


ggplot(doreco_rhythm_results_complete, aes(x = ioi_beat  , y = unbiased_cv)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm", color = "#0072B2")+
  #labs(x = "Median Age per Language",
  #     y = "Median IOI Beat [Hz] per Language") +
  theme_minimal()+
  my_custom_theme


## 04m beat precision plots ------

# histogram

hist_bp <- doreco_rhythm_results_complete %>% 
  ggplot(aes(x= ugof_ioi_sym))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 color = "white", fill = "darkblue", bins = 25)+
  
  theme_minimal()+
  #coord_cartesian(xlim = c(0,10))+
  #annotate("text", x = 5, y = 5, label = paste("n =", nrow(ioi_data_alternative)), size = 6)+
  ylab("Percentage [%]")+
  xlab("Beat Precision Values")+
  #annotate("text", x = 0.6, y = 9, label = "n = 1535")+
  my_custom_theme

# if we would want to add a vline at the median:
# geom_vline(aes(xintercept = median(ugof_ioi_sym)), 
# color = "red", linetype = "dashed", size = 1) +


# order languages by median beat precision
median_raw_bp_lan <- doreco_rhythm_results_complete %>%
  group_by(Language) %>%
  summarize(median_bp = median(ugof_ioi_sym, na.rm = TRUE))

doreco_rhythm_results_complete_bp_order <- doreco_rhythm_results_complete

doreco_rhythm_results_complete_bp_order$Language <- factor(doreco_rhythm_results_complete$Language,
                                                           levels = median_raw_bp_lan$Language[order(median_raw_bp_lan$median_bp)])

# boxplot comparison langugaes
boxplot_languages_bp <-doreco_rhythm_results_complete_bp_order %>% 
  dplyr::filter(is.na(ugof_ioi_sym) == FALSE) %>% 
  group_by(speaker, Language) %>% 
  ggplot(aes(y = ugof_ioi_sym, x = Language, fill = Family))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.2, size = 0.7)+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  coord_cartesian(ylim = c(0,0.55))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, 'cm'))+
  xlab("Languages")+
  ylab("Beat Precision Value")+
  my_custom_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))


# 05: revision - additional analysis ----

## 05a: analyse integer ratios -----

mean(all_results$ratio, na.rm = TRUE)
min(all_results$ratio, na.rm = TRUE)
max(all_results$ratio, na.rm = TRUE)
median(all_results$ratio, na.rm = TRUE)
sd(all_results$ratio, na.rm = TRUE)

# summarizing integer ratios per file
summary_ir_per_file <- all_results %>%
  group_by(filename) %>%
  summarise(mean_ratio = mean(ratio, na.rm = TRUE),
            median_ratio = median(ratio, na.rm = TRUE),
            min_ratio = min(ratio, na.rm = TRUE),
            max_ratio = max(ratio, na.rm = TRUE),
            n_intervals = n())


sequence_random <- runif(1000, min = 0.098, max = 12.274)
sequence_random <- as.data.frame(sequence_random)
random_ir <- expand.grid(i = 1:1000, j = 1:1000) %>%
  filter(i != j) %>%  # Entferne Fälle, wo i == j
  mutate(ratio = sequence_random$sequence_random[i] / (sequence_random$sequence_random[i] + sequence_random$sequence_random[j]))

# second version random seq

sequence_random_2 <- rnorm(1000, mean = 2.02, sd = 1.22)
# sd of ioi_data_meta$ioi, na.rm = TRUE
sequence_random_2 <- as.data.frame(sequence_random_2)
random_ir_2 <- expand.grid(i = 1:1000, j = 1:1000) %>%
  filter(i != j) %>%  # Entferne Fälle, wo i == j
  mutate(ratio = sequence_random_2$sequence_random_2[i] / (sequence_random_2$sequence_random_2[i] + sequence_random_2$sequence_random_2[j]))

# compare random seq to asimjeeg datoga (language looking most similar to "random distribution")
# and to cabecar (highest density peak)
# calculate 
## end second version random seq

density_plot_ir_raw <- all_results %>% 
  ggplot(aes(x = ratio))+
  geom_density(aes(color = Language))+
  #geom_density()+
  #geom_vline(xintercept = 1/2.25, linetype="dotted", linewidth = 2 )+  
  #geom_vline(xintercept = 0.56, linetype="dotted", linewidth = 2)+
  #geom_density(data = random_ir_2, aes(x = ratio),linewidth = 1.3, color = "black")+
  #coord_cartesian(xlim = c(0, 2))+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("Integer Ratios")+
  ylab(" Density")

count_below <- sum(all_results$ratio < 0.44, na.rm = TRUE)
count_above <- sum(all_results$ratio > 0.56, na.rm = TRUE)
total_count <- count_below + count_above
percentage <- 100- (total_count / nrow(all_results)) * 100

# which languages have lowest/highest density peaks?
# Build ggplot object to extract data
density_data_ir <- ggplot_build(density_plot_ir_raw)$data[[1]]

# Check available columns
print(colnames(density_data_ir))  

# Map group numbers back to languages
mapping <- all_results %>% 
  select(Language) %>% 
  distinct() %>% 
  mutate(group = as.numeric(as.factor(Language)))  # Ensure numeric group IDs match

# Find max density per language
peaks_density_ir <- density_data_ir %>%
  group_by(group) %>%
  summarise(max_density = max(y)) %>%
  arrange(max_density) %>%
  #head(3) %>%
  left_join(mapping, by = "group")  # Merge with language names

# new idea- density plot of ir per language but for iois 

density_plot_ioi_language_raw <- ioi_data_meta %>% 
  ggplot(aes(x = ioi))+
  geom_density(aes(color = Family))+
  #geom_density(data =sequence_random_2, aes(x = sequence_random_2), color = "black", size = 1.3)+
  #coord_cartesian(xlim = c(0, 2))+
  scale_color_manual(values = colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("IOI [sec]")+
  ylab(" Density")

# get information on peaks

# Build ggplot object to extract data
density_data_ioi <- ggplot_build(density_plot_ioi_language_raw)$data[[1]]

# Check available columns (useful for debugging)
print(colnames(density_data_ioi))  

# Map group numbers back to languages
mapping <- all_results %>% 
  select(Language) %>% 
  distinct() %>% 
  mutate(group = as.numeric(as.factor(Language)))  # Ensure numeric group IDs match

# Find max density per language
peaks_density_ioi <- density_data_ioi %>%
  group_by(group) %>%
  summarise(max_density = max(y),
            median= median(y)) %>%
  arrange(max_density) %>%
  #head(3) %>%
  left_join(mapping, by = "group")  # Merge with language names

# version 2, including x values for peak
peaks_density_ioi <- density_data_ioi %>%
  group_by(group) %>%
  slice_max(order_by = y, n = 1, with_ties = FALSE) %>%  # keep the x and y of the max density
  mutate(median_y = median(y)) %>%                       # optional: add median of y
  rename(max_density = y, max_density_x = x) %>%         # rename for clarity
  left_join(mapping, by = "group")


# 06: comparison to animal data ------

## 06a: load animal data ----
skylark_ioi <- read_delim("Burchardt_Briefer_Knoernschild_2021_skylarkSong_timepoints.csv", delim = ",")
skylark_ioi <- skylark_ioi %>% 
  mutate(ioi = `IOI[sec]`) %>% 
  filter(ioi > 0) %>% 
  select(ioi, ID)
ioi_cutoff_skylark <- quantile(skylark_ioi$ioi, .99)
skylark_ioi <- skylark_ioi %>% 
  filter(ioi <= ioi_cutoff_skylark)

brown_meagre_ioi <- read_delim("iois_BM_all.csv", delim = ",")
# some iois are 0.000 in a fashion that indicates values were doubled in calculation table
# therefore values of 0.000 are filtered out
brown_meagre_ioi <- brown_meagre_ioi %>% 
  filter(`unlist(ioi_all)` > 0)

#zebrafinch_ioi <- read_delim("Complete_ioi_Data_With_Meta_Information.csv", delim = ",")

capefur_seal <- read_delim("cape_fur_seal_ioi.csv", delim = ";")

frog_data <- read_delim("frog_results_appendix_2.csv", delim = ";")

##06b: rescale for comp -----

rescale_comparison <- data.frame()

rescale_humam <- rescale(ioi_data_meta$ioi)
rescale_skylark <- rescale(skylark_ioi$ioi)
rescale_BM <- rescale(brown_meagre_ioi$`unlist(ioi_all)`)
#rescale_ZF <- rescale(zebrafinch_ioi$ioi)
rescale_CFS <- rescale(capefur_seal$ioi)
rescale_random <- rescale(sequence_random_2$sequence_random_2)

ioi_comparison_rescaled <- bind_rows(
  data.frame(value = rescale_humam, source = "Human"),
  data.frame(value = rescale_skylark, source = "Skylark"),
  data.frame(value = rescale_BM, source = "Brown Meagre"),
  #data.frame(value = rescale_ZF, source = "Zebrafinch"),
  data.frame(value = rescale_CFS, source = "Cape Fur Seal"),
  data.frame(value = rescale_random, source = "Human Random")
)

## 06c: plot comp -----

#ioi_comparison_rescaled$source <- factor(ioi_comparison_rescaled$source, levels = c("Brown Meagre", "Cape Fur Seal", "Human", "Skylark"))

hist_ioi_rescaled_comparison_2 <- ioi_comparison_rescaled %>%
  #filter(source != "Skylark") %>% 
  ggplot(aes(x = value, fill = source)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = species_comp_color_2)+
  theme_minimal() +
  labs(
    x = "Rescaled IOI",
    y = "Density",
    fill = "Species",
    #title = "Density distribution of rescaled IOI across sources"
  ) +
  my_custom_theme+
  theme(
    text = element_text(size = 14),
    legend.position = "none"
  )

# run once, with legend.position in above plot "bottom" then comment out and change legend.position to "none"
#legend_comp <- get_legend(hist_ioi_rescaled_comparison_2)

legend_comp <- as.ggplot(legend_comp)

#ggsave("ioi_comparison_animal_rescaled.jpg", dpi = 300,
#       width = 12,
#       height = 8)

### comparison 2, beat precision ----

# Beat Precision positions and species labels
dots_bp <- data.frame(
  x = c(0.37, 0.38, 0.32, 0.48, 0.48),
  y = 0,
  species = factor(c("Human", "Brown Meagre", "Skylark", "Cape Fur Seal", "Human Random"), levels = c("Brown Meagre", "Cape Fur Seal", "Human", "Skylark", "Human Random"))
)


# line definition 
#line_bp <- data.frame(x = c(0, 0.5), y = 0)
#dots_df$label_y <- dots_df$y + rep(c(0.02, -0.02), length.out = nrow(dots_df))
#dots_df$vjust_pos <- rep(c(0, 1), length.out = nrow(dots_df))
library("ggrepel")

beat_precision_comp_2 <- ggplot() +
  geom_segment(aes(x = 0, xend = 0.5, y = 0, yend = 0), linewidth = 1) +
  geom_segment(aes(x = 0, xend = 0, y = -0.01, yend = 0.01), linewidth = 1) +
  geom_segment(aes(x = 0.5, xend = 0.5, y = -0.01, yend = 0.01), linewidth = 1) +
  geom_point(data = dots_bp, 
             aes(x = x, y = y, colour = species),
             size = 4) +
  geom_text(aes(x = 0, y = -0.03, label = "0"), size = 4, vjust = 1) +
  geom_text(aes(x = 0.5, y = -0.03, label = "0.5"), size = 4, vjust = 1) +
  geom_text_repel(
    data = dots_bp,
    aes(x = x, y = y, label = x, colour = species),
    nudge_y = 0.05,       # pushes labels slightly away from line
    size = 5,
    direction = "y",      # only move vertically (keeps x fixed)
    force = 8,
    min.segment.length = 0.05
  )+
  #geom_text(aes(x = 0.14, y = 0.1, label = "Average Beat Precision Across Species"), size = 4, vjust = 1) +
  scale_color_manual(values = c( "#1f77b4","#aec7e8", "#2ca02c","#ff7f0e","#d62728" )) +
  coord_fixed(ratio = 1/1.2) +
  scale_y_continuous(limits = c(-0.05, 0.15))+
  #theme_minimal(base_size = 14) +
  xlab("Average Beat Precision Across Species")+
  my_custom_theme+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",     # place legend below
    legend.box = "horizontal"
  ) +
  scale_x_continuous(limits = c(-0.05, 0.55))

# same figure for CV, but without end marker, because cv could be higher 
dots_cv <- data.frame(
  x = c(0.5, 0.18, 1.3, 0.23, 0.6),
  y = 0,
  species = factor(c("Human", "Brown Meagre", "Skylark", "Cape Fur Seal", "Human Random"), levels = c("Brown Meagre", "Cape Fur Seal", "Human", "Skylark", "Human Random"))
)
#line_cv <- data.frame(x = c(0, 1.5), y = 0)

cv_comp_2 <- ggplot() +
  geom_segment(aes(x = 0, xend = 1.5, y = 0, yend = 0), linewidth = 1) +
  geom_segment(aes(x = 0, xend = 0, y = -0.01, yend = 0.01), linewidth = 1) +
  #geom_segment(aes(x = 0.5, xend = 0.5, y = -0.01, yend = 0.01), linewidth = 1) +
  geom_point(data = dots_cv, 
             aes(x = x, y = y, colour = species),
             size = 4) +
  geom_text(aes(x = 0, y = -0.03, label = "0"), size = 4, vjust = 1) +
  #geom_text(aes(x = 0.5, y = -0.03, label = "0.5"), size = 4, vjust = 1) +
  geom_text_repel(
    data = dots_cv,
    aes(x = x, y = y, label = x, colour = species),
    nudge_y = 0.05,       # pushes labels slightly away from line
    size = 5,
    direction = "y",      # only move vertically (keeps x fixed)
    force = 5,
    min.segment.length = 0.02
  )+
  #geom_text(aes(x = 0.14, y = 0.1, label = "Average Beat Precision Across Species"), size = 4, vjust = 1) +
  scale_color_manual(values = c( "#1f77b4","#ff7f0e", "#2ca02c", "#aec7e8","#d62728" )) +
  #coord_fixed(ratio = 1/1.2) +
  scale_y_continuous(limits = c(-0.05, 0.15))+
  #theme_minimal(base_size = 14) +
  xlab("Average Coefficient of Variation Across Species")+
  my_custom_theme+
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",     
    legend.box = "horizontal"
  ) +
  scale_x_continuous(limits = c(-0.05, 1.55))


# plot grid, including legend
part_1_comp <- cowplot::plot_grid(hist_ioi_rescaled_comparison_2,cv_comp_2, beat_precision_comp_2, nrow= 3, labels = c("A", "B", "C"), rel_heights = c(0.5, 0.2, 0.3))

cowplot::plot_grid(part_1_comp, legend_comp, nrow = 2, align = "h", rel_heights = c(0.8, 0.2))

ggsave("ioi_comparison_animal_rescaled_bp_cv_random.jpg", dpi = 300,
       width = 10,
       height = 9)

## 06d: calculations ----

bm_mean <- mean(brown_meagre_ioi$`unlist(ioi_all)`, na.rm = TRUE) 
# 2.64 seconds
median(brown_meagre_ioi$`unlist(ioi_all)`, na.rm = TRUE)
# 2.48 seconds

sd(brown_meagre_ioi$`unlist(ioi_all)`, na.rm = TRUE)/mean(brown_meagre_ioi$`unlist(ioi_all)`, na.rm = TRUE) 

# frog cv

ESF_data_frog <- frog_data %>% 
  filter(species == "Eastern sedge frog") %>% 
  
WSF_data_frog <- frog_data %>% 
  filter(species == "Wallum sedge frog")
  
mean_cv_wsf <- mean(WSF_data_frog$`Coefficient of variation`, na.rm = TRUE)
median_cv_wsf <- median(WSF_data_frog$`Coefficient of variation`, na.rm = TRUE)

mean_cv_esf <- mean(ESF_data_frog$`Coefficient of variation`, na.rm = TRUE)
median_cv_esf <- median(ESF_data_frog$`Coefficient of variation`, na.rm = TRUE)


#skylark mean ioi

mean(skylark_ioi$`IOI[sec]`, na.rm = TRUE)
sd(skylark_ioi$`IOI[sec]`, na.rm = TRUE)/mean(skylark_ioi$`IOI[sec]`, na.rm = TRUE)

# 07: plot grids -----

## 07a: Figure 1 ----

# part 3 of Figure 1:  Histograms IOI distribution [sec]  

hist_plots <- cowplot::plot_grid(hist_ioi_raw, hist_cv,hist_n_element,
                                 labels = c("D", "E", "F"), ncol = 3)

ggsave("hist_plot_fig1_part3_v3.jpg", dpi = 300, hist_plots,
       width = 12,
       height = 6)

#part 1: annotations and part 2: map are generated and saved further up, figures where
# combined outside of R

## 07b: Figure 2-----

# Density Plot & Language Boxplots

cowplot::plot_grid(density_ioibeat, boxplot_languages,
                   nrow = 2,
                   rel_heights = c(0.4, 0.6),
                   labels = c("A", "B"))

ggsave("sup_figure_ioibeats.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

## 07c:old Figure 3 - now S5 ----

# original submission ioi beat
#cowplot::plot_grid(scatter_age, scatter_morph,
#                   box_gender, box_tone, ncol = 2,
#                   labels = c("A", "B", "C", "D"))


#ggsave("manuscript_figure3.jpg", dpi = 300,
#       width = 22,
#       height = 20,
#       units = "cm")

# revision: raw ioi 

cowplot::plot_grid(scatter_age_ioi, scatter_morph_ioi,
                   box_gender_ioi, box_tone_ioi, ncol = 2,
                   labels = c("A", "B", "C", "D"))

ggsave("supplements_figure_5_ioi.jpg", dpi = 300,
       width = 22,
       height = 20,
       units = "cm")

## 07d: Supp Figure 1 -----
# plot of IPU durations and pause durations histograms

cowplot::plot_grid(hist_sprachdauer, hist_pauseduration, ncol = 2,
                   labels = c("A", "B"))

ggsave("supplements_figure_ipu_pause_duration.jpg", dpi = 300,
       width= 12,
       height = 6)

## 07e: Supp Figure 2 ----
# main plot 2, but for cv values

cowplot::plot_grid(density_cv, boxplot_languages_cv,nrow = 2,
                   rel_heights = c(0.4, 0.6),
                   labels = c("A", "B"))

ggsave("supplements_figure_CV.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

## 07f: Figure 2 but raw IOI ----

cowplot::plot_grid(density_plot_ioi_language_raw, boxplot_languages_ioi,
                   nrow = 2,
                   rel_heights = c(0.4, 0.6),
                   labels = c("A", "B"))

ggsave("alt_figure_2_raw_ioi.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

##07g: beat precision supp figure ------

cowplot::plot_grid(hist_bp, boxplot_languages_bp,
                   nrow = 2,
                   rel_heights = c(0.3, 0.7),
                   labels = c("A", "B"))

ggsave("sup_figure_beatprecision.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

## 07h: density plots IOI & IR with random ----

cowplot::plot_grid(density_plot_ioi_language_raw, 
                   density_plot_ir_raw,
                   ncol = 2, labels = c("A", "B"))

ggsave("density_plots_per_language.jpg", dpi = 300,
       width = 14,
       height = 8)

