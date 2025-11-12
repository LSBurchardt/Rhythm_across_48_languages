# Rhythm Analysis of DoReCo languages
# Involved: Susanne Fuchs, Ludger Paschen, Lara S. Burchardt
# R Codes from: Lara S. Burchardt 

# Script  6 of 6

# Plots and Analyses for rhythm analysis


###############################################################################

# 00: preparations -----

## 00a: load packages -----

if (!require(install.load)) {
  install.packages("install.load")
}

library(install.load)

install_load("tidyverse", "psych", "tidygeocoder", "countrycode", "devtools",
             "lme4", "maps", "effsize","praatpicture", "grid", "ggplotify", "magick")

## 00b: prepare themes, color palettes, etc. ----

# color blind friendly palette with 28 colors for language families

colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
            "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#1f77b4", "#aec7e8", 
            "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728", "#ff9896", 
            "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2", 
            "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d")


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

doreco_rhythm_results_complete <- read_rds("rhythm_results_doreco_ioi_meta_complete.rds")

#ioi_data <- read_delim("iois_ipu_99quantilebreaks_means.csv", delim = ",")

ioi_data_alternative <- read_delim("unsplit_ioi_data_for_rhythm_analysis_including_meta_data_run_Oct24.csv", delim = ",")

# 02: data wrangling ----

## 02a: sort median ioi beat ----

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


## 02b: summarize by file ----

# after splitting file because of long silences, we summarize by file, to avoid pseudoreplication and inflation of n

doreco_rhythm_results_complete_summarized_file <- doreco_rhythm_results_complete %>%
  filter(!is.na(ioi_beat)) %>%
  group_by(Family, file) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), 
            across(where(is.character), first), 
            across(where(is.factor), first))


## 02c: z-scores ----

### 02c: raw scale ----
mean_ioi_beat <- mean(doreco_rhythm_results_complete$ioi_beat, na.rm = TRUE)
sd_ioi_beat <- sd(doreco_rhythm_results_complete$ioi_beat, na.rm = TRUE)

z_scores <- as.numeric((doreco_rhythm_results_complete$ioi_beat - mean_ioi_beat) / sd_ioi_beat)
z_scores <- as.data.frame(z_scores)

language_zscore <- cbind(z_scores, doreco_rhythm_results_complete$Language, as.numeric(doreco_rhythm_results_complete$ioi_beat))
colnames(language_zscore) <- c("zscore", "language", "ioi_beat")

outlier_indices_top <- language_zscore %>% 
  filter(z_scores >= 1.96)   # 1.96 corresponds to two tailed test, above/below will combine to be 5% --> significance level

min_outlier <- min(outlier_indices_top$ioi_beat) # vertical line x-position for density plot 
n_outliers_top <- nrow(outlier_indices_top)

outlier_indices_bottom <- language_zscore %>% 
  filter(z_scores <= -1.96)

max_outlier <- max(outlier_indices_bottom$ioi_beat) # vertical line x-position for density plot
n_outliers_bottom <- nrow(outlier_indices_bottom)

doreco_rhythm_results_complete <- cbind(doreco_rhythm_results_complete, z_scores)

all_zscore_outliers <- doreco_rhythm_results_complete %>% 
  filter(z_scores >= 1.96 | z_scores <= -1.96)

n_elements_95 <- nrow(doreco_rhythm_results_complete)-nrow(all_zscore_outliers)

###02c: log scale ----

mean_ioi_beat_log <- mean(log(doreco_rhythm_results_complete$ioi_beat), na.rm = TRUE)
sd_ioi_beat_log <- sd(log(doreco_rhythm_results_complete$ioi_beat), na.rm = TRUE)

z_scores_log <- as.numeric((log(doreco_rhythm_results_complete$ioi_beat) - mean_ioi_beat_log) / sd_ioi_beat_log)
z_scores_log <- as.data.frame(z_scores)
colnames(z_scores_log) <- c("z_score_log")

language_zscore_log <- cbind(z_scores, doreco_rhythm_results_complete$Language, as.numeric(log(doreco_rhythm_results_complete$ioi_beat)))
colnames(language_zscore_log) <- c("zscore_log", "language", "ioi_beat_log")

outlier_indices_top_log <- language_zscore_log %>% 
  filter(z_scores_log >= 1.96)   # 1.96 corresponds to two tailed test, above/below will combine to be 5% --> significance level

min_outlier_log <- min(outlier_indices_top_log$ioi_beat_log) # vertical line x-position for density plot 
n_outliers_top_log <- nrow(outlier_indices_top_log)

outlier_indices_bottom_log <- language_zscore_log %>% 
  filter(z_scores_log <= -1.96)

max_outlier_log <- max(outlier_indices_bottom_log$ioi_beat_log) # vertical line x-position for density plot
n_outliers_bottom_log <- nrow(outlier_indices_bottom_log)

doreco_rhythm_results_complete <- cbind(doreco_rhythm_results_complete, z_scores_log)

all_zscore_outliers_log <- doreco_rhythm_results_complete %>% 
  filter(z_scores_log >= 1.96 | z_scores_log <= -1.96)

n_elements_95_log <- nrow(doreco_rhythm_results_complete)-nrow(all_zscore_outliers_log)

# 02d: ioi data wrangling ----

meta_from_complete <- doreco_rhythm_results_complete %>% 
  select(filename, speaker, speaker_age, speaker_sex, synthesis, tone)

ioi_data_meta <- left_join(ioi_data_ir, meta_from_complete,
                           by = "filename")

# for plots: summarize ioi per file to get readable plots

ioi_data_meta_summarized_file <- ioi_data_meta %>% 
  group_by(filename) %>% 
  summarize(median_ioi = round(median(ioi, na.rm = TRUE), digits = 2))
            
ioi_data_meta_summarized_file <- left_join(ioi_data_meta_summarized_file, meta_from_complete, by = "filename")

# 03: statistics -----

# summarize

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
            filename_nr = length(unique(filename)))

summary_by_speaker <- doreco_rhythm_results_complete %>% 
  group_by(speaker) %>% 
  summarize(median_ioi_speaker = round(median(ioi_beat, na.rm = TRUE), digits = 2),
            filename_nr = length(unique(filename)),
            file_nr = length(unique(file)))

summary_by_language_iois <- ioi_data_alternative %>% 
  group_by(glottocode) %>% 
  summarize(nr_io_units = length(glottocode))
            
           
## 03a_1: correlations and effect sizes ioi beat----


### effect sizes ----
# (cohen's D, psych package)
# effect size gender

men <- doreco_rhythm_results_complete %>% 
  filter(speaker_sex == "m")
women <- doreco_rhythm_results_complete %>% 
  filter(speaker_sex == "f")

d_gender <- cohen.d(men$ioi_beat, women$ioi_beat)
# d = 0.10 --> negligible

t.test(men$ioi_beat, women$ioi_beat, var.equal = FALSE)
# p = 0.0002 --> significant but negligible

# effect size tone language

toneyes <- doreco_rhythm_results_complete %>% 
  filter(tone == "yes")
toneno <- doreco_rhythm_results_complete %>% 
  filter(tone == "no")

d_tone <- cohen.d(toneyes$ioi_beat, toneno$ioi_beat)
# d = 0.09 --> negligible

t.test(toneyes$ioi_beat, toneno$ioi_beat, var.equal = FALSE)
# p = 0.006 --> significant but negligible, change to former results most likely due to large sample size

### correlations ----

# correlation with median age

cor_age <- corr.test(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$speaker_age)
r_age <- cor_age$r
R_squared_age <- r_age^2

# ~ 0.013 % of variance are explained by age --> negligible

# correlation with morphological synthesis

cor_morph <- corr.test(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$synthesis)
r_morph <- cor_morph$r
R_squared_morph <- r_morph^2

# 1.08% of variance explained, so even though statistically significant, very small effect size


##03a_2: corr and eff size raw ioi -----

### effect sizes ----
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

### correlations ----

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

# plotting a histogram

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

## beat precision plots ------

#order languages by median beat precision
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

# 04: plots ----

## 04a: map plot ----

# map plot, where are languages spoken?
unique_coords <- unique(doreco_rhythm_results_complete[, c("Latitude", "Longitude")])

# Create a dataframe with unique coordinates
unique_coords_df <- data.frame(
  lon = unique_coords$Longitude,
  lat = unique_coords$Latitude)

unique_coords_df <- unique_coords_df %>%
  mutate(lon = ifelse(lon < 0, lon + 360, lon))

# Fortify the map data
mp1 <- fortify(map(fill=TRUE, plot=FALSE))
mp2 <- mp1
mp2$long <- mp2$long + 360
mp2$group <- mp2$group + max(mp2$group) + 1
mp <- rbind(mp1, mp2)

# Create the base plot with the world map

map <- ggplot() +  # No aesthetics here
  #geom_path(data = mp, aes(x = long, y = lat, group = group), fill = "grey75", color = "grey10") +  # World map outline
  geom_path(data = mp, aes(x = long, y = lat, group = group), fill = "#F2DDC1", color = "grey10") +  # World map outline
  geom_point(data = unique_coords_df, aes(x = lon, y = lat), color = "red", size = 2) +  # Your points
  scale_x_continuous(limits = c(0, 360)) +  # Set x limits
  coord_fixed(ratio = 1) +  # Ensure aspect ratio is fixed
  theme_minimal() + 
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 12, vjust = -0.5),
    axis.title.y = element_text(size = 12, vjust = 0.5)
  ) +
  labs(
    x = "Longitude [°]",
    y = "Latitude [°]"
  )


# version 2

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

## 04b: ioi distribution plots -----

## raw interval distribution

hist_ioi_raw <- ioi_data_alternative %>% 
  ggplot(aes(x= io_duration))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 color = "white", fill = "darkblue", bins = 100)+
  theme_minimal()+
  coord_cartesian(xlim = c(0,10))+
  annotate("text", x = 5, y = 5, label = paste("n =", nrow(ioi_data_alternative)), size = 6)+
  ylab("Percentage [%]")+
  xlab("IOI [sec]")+
  #annotate("text", x = 0.6, y = 9, label = "n = 1535")+
  my_custom_theme
print(hist_ioi_raw)


## ioi beat plots 
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
  annotate("text", x = 1.0, y = 15, label = paste("n =", nrow(doreco_rhythm_results_complete)), , size = 6)+
  my_custom_theme
print(hist_cv)

### revision version -----

hist_qcv <- qcv_file %>% 
  ggplot(aes(x= qcv))+
  geom_histogram(aes(y=stat((count)/sum(stat(count))*100)),
                 #color = "white", fill = "darkblue", binwidth = 0.1)+
                  color = "white", fill = "darkblue", bins = 100)+
  theme_minimal()+
  ylab("Percentage [%]")+
  xlab("QCV of IOI Durations per Sequence")+
  coord_cartesian(xlim = c(0,1))+
  annotate("text", x = 0.8, y = 15, label = paste("n =", nrow(doreco_rhythm_results_complete)), , size = 6)+
  my_custom_theme
print(hist_qcv)


### ioi in boxplot per language ----

ioi_data_ir <- left_join(ioi_data_ir, meta_languages, by = "Language", multiple = "any")

median_raw_ioi_lan <- ioi_data_ir %>%
  group_by(Language) %>%
  summarize(median_ioi = median(ioi, na.rm = TRUE))

ioi_data_ir$Language <- factor(ioi_data_ir$Language,
                               levels = median_raw_ioi_lan$Language[order(median_raw_ioi_lan$median_ioi)])

# boxplot - not sorted yet

boxplot_languages_ioi <- ioi_data_ir %>% 
  dplyr::filter(is.na(ioi) == FALSE) %>% 
  group_by(Language) %>% 
  ggplot(aes(y = ioi, x = Language, fill = Family))+
  #geom_boxplot(outliers = FALSE)+
  geom_boxplot()+
  #geom_jitter(alpha = 0.2, size = 0.5, shape = 01)+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  coord_cartesian(ylim = c(0,12.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, 'cm'))+
  xlab("Languages")+
  ylab("IOI [sec]")+
  my_custom_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))


## 04c: zscore plots ----

###04c: zscore plot raw ----
# density plot of all calculated ioi beats across all the dataset

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

#ggsave("plot_density_zscore.jpg", dpi = 300,
#       width = 12,
#       height = 4)


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

#additional plot for revision: qvc 
density_qcv <- 
  qcv_file %>% 
  ggplot(aes(x = qcv))+
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
  xlab("Quantile Coefficient of Variation")+
  my_custom_theme


cowplot::plot_grid(density_cv, density_qcv, labels = c("A", "B"))

ggsave("plot_density_variation.jpg", dpi = 300,
       width = 12,
       height = 4)


# additional plot for revision density of raw iois 

density_ioi <- 
  ioi_data_ir %>% 
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
  xlab("log- transformed ioi [sec]")+
  my_custom_theme

## what are the outliers? is there a trend here? age? 
doreco_rhythm_results_complete %>% 
  #ggplot(aes(x= unbiased_cv, y = z_scores$z_scores))+
  ggplot(aes(x= ugof_ioi, y = npvi))+
  geom_point()#+
  geom_smooth()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= npvi, y = z_scores))+
  geom_point()+
  geom_smooth()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= n_elements, y = zscores_log))+
  geom_point()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= synthesis, y = z_scores_log))+
  geom_point()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= speaker_age, y = z_scores))+
  geom_point()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x = speaker_sex, y = z_scores))+
  geom_boxplot()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x= speaker_age, y = avg_height))+
  geom_point()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x = Area, y = z_scores))+
  geom_boxplot()

doreco_rhythm_results_complete %>% 
  ggplot(aes(x = tone, y = z_scores))+
  geom_boxplot()

###04c: zscore plot log ----

density_ioibeat_log <- 
  doreco_rhythm_results_complete %>% 
  ggplot(aes(x = log(ioi_beat)))+
  geom_density()+
  geom_vline(xintercept = min_outlier_log, linetype="dotted", linewidth = 2 )+  # zscore > 2, corresponding ioi beat
  geom_vline(xintercept = max_outlier_log, linetype="dotted", linewidth = 2)+  # zscore < -2, corresponding ioi beat
  geom_jitter(aes(y = -1, color = Language), alpha = 0.5, size = 0.5)+
  annotate("text", x = max_outlier_log-0.3, y = -0.3, label = paste("n =", n_outliers_bottom_log), size = 5)+
  annotate("text", x = max_outlier_log-0.3, y = 1.5, label = expression("Z-Score "<="-1.96"), size = 5)+
  annotate("text", x = -0.8, y = -0.3, label = paste("n =",n_elements_95_log ), size = 5)+
  annotate("text", x = 0.5, y = -0.3, label = paste("n =", n_outliers_top_log), size = 5)+
  annotate("text", x = 0.5, y = 1.5, label = expression("Z-Score ">="1.96"), size = 5)+  
  #coord_cartesian(xlim = c(-0.5, 3))+
  ylab("Density")+
  xlab("log - IOI Beat [Hz]")+
  my_custom_theme


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

cor_morph <- cor(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$synthesis,
                 na.rm = TRUE)

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
  group_by(speaker) %>% 
  ggplot(aes(x= synthesis, y = median_ioi ))+
  geom_jitter(alpha = 0.3, size = 0.5)+
  geom_smooth(method = "lm", color = "#0072B2")+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab('Morphological Synthesis')+
  ylab("IOI [sec]")+
  my_custom_theme+
  annotate("text", x = 3, y = 6, label = paste("R² =", round(R_squared_morph_ioi, 2)), size = 5)



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
  annotate("text", x = 5, y = 4, label = expression("Cutoff Criterion, 99th Quantile"), size = 4, color = "darkred")+
  xlab("Pause Duration [sec]") +
  ylab("Percentage [%]") +
  coord_cartesian(ylim= c(0,10), xlim = c(0,7.5))+
  theme_minimal()

hist_pauseduration

## plot ioi beat vs cv


ggplot(doreco_rhythm_results_complete, aes(x = ioi_beat  , y = unbiased_cv)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm", color = "#0072B2")+
  #labs(x = "Median Age per Language",
  #     y = "Median IOI Beat [Hz] per Language") +
  theme_minimal()+
  my_custom_theme


##04m: fft plot -----


fourier_beat <- rhythm_analysis_rerun_fft %>% 
  dplyr::filter(is.na(fourier_beat) == FALSE) %>% 
  ggplot(aes(x = fourier_beat))+
  geom_density()
  
  geom_boxplot(outliers = FALSE)+
  geom_jitter(alpha = 0.2, size = 0.7)


boxplot_languages_fft <- rhythm_analysis_rerun_fft %>% 
  dplyr::filter(is.na(fourier_beat) == FALSE) %>% 
  #group_by(speaker, Language) %>% 
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


# 05: plot grids -----

## 05a: Figure 1 ----

# part 3 of Figure 1:  Histograms IOI distribution [sec]  

hist_plots <- cowplot::plot_grid(hist_ioi_raw, hist_cv,hist_n_element,
                   labels = c("D", "E", "F"), ncol = 3)

ggsave("hist_plot_fig1_part3_v2.jpg", dpi = 300, hist_plots,
       width = 12,
       height = 6)

#part 1: annotations and part 2: map are generated and saved further up, figures where
# combined outside of R

## 05b: Figure 2-----

# Density Plot & Language Boxplots

cowplot::plot_grid(density_ioibeat, boxplot_languages,
                   nrow = 2,
                   rel_heights = c(0.4, 0.6),
                   labels = c("A", "B"))

ggsave("sup_figure_ioibeats.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

## 05c: Figure 3 ----

# original submission ioi beat
cowplot::plot_grid(scatter_age, scatter_morph,
                   box_gender, box_tone, ncol = 2,
                   labels = c("A", "B", "C", "D"))


ggsave("manuscript_figure3.jpg", dpi = 300,
       width = 22,
       height = 20,
       units = "cm")

ggsave("manuscript_figure3.pdf", dpi = 300,
       width = 22,
       height = 20,
       units = "cm")

# revision: raw ioi 

cowplot::plot_grid(scatter_age_ioi, scatter_morph_ioi,
                   box_gender_ioi, box_tone_ioi, ncol = 2,
                   labels = c("A", "B", "C", "D"))

ggsave("manuscript_figure3_ioi.pdf", dpi = 300,
       width = 22,
       height = 20,
       units = "cm")

## 05d: Supp Figure 1 -----
# plot of IPU durations and pause durations histograms

cowplot::plot_grid(hist_sprachdauer, hist_pauseduration, ncol = 2,
                   labels = c("A", "B"))

ggsave("supplements_figure_ipu_pause_duration.jpg", dpi = 300,
       width= 12,
       height = 6)

## 05e: Supp Figure 2 ----
# main plot 2, but for cv values

cowplot::plot_grid(density_cv, boxplot_languages_cv,nrow = 2,
                   rel_heights = c(0.4, 0.6),
                   labels = c("A", "B"))

ggsave("supplements_figure_CV.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

## 05f: Figure 2 but raw IOI ----

cowplot::plot_grid(density_plot_ioi_language_raw, boxplot_languages_ioi,
                   nrow = 2,
                   rel_heights = c(0.4, 0.6),
                   labels = c("A", "B"))

ggsave("alt_figure_2_raw_ioi.pdf", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

##05g: beat precision supp figure ------

cowplot::plot_grid(hist_bp, boxplot_languages_bp,
                   nrow = 2,
                   rel_heights = c(0.3, 0.7),
                   labels = c("A", "B"))

ggsave("sup_figure_beatprecision.jpg", dpi = 300,
       width = 28,
       height = 26,
       units = "cm")

# 06: revision - additional analysis ----

## 06a: load integer ratios ----

ioi_data_ir <- read_delim("rhythm_analysis_results/revision/iois_doreco_rerun_ir.csv", delim = ",")
## 06b: integer ratios ------

# unique filename extraction
unique_files <- unique(ioi_data_ir$filename)

# folder for intermediate results
output_dir <- "rhythm_analysis_results/revision/intermediate_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# calculations per file

for (file in unique_files) {
  output_file <- file.path(output_dir, paste0("results_", file, ".csv"))
  
  # skip if output file already exists
  if (file.exists(output_file)) {
    next
  }
  
  # filter data
  data <- ioi_data_ir %>% filter(filename == file) %>% select(ioi)
  
  n <- nrow(data)
  
  if (n < 2) {
    next  # skip if not enough rows
  }
  
  # initialiese results table
  results <- expand.grid(i = 1:n, j = 1:n) %>%
    filter(i != j) %>%  # Entferne Fälle, wo i == j
    mutate(ratio = data$ioi[i] / (data$ioi[i] + data$ioi[j]))
  
  # save intermediate results per file 
  write_csv(results, output_file)
}

# combine all results
result_files <- list.files(output_dir, full.names = TRUE, pattern = "results_.*.csv")
all_results <- result_files %>% 
  map_dfr(~ read_csv(.x) %>% mutate(filename = basename(.x)))

write_csv(all_results, "doreco_interger_ratios_all_pairs_all_files.csv")
## 06c: analyse integer ratios -----

all_results <- read_delim("doreco_interger_ratios_all_pairs_all_files.csv")

mean(all_results$ratio, na.rm = TRUE)
min(all_results$ratio, na.rm = TRUE)
max(all_results$ratio, na.rm = TRUE)
median(all_results$ratio, na.rm = TRUE)
sd(all_results$ratio, na.rm = TRUE)



language_file <- doreco_rhythm_results_complete %>% 
  select(Language, filename)

#all_results$filename <- str_sub(all_results$filename, end = -5, start = 9)

all_results <- left_join(all_results, language_file, by = "filename")


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


density_plot_ir_raw <- all_results %>% 
  ggplot(aes(x = ratio))+
  geom_density(aes(color = Language))+
  #geom_density()+
  #geom_vline(xintercept = 1/2.25, linetype="dotted", linewidth = 2 )+  
  #geom_vline(xintercept = 0.56, linetype="dotted", linewidth = 2)+
  geom_density(data = random_ir, aes(x = ratio),size = 1.3, color = "black")+
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

# Check available columns (useful for debugging)
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

# new idea- this density plot of ir per language but for iois 

ioi_data_ir <- left_join(ioi_data_ir, language_file, by = "filename")

density_plot_ioi_language_raw <- ioi_data_ir %>% 
  ggplot(aes(x = ioi))+
  geom_density(aes(color = Family))+
  geom_density(data =sequence_random, aes(x = sequence_random), color = "black", size = 1.3)+
  #coord_cartesian(xlim = c(0, 2))+
  scale_color_manual(values = colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("IOI [sec]")+
  ylab(" Density")

density_plot_ioi <- ioi_data_ir %>% 
  ggplot(aes(x = ioi))+
  geom_density(aes(color = Family))+
  #geom_density(data =sequence_random, aes(x = sequence_random), color = "black", size = 1.3)+
  #coord_cartesian(xlim = c(0, 2))+
  #geom_jitter(aes(y = -1, color = Family), alpha = 0.2, size = 0.5, shape = 01)+
  theme_minimal()+
  coord_cartesian(xlim = c(0,7.5))+
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


#cowplot::plot_grid(density_plot_ioi_language_raw, density_plot_ioi_language_log,
#                   density_plot_ir_raw,density_plot_ir_log,
#                   ncol = 2, labels = c("A", "B", "C", "D"))

cowplot::plot_grid(density_plot_ioi_language_raw, 
                   density_plot_ir_raw,
                   ncol = 2, labels = c("A", "B"))

ggsave("density_plots_per_language.jpg", dpi = 300,
       width = 14,
       height = 8)

# transfered matlab code from roeske paper ----

# Histogrammfunktion für Ratios
plot_ratio_hist <- function(ratios, cycledurations, IS_CDSLOPASS = FALSE, IS_CDFASTPASS = FALSE, cd_cutoff = 0, IS_OPENNEWFIG = TRUE) {
  if (IS_OPENNEWFIG) {
    dev.new()
  }
  
  if (IS_CDSLOPASS) {
    allratio_limitedcd <- ratios[cycledurations >= cd_cutoff]
  } else if (IS_CDFASTPASS) {
    allratio_limitedcd <- ratios[cycledurations < cd_cutoff]
  } else {
    allratio_limitedcd <- ratios
  }
  
  nbins <- round(length(allratio_limitedcd) / 80)
  nbins <- max(100, min(nbins, 200))
  
  hist_color <- rgb(0.75, 0.75, 0.75, alpha = 1)
  
  hist(allratio_limitedcd, breaks = nbins, probability = TRUE, col = hist_color, border = hist_color, main = paste( "ratios pdf"), xlab = "ratio", ylab = "probability density", xlim = c(0, 1))
  
  abline(v = c(1/2, 1/3, 1/4, 2/3, 3/4), col = "green", lty = 2, lwd = 1.5)
  
  if (length(allratio_limitedcd) > 0) {
    density_data <- density(allratio_limitedcd, bw = 0.012)
    lines(density_data$x, density_data$y, col = "black", lwd = 2)
  } else {
    message("No ratio data, perhaps you excluded this tempo range?")
  }
}

plot_ratio_hist(all_results$ratio, cycledurations = NULL)

## rose plot

library(tidyverse)


# Function to compute angular deviations from the reference beat
compute_circular_deviation <- function(iois, reference_ioi) {
  # Convert IOIs to actual onset times
  actual_onsets <- cumsum(c(0, iois))  # Start at 0 and accumulate IOIs
  
  # Generate theoretical reference beats
  cumulative_beats <- seq(0, by = reference_ioi, length.out = length(actual_onsets))
  
  # Ensure closest_reference has the same length as actual_onsets
  closest_reference <- sapply(actual_onsets, function(t) {
    if (length(cumulative_beats) == 0) return(NA)  # Fallback safety check
    ref_index <- which.min(abs(cumulative_beats - t))
    if (length(ref_index) == 0) return(NA)  # Another safeguard
    return(cumulative_beats[ref_index])
  })
  
  # Compute deviations
  deviation <- actual_onsets - closest_reference
  angle_deviation <- (deviation / reference_ioi) * 360  # Convert to angular deviation
  angle_wrapped <- (closest_reference %% reference_ioi) / reference_ioi * 360  # Wrap angles
  
  # Final check: Ensure all lengths match before returning
  print(length(actual_onsets))  # Debug
  print(length(closest_reference))  # Debug
  stopifnot(length(actual_onsets) == length(closest_reference))
  
  # Return results
  tibble(
    interval = iois,
    actual_onset = actual_onsets[-1],  # Exclude the first (0) onset
    closest_beat = closest_reference[-1],
    deviation = deviation[-1],  
    angle_wrapped = angle_wrapped[-1],  
    angle_deviation = angle_deviation[-1]
  )
}

# Example Data
#ioi_data_ir <- tibble(
#  filename = rep("example", 6),
#  ioi = c(0.48, 0.52, 0.47, 0.49, 0.51, 0.50)  # Sample inter-onset intervals
#)

# Create a list of IOIs per filename
ioi_list <- ioi_data_ir %>%
  group_by(filename) %>%
  summarise(iois = list(ioi)) %>%
  pull(iois)

# Create a dataframe with references
reference_df <- ioi_data_ir %>%
  group_by(filename) %>%
  summarise(reference_ioi = mean(ioi, na.rm = TRUE))

# Initialize a list to store results
deviation_results_list <- list()

# Loop through the list of IOIs and corresponding references
for (i in seq_along(ioi_list)) {
  iois <- ioi_list[[i]]
  reference_ioi <- reference_df$reference_ioi[i]
  
  # Compute deviations for the current filename
  deviation_results <- compute_circular_deviation(iois, reference_ioi)
  
  # Add filename to the results
  deviation_results <- deviation_results %>%
    mutate(filename = reference_df$filename[i])
  
  # Store results in the list
  deviation_results_list[[i]] <- deviation_results
}

# Combine all results into a single dataframe
circular_deviation_results <- bind_rows(deviation_results_list)

# Plot as a rose plot
#ggplot(circular_deviation_results, aes(x = (angle_wrapped + angle_deviation) %% 360)) +

ggplot(circular_deviation_results,aes(x = ((closest_reference - min(closest_reference) + deviation) / reference_ioi * 360) %% 360)) +
  #aes(x = ((closest_reference - min(closest_reference) + deviation) / reference_ioi * 360) %% 360)
    geom_histogram(aes(y = ..count..), bins = 36, fill = "blue", alpha = 0.7, color = "black") +
  coord_polar(start = 0) +  
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45), labels = seq(0, 360, by = 45)) +
  scale_y_continuous(limits = c(0,10000))+
  labs(
    title = "Deviation from Reference Beat (Rose Plot)",
    x = "Phase (°)",
    y = "Count"
  ) +
  theme_minimal()
### deviation plot

ggplot(circular_deviation_results,aes(x = deviation)) +
  #aes(x = ((closest_reference - min(closest_reference) + deviation) / reference_ioi * 360) %% 360)
  geom_histogram(aes(y = ..count..), bins = 100, fill = "blue", alpha = 0.7, color = "black")

## 06d: quartile cv ----

#cv overall
ioi_cv <- sd(ioi_data_ir$ioi, na.rm = TRUE)/mean(ioi_data_ir$ioi, na.rm = TRUE)

# calcualte quartile cv values for all sequences
# source: https://doi.org/10.1038/s41598-023-31711-8
# function: "cv" 

qcv_lan <- ioi_data_ir %>%
  drop_na() %>% 
  group_by(Language) %>%
  summarize(qcv = cv(ioi, method = "quartile"))

qcv_file <- ioi_data_ir %>% 
  drop_na() %>% 
  group_by(filename) %>% 
  summarize(qcv = cv(ioi, method = "quartile"))

qcv_file <- left_join(qcv_file, language_file, by = "filename")

for(i in 1:length(qcv_file$filename)){
  names <- str_split(qcv_file$filename[i], pattern = "_")
  names <- names[[1]]
  qcv_file$speaker[i] <- names[[3]]   #filename from raw file, to use to join rhythm results with meta data file
}


qcv_file %>% 
  group_by(Language) %>% 
  ggplot(aes(y = log(qcv), group = Language))+
  geom_boxplot()+
  theme_minimal()

qcv_file %>% 
  ggplot(aes(x = qcv))+
  geom_density()+
  theme_minimal()

density_qcv <- 
  qcv_file %>% 
  ggplot(aes(x = sqrt(qcv)))+
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
  xlab("sqrt Quartile Coefficient of Variation")+
  my_custom_theme


library(lme4)
#using lmer

model_ioi_beat_lmer <- lmer(log(ioi_beat) ~ (1 | Language), 
                            data = doreco_rhythm_results_complete)
icc(model_ioi_beat_lmer)

# Fit Gamma GLMM with log link
model <- glmer(qcv ~ (1 | Language), 
               data = qcv_file, 
               family = Gamma(link = "log"))

# View model summary
summary(model)


# Model 1: qcv explained by Language (Gamma GLMM)
model_qcv <- glmer(qcv ~ (1 | Language), 
                   data = qcv_file, 
                   family = Gamma(link = "log"))

library(performance)
icc(model_qcv)

model_qcv_speaker <- glmer(qcv ~ (1| speaker),
                           data = qcv_file, 
                           family = Gamma(link = "log"))
icc(model_qcv_speaker)

model_ioi <- glmer(ioi_beat ~  (1 | speaker), 
                   data = doreco_rhythm_results_complete, 
                   family = Gamma(link = "log"))

summary(model_ioi)
icc(model_ioi)

model_ioi_raw <- glmer(io_duration ~ (1| glottocode),
                       data = data_ipu,
                       family = Gamma(link = "log"))
summary(model_ioi_raw)
icc(model_ioi_raw)


model_cv <- glmer(unbiased_cv ~ (1 | Language), 
                   data = doreco_rhythm_results_complete, 
                   family = Gamma(link = "log"))

icc(model_cv)

model_ioi_speaker <- glmer(ioi_beat ~(1| speaker),
                           data = doreco_rhythm_results_complete,
                           family = Gamma(link = "log"))
icc(model_ioi_speaker)

# Model 2: integer_ratio explained by Language (LMM)
model_integer_ratio <- lmer(ratio ~ (1 | Language), 
                            data = all_results)

summary(model_integer_ratio)

# model 3: ioi beat explained by language

model_ioi_beat <- glmer(unbiased_cv ~ (1 | Language), 
                        data = doreco_rhythm_results_complete, 
                        family = Gamma(link = "log"))
summary(model_ioi_beat)

# model 4: qcv by language, speaker as random factor

full_model <- glmer(qcv ~ Language + (1 | speaker), data = qcv_file, family = Gamma(link = "log"))
summary(model_qcv_4)

null_model <- glmer(qcv ~ 1 + (1 | speaker), data = qcv_file, family = Gamma(link = "log"))
anova(null_model, full_model, test = "Chisq")

full_model_fixed <- glm(qcv ~ Language, data = qcv_file, family = Gamma(link = "log"))
summary(full_model_fixed)

AIC(full_model, full_model_fixed)
BIC(full_model, full_model_fixed)

anova(full_model, full_model_fixed, test = "Chisq")

library(MuMIn)
r.squaredGLMM(full_model)  # With speaker effect
r.squaredGLMM(full_model_fixed)  # Without speaker effect

### 06d-1 is language a good predictor of qcv? ------
# step 1) fit fixed effect model using language

full_model_fixed <- glm(qcv ~ Language, data = qcv_file, family = Gamma(link = "log"))

# step 2) fit null-model
null_model_fixed <- glm(qcv ~ 1, data = qcv_file, family = Gamma(link = "log"))

# step 3) likelihood ratio test

anova(null_model_fixed, full_model_fixed, test = "Chisq")

# step 4: interpretation and effect sizes

summary(full_model_fixed)
confint(full_model_fixed)
r.squaredGLMM(full_model_fixed)


# estimates of full_model_fixed are effect sizes (log transformed)
# we exponentiate to get the raw value and add a column whether the effect is low, moderate. strong

estimates <- coef(full_model_fixed)

# Create a dataframe
estimates_df <- data.frame(
  Predictor = names(estimates),
  Estimate = estimates
)

estimates_df$raw_estimates <- exp(estimates_df$Estimate)
estimates_df$diff_reference <- abs(1-estimates_df$raw_estimates)

estimates_df <- estimates_df %>%
  mutate(interpretation = case_when(
    abs(diff_reference) < 0.10 ~ "low",
    abs(diff_reference) >= 0.10 & abs(diff_reference) < 0.30 ~ "moderate",
    abs(diff_reference) >= 0.30 ~ "high"
  ))

## 06e: is ioi beat a good predictor for language 
# same appraoch as above, but with ioi beat instead of qcv

full_model_ioi <- glmer(ioi_beat ~ Language + (1 | speaker), data = doreco_rhythm_results_complete, family = Gamma(link = "log"))
summary(full_model_ioi)

null_model_ioi <- glmer(ioi_beat ~ 1 + (1 | speaker), data = doreco_rhythm_results_complete, family = Gamma(link = "log"))
anova(null_model_ioi, full_model_ioi, test = "Chisq")

full_model_fixed_ioi <- glm(ioi_beat ~ Language, data = doreco_rhythm_results_complete, family = Gamma(link = "log"))
summary(full_model_fixed_ioi)

AIC(full_model_ioi, full_model_fixed_ioi)
BIC(full_model_ioi, full_model_fixed_ioi)

anova(full_model_ioi, full_model_fixed_ioi, test = "Chisq")

r.squaredGLMM(full_model_ioi)  # With speaker effect
r.squaredGLMM(full_model_fixed_ioi)  # Without speaker effect

### 06e-1 is language a good predictor of ioi beat? ------
# step 1) fit fixed effect model using language

full_model_fixed_ioi <- glm(ioi_beat ~ Language, data = doreco_rhythm_results_complete, family = Gamma(link = "log"))

# step 2) fit null-model
null_model_fixed_ioi <- glm(ioi_beat ~ 1, data = doreco_rhythm_results_complete, family = Gamma(link = "log"))

# step 3) likelihood ratio test

anova(null_model_fixed_ioi, full_model_fixed_ioi, test = "Chisq")

# step 4: interpretation and effect sizes

summary(full_model_fixed_ioi)
confint(full_model_fixed_ioi)
r.squaredGLMM(full_model_fixed_ioi)


# estimates of full_model_fixed are effect sizes (log transformed)
# we exponentiate to get the raw value and add a column whether the effect is low, moderate. strong

estimates <- coef(full_model_fixed_ioi)

# Create a dataframe
estimates_df_ioi <- data.frame(
  Predictor = names(estimates),
  Estimate = estimates
)

estimates_df_ioi$raw_estimates <- exp(estimates_df$Estimate)
estimates_df_ioi$diff_reference <- abs(1-estimates_df$raw_estimates)

estimates_df_ioi <- estimates_df %>%
  mutate(interpretation = case_when(
    abs(diff_reference) <= 0.05 ~ "very low",
    abs(diff_reference) < 0.10 ~ "low",
    abs(diff_reference) >= 0.10 & abs(diff_reference) < 0.30 ~ "moderate",
    abs(diff_reference) >= 0.30 ~ "high"
  ))

estimates_df_ioi %>% 
  count(interpretation)

### 06f-1 is language a good predictor of unbiased cv? ------

full_model_cv <- glmer(unbiased_cv ~ Language + (1 | speaker), data = doreco_rhythm_results_complete, family = Gamma(link = "log"))
summary(full_model_cv)

null_model_cv <- glmer(unbiased_cv ~ 1 + (1 | speaker), data = doreco_rhythm_results_complete, family = Gamma(link = "log"))
anova(null_model_cv, full_model_cv, test = "Chisq")

full_model_fixed_cv <- glm(unbiased_cv ~ Language, data = doreco_rhythm_results_complete, family = Gamma(link = "log"))
summary(full_model_fixed_cv)

AIC(full_model_ioi, full_model_fixed_cv)
BIC(full_model_ioi, full_model_fixed_cv)

anova(full_model_cv, full_model_fixed_cv, test = "Chisq")

r.squaredGLMM(full_model_cv)  # With speaker effect
r.squaredGLMM(full_model_fixed_cv)  # Without speaker effect

### 06e-1 is language a good predictor of qcv? ------
# step 1) fit fixed effect model using language

full_model_fixed_cv <- glm(unbiased_cv ~ Language, data = doreco_rhythm_results_complete, family = Gamma(link = "log"))

# step 2) fit null-model
null_model_fixed_cv <- glm(unbiased_cv ~ 1, data = doreco_rhythm_results_complete, family = Gamma(link = "log"))

# step 3) likelihood ratio test

anova(null_model_fixed_cv, full_model_fixed_cv, test = "Chisq")

# step 4: interpretation and effect sizes

summary(full_model_fixed_cv)
confint(full_model_fixed_cv)
r.squaredGLMM(full_model_fixed_cv)


# estimates of full_model_fixed are effect sizes (log transformed)
# we exponentiate to get the raw value and add a column whether the effect is low, moderate. strong

estimates <- coef(full_model_fixed_cv)

# Create a dataframe
estimates_df_cv <- data.frame(
  Predictor = names(estimates),
  Estimate = estimates
)

estimates_df_cv$raw_estimates <- exp(estimates_df$Estimate)
estimates_df_cv$diff_reference <- abs(1-estimates_df$raw_estimates)

estimates_df_cv <- estimates_df %>%
  mutate(interpretation = case_when(
    abs(diff_reference) <= 0.05 ~ "very low",
    abs(diff_reference) < 0.10 ~ "low",
    abs(diff_reference) >= 0.10 & abs(diff_reference) < 0.30 ~ "moderate",
    abs(diff_reference) >= 0.30 ~ "high"
  ))

estimates_df_cv %>% 
  count(interpretation)

# 07: bar plot, explanation sequences ----

library(jpeg)
library(grid)

img <- readJPEG("part_figure1_stilizied_sequence.jpg")  # Replace with your actual file path
img_grob <- rasterGrob(img, interpolate = TRUE)




# Create a dataframe with event times, widths, and sequences
data <- data.frame(
  time = c(1, 2.5, 4, 5.5, 7, 12, 13.5, 15, 16.8, 18.5, 20.2),
  width = c(0.8, 1, 0.6, 1.2, 0.9, 1, 0.7, 1.1, 0.8, 1, 0.9),
  sequence = c(rep("Sequence 1", 5), rep("Sequence 2", 6))
)

# Plot
p <- ggplot(data, aes(x = time, y = 1, fill = sequence, width = width)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("Sequence 1" = "#1B4F72", "Sequence 2" = "#1B4F72")) +
  
  # Add annotations
  annotate("text", x = c(1.75, 3.25, 4.75, 6.25), y = 1.3, label = "IOI", size = 5) +
  annotate("text", x = c(13, 15, 17), y = 1.3, label = "IPU", size = 5) +
  annotate("text", x = 7, y = 0.5, label = "Silent Pause > Cut-Off\nCriterion", size = 4, fontface = "bold", color = "red") +
  
  # Add sequence labels
  annotate("text", x = 3, y = -0.2, label = "Sequence 1", size = 5, fontface = "bold", color = "red") +
  annotate("text", x = 15, y = -0.2, label = "Sequence 2", size = 5, fontface = "bold", color = "red") +
  
  # Add IOI arrows from bar start to next bar start
  geom_segment(aes(x = 1, xend = 2.5, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 2.5, xend = 4, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 4, xend = 5.5, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 5.5, xend = 7, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  
  # Add IPU arrows from bar start to bar end
  geom_segment(aes(x = 12, xend = 12 + 1, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 13.5, xend = 13.5 + 0.7, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 15, xend = 15 + 1.1, y = 1.4, yend = 1.4), arrow = arrow(length = unit(0.2, "cm"))) +
  
  # Adjust theme
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  labs(x = "Time", y = "")

print(p)

# 08: distribution ----

library(fitdistrplus)
fit_norm <- fitdist(doreco_rhythm_results_complete$ioi_beat, "norm")
fit_lognorm <- fitdist(doreco_rhythm_results_complete$ioi_beat, "lnorm")
fit_gamma <- fitdist(doreco_rhythm_results_complete$ioi_beat, "gamma")
plot(fit_norm)
plot(fit_lognorm)
plot(fit_gamma)
gofstat(list(fit_norm, fit_lognorm, fit_gamma))

# results: log normal distribution is the best fit for all statistics; KS, Cramer von Mises, Anders Darling, AIC and BIC
# therefore the best modelling approach would be a Gamma model with log link

# 09: Modelling -----

## 09a: ANOVA ----
library(effectsize)
anova_model <- aov(log(ioi_beat) ~ Language, data = doreco_rhythm_results_complete)
summary(anova_model)

eta_squared(anova_model, generalized = TRUE)  # ηg²
TukeyHSD(anova_model)


library(effsize)
cohen.d(doreco_rhythm_results_complete$ioi_beat, doreco_rhythm_results_complete$Language)
