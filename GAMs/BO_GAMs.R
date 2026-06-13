##--------------------------------
#GAM of in situ transects

# Load necessary library
library(tidyverse)
library(ggplot2)
library(ggdist)


# Read the file (adjust separator if needed: "\t" for tab, "" for whitespace)
raw_data <- read.delim("../in situ radiometry/Bolonia - okamurae/20250613_okamurae_bolonia.txt", check.names = FALSE, stringsAsFactors = FALSE)

# Convert all comma decimals to dots (except the first column with 'Wavelength')
raw_data[-1] <- lapply(raw_data[-1], function(x) as.numeric(gsub(",", ".", x)))

# Make sure Wavelength is numeric
raw_data$Wavelength <- as.numeric(raw_data$Wavelength)

# raw_data$MeanReflectance <- rowMeans(raw_data[-1], na.rm = TRUE)
raw_data_long <- raw_data %>% 
  pivot_longer(-1,names_to = 'old_name', values_to = 'reflec') %>% 
  group_by(old_name) %>% 
  mutate(id=cur_group_id(),
         #reflec = reflec*1000,
         group=as.factor(case_when(id%in%c(22:41)~'RB1',
                                   id%in%c(42:61)~'RB2',
                                   id%in%c(62:81)~'RB3',
                                   id%in%c(82:101)~'RB4',
                                   id%in%c(122:181)~'RB5',
                                   id%in%c(182:241)~'RB6',
                                   id%in%c(262:281)~'RB7',
                                   id%in%c(282:321)~'RB8',
                                   id%in%c(322:361)~'RB9'))) %>% 
  filter(Wavelength>399&Wavelength<901) %>% 
  filter(!is.na(group))


ggplot(raw_data_long, aes(x=Wavelength, y= reflec, color=group)) + 
  geom_smooth(se = TRUE)


library(mgcv)
#gam_1 <- gam(reflec~s(Wavelength, by=group)+group, family=Gamma(),data=raw_data_long) #gaussian was the other option.
library(mgcv)
#gam_1 <- gam(reflec~s(Wavelength, by=group)+group, family=Gamma(),data=raw_data_long) #gaussian() was the other option.
gam_1 <- bam(
  reflec ~ s(Wavelength, by = group, k = 50) + group,
  family = Gamma(),
  data = raw_data_long,
  method = "REML",
  discrete = TRUE,
  nthreads = parallel::detectCores() - 1  # use multiple cores if available
)


ModelOutputs<-data.frame(Fitted=fitted(gam_1),
                         Residuals=resid(gam_1))

p3<-ggplot(ModelOutputs)+
  geom_point(aes(x=Fitted,y=Residuals))+
  theme_classic()+
  labs(y="Residuals",x="Fitted Values")
#p3

p4<-ggplot(ModelOutputs) +
  stat_qq(aes(sample=Residuals))+
  stat_qq_line(aes(sample=Residuals))+
  theme_classic()+
  labs(y="Sample Quartiles",x="Theoretical Quartiles")
#p4

#install.packages("patchwork")
library(patchwork)
p3+p4

summary(gam_1)


NewData_1<-expand_grid(Wavelength=seq(min(raw_data_long$Wavelength),max(raw_data_long$Wavelength),length.out=500),
                       id=min(raw_data_long$id):max(raw_data_long$id)
) %>% 
  mutate(group=as.factor(case_when(id%in%c(22:41)~'RB1',
                                   id%in%c(42:61)~'RB2',
                                   id%in%c(62:81)~'RB3',
                                   id%in%c(82:101)~'RB4',
                                   id%in%c(122:181)~'RB5',
                                   id%in%c(182:241)~'RB6',
                                   id%in%c(262:281)~'RB7',
                                   id%in%c(282:321)~'RB8',
                                   id%in%c(322:361)~'RB9')))


Pred<-predict(gam_1,NewData_1,se.fit=T,type="response") # creates a list of predicted mean value y'

NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

# Plot using ggplot2 with custom palette
p <- ggplot(NewData) +
  geom_ribbon(aes(x = Wavelength,
                  ymax = Upr,
                  ymin = Lwr,
                  fill = as.factor(group)), alpha = 0.2) +
  geom_line(aes(x = Wavelength,
                y = response,
                colour = as.factor(group))) +
  scale_color_manual(values = c(
    "#F4DBB5",
    "#E8B878",
    "#CF8552",  
    "#B45035",  
    "#8C6300",  
    "#654500",  
    "#8DBBDC",  
    "#1f78b4",   
    "#376491"),
    na.translate = FALSE) +   # <— removes NA from color legend)
  scale_fill_manual(values = c(
    "#F4DBB5",
    "#E8B878",
    "#CF8552",  
    "#B45035",  
    "#8C6300",  
    "#654500",  
    "#8DBBDC",  
    "#1f78b4",   
    "#376491"),
    na.translate = FALSE) +
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # only gridlines every 100 nm
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1)) +    # only gridlines every 0.1
  labs(
    title = "Rugulopteryx okamurae - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.5)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

p

## Export graph
ggsave(
  filename = "GAM_Rugulopteryx_insitu_BAM_50.tiff",
  plot = p,
  width = 730/100, height = 630/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)


##-----------------------
# GAMs in PlanetScope windows

# raw_data$MeanReflectance <- rowMeans(raw_data[-1], na.rm = TRUE)

planetScope <- NewData %>%
  mutate(
    wv_planet = as.factor(case_when(
      (Wavelength > 430 & Wavelength < 453) ~ '443',
      (Wavelength > 464 & Wavelength < 516) ~ '490',
      (Wavelength > 512 & Wavelength < 550) ~ '531',
      (Wavelength > 546 & Wavelength < 584) ~ '565',
      (Wavelength > 599 & Wavelength < 621) ~ '610',
      (Wavelength > 649 & Wavelength < 681) ~ '665',
      (Wavelength > 696 & Wavelength < 712) ~ '705',
      (Wavelength > 844 & Wavelength < 886) ~ '865'
    ))
  ) %>% 
  filter(!is.na(wv_planet))


# Load necessary libraries
library(dplyr)

# Example of complete set of wavelengths (modify this as needed)
full_wavelengths <- seq(400, 900, by = 1)  # Assuming the range is from 400 to 900 nm

# Group by 'group' and 'wv_planet' and calculate mean_reflec and mean_se.fit
summarized_data <- planetScope %>%
  group_by(group, wv_planet) %>%
  summarise(
    mean_reflec = mean(response, na.rm = TRUE),  
    mean_se.fit = mean(se.fit, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  # Create a complete set of all combinations of 'group' and 'wv_planet'
  complete(group, wv_planet = as.factor(full_wavelengths), fill = list(mean_reflec = NA, mean_se.fit = NA)) %>%
  filter(!is.na(group))  # Remove rows where group is NA


## plot lines and violins
# Plot using ggplot2 with custom palette
ps <- ggplot(summarized_data, aes( x=as.numeric(as.character(wv_planet)), y = mean_reflec)) +
  geom_violin(aes(fill = as.factor(wv_planet)), alpha = 0.6) +
  scale_fill_manual(
    values = c(
    '443' = '#98CDD5',   # Coastal-Blue
    '490' = '#39909E',   # Blue
    '531' = '#C2D7AD',   # Green_i
    '565' = '#6C9E39',   # Green
    '610' = '#e2c800ff', # Yellow
    '665' = '#AD4C3E',   # Red
    '705' = '#9C9263',   # Red-edge
    '865' = '#CC7349'    # NIR
  ),
  labels = c(
    '443' = 'Coastal-Blue',
    '490' = 'Blue',
    '531' = 'Green_i',
    '565' = 'Green',
    '610' = 'Yellow',
    '665' = 'Red',
    '705' = 'Red-edge',
    '865' = 'NIR'
  ),
  na.translate = FALSE
  ) + # Removes NA from color legendna.translate = FALSE) +
  geom_jitter(width = 0.2, size = 2, color = 'grey', alpha = 0.7) +
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # X-axis range and breaks
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1)) +    # Y-axis range and breaks
  labs(
    title = "Rugulopteryx okamurae - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.5)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

ps

ggsave(
  filename = "PS_Rugulopteryx_BAM.tiff",
  plot = ps,
  width = 730/100, height = 630/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)


## halfeye
halfeye_ps <- ggplot(summarized_data, aes( x=as.numeric(as.character(wv_planet)), y = mean_reflec)) +
  stat_halfeye(aes(fill = as.factor(wv_planet)), alpha = 0.9, justification = -0.12,
              position = position_dodge(0.7)
  ) +
  scale_fill_manual(
    values = c(
      '443' = '#98CDD5',   # Coastal-Blue
      '490' = '#39909E',   # Blue
      '531' = '#C2D7AD',   # Green_i
      '565' = '#6C9E39',   # Green
      '610' = '#e2c800ff', # Yellow
      '665' = '#AD4C3E',   # Red
      '705' = '#9C9263',   # Red-edge
      '865' = '#CC7349'    # NIR
    ),
    labels = c(
      '443' = 'Coastal-Blue',
      '490' = 'Blue',
      '531' = 'Green_i',
      '565' = 'Green',
      '610' = 'Yellow',
      '665' = 'Red',
      '705' = 'Red-edge',
      '865' = 'NIR'
    ),
    na.translate = FALSE
  ) + # Removes NA from color legendna.translate = FALSE) +
  geom_jitter(width = 0.2, size = 1, color = 'grey', alpha = 0.7) +
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # X-axis range and breaks
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1)) +    # Y-axis range and breaks
  labs(
    title = "Rugulopteryx okamurae - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    #legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    #legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

halfeye_ps

ggsave(
  filename = "PS_Rugulopteryx_GAM.tiff",
  plot = halfeye_ps,
  width = 730/100, height = 450/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)


##----------------------------
# GAMs in Sentinel-2 windows

## To include the 842 overlapping band:
# First, filter and mutate for the '842' range (Wavelength > 784 & Wavelength < 901)
s2_842 <- NewData %>%
  filter(Wavelength > 784 & Wavelength < 901) %>%
  mutate(wv_s2 = '842')  # Assign '842' directly

# Then, filter and mutate for the '865' range (Wavelength > 854 & Wavelength < 875)
s2_865 <- NewData %>%
  filter(Wavelength > 854 & Wavelength < 875) %>%
  mutate(wv_s2 = '865')  # Assign '865' directly

# Combine both datasets
s2_dataset <- bind_rows(s2_842, s2_865)

# Apply the remaining 'case_when' for the other wavelength ranges
s2_dataset <- s2_dataset %>%
  bind_rows(
    NewData %>%
      filter(!(Wavelength > 784 & Wavelength < 901) & !(Wavelength > 854 & Wavelength < 875)) %>%
      mutate(wv_s2 = case_when(
        (Wavelength > 432 & Wavelength < 454) ~ '443',
        (Wavelength > 457 & Wavelength < 524) ~ '490',
        (Wavelength > 542 & Wavelength < 579) ~ '560',
        (Wavelength > 649 & Wavelength < 681) ~ '665',
        (Wavelength > 697 & Wavelength < 714) ~ '705',
        (Wavelength > 732 & Wavelength < 749) ~ '740',
        (Wavelength > 772 & Wavelength < 794) ~ '783'
      )) %>% 
      filter(!is.na(wv_s2))
  )

library(dplyr)

# Example of complete set of wavelengths (modify this as needed)
full_wavelengths <- seq(400, 900, by = 1)  # Assuming the range is from 400 to 900 nm

# Group by 'group' and 'wv_planet' and calculate mean_reflec and mean_se.fit
summarized_data <- s2_dataset %>%
  group_by(group, wv_s2) %>%
  summarise(
    mean_reflec = mean(response, na.rm = TRUE),  
    mean_se.fit = mean(se.fit, na.rm = TRUE), 
    .groups = "drop"
  ) %>%
  # Create a complete set of all combinations of 'group' and 'wv_planet'
  complete(group, wv_s2 = as.factor(full_wavelengths), fill = list(mean_reflec = NA, mean_se.fit = NA)) %>%
  filter(!is.na(group))  # Remove rows where group is NA


ggplot(summarized_data, aes(x = as.numeric(wv_s2), y = mean_reflec)) +
  geom_smooth(aes(color = group)) + # Jittered points for better visibility
  labs(x = "Wavelength (nm)", y = "Reflectance", title = "Violin Plots for Different Wavelengths") +
  theme_minimal() +
  #geom_violin(inherit.aes = F, data = x, aes(x = as.numeric(wv_planet), y = mean_reflec, group = group))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for readability


## plot lines and violins
# Plot using ggplot2 with custom palette
s2 <- ggplot(summarized_data, aes( x=as.numeric(as.character(wv_s2)), y = mean_reflec)) +
  geom_violin(aes(fill = as.factor(wv_s2)), alpha = 0.6) +
  scale_fill_manual(
    values = c(
      '443'= '#98CDD5',
      '490'= '#39909E',
      '560'= '#6C9E39',
      '665'= '#AD4C3E',
      '705'= '#9C9263',
      '740'= '#E3CA4B',
      '783'= '#E6A941',
      '842'= '#E6CE94',
      '865'= '#CC7349'
    ),

    labels = c(
      '443'= 'Coastal-Blue',
      '490'= 'Blue',
      '560'= 'Green',
      '665'= 'Red',
      '705'= 'Red-Edge 1',
      '740'= 'Red-Edge 2',
      '783'= 'Red-Edge 3',
      '842'= 'NIR',
      '865'= 'narrow NIR'
    ),
    na.translate = FALSE
  ) + # Removes NA from color legendna.translate = FALSE) +
  geom_jitter(width = 0.2, size = 2, color = 'grey', alpha = 0.7) +
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # X-axis range and breaks
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1)) +    # Y-axis range and breaks
  labs(
    title = "Rugulopteryx okamurae - S2",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.5)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

s2

ggsave(
  filename = "S2_Rugulopteryx_continuous.tiff",
  plot = ps,
  width = 730/100, height = 630/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)

library(ggdist)
## halfeye
halfeye_s2 <- ggplot(summarized_data, aes( x=as.numeric(as.character(wv_s2)), y = mean_reflec)) +
  stat_halfeye(aes(fill = as.factor(wv_s2)), alpha = 0.9, justification = -0.12,
               position = position_dodge(0.7)
  ) +
  scale_fill_manual(
    values = c(
      '443'= '#98CDD5',
      '490'= '#39909E',
      '560'= '#6C9E39',
      '665'= '#AD4C3E',
      '705'= '#9C9263',
      '740'= '#E3CA4B',
      '783'= '#E6A941',
      '842'= '#E6CE94',
      '865'= '#CC7349'
    ),
    
    labels = c(
      '443'= 'Coastal-Blue',
      '490'= 'Blue',
      '560'= 'Green',
      '665'= 'Red',
      '705'= 'Red-Edge 1',
      '740'= 'Red-Edge 2',
      '783'= 'Red-Edge',
      '842'= 'NIR',
      '865'= 'narrow NIR'
    ),
    na.translate = FALSE
  ) + # Removes NA from color legenda.translate = FALSE) +
  geom_jitter(width = 0.2, size = 1, color = 'grey', alpha = 0.7) +
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # X-axis range and breaks
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1)) +    # Y-axis range and breaks
  labs(
    title = "Rugulopteryx okamurae - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.5)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    #legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    #legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

halfeye_s2

ggsave(
  filename = "S2_Rugulopteryx_BAM.tiff",
  plot = halfeye_s2,
  width = 730/100, height = 450/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)








##-------------------
# in situ signatures

# Initialize result data frame
summary_df <- data.frame(
  Wavelength = raw_data$Wavelength
)


# Compute mean and sd for each DZ group
for (group in names(dz_ranges)) {
  cols <- dz_ranges[[group]]
  summary_df[[paste0(group, "_mean")]] <- rowMeans(raw_data[, cols], na.rm = TRUE)
  summary_df[[paste0(group, "_sd")]] <- apply(raw_data[, cols], 1, sd, na.rm = TRUE)
}

library(ggplot2)
library(tidyr)
library(dplyr)


# Pivot data to long format for ggplot
long_data <- summary_df %>%
  pivot_longer(
    cols = -Wavelength,
    names_to = c("Group", "Stat"),
    names_sep = "_",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Stat,
    values_from = Value
  ) %>%
  mutate(
    ymin = mean - sd,
    ymax = mean + sd
  )


# library(paletteer)

# Extract 10 colors from BrwnYl palette
# brwnyl_palette <- rev(paletteer_c("grDevices::BrwnYl", 10))

manual_colors <- c(
  "#F8EDBF",
  "#d7baa3",
  "#A36B2B",  
  "#E0C9AB",  
  "#DAA66D",  
  "#C2902A",  
  "#8C6300",  
  "#654500",   
  "#4981BF"  
  )


# Ensure Group is a factor with consistent order
#long_data$Group <- factor(long_data$Group, levels = unique(long_data$Group))

# Plot using ggplot2 with custom palette
ggplot(long_data, aes(x = Wavelength, y = mean, color = Group, fill = Group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  labs(
    title = "Rugulopteryx okamurae - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  xlim(400, 900) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = manual_colors) +
  scale_fill_manual(values = manual_colors)

