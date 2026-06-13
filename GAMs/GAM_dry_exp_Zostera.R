# Load necessary library
library(tidyverse)
library(ggplot2)

# Read the file (adjust separator if needed: "\t" for tab, "" for whitespace)
raw_data <- read.delim("../in situ radiometry/Experiment/zostera/20250612_dry_exp_zostera.txt", check.names = FALSE, stringsAsFactors = FALSE)

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
         group=as.factor(case_when(id%in%c(2:61)~'0',
                                   id%in%c(62:121)~'1',
                                   id%in%c(122:181)~'2',
                                   id%in%c(182:241)~'3',
                                   id%in%c(242:301)~'4',
                                   id%in%c(302:361)~'5',
                                   id%in%c(362:421)~'6',
                                   id%in%c(422:481)~'7',
                                   id%in%c(482:541)~'21',
                                   id%in%c(542:601)~'32'))) %>% 
  filter(Wavelength>399&Wavelength<901) %>% 
  filter(!is.na(group))


ggplot(raw_data_long, aes(x=Wavelength, y= reflec, color=group)) + 
  geom_smooth(se = TRUE)


library(mgcv)
##gam_1 <- gam(reflec~s(Wavelength, by=group, k=30)+group, family=Gamma(),data=raw_data_long, method = "REML", control = gam.control(trace = TRUE)) #gaussian() was the other option.
library(mgcv)

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
  mutate(group=as.factor(case_when(id%in%c(2:61)~'0',
                                   id%in%c(62:121)~'1',
                                   id%in%c(122:181)~'2',
                                   id%in%c(182:241)~'3',
                                   id%in%c(242:301)~'4',
                                   id%in%c(302:361)~'5',
                                   id%in%c(362:421)~'6',
                                   id%in%c(422:481)~'7',
                                   id%in%c(482:541)~'21',
                                   id%in%c(542:601)~'32')))


Pred<-predict(gam_1,NewData_1,se.fit=T,type="response") # creates a list of predicted mean value y'


NewData<-NewData_1 %>% 
  mutate(response=Pred$fit,
         se.fit=Pred$se.fit,
         Upr=response+(se.fit*1.96),
         Lwr=response-(se.fit*1.96))

library(paletteer)

NewData$group <- factor(NewData$group, levels = 0:32)
# Extract 10 colors from BrwnYl palette (since you have 10 groups)
brwnyl_palette <- (paletteer_c("grDevices::Green-Yellow", 10))

# Plot using ggplot2 with custom palette
p <- ggplot(NewData) +
  geom_ribbon(aes(x = Wavelength,
                  ymax = Upr,
                  ymin = Lwr,
                  fill = as.factor(group)), alpha = 0.2) +
  geom_line(aes(x = Wavelength,
                y = response,
                colour = as.factor(group))) + 
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # only gridlines every 100 nm
  scale_y_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, by = 0.1)) +    # only gridlines every 0.1
  labs(
    title = "Zostera noltii - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.4)) +
  scale_fill_manual(values = brwnyl_palette) +  # Apply palette to the fill aesthetic
  scale_color_manual(values = brwnyl_palette) +# Apply palette to the color aesthetic
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Print the plot
print(p)

## Export graph
ggsave(
  filename = "GAM_Zostera_exp_BAM_50.tiff",
  plot = p,
  width = 730/100, height = 630/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)



-------------
## Summarise by wavelength groups
## By groups

wv_groups <- NewData %>%
  mutate(
    wv_group = as.factor(case_when(
      (Wavelength > 399 & Wavelength < 449) ~ 'wv_450',
      (Wavelength > 449 & Wavelength < 499) ~ 'wv_500',
      (Wavelength > 499 & Wavelength < 549) ~ 'wv_550',
      (Wavelength > 549 & Wavelength < 599) ~ 'wv_600',
      (Wavelength > 599 & Wavelength < 649) ~ 'wv_650',
      (Wavelength > 649 & Wavelength < 699) ~ 'wv_700',
      (Wavelength > 699 & Wavelength < 749) ~ 'wv_750',
      (Wavelength > 749 & Wavelength < 799) ~ 'wv_800',
      (Wavelength > 799 & Wavelength < 849) ~ 'wv_850',
      (Wavelength > 849 & Wavelength < 901) ~ 'wv_900'
    ))
  )

## Group by times and wavelength groups
summarized_data <- wv_groups %>%
  group_by(group, wv_group) %>%
  summarise(
    mean_reflec = mean(response, na.rm = TRUE),  
    mean_se.fit = mean(se.fit, na.rm = TRUE), 
    .groups = "drop" 
  )


# Pivot the summarized data to wide format
wide_data <- summarized_data %>%
  pivot_wider(
    names_from = group, 
    values_from = c(mean_reflec, mean_se.fit),
    names_glue = "{group}_{.value}"
  )


grouped_50<- wide_data %>%
  pivot_longer(
    cols = -c(wv_group,`0_mean_reflec`,`0_mean_se.fit`), 
    names_to = c("group", ".value"),
    names_pattern = "(\\d+)_(.*)" 
  ) %>% 
  mutate(change_reflec = `0_mean_reflec`- mean_reflec,
         change_se = sqrt((`0_mean_se.fit`)^2+(mean_se.fit)^2)
  )


# Plot the change over time with error bars
ggplot(grouped_50, aes(x = as.numeric(group), y = change_reflec, color = wv_group)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = change_reflec - change_se, ymax = change_reflec + change_se), width = 0.1) +
  labs(title = "Change in Reflectance Over Time Compared to Time 0",
       x = "Time (hours)", y = "Reflectance Change (Time - 0)") +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "none")


-------------
  ## Summarise by continuous wavelength
  ## continuous, smooth

## Group by times and wavelength groups
summarized_data <- NewData %>%
  group_by(group, Wavelength) %>%
  summarise(
    mean_reflec = mean(response, na.rm = TRUE),  
    mean_se.fit = mean(se.fit, na.rm = TRUE), 
    .groups = "drop" 
  )


# Pivot the summarized data to wide format
wide_data <- summarized_data %>%
  pivot_wider(
    names_from = group, 
    values_from = c(mean_reflec, mean_se.fit),
    names_glue = "{group}_{.value}"
  )


pivot_continuous<- wide_data %>%
  pivot_longer(
    cols = -c(Wavelength,`0_mean_reflec`,`0_mean_se.fit`), 
    names_to = c("group", ".value"),
    names_pattern = "(\\d+)_(.*)" 
  ) %>% 
  mutate(change_reflec = `0_mean_reflec`- mean_reflec,
         change_se = sqrt((`0_mean_se.fit`)^2+(mean_se.fit)^2)
  )


---------

library(ggplot2)
library(paletteer)

# Ensure 'group' is a factor and set its levels in the desired order
pivot_continuous$group <- factor(pivot_continuous$group, levels = as.character(0:32))

# Check the unique groups and their order
unique_groups <- unique(pivot_continuous$group)
print(unique_groups)

# Generate the color palette dynamically based on the number of unique groups
brwnyl_palette <- paletteer_c("grDevices::Green-Yellow", length(unique_groups))

# Multiply bu -1 as change_reflec has values flipped
pivot_continuous$change_reflec <- pivot_continuous$change_reflec*(-1)
# Convert `0_mean_reflec` to numeric if needed
pivot_continuous$`0_mean_reflec` <- as.numeric(pivot_continuous$`0_mean_reflec`)


# Plot using ggplot2 with the custom color palette and ordered legend
p <- ggplot(pivot_continuous) +
  geom_ribbon(aes(x = Wavelength,
                  ymax = change_reflec + change_se,
                  ymin = change_reflec - change_se,
                  fill = factor(group))) +  # Use factor(group) for fill
  geom_line(aes(x = Wavelength,
                y = change_reflec,
                colour = factor(group))) +
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # only gridlines every 100 nm
  scale_y_continuous(
    limits = c(-0.06, 0.10),
    breaks = seq(-0.06, 0.10, by = 0.02)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # only gridlines every 0.005
  labs(
    title = "Zostera noltii - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1})~"changes over time")
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(-0.06, 0.10))+
  scale_fill_manual(values = brwnyl_palette) +  # Apply the custom palette to the fill aesthetic
  scale_color_manual(values = brwnyl_palette) +  # Apply the custom palette to the color aesthetic
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Print the plot
print(p)

## Export graph
ggsave(
  filename = "GAM_Zostera_exp_BAM_50.tiff",
  plot = p,
  width = 730/100, height = 630/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)

  
  
  
-------
  
library(ggplot2)
library(paletteer)

# Check how many unique groups you have
unique_groups <- unique(pivot_continuous$group)

# Generate the color palette dynamically based on the number of unique groups
brwnyl_palette <- paletteer_c("grDevices::Green-Yellow", length(unique_groups))

# Plot using ggplot2 with the custom color palette
p <- ggplot(pivot_continuous) +
  # Use `group` as a factor for both fill and color
  geom_ribbon(aes(x = Wavelength,
                  ymax = change_reflec + change_se,
                  ymin = change_reflec - change_se,
                  fill = as.factor(group)), alpha = 0.2) +  # Use factor(group) for fill
  geom_line(aes(x = Wavelength,
                y = change_reflec,
                colour = as.factor(group))) +
  geom_ribbon(aes(x = Wavelength,
                  ymax = `0_mean_reflec`+ `0_mean_se.fit`,
                  ymin = `0_mean_reflec` - `0_mean_se.fit`)) +
  geom_line(aes(x = Wavelength,
                y = `0_mean_reflec`))
  scale_x_continuous(
    limits = c(400, 900),
    breaks = seq(400, 900, by = 100)) +  # only gridlines every 100 nm
  scale_y_continuous(
    limits = c(-0.1, 0.1),
    breaks = seq(-0.1, 0.1, by = 0.05)) +    # only gridlines every 0.005
  labs(
    title = "Zostera noltei - Water Content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(-0.1, 0.1)) +
  scale_fill_manual(values = brwnyl_palette) +  # Apply the custom palette to the fill aesthetic
  scale_color_manual(values = brwnyl_palette) +  # Apply the custom palette to the color aesthetic
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

# Print the plot
print(p)
