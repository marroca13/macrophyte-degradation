# Load necessary library
library(tidyverse)

# Read the file (adjust separator if needed: "\t" for tab, "" for whitespace)
raw_data <- read.delim("../in situ radiometry/Experiment/rugulopteryx/20250612_dry_exp_okamurae.txt", check.names = FALSE, stringsAsFactors = FALSE)

# Convert all comma decimals to dots (except the first column with 'Wavelength')
raw_data[-1] <- lapply(raw_data[-1], function(x) as.numeric(gsub(",", ".", x)))

# Make sure Wavelength is numeric
raw_data$Wavelength <- as.numeric(raw_data$Wavelength)

raw_data$MeanReflectance <- rowMeans(raw_data[-1], na.rm = TRUE)

## As there are additional data that shouldn't be in the txt file, we clean first
library(dplyr)

raw_data <- raw_data %>%
  select(-starts_with("FR"), -starts_with("RA"))


# Define column ranges for each RO group
dz_ranges <- list(
  RO0 = 62:121, # here we have 20 measurements more, 80 in total
  RO1 = 142:201, # this one is doubled but consistent, 120 measurements
  RO2 = 202:261,
  RO3 = 262:321,
  RO4 = 2:61, # first rows are RO4
  RO5 = 322:381, # skipped when we realized R1 was doubled
  RO6 = 382:441,
  RO7 = 442:501,
  RO8 = 502:561,
  RO9 = 562:621
)

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

library(paletteer)

# Extract 10 colors from BrwnYl palette (since you have 10 groups)
brwnyl_palette <- (paletteer_c("grDevices::BrwnYl", 10))

p <- ggplot(long_data, aes(x = Wavelength, y = mean, color = Group, fill = Group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.05, color = NA) +
  labs(
    title = "Rugulopteryx okamurae - radiometry/water content",
    x = "Wavelength (nm)",
    y = expression(Rrs~(sr^{-1}))
  ) +
  coord_cartesian(xlim = c(400, 900), ylim = c(0, 0.4)) +
  theme_minimal(base_family = "roboto") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),      
    legend.text = element_text(size = 14),     
    legend.title = element_blank()
  ) +
  scale_color_manual(values = brwnyl_palette) +
  scale_fill_manual(values = brwnyl_palette) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.4)))

p

## Export graph
ggsave(
  filename = "SM2_Rugu_exp_raw.tiff",
  plot = p,
  width = 730/100, height = 630/100,   # Convert pixels to inches (1 inch = 100px for 300dpi)
  dpi = 300
)