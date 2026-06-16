library(data.table)
library(dplyr)
library(ggplot2)
library(ggtext)
library(patchwork)
library(scales)
# set the working directory
setwd("D:/MaxTreeHeightProject")

# 1. Load megatree node information and extract family-level divergence time.
mega_divergence <- fread("D:/R code/U.PhyloMaker建系统发育树/Datafiles_Megatrees&Genus-family relationship/plant_20221214-main/nodes.info.1.WCVP.csv")

family_divtime <- mega_divergence %>%
  filter(level == "F") %>%
  dplyr::select(family, divergence = rn.bl)

# 2. Load species-level maximum tree height data and merge family-level divergence time.
wood_height_species <- fread("Data/MaxTreeheight_20260410_for_single_trees_filter_by_meanmaxHeight_lonlatrange.csv") %>%
  dplyr::select(species, q975Height, Family, Genus, ORDER, GROUP) %>%
  rename(family = Family, genus = Genus, maxHeight = q975Height, group = GROUP) %>%
  left_join(family_divtime, by = "family") %>%
  filter(divergence <= 200)


# 3. Summarize maximum tree height and species richness at the family level.
wood_height_family <- wood_height_species %>%
  group_by(family) %>%
  summarise(
    speciesNum = n_distinct(species),
    group = first(na.omit(group)),
    max_height = max(maxHeight, na.rm = TRUE),
    mean_height = mean(maxHeight, na.rm = TRUE),
    divergence = first(na.omit(divergence)),
    .groups = "drop"
  ) 

# 4. Fit linear models between family-level divergence time and maximum tree height.
model_all <- lm(mean_height ~ divergence, data = wood_height_family)
summary(model_all)

model_angiosperms <- lm(mean_height ~ divergence, data = filter(wood_height_family, group == "Angiosperms"))
summary(model_angiosperms)

model_gymnosperms <- lm(mean_height ~ divergence, data = filter(wood_height_family, group == "Gymnosperms"))
summary(model_gymnosperms)

# 5. Extract regression statistics for each plant group.
stats_df <- wood_height_family %>%
  filter(group %in% c("Angiosperms", "Gymnosperms")) %>%
  group_by(group) %>%
  do({
    fit <- lm(mean_height ~ divergence, data = .)
    fit_sum <- summary(fit)
    data.frame(slope = coef(fit)[2], r2 = fit_sum$r.squared, p = coef(fit_sum)[2, 4], n = nobs(fit))
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("<i>Slope</i> = ", round(slope, 2), "<br><i>P</i> ", ifelse(p < 0.001, "&lt; 0.001", paste0("= ", signif(p, 2)))),
    x = max(wood_height_family$divergence, na.rm = TRUE) * 0.75,
    y = ifelse(group == "Angiosperms", max(wood_height_family$mean_height, na.rm = TRUE) * 0.92, max(wood_height_family$mean_height, na.rm = TRUE) * 1.05)
  )

# 6. Plot the relationship between family divergence time and maximum tree height.
group_cols <- c("Angiosperms" = "#6baed6", "Gymnosperms" = "#fdae6b")

group_cols <- c(
  "Angiosperms" = "#4C78A8",
  "Gymnosperms" = "#E69542"
)

p1 <- ggplot(wood_height_family, aes(x = divergence, y = mean_height, group = group)) +
  geom_point(aes(size = speciesNum, fill = group), shape = 21, color = "black", alpha = 0.8) +
  geom_smooth(aes(color = group, linetype = group), method = "lm", se = FALSE) +
  geom_richtext(data = stats_df, aes(x = x, y = y, label = label, color = group), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4.5, fill = NA, label.color = NA, show.legend = FALSE) +
  scale_size_continuous(range = c(3, 8)) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  scale_linetype_manual(values = c("Angiosperms" = 1, "Gymnosperms" = 1)) +
  labs(x = "Divergence time (Ma)", y = "Mean maximum tree height (m)") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    legend.justification = c("left", "top"),
    legend.position = "inside",
    legend.position.inside = c(0.01, 0.99),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  ) +
  guides(size = "none")
p1

# 7. Calculate the percentage of angiosperms and gymnosperms across divergence-time bins.
woodHeight_family1 <- fread("Data/MaxTreeheight_20260410_for_single_trees_filter_by_meanmaxHeight_lonlatrange_gt5.csv") %>%
  dplyr::select(species, q975Height, Family, Genus, ORDER, GROUP) %>%
  rename(family = Family, genus = Genus, maxHeight = q975Height, group = GROUP) %>%
  left_join(family_divtime, by = "family")%>%
  filter(divergence <= 200)
max(woodHeight_family1$divergence, na.rm = TRUE)
min(woodHeight_family1$divergence, na.rm = TRUE)
woodHeight_family2 <- woodHeight_family1 %>%
  filter(!is.na(divergence)) %>% 
  mutate(
    div_bin = cut(
      divergence,
      breaks = c(20, 30, 50,70, 100, 150),
      labels = c("20-30", "30-50", "50-70","70-100", ">100"),
      right = FALSE,       
      include.lowest = TRUE
    )
  )
div_summary <- woodHeight_family2 %>%
  group_by(div_bin, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(div_bin) %>%
  mutate(percentage = n / sum(n) * 100)

# 8. Plot the proportional composition of angiosperms and gymnosperms across divergence-time bins.
p2<-ggplot(div_summary, aes(x = div_bin, y = percentage, fill = group)) +
  geom_bar(stat = "identity", position = "fill") +   # position="fill" 自动百分比
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(
    values = c("Angiosperms" = "#6baed6","Gymnosperms" = "#fdae6b")) +
  labs(x = "Divergence time (Ma)",#y = "Percentage of angiosperms and gymnosperms (%)",
       y = "Percentage of angiosperms and gymnosperms (%)",
       fill = "Group") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        # legend.position = "none",
        axis.text.y = element_text(angle = 90),
        legend.justification= c("left", "top"), # 图例锚点对齐方式
        legend.position = "inside",
        legend.position.inside = c(0.026, 0.95),        # 图内位置 (x, y)，范围 0-1
        legend.title = element_blank())
p2

# 9. Combine and export the final figure.
combined_plot <- (p1 + p2) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22))
combined_plot
ggsave("Plots/Figure_02_Divergence_Famliy_TreeHeight_combined_PercentageAG.pdf", combined_plot, width = 12, height = 6, units = "in")
ggsave("Plots/Figure_02_Divergence_Famliy_TreeHeight_combined_PercentageAG.jpeg",combined_plot,width = 12,height = 6,units = "in",dpi = 300)


