# Plot the phenotypic distributions and correlations of AFB and AFS in the UKB, by cohort

library(haven)
library(dplyr)
library(cowplot)
library(gganimate)
library(ggridges)
library(viridis)
library(forcats)

df <- read_dta("nicola_ukb_copy/ukbiobank.dta")
df_sel <- select(df, eid = n_eid, sex = n_31_0_0, birth_year = n_34_0_0, afb = n_2754_0_0, afs = n_2139_0_0)

# If missing codes are present, set to missing
df_sel$afb[df_sel$afb < 0] <- NA
df_sel$afs[df_sel$afs < 0] <- NA

df_sel <- mutate(
    df_sel,
    old = birth_year < 1950,
    cohort = cut(
        birth_year,
        breaks = c(1933, 1940, 1945, 1950, 1955, 1960, 1971),
        labels = c("<1941", "1941-1945", "1946-1950", "1951-1955", "1956-1960", ">1960")
        )
    )

# Get rid of Stata labels
df_sel <- zap_labels(df_sel)

# Select only women because afb was asked only of women
df_sel <- filter(df_sel, sex == 0)

df_sel_no_na <- na.omit(df_sel)

afb_anim <- ggplot(df_sel_no_na, aes(afb)) +
    geom_density(bw = 0.55) +
    transition_states(cohort) +
    labs(
        title = "Birth year: {closest_state}"
    ) +
    xlab("Age at first birth") +
    shadow_trail(color = "grey")

anim_save("figs/afb_cohort_anim.gif", afb_anim)

afs_anim <- ggplot(df_sel_no_na, aes(afs)) +
    geom_density(bw = 0.55) +
    transition_states(cohort) +
    labs(
        title = "Birth year: {closest_state}"
    ) +
    xlab("Age at first sexual intercourse") +
    shadow_trail(color = "grey")

anim_save("figs/afs_cohort_anim.gif", afs_anim)

afb_static <- ggplot(df_sel, aes(afb, cohort, fill = ..x..)) +
    geom_density_ridges_gradient(rel_min_height = 0.01, bandwidth = 0.55) +
    xlab("Age at first birth") +
    ylab("Birth year") +
    scale_fill_viridis(limits = c(10, 40), guide = F) +
    coord_cartesian(xlim = c(10, 45))

afs_static <- ggplot(df_sel, aes(afs, cohort, fill = ..x..)) +
    geom_density_ridges_gradient(rel_min_height = 0.01, bandwidth = 0.55) +
    xlab("Age at first sexual intercourse") +
    ylab("Birth year") +
    scale_fill_viridis("Age", limits = c(10, 40), guide = guide_colorbar(barheight = unit(3, "inches"))) +
    coord_cartesian(xlim = c(10, 45)) +
    theme(legend.position = c(0.8, 0.6))

static_plts <- plot_grid(afb_static, afs_static, labels = c("A", "B"))

save_plot("figs/afb_afs_ridge_plts.pdf",
    static_plts,
    nrow = 2,
    base_aspect_ratio = 2.3)

cors_df <- mutate(df_sel_no_na, cohort = fct_relabel(cohort, ~ paste0("Birth year: ", .))) %>%
    group_by(cohort) %>%
    summarize(corr = round(cor(afs, afb), digits = 2))

cors_plt <- mutate(df_sel_no_na, cohort = fct_relabel(cohort, ~ paste0("Birth year: ", .))) %>%
    ggplot(aes(afs, afb)) +
    geom_abline(color = "grey") +
    geom_count(aes(color = ..prop.., size = ..prop..)) +
    geom_text(data = cors_df, aes(30, 12, label = paste0("Pearson's r: ", corr))) +
    facet_wrap(~cohort) +
    scale_color_viridis("Proportion", option = "C", guide = guide_colorbar(barheight = unit(2, "inches"))) +
    scale_size_area("Proportion", max_size = 4) +
    coord_cartesian(xlim = c(12, 40), ylim = c(12, 45)) +
    xlab("Age at first sexual intercourse") +
    ylab("Age at first birth") +
    theme(strip.background = element_rect(fill = "#DCDCDC"))

ggsave("figs/afb_afs_corr_plt.pdf", cors_plt, width = 15, height = 10)

cors_anim <- mutate(df_sel_no_na, cohort = fct_relabel(cohort, ~ paste0("Birth year: ", .))) %>%
    ggplot(aes(afs, afb)) +
    geom_abline(color = "grey") +
    geom_count(aes(color = ..prop.., size = ..prop.., group = cohort)) +
    geom_text(data = cors_df, aes(30, 12, label = paste0("Pearson's r: ", corr))) +
    scale_color_viridis("Proportion", option = "C", guide = guide_colorbar(barheight = unit(2, "inches"))) +
    scale_size_area("Proportion", max_size = 4) +
    coord_cartesian(xlim = c(12, 40), ylim = c(12, 45)) +
    xlab("Age at first sexual intercourse") +
    ylab("Age at first birth") +
    labs(title = "{closest_state}") +
    transition_states(cohort, state_length = 3) +
    enter_fade() +
    exit_fade()

anim_save("figs/afb_afs_corr_anim.gif", cors_anim, nframes = 200)

# Animate year by year
cors_by_df <- filter(df_sel_no_na, !birth_year %in% c(1936, 1937, 1938, 1970)) %>%
    mutate(birth_year = fct_drop(as.character(birth_year))) %>%
    group_by(birth_year) %>%
    summarize(corr = round(cor(afs, afb), digits = 2))

# Remove years with few subjects because they mess up the scale for the proportions
# This is, admittedly, quick and dirty
cors_by_anim <- filter(df_sel_no_na, !birth_year %in% c(1936, 1937, 1938, 1970)) %>%
    mutate(birth_year = fct_drop(as.character(birth_year))) %>%
    ggplot(aes(afs, afb)) +
    geom_abline(color = "grey") +
    geom_count(aes(color = ..prop.., size = ..prop.., group = birth_year)) +
    geom_text(data = cors_by_df, aes(30, 12, label = paste0("Pearson's r: ", corr))) +
    scale_color_viridis("Proportion", option = "C", guide = guide_colorbar(barheight = unit(2, "inches"))) +
    scale_size_area("Proportion", max_size = 4) +
    coord_cartesian(xlim = c(12, 40), ylim = c(12, 45)) +
    xlab("Age at first sexual intercourse") +
    ylab("Age at first birth") +
    labs(title = "Birth year: {closest_state}") +
    transition_states(birth_year)

anim_save("figs/afb_afs_corr_anim_by_year.gif", cors_by_anim, nframes = 100)
