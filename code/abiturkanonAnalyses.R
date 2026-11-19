####################### CHR2027: German School Canon ###########################
##
## Analysis code for "Putting the Canon under a Microscope: Canonization
## Processes in German School Curricula" (CHR 2027).
##
################################################################################

#------------------------------------------------------------------------------#
## 1 Setup --------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# Set this to the folder where you cloned the repository
project_dir <- "path/to/abiturkanon"

#------------------------------------------------------------------------------#
## 1.1 Packages ---------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(tidyverse)
library(ggridges)
library(ggthemes)
library(igraph)
library(geodata)
library(sf)
library(survival)

#------------------------------------------------------------------------------#
## 1.2 Visual language --------------------------------------------------------#
#------------------------------------------------------------------------------#

theme_chr <- function(base_size = 13, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      line =              element_line(colour = '#DADADA', size = 1,
                                       linetype = 1, lineend = "butt"),
      rect =              element_rect(fill = "transparent", colour = NA,
                                       size = 0.5, linetype = 1),
      text =              element_text(family = base_family, face = "plain",
                                       colour = "#656565", size = base_size,
                                       hjust = 0.5, vjust = 0.5, angle = 0,
                                       lineheight = 0.9),
      plot.title =        element_text(size = rel(1.5), family = '',
                                       face = 'bold', hjust = -0.05,
                                       vjust = 1.5, colour = '#3B3B3B'),
      axis.text =         element_text(),
      axis.ticks =        element_line(),
      panel.grid.major =  element_line(),
      panel.grid.minor =  element_blank(),
      plot.background =   element_rect(),
      panel.background =  element_rect(),
      legend.key =        element_rect(colour = '#DADADA'),
      legend.background = element_rect(fill = "transparent", colour = NA),
      panel.border =      element_blank(),
      axis.line =         element_blank(),
      strip.background =  element_rect(fill = "transparent", colour = NA),
      legend.position = 'none',
      complete = TRUE
    )
}

#------------------------------------------------------------------------------#
## 1.3 Paths, data import and constants ---------------------------------------#
#------------------------------------------------------------------------------#

setwd(project_dir)

dnb <- read.csv(
  "data\\260727_dnb_results_manuallyFiltered.csv",
  sep = ";", stringsAsFactors = FALSE, encoding = "UTF-8"
)

df <- read.csv(
  "data\\metadata_abitur.csv",
  sep = ";", encoding = "UTF-8"
)

dirs <- c("img", "analyses")
for (dir in dirs) {
  if (dir.exists(dir) == FALSE) dir.create(dir)
}

state_cols <- c("BW", "BY", "BE", "BB", "HB", "HH", "HE", "MV",
                "NI", "NW", "RP", "SL", "SN", "ST", "SH", "TH")

state_abbr_map <- c(
  "Baden-Württemberg" = "BW", "Bayern" = "BY", "Berlin" = "BE",
  "Brandenburg" = "BB", "Bremen" = "HB", "Hamburg" = "HH",
  "Hessen" = "HE", "Mecklenburg-Vorpommern" = "MV",
  "Niedersachsen" = "NI", "Nordrhein-Westfalen" = "NW", "Rheinland-Pfalz" = "RP",
  "Saarland" = "SL", "Sachsen" = "SN", "Sachsen-Anhalt" = "ST",
  "Schleswig-Holstein" = "SH", "Thüringen" = "TH"
)

foreign_origin_authors <- c("Shakespeare", "Sophokles", "Euripides", "Anouilh",
                            "Moli?re", "Ibsen", "Camus", "Barbusse")

unit_key <- c("textID", "title", "author_lastname")

#------------------------------------------------------------------------------#
## 2 Reading-list data and coverage -------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## 2.1 Map of list coverage by state ------------------------------------------#
#------------------------------------------------------------------------------#

data <- read.csv(".\\data\\overview_abitur.csv", sep = ";",
                 stringsAsFactors = FALSE, check.names = FALSE)

data_summary <- data %>%
  rowwise() %>%
  mutate(total = sum(c_across(-state), na.rm = TRUE)) %>%
  ungroup()

east_states <- c("Brandenburg", "Mecklenburg-Vorpommern", "Sachsen-Anhalt",
                 "Sachsen", "Thüringen")

data_summary <- data_summary %>%
  mutate(region = if_else(state %in% east_states, "East", "West"))

germany_spat <- geodata::gadm("Germany", level = 1, path = tempdir())
germany_sf <- sf::st_as_sf(germany_spat)

for (col in names(germany_sf)) {
  if (is.character(germany_sf[[col]])) {
    germany_sf[[col]] <- iconv(germany_sf[[col]], from = "UTF-8", to = "latin1")
  }
}

germany_merged <- germany_sf %>%
  left_join(data_summary, by = c("NAME_1" = "state"))

east_union <- germany_merged %>%
  filter(region == "East") %>%
  st_union()

make_hatch_lines <- function(bbox, spacing = 0.3, angle = 45) {
  angle_rad <- angle * pi / 180
  diag_len <- sqrt((bbox["xmax"] - bbox["xmin"])^2 + (bbox["ymax"] - bbox["ymin"])^2)
  centre <- c(mean(bbox[c("xmin", "xmax")]), mean(bbox[c("ymin", "ymax")]))
  offsets <- seq(-diag_len, diag_len, by = spacing)
  
  lines <- map(offsets, function(off) {
    p1 <- c(-diag_len, off)
    p2 <- c(diag_len, off)
    rotate_pt <- function(p) {
      c(p[1] * cos(angle_rad) - p[2] * sin(angle_rad),
        p[1] * sin(angle_rad) + p[2] * cos(angle_rad)) + centre
    }
    st_linestring(rbind(rotate_pt(p1), rotate_pt(p2)))
  })
  
  st_sfc(lines, crs = st_crs(germany_merged))
}

hatch_lines <- make_hatch_lines(st_bbox(germany_merged), spacing = 0.3, angle = 45)
east_stripes <- st_intersection(st_sf(geometry = hatch_lines), east_union)

ggplot(germany_merged) +
  geom_sf(aes(fill = total), colour = "grey50", linewidth = 0.25) +
  geom_sf(data = east_stripes, colour = "grey80", linewidth = 0.65) +
  scale_fill_gradient(low = "white", high = "black", guide = "none") +
  theme_chr(base_size = 11) +
  theme(
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 4, r = 20, b = 4, l = 4),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

ggsave("img//map_states.png", dpi = 360, width = 20, height = 22, units = "cm")

#------------------------------------------------------------------------------#
## 2.2 Coverage heatmap: states x years ---------------------------------------#
#------------------------------------------------------------------------------#

df_long_states <- data %>%
  pivot_longer(-state, names_to = "year", values_to = "present") %>%
  mutate(
    year = as.integer(year),
    present = factor(if_else(is.na(present), "absent", "present"),
                     levels = c("absent", "present")),
    abbreviation = state_abbr_map[state]
  )

state_order <- df_long_states %>%
  group_by(abbreviation) %>%
  summarise(n_present = sum(present == "present"), .groups = "drop") %>%
  arrange(n_present) %>%
  pull(abbreviation)

df_long_states <- df_long_states %>% mutate(state = factor(abbreviation, levels = state_order))

ggplot(df_long_states, aes(x = year, y = state, fill = present)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("absent" = "white", "present" = "black"), name = NULL) +
  scale_x_continuous(breaks = unique(df_long_states$year), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  guides(fill = "none") +
  theme_chr(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid = element_blank())

ggsave("img//heatmap_data.png", dpi = 360, width = 16, height = 8, units = "cm")

#------------------------------------------------------------------------------#
## 2.3 Title-level summary of reading lists -----------------------------------#
#------------------------------------------------------------------------------#

abitur_summary <- df %>%
  mutate(
    states = rowSums(select(., all_of(state_cols)) == 1, na.rm = TRUE),
    title = str_squish(title)
  ) %>%
  group_by(textID, title, author_lastname, author_firstname, pub_year_clean, gender) %>%
  summarise(
    total_states = sum(states, na.rm = TRUE),
    first_abitur_year = min(year, na.rm = TRUE),
    n_years = n_distinct(year),
    .groups = "drop"
  )

#------------------------------------------------------------------------------#
## 2.4 Long format and most frequent authors and texts ------------------------#
#------------------------------------------------------------------------------#

df_long <- df %>%
  pivot_longer(cols = all_of(state_cols), names_to = "state", values_to = "included") %>%
  filter(!is.na(title) & title != "", included == 1) %>%
  select(year, state, title, textID, author_firstname, author_lastname)

all_unique_titles <- df_long %>% distinct(title, textID) %>% arrange(title)

# df_long carries no gender column; join it from df
author_gender_lookup <- df %>% distinct(author_lastname, author_firstname, gender)

author_counts_overall <- df_long %>%
  left_join(author_gender_lookup, by = c("author_lastname", "author_firstname")) %>%
  count(author_lastname, author_firstname, gender, sort = TRUE)

top_author_overall <- author_counts_overall %>% slice_max(n, n = 1, with_ties = TRUE)

top_author_by_gender <- df_long %>%
  left_join(author_gender_lookup, by = c("author_lastname", "author_firstname")) %>%
  count(gender, author_lastname, author_firstname, sort = TRUE) %>%
  group_by(gender) %>%
  slice_max(n, n = 1, with_ties = TRUE) %>%
  ungroup()

text_counts <- df_long %>%
  count(textID, title, sort = TRUE)

top_text_overall <- text_counts %>% slice_max(n, n = 1, with_ties = TRUE)

#------------------------------------------------------------------------------#
## 2.5 Coverage mask and balanced panel ---------------------------------------#
#------------------------------------------------------------------------------#

# Match state names on ASCII letters only, so umlaut encoding issues
# (Wuerttemberg, Thueringen) cannot break the lookup.
state_key <- function(x) tolower(gsub("[^A-Za-z]", "", x, useBytes = TRUE))
state_lookup <- setNames(unname(state_abbr_map), state_key(names(state_abbr_map)))

coverage <- data %>%
  pivot_longer(-state, names_to = "year", values_to = "v") %>%
  transmute(
    state = unname(state_lookup[state_key(state)]),
    year = as.integer(year),
    has_list = !is.na(v) & as.character(v) == "1"
  )

stopifnot(!anyNA(coverage$state), n_distinct(coverage$state) == 16)
message("State-years with a list: ", sum(coverage$has_list), " of ", nrow(coverage))

# Sanity check: titles should never appear in a state-year without a list
uncovered <- df_long %>%
  distinct(state, year) %>%
  left_join(coverage, by = c("state", "year")) %>%
  filter(is.na(has_list) | !has_list)
if (nrow(uncovered) > 0) {
  warning("Titles listed in state-years marked as uncovered in overview_abitur.csv: ",
          paste(uncovered$state, uncovered$year, collapse = ", "))
}

# Balanced panel: states with a list in every year of the window
balanced_window <- 2013:2026
balanced_states <- coverage %>%
  filter(year %in% balanced_window) %>%
  group_by(state) %>%
  summarise(full = all(has_list), .groups = "drop") %>%
  filter(full) %>%
  pull(state)
message("Balanced panel ", min(balanced_window), "-", max(balanced_window), ": ",
        paste(balanced_states, collapse = ", "))

df_long_bal  <- df_long %>% filter(state %in% balanced_states, year %in% balanced_window)
coverage_bal <- coverage %>% filter(state %in% balanced_states, year %in% balanced_window)

# Earliest year a title can be listed (publication year)
title_pub_year <- df %>%
  distinct(textID, pub_year_clean) %>%
  group_by(textID) %>%
  summarise(pub_year = suppressWarnings(min(as.numeric(pub_year_clean), na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(pub_year = if_else(is.finite(pub_year), pub_year, -Inf))  # unknown: always eligible

#------------------------------------------------------------------------------#
## 2.6 Panel builder ----------------------------------------------------------#
#------------------------------------------------------------------------------#
## Title x state x year grid with presence, list coverage and eligibility
## (a title can only be listed from its publication year onwards). Shared by
## all dynamic measures below.

build_panel <- function(df_long, coverage, unit_col = unit_key,
                        states = state_cols, years = NULL) {
  if (is.null(years)) years <- full_seq(range(df_long$year, na.rm = TRUE), 1)
  units <- df_long %>% distinct(across(all_of(unit_col)))
  presence <- df_long %>%
    filter(state %in% states, year %in% years) %>%
    distinct(across(all_of(c(unit_col, "state", "year")))) %>%
    mutate(present = 1L)

  units %>%
    expand_grid(state = states, year = years) %>%
    left_join(presence, by = c(unit_col, "state", "year")) %>%
    left_join(coverage, by = c("state", "year")) %>%
    left_join(title_pub_year, by = "textID") %>%
    mutate(
      present  = replace_na(present, 0L),
      has_list = replace_na(has_list, FALSE),
      eligible = year >= replace_na(pub_year, -Inf)
    ) %>%
    arrange(across(all_of(unit_col)), state, year) %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    mutate(
      present_next  = lead(present),
      has_list_next = lead(has_list),
      eligible_next = lead(eligible),
      has_list_prev = lag(has_list, default = FALSE)
    ) %>%
    ungroup() %>%
    select(-pub_year)
}

# Rescue effect: is the title on at least one OTHER state's list this year?
add_presence_elsewhere <- function(panel, unit_col = unit_key) {
  panel %>%
    group_by(across(all_of(c(unit_col, "year")))) %>%
    mutate(n_present = sum(present), other_present = (n_present - present) > 0) %>%
    ungroup()
}

# Only transitions where both years have a list and the title exists at t+1
valid_transitions <- function(panel) {
  panel %>% filter(!is.na(present_next), has_list, has_list_next, eligible_next)
}

#------------------------------------------------------------------------------#
## 3 Bibliographic data (German National Library) -----------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## 3.1 Clean DNB records ------------------------------------------------------#
#------------------------------------------------------------------------------#

dnb_detail <- dnb %>%
  mutate(
    pub_year = str_extract(pub_year, "\\d{4}"),
    pub_year = as.numeric(pub_year),
    language_codes = str_remove_all(language_codes, "1\\\\p"),
    language_codes = str_squish(str_remove_all(language_codes, "^\\|\\||\\|\\|$")),
    language_codes = str_replace_all(language_codes, "\\|\\|\\s*\\|\\|", "||"),
    language_codes = na_if(language_codes, ""),
    series_first = str_split(series, "\\|\\|", simplify = TRUE)[, 1],
    title = str_squish(title)
  ) %>%
  select(title, author_lastname, author_firstname,
         pub_year, language_codes, series_first)

write.csv(dnb_detail, "data\\dnb_detail.csv")

#------------------------------------------------------------------------------#
## 3.2 Merge DNB and reading-list data; classify entries ----------------------#
#------------------------------------------------------------------------------#

## Merge key: title + author (DNB records carry no textID). Each DNB record
## is classified as edition, translation, or unresolved.

dnb_with_abitur <- merge(
  dnb_detail, abitur_summary,
  by.x = c("title", "author_lastname", "author_firstname"),
  by.y = c("title", "author_lastname", "author_firstname"),
  all = TRUE
)

dnb_with_abitur <- dnb_with_abitur %>%
  mutate(literary_period = ifelse(pub_year_clean < 1950, "Classic", "Contemporary"))

dnb_with_abitur <- dnb_with_abitur %>%
  mutate(
    language_codes = replace_na(language_codes, ""),
    lang_list = str_split(language_codes, "\\|\\|") %>%
      map(str_trim) %>%
      map(~ .x[.x != ""]),
    n_languages = map_int(lang_list, length),
    is_foreign_origin = author_lastname %in% foreign_origin_authors,
    has_no_language_data = map_lgl(lang_list, ~ length(.x) == 0),
    entry_type = case_when(
      has_no_language_data ~ NA_character_,
      is_foreign_origin ~ "unclear - foreign-origin text, check original_language manually",
      map_lgl(lang_list, ~ all(.x == "ger")) ~ "edition",
      TRUE ~ "translation"
    )
  ) %>%
  select(-lang_list, -is_foreign_origin, -has_no_language_data)

write.csv(dnb_with_abitur, "data\\dnb_with_abitur.csv")

#------------------------------------------------------------------------------#
## 3.3 Title-level editions, translations and coverage ------------------------#
#------------------------------------------------------------------------------#

## Single source for n_editions, n_translations and dnb_coverage, used by the
## histograms, correlations and survival models.

title_summary <- dnb_with_abitur %>%
  group_by(textID, title, author_lastname) %>%
  summarise(
    abi_coverage = first(total_states) * first(n_years),
    n_editions = sum(entry_type == "edition", na.rm = TRUE),
    n_translations = sum(entry_type == "translation", na.rm = TRUE),
    n_unclear = sum(entry_type == "unclear - foreign-origin text, check original_language manually", na.rm = TRUE),
    n_no_language_data = sum(is.na(entry_type)),
    dnb_coverage = n_editions + n_translations,
    literary_period = first(literary_period),
    .groups = "drop"
  ) %>%
  mutate(coverage_status = case_when(
    dnb_coverage == 0 & (n_unclear > 0 | n_no_language_data > 0) ~ "unresolved - needs manual check",
    dnb_coverage == 0 ~ "confirmed zero",
    TRUE ~ "resolved"
  ))

write.csv(title_summary, "data\\title_summary.csv")

# Exclude unresolved (dependent publications + non-German-origin texts
# pending manual language check) -- see methods section for exact counts
title_summary_resolved <- title_summary %>%
  filter(coverage_status != "unresolved - needs manual check")

#------------------------------------------------------------------------------#
## 3.4 Editions and translations: distributions and correlation ---------------#
#------------------------------------------------------------------------------#

ggplot(title_summary_resolved, aes(x = n_editions)) +
  geom_histogram(fill = "black", colour = "white", alpha = 0.8, binwidth = 10) +
  labs(x = "Count", y = "Number of texts") +
  theme_chr(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("img//hist_DNB_editions.png", dpi = 360, width = 16, height = 12, units = "cm")

ggplot(title_summary_resolved, aes(x = n_translations)) +
  geom_histogram(fill = "black", colour = "white", alpha = 0.8, binwidth = 10) +
  labs(x = "Count", y = "Number of texts") +
  theme_chr(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("img//hist_DNB_translations.png", dpi = 360, width = 16, height = 12, units = "cm")

cor_test_ed_tr <- cor.test(title_summary_resolved$n_editions,
                           title_summary_resolved$n_translations,
                           method = "spearman")

cor_label_ed_tr <- paste0("Spearman's ?? = ", round(cor_test_ed_tr$estimate, 2),
                          ifelse(cor_test_ed_tr$p.value < 0.001, ", p < 0.001",
                                 paste0(", p = ", round(cor_test_ed_tr$p.value, 3))))

ggplot(title_summary_resolved, aes(x = n_editions, y = n_translations)) +
  geom_point() +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  labs(x = "Number of editions", y = "Number of translations") +
  theme_chr(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("img//scatter_DNB.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 3.5 Reading-list coverage x DNB coverage -----------------------------------#
#------------------------------------------------------------------------------#

cor_test <- cor.test(title_summary_resolved$abi_coverage, title_summary_resolved$dnb_coverage,
                     method = "spearman")

only_classics <- subset(title_summary_resolved, literary_period == "Classic")
cor_test_only_classics <- cor.test(only_classics$abi_coverage, only_classics$dnb_coverage,
                                   method = "spearman")

only_contemp <- subset(title_summary_resolved, literary_period == "Contemporary")
cor_test_only_contemp <- cor.test(only_contemp$abi_coverage, only_contemp$dnb_coverage,
                                  method = "spearman")

ggplot(title_summary_resolved, aes(x = abi_coverage, y = dnb_coverage)) +
  geom_jitter(aes(fill = literary_period), shape = 21, colour = "black",
              width = 0.05, height = 0, size = 2.5, alpha = 0.8, stroke = 0.3) +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  scale_fill_manual(values = c("Classic" = "white", "Contemporary" = "black"), name = "Period") +
  scale_x_log10() +
  labs(x = "Abitur coverage (states ? years, log scale)",
       y = "DNB coverage (editions + translations)") +
  theme_chr(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("img//scatter_DNBxAbi.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 4 Descriptive overview -----------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## 4.1 Text- and author-level summaries ---------------------------------------#
#------------------------------------------------------------------------------#

df <- df %>%
  mutate(states = rowSums(select(., all_of(state_cols)) == 1, na.rm = TRUE))

text_summary <- df %>%
  group_by(textID, title, author_lastname, author_firstname, pub_year_clean, gender) %>%
  summarise(
    total_states = sum(states, na.rm = TRUE),
    n_years = n_distinct(year),
    first_year_listed = min(year, na.rm = TRUE),
    last_year_listed  = max(year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_states))

write.csv(text_summary, "analyses\\text_summary.csv")

author_summary <- df %>%
  group_by(author_lastname, author_firstname) %>%
  summarise(
    total_states = sum(states, na.rm = TRUE),
    n_titles = n_distinct(title),
    n_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  arrange(desc(total_states))

write.csv(author_summary, "analyses\\author_summary.csv")

#------------------------------------------------------------------------------#
## 4.2 Author heatmap ---------------------------------------------------------#
#------------------------------------------------------------------------------#

author_yearly <- df %>%
  group_by(year, author_lastname, author_firstname) %>%
  summarise(states = sum(states, na.rm = TRUE), .groups = "drop")

authors_multi <- author_yearly %>%
  group_by(author_lastname, author_firstname) %>%
  summarise(n_years = n_distinct(year), .groups = "drop") %>%
  filter(n_years > 2)

author_heatmap <- author_yearly %>%
  semi_join(authors_multi, by = c("author_lastname", "author_firstname")) %>%
  group_by(year, author_lastname, author_firstname) %>%
  summarise(total_states = sum(states, na.rm = TRUE), .groups = "drop")

peak_years_author <- author_heatmap %>%
  group_by(author_lastname, author_firstname) %>%
  slice_max(order_by = total_states, n = 1, with_ties = FALSE) %>%
  select(author_lastname, author_firstname, peak_year = year)

author_heatmap <- author_heatmap %>%
  left_join(peak_years_author, by = c("author_lastname", "author_firstname")) %>%
  mutate(author_label = paste(author_firstname, author_lastname)) %>%
  mutate(author_label = gsub("^nan ", "", author_label))

author_order <- author_heatmap %>%
  distinct(author_label, peak_year) %>%
  arrange(peak_year) %>%
  pull(author_label)

author_heatmap$author_label <- factor(author_heatmap$author_label, levels = author_order)

ggplot(author_heatmap, aes(x = year, y = author_label, fill = total_states)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "lightgrey", high = "black", guide = "none") +
  theme_chr(base_size = 11) +
  theme(axis.text.y = element_text(size = 7),
        plot.margin = margin(t = 4, r = 20, b = 4, l = 4),
        axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave("img/author_heatmap.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 4.3 Text heatmap -----------------------------------------------------------#
#------------------------------------------------------------------------------#

text_yearly <- df %>%
  group_by(year, textID, title, author_lastname, author_firstname) %>%
  summarise(states = sum(states, na.rm = TRUE), .groups = "drop") %>%
  mutate(title_label = paste0(title, " (", author_lastname, ")"))

texts_multi <- text_yearly %>%
  group_by(textID, title_label) %>%
  summarise(n_years = n_distinct(year), .groups = "drop") %>%
  filter(n_years > 4)

text_multi <- text_yearly %>%
  semi_join(texts_multi, by = c("textID", "title_label"))

title_order <- text_multi %>%
  group_by(title_label, textID) %>%
  summarise(total_states = sum(states, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_states)) %>%
  pull(title_label)

text_multi <- text_multi %>%
  mutate(title_label = factor(title_label, levels = title_order))

peak_years_text <- text_multi %>%
  group_by(title) %>%
  slice_max(states, n = 1, with_ties = FALSE) %>%
  select(title, peak_year = year)

text_multi <- text_multi %>%
  left_join(peak_years_text, by = "title") %>%
  mutate(title = fct_reorder(title, peak_year))

ggplot(text_multi, aes(x = year, y = title, fill = states)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "lightgrey", high = "black", guide = "none") +
  theme_chr(base_size = 11) +
  theme(axis.text.y = element_text(size = 7),
        plot.margin = margin(t = 4, r = 20, b = 4, l = 4),
        axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave("img//text_heatmap.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 4.4 Share of female authors ------------------------------------------------#
#------------------------------------------------------------------------------#

gender_share <- df %>%
  group_by(year, gender) %>%
  summarise(n_authors = n_distinct(paste(author_firstname, author_lastname)), .groups = "drop") %>%
  group_by(year) %>%
  mutate(share = n_authors / sum(n_authors)) %>%
  filter(gender == "f") %>%
  ungroup() %>%
  complete(year = full_seq(min(df$year, na.rm = TRUE):max(df$year, na.rm = TRUE), 1),
           fill = list(share = 0))
# Years with no female authors are filled with 0

ggplot(gender_share, aes(x = year, y = share)) +
  geom_line(colour = "black", linewidth = 1) +
  geom_point(colour = "black", size = 2) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = "Year", y = "Proportion of female authors") +
  theme_chr(base_size = 11)

ggsave("img//gender_proportion.png", dpi = 360, width = 16, height = 8, units = "cm")

#------------------------------------------------------------------------------#
## 5 Canon dynamics -----------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## 5.1 Shannon diversity and its partition ------------------------------------#
#------------------------------------------------------------------------------#

## gamma = Shannon of pooled abundances (number of states listing each title)
## alpha = list-size-weighted mean within-list Shannon (= ln list size)
## beta  = gamma - alpha >= 0: differentiation between states;
##         exp(beta) = effective number of distinct lists

shannon_partition <- function(df_long) {
  occ <- df_long %>% distinct(textID, state, year)

  gamma <- occ %>%
    count(year, textID, name = "n") %>%
    group_by(year) %>%
    summarise(shannon_index = -sum((n / sum(n)) * log(n / sum(n))),
              richness = n(), .groups = "drop")

  alpha <- occ %>%
    count(year, state, name = "list_size") %>%
    group_by(year) %>%
    summarise(n_lists = n(),
              mean_list_size = mean(list_size),
              alpha = sum(list_size / sum(list_size) * log(list_size)),
              .groups = "drop")

  gamma %>%
    left_join(alpha, by = "year") %>%
    mutate(beta = shannon_index - alpha,
           eff_n_distinct_lists = exp(beta)) %>%
    arrange(year)
}

shannon_by_year     <- shannon_partition(df_long)
shannon_by_year_bal <- shannon_partition(df_long_bal)
write.csv(shannon_by_year, "analyses\\shannon_diversity_by_year.csv")
write.csv(shannon_by_year_bal, "analyses\\shannon_diversity_by_year_balanced.csv")

pre_2026  <- shannon_by_year %>% filter(year != 2026)
band_min  <- min(pre_2026$shannon_index, na.rm = TRUE)
band_max  <- max(pre_2026$shannon_index, na.rm = TRUE)
band_mean <- mean(pre_2026$shannon_index, na.rm = TRUE)  # excludes 2026

ggplot(shannon_by_year, aes(x = year, y = shannon_index)) +
  geom_line(colour = "black") +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(colour = "black", size = 2.5) +
  labs(x = "Year", y = "Shannon index") +
  theme_chr(base_size = 11)

ggsave("img//shannon_diversity_by_year.png", dpi = 360, width = 16, height = 12, units = "cm")

shannon_long <- shannon_by_year %>%
  select(year, `Pooled (gamma)` = shannon_index,
         `Within lists (alpha)` = alpha, `Between states (beta)` = beta) %>%
  pivot_longer(-year, names_to = "component", values_to = "value")

ggplot(shannon_long, aes(x = year, y = value, linetype = component, shape = component)) +
  geom_line(colour = "black", linewidth = 0.6) +
  geom_point(colour = "black", size = 2) +
  scale_linetype_manual(
    values = c("Pooled (gamma)"        = "solid",
               "Within lists (alpha)"  = "dashed",
               "Between states (beta)" = "dotted"),
    name = NULL) +
  scale_shape_manual(
    values = c("Pooled (gamma)"        = 16,  
               "Within lists (alpha)"  = 1, 
               "Between states (beta)" = 17), 
    name = NULL) +
  labs(x = "Year", y = "Shannon index") +
  theme_chr(base_size = 11) +
  theme(legend.position = "top",
        legend.key.width = unit(1.5, "cm"))

ggsave("img//shannon_partition_by_year.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 5.2 Yearly additions, removals and within-state turnover -------------------#
#------------------------------------------------------------------------------#

## Titles added to and removed from lists, normalised by the number of states
## with a list in both years. Within-state turnover is the pooled Jaccard
## dissimilarity between consecutive lists of the same state and does not
## depend on which states are observed.

yearly_dynamics <- function(df_long, coverage, unit_col = unit_key,
                            states = state_cols, years = NULL) {
  build_panel(df_long, coverage, unit_col, states, years) %>%
    valid_transitions() %>%
    group_by(year_from = year) %>%
    summarise(
      n_lists_both_years = n_distinct(state),
      n_added   = sum(present == 0 & present_next == 1),
      n_removed = sum(present == 1 & present_next == 0),
      n_kept    = sum(present == 1 & present_next == 1),
      list_size_from = sum(present) / n_lists_both_years,
      list_size_to   = sum(present_next) / n_lists_both_years,
      colonisation_rate = n_added / sum(present == 0),   # reported, not compared
      extinction_rate   = n_removed / sum(present == 1),
      .groups = "drop"
    ) %>%
    mutate(
      year_to = year_from + 1L,
      added_per_list   = n_added / n_lists_both_years,
      removed_per_list = n_removed / n_lists_both_years,
      net_change_per_list = added_per_list - removed_per_list,
      within_state_turnover = (n_added + n_removed) / (n_added + n_removed + n_kept)
    )
}

period_rates_titles <- yearly_dynamics(df_long, coverage)
write.csv(period_rates_titles, "analyses\\period_rates_titles.csv")

period_rates_long <- period_rates_titles %>%
  select(year_to, Added = added_per_list, Removed = removed_per_list) %>%
  pivot_longer(-year_to, names_to = "change", values_to = "titles_per_list")

ggplot(period_rates_long, aes(x = year_to, y = titles_per_list, linetype = change)) +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85",
              linewidth = 0.6, alpha = 0.5) +
  geom_line(colour = "black") +
  geom_point(aes(shape = change), colour = "black", size = 2.5) +
  scale_shape_manual(values = c("Added" = 16, "Removed" = 1), name = NULL) +
  scale_linetype_manual(values = c("Added" = "solid", "Removed" = "dashed"), name = NULL) +
  labs(x = "Year", y = "Titles added / removed per list") +
  theme_chr(base_size = 11) +
  theme(legend.position = "top")

ggsave("img//period_rates.png", dpi = 360, width = 16, height = 12, units = "cm")

ggplot(period_rates_titles, aes(x = year_to, y = within_state_turnover)) +
  geom_line(colour = "black") +
  geom_smooth(data = filter(period_rates_titles, year_to < 2026),
              method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(colour = "black", size = 2) +
  labs(x = "Year", y = "Within-state turnover (1 - Jaccard)") +
  theme_chr(base_size = 11)

ggsave("img//turnover_within_state.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 5.3 Canon-wide turnover (Jaccard/Soerensen) --------------------------------#
#------------------------------------------------------------------------------#

## Compares the union of titles across all observed states in consecutive
## years. Because the set of observed states changes, the balanced-panel
## series is the comparable one.

turnover_time_series <- function(df_long) {
  years <- sort(unique(df_long$year))

  yearly_titles <- df_long %>%
    group_by(year) %>%
    summarise(titles = list(unique(textID)), .groups = "drop")

  map_dfr(seq_len(length(years) - 1), function(i) {
    y1 <- years[i]; y2 <- years[i + 1]
    t1 <- yearly_titles$titles[[which(yearly_titles$year == y1)]]
    t2 <- yearly_titles$titles[[which(yearly_titles$year == y2)]]

    intersection_n <- length(intersect(t1, t2))
    union_n <- length(union(t1, t2))

    tibble(
      year_from = y1, year_to = y2,
      n_titles_from = length(t1), n_titles_to = length(t2),
      n_added = length(setdiff(t2, t1)), n_removed = length(setdiff(t1, t2)),
      jaccard_similarity = intersection_n / union_n,
      sorensen_similarity = 2 * intersection_n / (length(t1) + length(t2)),
      turnover = 1 - jaccard_similarity
    )
  })
}

turnover_titles     <- turnover_time_series(df_long)
turnover_titles_bal <- turnover_time_series(df_long_bal)
write.csv(turnover_titles, "analyses\\turnover_time_series.csv")
write.csv(turnover_titles_bal, "analyses\\turnover_time_series_balanced.csv")

ggplot(turnover_titles, aes(x = year_to, y = turnover)) +
  geom_line(colour = "black") +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(colour = "black", size = 2) +
  xlim(2005, 2026) +
  labs(x = "Year", y = "Turnover (1 - Jaccard similarity to previous year)") +
  theme_chr(base_size = 11)

ggsave("img//turnover_time_series.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 5.4 Colonisation and extinction rates per title ----------------------------#
#------------------------------------------------------------------------------#

## colonisation = P(present at t+1 | absent at t, present elsewhere at t)
## extinction   = P(absent at t+1 | present at t)
## Whole-period rates are descriptive; the survival models use a lagged rate
## (section 7.3).

metapopulation_rates <- function(df_long, coverage, unit_col = unit_key,
                                 rescue_effect = TRUE, states = state_cols,
                                 years = NULL) {
  panel <- build_panel(df_long, coverage, unit_col, states, years) %>%
    add_presence_elsewhere(unit_col)

  footprint <- panel %>%
    group_by(across(all_of(unit_col))) %>%
    summarise(n_states_ever = n_distinct(state[present == 1]),
              years_present_total = sum(present), .groups = "drop")

  valid_transitions(panel) %>%
    group_by(across(all_of(unit_col))) %>%
    summarise(
      colon_opp  = if (rescue_effect) sum(present == 0 & other_present) else sum(present == 0),
      colon_succ = if (rescue_effect) sum(present == 0 & other_present & present_next == 1)
                   else sum(present == 0 & present_next == 1),
      ext_opp    = sum(present == 1),
      ext_succ   = sum(present == 1 & present_next == 0),
      .groups = "drop"
    ) %>%
    left_join(footprint, by = unit_col) %>%
    mutate(colonisation_rate = colon_succ / colon_opp,
           extinction_rate   = ext_succ / ext_opp) %>%
    arrange(desc(colonisation_rate))
}

metapop_titles <- metapopulation_rates(df_long, coverage)

expectation_model <- lm(colonisation_rate ~ years_present_total,
                        data = metapop_titles, na.action = na.exclude)

metapop_titles <- metapop_titles %>%
  mutate(
    expected_colonisation_rate = predict(expectation_model, newdata = .),
    colonisation_residual = colonisation_rate - expected_colonisation_rate
  ) %>%
  arrange(desc(colonisation_residual))

write.csv(metapop_titles, "analyses\\metapop_rates_titles.csv")

#------------------------------------------------------------------------------#
## 5.5 Contagion-dependence by year -------------------------------------------#
#------------------------------------------------------------------------------#

pooled_contagion_by_year <- function(df_long, coverage, unit_col = unit_key,
                                     gap_threshold = 0.15, low_rate_threshold = 0.15) {
  opportunities <- build_panel(df_long, coverage, unit_col) %>%
    add_presence_elsewhere(unit_col) %>%
    valid_transitions() %>%
    filter(present == 0)

  rates <- opportunities %>%
    group_by(year, other_present) %>%
    summarise(opportunities = n(), rate = mean(present_next == 1), .groups = "drop") %>%
    mutate(condition = if_else(other_present, "elsewhere_present", "elsewhere_absent")) %>%
    select(-other_present) %>%
    pivot_wider(names_from = condition, values_from = c(opportunities, rate),
                names_glue = "{condition}_{.value}")

  required_cols <- c("elsewhere_present_opportunities", "elsewhere_absent_opportunities",
                     "elsewhere_present_rate", "elsewhere_absent_rate")
  for (col in setdiff(required_cols, names(rates))) rates[[col]] <- NA_real_

  rates %>%
    mutate(
      contagion_gap = elsewhere_present_rate - elsewhere_absent_rate,
      classification = case_when(
        is.na(elsewhere_present_rate) | is.na(elsewhere_absent_rate) ~ "insufficient data",
        contagion_gap >= gap_threshold & elsewhere_absent_rate < low_rate_threshold ~ "contagion-dependent",
        abs(contagion_gap) < gap_threshold ~ "independently adopted",
        TRUE ~ "mixed / unclear"
      )
    ) %>%
    arrange(year)
}

period_contagion_titles <- pooled_contagion_by_year(df_long, coverage)
write.csv(period_contagion_titles, "analyses\\period_contagion_titles.csv")

ggplot(period_contagion_titles, aes(x = year, y = contagion_gap, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(colour = "black") +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(colour = "black", size = 2.5) +
  labs(x = "Year", y = "Contagion gap (elsewhere-present rate - elsewhere-absent rate)") +
  theme_chr(base_size = 11)

ggsave("img//period_contagion.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 5.6 Balanced panel: additional outputs -------------------------------------#
#------------------------------------------------------------------------------#

period_rates_bal <- yearly_dynamics(df_long_bal, coverage_bal, states = balanced_states)
metapop_bal      <- metapopulation_rates(df_long_bal, coverage_bal, states = balanced_states)

write.csv(period_rates_bal, "analyses\\period_rates_titles_balanced.csv")
write.csv(metapop_bal, "analyses\\metapop_rates_titles_balanced.csv")

#------------------------------------------------------------------------------#
## 6 Change-point tests -------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## 6.1 2026 endpoint test and 2020 segmented test -----------------------------#
#------------------------------------------------------------------------------#

## (a) Fit the 2006-2025 linear trend and test whether 2026 lies outside the
##     95% prediction interval (equivalent to a t-test on a 2026 dummy).
## (b) Broken-stick model with a change in level and slope from 2020 (shared
##     national pool), compared with a straight line by F-test, on years
##     before 2026. Exploratory, given the small number of years.

test_endpoint_break <- function(dat, x, y, label, break_year = 2026) {
  pre  <- dat %>% filter(.data[[x]] < break_year, !is.na(.data[[y]]))
  post <- dat %>% filter(.data[[x]] == break_year, !is.na(.data[[y]]))
  if (nrow(post) == 0 || nrow(pre) < 4) return(NULL)

  fit <- lm(reformulate(x, y), data = pre)
  pr  <- predict(fit, newdata = post[1, ], interval = "prediction", level = 0.95, se.fit = TRUE)
  obs <- post[[y]][1]
  pred <- pr$fit[1, "fit"]
  se_pred <- sqrt(pr$se.fit[1]^2 + pr$residual.scale^2)
  t_stat <- (obs - pred) / se_pred

  tibble(
    measure = label, n_pre = nrow(pre),
    observed_2026 = obs, predicted_2026 = pred,
    pi95_lower = pr$fit[1, "lwr"], pi95_upper = pr$fit[1, "upr"],
    t = t_stat, df = pr$df, p_two_sided = 2 * pt(-abs(t_stat), pr$df),
    pre_slope = unname(coef(fit)[x]),
    pre_slope_p = summary(fit)$coefficients[x, 4]
  )
}

test_segmented <- function(dat, x, y, label, cut = 2020, end_before = 2026) {
  d <- dat %>%
    filter(.data[[x]] < end_before, !is.na(.data[[y]])) %>%
    mutate(xx = .data[[x]], yy = .data[[y]],
           post = as.numeric(xx >= cut), post_slope = pmax(xx - cut, 0))
  if (nrow(d) < 6 || length(unique(d$post)) < 2) return(NULL)

  m0 <- lm(yy ~ xx, data = d)
  m1 <- lm(yy ~ xx + post + post_slope, data = d)
  a  <- anova(m0, m1)
  co <- summary(m1)$coefficients

  tibble(
    measure = label, cut = cut, n = nrow(d),
    level_shift = co["post", 1], level_shift_p = co["post", 4],
    slope_change = co["post_slope", 1], slope_change_p = co["post_slope", 4],
    F = a$F[2], F_p = a$`Pr(>F)`[2]
  )
}

break_inputs <- list(
  list(shannon_by_year,     "year",    "shannon_index",         "Shannon gamma (all states)"),
  list(shannon_by_year,     "year",    "alpha",                 "Shannon alpha (all states)"),
  list(shannon_by_year,     "year",    "beta",                  "Shannon beta (all states)"),
  list(shannon_by_year_bal, "year",    "shannon_index",         "Shannon gamma (balanced)"),
  list(shannon_by_year_bal, "year",    "alpha",                 "Shannon alpha (balanced)"),
  list(shannon_by_year_bal, "year",    "beta",                  "Shannon beta (balanced)"),
  list(turnover_titles,     "year_to", "turnover",              "Canon-wide turnover (all states)"),
  list(turnover_titles_bal, "year_to", "turnover",              "Canon-wide turnover (balanced)"),
  list(period_rates_titles, "year_to", "within_state_turnover", "Within-state turnover"),
  list(period_rates_titles, "year_to", "added_per_list",        "Titles added per list"),
  list(period_rates_titles, "year_to", "removed_per_list",      "Titles removed per list")
)

break_2026_tests <- map_dfr(break_inputs, ~ test_endpoint_break(.x[[1]], .x[[2]], .x[[3]], .x[[4]]))
break_2020_tests <- map_dfr(break_inputs, ~ test_segmented(.x[[1]], .x[[2]], .x[[3]], .x[[4]]))

write.csv(break_2026_tests, "analyses\\break_2026_tests.csv")
write.csv(break_2020_tests, "analyses\\break_2020_segmented_tests.csv")
break_2026_tests
break_2020_tests

#------------------------------------------------------------------------------#
## 7 Persistence --------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
## 7.1 Survival spells --------------------------------------------------------#
#------------------------------------------------------------------------------#

## A spell is a run of consecutive years a title is on one state's list.
##   event = 1      the following year has a list and the title is absent
##   event = 0      the following year has no list, or lies beyond 2026
##   left_censored  the year before the spell has no list (start unknown)
## Spells of the same title are clustered (text_id) in the Cox models.

build_spells <- function(df_long, coverage, unit_col = unit_key,
                         states = state_cols, years = NULL) {
  build_panel(df_long, coverage, unit_col, states, years) %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    mutate(spell_no = cumsum(present == 1L & lag(present, default = 0L) == 0L)) %>%
    filter(present == 1L) %>%
    group_by(across(all_of(c(unit_col, "state", "spell_no")))) %>%
    summarise(
      start_year = min(year),
      end_year = max(year),
      left_censored = !first(has_list_prev),
      next_has_list = last(has_list_next),
      .groups = "drop"
    ) %>%
    mutate(
      duration = end_year - start_year + 1L,
      event = if_else(!is.na(next_has_list) & next_has_list, 1L, 0L),
      text_id = paste(title, author_lastname, sep = " - ")
    ) %>%
    select(-next_has_list, -spell_no)
}

spells_titles <- build_spells(df_long, coverage)
write.csv(spells_titles, "analyses\\spells_titles.csv")

spells_titles %>%
  summarise(n_spells = n(), n_removals = sum(event),
            n_censored = sum(event == 0), n_left_censored = sum(left_censored))

## Persistence: core vs periphery (descriptive) -------------------------------#

title_persistence <- spells_titles %>%
  group_by(across(all_of(unit_key))) %>%
  summarise(
    total_presence = sum(duration),
    n_spells = n(),
    median_duration = median(duration),
    .groups = "drop"
  ) %>%
  left_join(metapop_titles %>% select(all_of(unit_key), colonisation_rate, n_states_ever),
            by = unit_key) %>%
  filter(!is.na(colonisation_rate))

presence_cutoff <- median(title_persistence$total_presence, na.rm = TRUE)

title_persistence <- title_persistence %>%
  mutate(tier = if_else(total_presence > presence_cutoff, "core", "periphery"))

title_persistence %>%
  group_by(tier) %>%
  summarise(n_titles = n(),
            mean_colonisation_rate = mean(colonisation_rate, na.rm = TRUE),
            median_colonisation_rate = median(colonisation_rate, na.rm = TRUE))

core_examples <- title_persistence %>%
  filter(tier == "core", n_states_ever >= 3) %>%
  arrange(desc(total_presence)) %>%
  select(title, author_lastname, total_presence, n_spells, median_duration,
         colonisation_rate, n_states_ever) %>%
  head(10)

periphery_examples <- title_persistence %>%
  filter(tier == "periphery") %>%
  arrange(total_presence) %>%
  select(title, author_lastname, total_presence, n_spells, median_duration, colonisation_rate) %>%
  head(10)

core_examples
periphery_examples

## Kaplan-Meier ---------------------------------------------------------------#

km_fit     <- survfit(Surv(duration, event) ~ 1, data = spells_titles)
km_fit_obs <- survfit(Surv(duration, event) ~ 1, data = filter(spells_titles, !left_censored))
km_fit
km_fit_obs   # sensitivity: spells with an observed start only

km_curve <- bind_rows(
  tibble(time = 0, surv = 1, n_censor = 0),
  tibble(time = km_fit$time, surv = km_fit$surv, n_censor = km_fit$n.censor)
)

median_survival <- km_curve %>% filter(surv <= 0.5) %>% slice(1) %>% pull(time)
if (length(median_survival) == 0) {
  message("Median spell survival not reached within the observation window.")
}

ggplot(km_curve, aes(x = time, y = surv)) +
  geom_step(colour = "black") +
  geom_point(data = km_curve %>% filter(n_censor > 0), shape = 3, size = 2, colour = "black") +
  labs(x = "Years since colonisation", y = "Proportion of spells still present") +
  theme_chr(base_size = 11)

ggsave("img//km_survival_curve.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## 7.2 Cumulative DNB counts by year ------------------------------------------#
#------------------------------------------------------------------------------#

dnb_counts_by_year <- dnb_with_abitur %>%
  semi_join(title_summary_resolved, by = unit_key) %>%
  filter(entry_type %in% c("edition", "translation"), !is.na(pub_year)) %>%
  count(across(all_of(unit_key)), entry_type, pub_year, name = "n_pub")

n_dropped_pub <- dnb_with_abitur %>%
  filter(entry_type %in% c("edition", "translation"), is.na(pub_year)) %>%
  nrow()
if (n_dropped_pub > 0) message("DNB entries without publication year (dropped): ", n_dropped_pub)

add_cumulative_dnb <- function(ap) {
  keys <- ap %>% distinct(across(all_of(unit_key)), calendar_year)
  cum <- keys %>%
    left_join(dnb_counts_by_year,
              by = join_by(textID, title, author_lastname, calendar_year >= pub_year)) %>%
    group_by(across(all_of(unit_key)), calendar_year) %>%
    summarise(
      cumulative_editions     = sum(n_pub[entry_type %in% "edition"]),
      cumulative_translations = sum(n_pub[entry_type %in% "translation"]),
      .groups = "drop"
    )
  ap %>%
    left_join(cum, by = c(unit_key, "calendar_year")) %>%
    mutate(cumulative_editions_10     = cumulative_editions / 10,
           cumulative_translations_10 = cumulative_translations / 10)
}

#------------------------------------------------------------------------------#
## 7.3 Lagged colonisation rate -----------------------------------------------#
#------------------------------------------------------------------------------#

## For a title in state s in year t: colonisation successes / opportunities in
## other states, using transitions up to t-1 only. Shrunk towards the pooled
## rate p0 with k pseudo-opportunities.

lagged_colonisation <- function(df_long, coverage, unit_col = unit_key,
                                states = state_cols, years = NULL, k = 10) {
  panel <- build_panel(df_long, coverage, unit_col, states, years) %>%
    add_presence_elsewhere(unit_col)

  opp <- valid_transitions(panel) %>%
    filter(present == 0, other_present) %>%
    mutate(succ = as.integer(present_next == 1L))

  p0 <- mean(opp$succ)

  by_title <- opp %>%
    group_by(across(all_of(c(unit_col, "year")))) %>%
    summarise(n_opp = n(), n_succ = sum(succ), .groups = "drop")
  by_title_state <- opp %>%
    group_by(across(all_of(c(unit_col, "state", "year")))) %>%
    summarise(n_opp_s = n(), n_succ_s = sum(succ), .groups = "drop")

  panel %>%
    select(all_of(c(unit_col, "state", "year"))) %>%
    left_join(by_title, by = c(unit_col, "year")) %>%
    left_join(by_title_state, by = c(unit_col, "state", "year")) %>%
    mutate(across(c(n_opp, n_succ, n_opp_s, n_succ_s), ~ replace_na(.x, 0L))) %>%
    arrange(across(all_of(c(unit_col, "state"))), year) %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    mutate(cum_opp  = lag(cumsum(n_opp - n_opp_s),  default = 0L),
           cum_succ = lag(cumsum(n_succ - n_succ_s), default = 0L)) %>%
    ungroup() %>%
    transmute(across(all_of(unit_col)), state, calendar_year = year,
              colon_opp_lag = cum_opp,
              colonisation_rate_lag = (cum_succ + k * p0) / (cum_opp + k))
}

#------------------------------------------------------------------------------#
## 7.4 Cox proportional hazards models ----------------------------------------#
#------------------------------------------------------------------------------#

## Main model: time-varying lagged colonisation rate and cumulative DNB counts.
## Sensitivity: 2026-list control, window ending 2025, observed starts only,
## balanced states, literary-period interaction, split period, and the
## whole-period (non-lagged) rate for comparison.

titles_on_2026_list <- df_long %>% filter(year == 2026) %>% distinct(textID) %>% pull(textID)

title_period <- dnb_with_abitur %>%
  filter(!is.na(textID)) %>%
  distinct(across(all_of(unit_key)), literary_period)

build_age_period_data <- function(spells, lagged_colon) {
  spells %>%
    mutate(spell_id = row_number()) %>%
    rowwise() %>%
    reframe(
      spell_id, textID, title, author_lastname, state, left_censored, text_id,
      calendar_year = seq(start_year, end_year),
      tstart = seq(0, duration - 1),
      tstop  = seq(1, duration),
      event  = c(rep(0L, duration - 1), event)
    ) %>%
    left_join(lagged_colon, by = c(unit_key, "state", "calendar_year")) %>%
    mutate(colonisation_rate_lag_x10 = colonisation_rate_lag * 10,
           on_2026_list = as.integer(textID %in% titles_on_2026_list)) %>%
    left_join(title_period, by = unit_key) %>%
    add_cumulative_dnb()
}

tidy_cox <- function(fit, label) {
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  se_col <- if ("robust se" %in% colnames(s)) "robust se" else "se(coef)"
  tibble(
    model = label, term = rownames(s),
    hr = exp(s[, "coef"]),
    ci95_low  = exp(s[, "coef"] - 1.96 * s[, se_col]),
    ci95_high = exp(s[, "coef"] + 1.96 * s[, se_col]),
    p = s[, ncol(s)],
    n_rows = fit$n, n_events = fit$nevent
  )
}

fit_cox <- function(formula, data, label) {
  tryCatch(coxph(formula, data = data),
           error = function(e) { message("Cox model '", label, "' failed: ", conditionMessage(e)); NULL })
}

f_main <- Surv(tstart, tstop, event) ~ colonisation_rate_lag_x10 +
  cumulative_editions_10 + cumulative_translations_10 + cluster(text_id)

# Main model ------------------------------------------------------------------#

lagged_colon    <- lagged_colonisation(df_long, coverage)
age_period_data <- build_age_period_data(spells_titles, lagged_colon)

cox_main <- fit_cox(f_main, age_period_data, "main")

# Sensitivity 1: control for titles on the shared 2026 list -------------------#

cox_2026_ctrl <- fit_cox(update(f_main, . ~ . + on_2026_list), age_period_data,
                         "control: on 2026 list")

# Sensitivity 2: observation window ends 2025 (no shared list) ----------------#

years_to_2025   <- 2006:2025
spells_2025     <- build_spells(df_long, coverage, years = years_to_2025)
lagged_2025     <- lagged_colonisation(df_long, coverage, years = years_to_2025)
age_period_2025 <- build_age_period_data(spells_2025, lagged_2025)
cox_to_2025     <- fit_cox(f_main, age_period_2025, "window 2006-2025")

# Sensitivity 3: spells with an observed start only ---------------------------#

cox_obs_start <- fit_cox(f_main, filter(age_period_data, !left_censored),
                         "excl. left-censored spells")

# Sensitivity 4: balanced panel ------------------------------------------------#

spells_bal     <- build_spells(df_long_bal, coverage_bal, states = balanced_states)
lagged_bal     <- lagged_colonisation(df_long_bal, coverage_bal, states = balanced_states)
age_period_bal <- build_age_period_data(spells_bal, lagged_bal)
cox_bal        <- fit_cox(f_main, age_period_bal, "balanced states 2013-2026")

# Sensitivity 5: literary period interaction (former interaction model) --------#

cox_period <- fit_cox(
  Surv(tstart, tstop, event) ~ literary_period * (cumulative_editions_10 + cumulative_translations_10) +
    colonisation_rate_lag_x10 + cluster(text_id),
  age_period_data, "period interaction"
)

# Sensitivity 6: split period -- colonisation 2006-2015 predicts spells
# starting 2016 or later (static covariate, strictly prior information) -------#

colon_early <- metapopulation_rates(df_long, coverage, years = 2006:2015) %>%
  select(all_of(unit_key), colonisation_rate_early = colonisation_rate)

cox_split_data <- spells_titles %>%
  filter(start_year >= 2016) %>%
  left_join(colon_early, by = unit_key) %>%
  mutate(colonisation_rate_early_x10 = colonisation_rate_early * 10)

cox_split <- fit_cox(Surv(duration, event) ~ colonisation_rate_early_x10 + cluster(text_id),
                     cox_split_data, "split: 2006-15 rate -> spells from 2016")

# Comparison: whole-period (non-lagged) colonisation rate ---------------------#

cox_whole_period_data <- spells_titles %>%
  left_join(metapop_titles %>% select(all_of(unit_key), colonisation_rate), by = unit_key) %>%
  mutate(colonisation_rate_x10 = colonisation_rate * 10)

cox_whole_period <- fit_cox(Surv(duration, event) ~ colonisation_rate_x10 + cluster(text_id),
                   cox_whole_period_data, "whole-period rate (not lagged)")

# Collect results -------------------------------------------------------------#

cox_results <- bind_rows(
  tidy_cox(cox_main,      "main"),
  tidy_cox(cox_2026_ctrl, "control: on 2026 list"),
  tidy_cox(cox_to_2025,   "window 2006-2025"),
  tidy_cox(cox_obs_start, "excl. left-censored spells"),
  tidy_cox(cox_bal,       "balanced states 2013-2026"),
  tidy_cox(cox_period,    "period interaction"),
  tidy_cox(cox_split,     "split: 2006-15 rate -> spells from 2016"),
  tidy_cox(cox_whole_period,       "whole-period rate (not lagged)")
)

write.csv(cox_results, "analyses\\cox_models_all.csv")
write.csv(filter(cox_results, model == "main"), "analyses\\cox_model_coefficients.csv")
cox_results

# Proportional-hazards check for the main model
if (!is.null(cox_main)) {
  zph_test <- tryCatch(cox.zph(cox_main), error = function(e) NULL)
  if (!is.null(zph_test)) {
    write.csv(as.data.frame(zph_test$table) %>% rownames_to_column("term"),
              "analyses\\cox_ph_test.csv")
  }
}

cox.zph(cox_main)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#