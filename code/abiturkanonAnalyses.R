####################### CHR2027: German School Canon ###########################

#------------------------------------------------------------------------------#
## Required packages ----------------------------------------------------------#
#------------------------------------------------------------------------------#

library(tidyverse)
library(broom)
library(ggridges)
library(ggthemes)
library(igraph)
library(geodata)
library(sf)
library(survival)
library(vegan)

#------------------------------------------------------------------------------#
## Visual language ------------------------------------------------------------#
#------------------------------------------------------------------------------#

theme_chr <- function(base_size = 13, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      line =              element_line(colour = '#DADADA', linewidth = 1,
                                       linetype = 1, lineend = "butt"),
      rect =              element_rect(fill = "transparent", colour = NA,
                                       linewidth = 0.5, linetype = 1),
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
## Settings / paths -----------------------------------------------------------#
#------------------------------------------------------------------------------#

setwd("...\\abiturkanon")

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
                            "Molière", "Ibsen", "Camus", "Barbusse")

#------------------------------------------------------------------------------#
## Map ------------------------------------------------------------------------#
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
## Heatmap: States x time -----------------------------------------------------#
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
## Clean DNB data -------------------------------------------------------------#
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
## Summarise school reading data ----------------------------------------------#
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
## Counts: most frequent author / text, overall and by gender -----------------#
#------------------------------------------------------------------------------#

df_long <- df %>%
  pivot_longer(cols = all_of(state_cols), names_to = "state", values_to = "included") %>%
  filter(!is.na(title) & title != "", included == 1) %>%
  select(year, state, title, textID, author_firstname, author_lastname)

all_unique_titles <- df_long %>% distinct(title, textID) %>% arrange(title)

author_counts_overall <- df_long %>%
  count(author_lastname, author_firstname, gender <- NA, sort = TRUE) 

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
## DNB x Abitur ---------------------------------------------------------------#
#------------------------------------------------------------------------------#

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
## Title-level summary: editions, translations, coverage ----------------------#
#------------------------------------------------------------------------------#

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

title_summary_resolved <- title_summary %>%
  filter(coverage_status != "unresolved - needs manual check")

#------------------------------------------------------------------------------#
## DNB histograms & editions x translations correlation -----------------------#
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
## Abi x DNB correlations (overall + by period) -------------------------------#
#------------------------------------------------------------------------------#

cor_test <- cor.test(title_summary_resolved$abi_coverage, title_summary_resolved$dnb_coverage,
                     method = "spearman")
print(cor_test)

only_classics <- subset(title_summary_resolved, literary_period == "Classic")
cor_test_only_classics <- cor.test(only_classics$abi_coverage, only_classics$dnb_coverage,
                                   method = "spearman")
print(cor_test_only_classics)

only_contemp <- subset(title_summary_resolved, literary_period == "Contemporary")
cor_test_only_contemp <- cor.test(only_contemp$abi_coverage, only_contemp$dnb_coverage,
                                  method = "spearman")
print(cor_test_only_contemp)

ggplot(title_summary_resolved, aes(x = abi_coverage, y = dnb_coverage)) +
  geom_jitter(aes(fill = literary_period), shape = 21, colour = "black",
              width = 0.05, height = 0, size = 2.5, alpha = 0.8, stroke = 0.3) +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  scale_fill_manual(values = c("Classic" = "white", "Contemporary" = "black"), name = "Period") +
  scale_x_log10() +
  labs(x = "Abitur coverage (states × years, log scale)",
       y = "DNB coverage (editions + translations)") +
  theme_chr(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave("img//scatter_DNBxAbi.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## Text-level and author-level summaries --------------------------------------#
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
## Author heatmap -------------------------------------------------------------#
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
## Text heatmap ---------------------------------------------------------------#
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
## Gender proportion ----------------------------------------------------------#
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

ggplot(gender_share, aes(x = year, y = share)) +
  geom_line(colour = "black", linewidth = 1) +
  geom_point(colour = "black", size = 2) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(x = "Year", y = "Proportion of female authors") +
  theme_chr(base_size = 11)

ggsave("img//gender_proportion.png", dpi = 360, width = 16, height = 8, units = "cm")

#------------------------------------------------------------------------------#
## Metapopulation colonisation / extinction rates -----------------------------#
#------------------------------------------------------------------------------#

metapopulation_rates <- function(df_long, unit_col = c("textID", "title", "author_lastname"),
                                 rescue_effect = FALSE) {
  
  years <- full_seq(range(df_long$year, na.rm = TRUE), 1)
  units <- df_long %>% distinct(across(all_of(unit_col)))
  
  panel <- units %>%
    expand_grid(state = state_cols, year = years) %>%
    left_join(
      df_long %>%
        distinct(across(all_of(c(unit_col, "state", "year")))) %>%
        mutate(present = 1L),
      by = c(unit_col, "state", "year")
    ) %>%
    mutate(present = replace_na(present, 0L)) %>%
    arrange(across(all_of(unit_col)), state, year)
  
  panel <- panel %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    mutate(present_next = lead(present)) %>%
    ungroup()
  
  if (rescue_effect) {
    presence_elsewhere <- panel %>%
      group_by(across(all_of(c(unit_col, "year")))) %>%
      summarise(n_present = sum(present), .groups = "drop")
    panel <- panel %>% left_join(presence_elsewhere, by = c(unit_col, "year"))
  }
  
  transitions <- panel %>% filter(!is.na(present_next))
  
  if (rescue_effect) {
    transitions <- transitions %>% mutate(other_present = (n_present - present) > 0)
  }
  
  rates <- transitions %>%
    group_by(across(all_of(unit_col))) %>%
    summarise(
      colon_opp = if (rescue_effect) sum(present == 0 & other_present) else sum(present == 0),
      colon_succ = if (rescue_effect) sum(present == 0 & other_present & present_next == 1)
      else sum(present == 0 & present_next == 1),
      ext_opp  = sum(present == 1),
      ext_succ = sum(present == 1 & present_next == 0),
      n_states_ever = n_distinct(state[present == 1]),
      years_present_total = sum(present),
      .groups = "drop"
    ) %>%
    mutate(
      colonisation_rate = colon_succ / colon_opp,
      extinction_rate   = ext_succ / ext_opp,
      turnover_ratio     = colonisation_rate / extinction_rate
    ) %>%
    select(-colon_opp, -colon_succ, -ext_opp, -ext_succ) %>%
    arrange(desc(colonisation_rate))
  
  rates
}

metapop_titles <- metapopulation_rates(df_long, unit_col = c("textID", "title", "author_lastname"),
                                       rescue_effect = TRUE)

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
## Corpus-wide turnover over time (Jaccard/Sørensen) --------------------------#
#------------------------------------------------------------------------------#

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

turnover_titles <- turnover_time_series(df_long)
write.csv(turnover_titles, "analyses\\turnover_time_series.csv")

ggplot(turnover_titles, aes(x = year_to, y = turnover)) +
  geom_line(colour = "black") +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(colour = "black", size = 2) +
  #xlim(2006, 2026) +
  labs(x = "Year", y = "Turnover (1 - Jaccard similarity to previous year)") +
  theme_chr(base_size = 11)

ggsave("img//turnover_time_series.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## Shannon diversity per year (pooled abundance) ------------------------------#
#------------------------------------------------------------------------------#

year_text_states <- df_long %>%
  distinct(textID, state, year) %>%
  count(year, textID, name = "n_states_present")

year_community <- year_text_states %>%
  pivot_wider(names_from = textID, values_from = n_states_present, values_fill = 0) %>%
  column_to_rownames("year") %>%
  as.matrix()

shannon_by_year <- tibble(
  year = as.integer(rownames(year_community)),
  shannon_index = vegan::diversity(year_community, index = "shannon"),
  richness = rowSums(year_community > 0)
) %>%
  arrange(year)

write.csv(shannon_by_year, "analyses\\shannon_diversity_by_year.csv")

pre_2026 <- shannon_by_year %>% filter(year != 2026)
band_min  <- min(pre_2026$shannon_index, na.rm = TRUE)
band_max  <- max(pre_2026$shannon_index, na.rm = TRUE)
band_mean <- mean(pre_2026$shannon_index, na.rm = TRUE)  

ggplot(shannon_by_year, aes(x = year, y = shannon_index)) +
  geom_line(colour = "black") +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(aes(colour = year == 2026), size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "black")) +
  labs(x = "Year", y = "Shannon index") +
  theme_chr(base_size = 11)

ggsave("img//shannon_diversity_by_year.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## Colonisation/extinction rates pooled across all texts, by year -------------#
#------------------------------------------------------------------------------#

pooled_rates_by_year <- function(df_long, unit_col = "textID") {
  years <- full_seq(range(df_long$year, na.rm = TRUE), 1)
  units <- unique(df_long[[unit_col]])
  
  panel <- expand_grid(unit = units, state = state_cols, year = years) %>%
    rename(!!unit_col := unit) %>%
    left_join(
      df_long %>%
        distinct(across(all_of(c(unit_col, "state", "year")))) %>%
        mutate(present = 1L),
      by = c(unit_col, "state", "year")
    ) %>%
    mutate(present = replace_na(present, 0L)) %>%
    arrange(across(all_of(unit_col)), state, year)
  
  panel <- panel %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    mutate(present_next = lead(present), year_next = lead(year)) %>%
    ungroup() %>%
    filter(!is.na(present_next), year_next == year + 1)
  
  panel %>%
    group_by(year) %>%
    summarise(
      colonisation_rate = sum(present == 0 & present_next == 1) / sum(present == 0),
      extinction_rate    = sum(present == 1 & present_next == 0) / sum(present == 1),
      n_transitions      = n(),
      .groups = "drop"
    ) %>%
    mutate(turnover_ratio = colonisation_rate / extinction_rate) %>%
    arrange(year)
}

period_rates_titles <- pooled_rates_by_year(df_long, unit_col = "textID")
write.csv(period_rates_titles, "analyses\\period_rates_titles.csv")

ggplot(period_rates_titles, aes(x = year, y = turnover_ratio, group = 1)) +
  geom_line(colour = "black") +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey85", linewidth = 0.6) +
  geom_point(colour = "black", size = 2.5) +
  labs(x = "Year", y = "Colonization extinction ratio") +
  theme_chr(base_size = 11)

ggsave("img//period_rates.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## Contagion-dependence pooled across all texts, by year ----------------------#
#------------------------------------------------------------------------------#

pooled_contagion_by_year <- function(df_long, unit_col = "textID",
                                     gap_threshold = 0.15, low_rate_threshold = 0.15) {
  years <- full_seq(range(df_long$year, na.rm = TRUE), 1)
  units <- unique(df_long[[unit_col]])
  
  panel <- expand_grid(unit = units, state = state_cols, year = years) %>%
    rename(!!unit_col := unit) %>%
    left_join(
      df_long %>%
        distinct(across(all_of(c(unit_col, "state", "year")))) %>%
        mutate(present = 1L),
      by = c(unit_col, "state", "year")
    ) %>%
    mutate(present = replace_na(present, 0L)) %>%
    arrange(across(all_of(unit_col)), state, year)
  
  panel <- panel %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    mutate(present_next = lead(present), year_next = lead(year)) %>%
    ungroup()
  
  presence_elsewhere <- panel %>%
    group_by(across(all_of(c(unit_col, "year")))) %>%
    summarise(n_present = sum(present), .groups = "drop")
  
  opportunities <- panel %>%
    left_join(presence_elsewhere, by = c(unit_col, "year")) %>%
    mutate(other_present = (n_present - present) > 0) %>%
    filter(!is.na(present_next), year_next == year + 1, present == 0)
  
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

period_contagion_titles <- pooled_contagion_by_year(df_long, unit_col = "textID")
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
## Stickiness: survival analysis of time-to-extinction ("spells") -------------#
#------------------------------------------------------------------------------#

build_spells <- function(df_long, unit_col = c("textID", "title", "author_lastname")) {
  years <- full_seq(range(df_long$year, na.rm = TRUE), 1)
  max_year <- max(years)
  
  units <- df_long %>% distinct(across(all_of(unit_col)))
  
  panel <- units %>%
    expand_grid(state = state_cols, year = years) %>%
    left_join(
      df_long %>%
        distinct(across(all_of(c(unit_col, "state", "year")))) %>%
        mutate(present = 1L),
      by = c(unit_col, "state", "year")
    ) %>%
    mutate(present = replace_na(present, 0L)) %>%
    arrange(across(all_of(unit_col)), state, year)
  
  panel %>%
    group_by(across(all_of(c(unit_col, "state")))) %>%
    group_modify(function(.x, .y) {
      r <- rle(.x$present)
      ends <- cumsum(r$lengths)
      starts <- ends - r$lengths + 1
      
      tibble(run_value = r$values, start_idx = starts, end_idx = ends) %>%
        filter(run_value == 1) %>%
        mutate(
          start_year = .x$year[start_idx],
          end_year   = .x$year[end_idx],
          duration   = end_year - start_year + 1,
          event      = if_else(end_year < max_year, 1L, 0L)
        ) %>%
        select(start_year, end_year, duration, event)
    }) %>%
    ungroup()
}

## Persistence

spells_titles <- build_spells(df_long, unit_col = c("textID", "title", "author_lastname")) %>%
  mutate(text_id = paste(title, author_lastname, sep = " - "))

write.csv(spells_titles, "analyses\\spells_titles.csv")

title_persistence <- spells_titles %>%
  group_by(textID, title, author_lastname) %>%
  summarise(
    total_presence = sum(duration),   # combines n_spells AND spell length
    n_spells = n(),
    median_duration = median(duration),
    .groups = "drop"
  ) %>%
  left_join(
    metapop_titles %>% select(textID, title, author_lastname, colonisation_rate, n_states_ever),
    by = c("textID", "title", "author_lastname")
  ) %>%
  filter(!is.na(colonisation_rate))

presence_cutoff <- median(title_persistence$total_presence, na.rm = TRUE)

title_persistence <- title_persistence %>%
  mutate(tier = if_else(total_presence > presence_cutoff, "core", "periphery"))

# check: does colonisation rate differ meaningfully between the two tiers?
title_persistence %>%
  group_by(tier) %>%
  summarise(
    n_titles = n(),
    mean_colonisation_rate = mean(colonisation_rate, na.rm = TRUE),
    median_colonisation_rate = median(colonisation_rate, na.rm = TRUE)
  )

## Core examples

core_examples <- title_persistence %>%
  filter(tier == "core", n_states_ever >= 3) %>%
  arrange(desc(total_presence)) %>%
  select(title, author_lastname, total_presence, n_spells, median_duration, colonisation_rate, n_states_ever) %>%
  head(10)

periphery_examples <- title_persistence %>%
  filter(tier == "periphery") %>%
  arrange(total_presence) %>%
  select(title, author_lastname, total_presence, n_spells, median_duration, colonisation_rate) %>%
  head(10)

print(core_examples)
print(periphery_examples)

km_fit <- survfit(Surv(duration, event) ~ 1, data = spells_titles)

km_curve <- bind_rows(
  tibble(time = 0, surv = 1, n_censor = 0),
  tibble(time = km_fit$time, surv = km_fit$surv, n_censor = km_fit$n.censor)
)

median_survival <- km_curve %>% filter(surv <= 0.5) %>% slice(1) %>% pull(time)
if (length(median_survival) == 0) {
  message("Median spell survival time not reached within the observation window ",
          "(fewer than half of all spells have ended) - the corpus overall is ",
          "still more durable than not, over the years actually observed.")
}

ggplot(km_curve, aes(x = time, y = surv)) +
  geom_step(colour = "black") +
  geom_point(data = km_curve %>% filter(n_censor > 0), shape = 3, size = 2, colour = "black") +
  labs(x = "Years since colonization", y = "Proportion of spells still present") +
  theme_chr(base_size = 11)

ggsave("img//km_survival_curve.png", dpi = 360, width = 16, height = 12, units = "cm")

#------------------------------------------------------------------------------#
## Stickiness: Cox model -- static bibliographic covariates -------------------#
#------------------------------------------------------------------------------#

covariate_terms <- c("colonisation_rate_x10", "n_editions", "n_translations")

cox_data <- spells_titles %>%
  left_join(
    dnb_with_abitur %>% distinct(textID, title, author_lastname, literary_period),
    by = c("textID", "title", "author_lastname")
  ) %>%
  left_join(
    metapop_titles %>% select(textID, title, author_lastname, colonisation_rate),
    by = c("textID", "title", "author_lastname")
  ) %>%
  left_join(
    title_summary_resolved %>% select(textID, title, author_lastname, n_editions, n_translations),
    by = c("textID", "title", "author_lastname")
  ) %>%
  mutate(colonisation_rate_x10 = colonisation_rate * 10) %>%
  filter(complete.cases(across(all_of(covariate_terms))))
  
cox_model <- tryCatch({
  cox_formula <- as.formula(
    paste0("Surv(duration, event) ~ ", paste(covariate_terms, collapse = " + "), " + cluster(text_id)")
  )
  coxph(cox_formula, data = cox_data)
}, error = function(e) {
  message("Cox model failed to fit: ", conditionMessage(e))
  NULL
})

if (!is.null(cox_model)) {
  cox_coefficients <- as.data.frame(summary(cox_model)$coefficients) %>%
    rownames_to_column("term")
  write.csv(cox_coefficients, "analyses\\cox_model_coefficients.csv")
  
  zph_test <- tryCatch(cox.zph(cox_model), error = function(e) NULL)
  if (!is.null(zph_test)) {
    zph_table <- as.data.frame(zph_test$table) %>% rownames_to_column("term")
    write.csv(zph_table, "analyses\\cox_ph_test.csv")
  }
} else {
  message("Skipping cox_model_coefficients.csv / cox_ph_test.csv (model did not fit).")
}

#------------------------------------------------------------------------------#
## Stickiness: age-period Cox model -- cumulative bibliographic covariates ----#
#------------------------------------------------------------------------------#

build_cumulative_count <- function(data, type) {
  data_filtered <- data %>%
    semi_join(title_summary_resolved, by = c("textID", "title", "author_lastname")) %>%
    filter(entry_type == type, !is.na(pub_year), pub_year >= 1450, pub_year <= 2026)
  
  n_dropped <- sum(data$entry_type == type & is.na(data$pub_year), na.rm = TRUE)
  if (n_dropped > 0) {
    message("build_cumulative_count(", type, "): dropped ", n_dropped,
            " row(s) with no clear publication year.")
  }
  
  data_filtered %>%
    count(textID, title, author_lastname, pub_year, name = "n_new") %>%
    arrange(textID, title, author_lastname, pub_year) %>%
    group_by(textID, title, author_lastname) %>%
    mutate(cumulative = cumsum(n_new)) %>%
    ungroup() %>%
    select(textID, title, author_lastname, calendar_year = pub_year, cumulative)
}

cumulative_editions <- build_cumulative_count(dnb_with_abitur, "edition") %>%
  rename(cumulative_editions = cumulative)

cumulative_translations <- build_cumulative_count(dnb_with_abitur, "translation") %>%
  rename(cumulative_translations = cumulative)

all_years <- full_seq(c(1450, 2026), 1)

build_age_period_data <- function(spells_titles) {
  spells_titles %>%
    mutate(spell_id = row_number()) %>%
    pmap_dfr(function(textID, title, author_lastname, state, start_year, end_year,
                      duration, event, spell_id, text_id) {
      tibble(
        spell_id = spell_id, textID = textID, title = title,
        author_lastname = author_lastname, state = state,
        calendar_year = seq(start_year, end_year),
        tstart = 0:(duration - 1), tstop = 1:duration,
        event = c(rep(0L, duration - 1), event)
      )
    })
}

age_period_data <- build_age_period_data(spells_titles) %>%
  left_join(
    dnb_with_abitur %>% distinct(textID, title, author_lastname, literary_period),
    by = c("textID", "title", "author_lastname")
  ) %>%
  left_join(
    metapop_titles %>% select(textID, title, author_lastname, colonisation_rate),
    by = c("textID", "title", "author_lastname")
  ) %>%
  mutate(colonisation_rate_x10 = colonisation_rate * 10)

title_year_scaffold <- age_period_data %>%
  distinct(textID, title, author_lastname) %>%
  expand_grid(calendar_year = all_years)

cumulative_editions_full <- title_year_scaffold %>%
  left_join(cumulative_editions, by = c("textID", "title", "author_lastname", "calendar_year")) %>%
  arrange(textID, title, author_lastname, calendar_year) %>%
  group_by(textID, title, author_lastname) %>%
  fill(cumulative_editions, .direction = "down") %>%
  mutate(cumulative_editions = replace_na(cumulative_editions, 0)) %>%
  ungroup()

cumulative_translations_full <- title_year_scaffold %>%
  left_join(cumulative_translations, by = c("textID", "title", "author_lastname", "calendar_year")) %>%
  arrange(textID, title, author_lastname, calendar_year) %>%
  group_by(textID, title, author_lastname) %>%
  fill(cumulative_translations, .direction = "down") %>%
  mutate(cumulative_translations = replace_na(cumulative_translations, 0)) %>%
  ungroup()

age_period_data <- age_period_data %>%
  semi_join(title_summary_resolved, by = c("textID", "title", "author_lastname")) %>%
  left_join(cumulative_editions_full, by = c("textID", "title", "author_lastname", "calendar_year")) %>%
  left_join(cumulative_translations_full, by = c("textID", "title", "author_lastname", "calendar_year")) %>%
  mutate(text_id = paste(title, author_lastname, sep = " - ")) %>%
  filter(complete.cases(across(all_of(c("cumulative_editions", "cumulative_translations")))))

cox_model_interaction <- tryCatch({
  coxph(
    Surv(tstart, tstop, event) ~ cumulative_editions + cumulative_translations +
      colonisation_rate_x10 + cluster(text_id),
    data = age_period_data,
    id = spell_id
  )
}, error = function(e) {
  message("Interaction Cox model failed to fit: ", conditionMessage(e))
  NULL
})

if (!is.null(cox_model_interaction)) {
  interaction_coefficients <- as.data.frame(summary(cox_model_interaction)$coefficients) %>%
    rownames_to_column("term")
  write.csv(interaction_coefficients, "analyses\\cox_model_interaction_coefficients.csv")
}

exp(confint(cox_model_interaction))

#------------------------------------------------------------------------------#
## Robustness check without final years ---------------------------------------#
#------------------------------------------------------------------------------#

fit_trend_with_influence <- function(data, y_var, label, flag_year) {
  formula_obj <- as.formula(paste(y_var, "~ year"))
  model <- lm(formula_obj, data = data)
  
  tidy_model <- tidy(model) %>% filter(term == "year")
  glance_model <- glance(model)
  
  flag_idx <- which(data$year == flag_year)
  
  if (length(flag_idx) != 1) {
    warning("fit_trend_with_influence(): could not uniquely locate year ",
            flag_year, " in ", label, " -- skipping influence diagnostics.")
    dfbeta_year <- NA_real_
    dfbetas_year <- NA_real_
    cooks_d <- NA_real_
    threshold_dfbetas <- NA_real_
    threshold_cooks <- NA_real_
  } else {
    n <- nrow(data)
    dfbeta_mat <- dfbeta(model)
    dfbetas_mat <- dfbetas(model)
    cooks <- cooks.distance(model)
    
    dfbeta_year <- dfbeta_mat[flag_idx, "year"]
    dfbetas_year <- dfbetas_mat[flag_idx, "year"]
    cooks_d <- cooks[flag_idx]
    threshold_dfbetas <- 2 / sqrt(n)
    threshold_cooks <- 4 / n
  }
  
  tibble(
    dataset = label,
    flagged_year = flag_year,
    n_obs = nrow(data),
    slope_full = tidy_model$estimate,
    se_full = tidy_model$std.error,
    p_value_full = tidy_model$p.value,
    r_squared_full = glance_model$r.squared,
    dfbeta_year = dfbeta_year,
    dfbetas_year = dfbetas_year,
    dfbetas_threshold = threshold_dfbetas,
    dfbetas_flag = abs(dfbetas_year) > threshold_dfbetas,
    cooks_d = cooks_d,
    cooks_threshold = threshold_cooks,
    cooks_flag = cooks_d > threshold_cooks
  )
}

fit_trend <- function(data, y_var, label) {
  formula_obj <- as.formula(paste(y_var, "~ year"))
  model <- lm(formula_obj, data = data)
  tidy_model <- tidy(model) %>% filter(term == "year")
  glance_model <- glance(model)
  
  tibble(
    dataset = label,
    n_obs = nrow(data),
    slope = tidy_model$estimate,
    se = tidy_model$std.error,
    p_value = tidy_model$p.value,
    r_squared = glance_model$r.squared
  )
}

trend_summary <- list()
influence_summary <- list()

## 1. Shannon diversity -------------------------------------------------------#

shannon_full <- shannon_by_year
shannon_pre_endpoint <- shannon_by_year %>% filter(year != 2026)

trend_summary$shannon_full <- fit_trend(shannon_full, "shannon_index", "Shannon (full)")
trend_summary$shannon_pre_endpoint <- fit_trend(shannon_pre_endpoint, "shannon_index", "Shannon (excl. 2026)")
influence_summary$shannon_2026 <- fit_trend_with_influence(shannon_full, "shannon_index", "Shannon", flag_year = 2026)

## 2. Turnover (Jaccard dissimilarity) ----------------------------------------#

turnover_full <- turnover_titles %>% rename(year = year_to)
turnover_pre_endpoint <- turnover_full %>% filter(year != 2026)

trend_summary$turnover_full <- fit_trend(turnover_full, "turnover", "Turnover (full)")
trend_summary$turnover_pre_endpoint <- fit_trend(turnover_pre_endpoint, "turnover", "Turnover (excl. 2026)")
influence_summary$turnover_2026 <- fit_trend_with_influence(turnover_full, "turnover", "Turnover", flag_year = 2026)

## 3. Colonization-extinction ratio -------------------------------------------#

ratio_full <- period_rates_titles
ratio_pre_endpoint <- ratio_full %>% filter(year != 2025)

trend_summary$ratio_full <- fit_trend(ratio_full, "turnover_ratio", "Colon/ext ratio (full)")
trend_summary$ratio_pre_endpoint <- fit_trend(ratio_pre_endpoint, "turnover_ratio", "Colon/ext ratio (excl. 2025)")
influence_summary$ratio_2025 <- fit_trend_with_influence(ratio_full, "turnover_ratio", "Colonization-extinction ratio", flag_year = 2025)

## 4. Contagion gap -----------------------------------------------------------#

contagion_full <- period_contagion_titles
contagion_pre_endpoint <- contagion_full %>% filter(year != 2025)

trend_summary$contagion_full <- fit_trend(contagion_full, "contagion_gap", "Contagion gap (full)")
trend_summary$contagion_pre_endpoint <- fit_trend(contagion_pre_endpoint, "contagion_gap", "Contagion gap (excl. 2025)")
influence_summary$contagion_2025 <- fit_trend_with_influence(contagion_full, "contagion_gap", "Contagion gap", flag_year = 2025)

trend_summary_table <- bind_rows(trend_summary)
write.csv(trend_summary_table, "analyses\\pre_endpoint_trend_comparison.csv", row.names = FALSE)

influence_summary_table <- bind_rows(influence_summary)
write.csv(influence_summary_table, "analyses\\pre_endpoint_influence_diagnostics.csv", row.names = FALSE)

print(trend_summary_table)
print(influence_summary_table)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
