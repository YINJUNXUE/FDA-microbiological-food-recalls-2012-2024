# Validation pipeline from openFDA JSON export
# Input: food-enforcement-0001-of-0001.json
# Output: all manuscript-ready datasets, tables, figures, supplementary tables, and dictionaries

options(stringsAsFactors = FALSE)

wd <- "D:/桌面临时文件/已经发表论文/慧玲论文/FDA微生物/分析数据/分析验证"
setwd(wd)

input_json <- "food-enforcement-0001-of-0001.json"
out_dir <- "validation_outputs_from_json"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_json)) stop("Missing input JSON file.")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package jsonlite is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package ggplot2 is required.")

library(jsonlite)
library(ggplot2)

message("Reading JSON: ", input_json)
json_obj <- fromJSON(input_json, flatten = TRUE)
df_raw <- json_obj$results
if (!is.data.frame(df_raw)) stop("JSON results could not be converted to a data frame.")

# Drop list columns such as openfda when present, so CSV export is stable.
is_list_col <- vapply(df_raw, is.list, logical(1))
if (any(is_list_col)) df_raw <- df_raw[, !is_list_col, drop = FALSE]

write.csv(df_raw, file.path(out_dir, "food_enforcement_from_json.csv"), row.names = FALSE, fileEncoding = "UTF-8")

key_vars <- c(
  "classification", "product_type", "reason_for_recall", "product_description",
  "distribution_pattern", "recall_initiation_date", "recalling_firm", "event_id"
)
missing_vars <- setdiff(key_vars, names(df_raw))
if (length(missing_vars) > 0) stop(paste("Missing key variables:", paste(missing_vars, collapse = ", ")))

norm_text <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

has <- function(text, pattern) grepl(pattern, text, perl = TRUE)

# 1. Rebuild microbial analysis dataset ---------------------------------------

n_raw <- nrow(df_raw)
df_food <- df_raw[tolower(trimws(as.character(df_raw$product_type))) == "food", , drop = FALSE]
n_food <- nrow(df_food)

reason_food_lc <- tolower(norm_text(df_food$reason_for_recall))
microbial_pattern <- paste(
  c(
    "listeria", "monocytogenes", "salmonella",
    "\\be\\.?\\s*coli\\b", "escherichia\\s+coli", "\\bstec\\b", "shiga",
    "\\bo157\\b", "\\bo26\\b", "\\bo45\\b", "\\bo103\\b", "\\bo111\\b", "\\bo121\\b", "\\bo145\\b",
    "clostridium", "botulinum", "botulism", "hepatitis\\s*a", "\\bhav\\b",
    "norovirus", "noro\\s*virus", "campylobacter", "staphylococcus", "staph\\.?\\s*aureus",
    "cyclospora", "cronobacter", "sakazakii", "bacillus\\s+cereus", "\\bvibrio\\b",
    "pathogen", "pathogenic", "microbial", "microbiological", "micro-organism", "microorganism",
    "bacteria", "bacterial", "bacterium", "mold", "mould", "yeast", "fungal", "fungus",
    "spoilage", "spoiled", "spoil", "temperature abuse", "temperature control",
    "temperature excursion", "time/temperature", "time and temperature", "refrigerat",
    "cold chain", "held at improper temperature", "commercial sterility",
    "underprocess", "under-process", "under processed", "inadequate processing",
    "contamination.*organism", "organism.*contamination"
  ),
  collapse = "|"
)

df_micro <- df_food[has(reason_food_lc, microbial_pattern), , drop = FALSE]
n_micro <- nrow(df_micro)

reason_micro_lc <- tolower(norm_text(df_micro$reason_for_recall))
desc_micro_lc <- tolower(norm_text(df_micro$product_description))
combined_lc <- paste(reason_micro_lc, desc_micro_lc)

pathogen_cat <- rep("Unspecified microbial contamination", nrow(df_micro))
pathogen_cat[has(combined_lc, "listeria\\s+monocytogenes|l\\.?\\s*monocytogenes")] <- "Listeria monocytogenes"
pathogen_cat[has(combined_lc, "salmonella")] <- "Salmonella"
pathogen_cat[has(combined_lc, "escherichia\\s+coli|\\be\\.?\\s*coli\\b|\\bstec\\b|shiga\\s*toxin|o157|o26|o45|o103|o111|o121|o145")] <- "Escherichia coli"
pathogen_cat[has(combined_lc, "clostridium\\s+botulinum|botulinum|botulism")] <- "Clostridium botulinum"
pathogen_cat[has(combined_lc, "hepatitis\\s*a|\\bhav\\b")] <- "Hepatitis A"
pathogen_cat[has(combined_lc, "norovirus|noro\\s*virus")] <- "Norovirus"
pathogen_cat[has(combined_lc, "campylobacter")] <- "Campylobacter"
pathogen_cat[has(combined_lc, "staphylococcus\\s+aureus|staph\\.?\\s*aureus")] <- "Staphylococcus aureus"
pathogen_cat[has(combined_lc, "cyclospora")] <- "Cyclospora"
pathogen_cat[has(combined_lc, "cronobacter|sakazakii")] <- "Cronobacter"
pathogen_cat[has(combined_lc, "bacillus\\s+cereus")] <- "Bacillus cereus"
pathogen_cat[has(combined_lc, "\\bvibrio\\b|vibrio\\s+(vulnificus|parahaemolyticus|cholerae)")] <- "Vibrio"
other_specified <- "yersinia|shigella|brucella|parasite|parasitic|cryptosporidium|giardia|toxoplasma|mycobacterium|pseudomonas|enterobacter|klebsiella|proteus|scombrotoxin|histamine"
pathogen_cat[pathogen_cat == "Unspecified microbial contamination" & has(combined_lc, other_specified)] <- "Other specified pathogens"
df_micro$pathogen_cat <- pathogen_cat

df_micro$class1 <- ifelse(
  trimws(df_micro$classification) == "Class I", 1,
  ifelse(trimws(df_micro$classification) %in% c("Class II", "Class III"), 0, NA)
)
df_micro$nationwide <- ifelse(has(tolower(norm_text(df_micro$distribution_pattern)), "\\bnationwide\\b"), 1, 0)

date_chr <- norm_text(df_micro$recall_initiation_date)
year_from_8digit <- ifelse(grepl("^\\d{8}$", date_chr), substr(date_chr, 1, 4), NA)
year_from_iso <- ifelse(grepl("^\\d{4}[-/]\\d{1,2}[-/]\\d{1,2}$", date_chr), substr(date_chr, 1, 4), NA)
year_from_us <- ifelse(grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", date_chr), sub("^\\d{1,2}/\\d{1,2}/(\\d{4})$", "\\1", date_chr), NA)
df_micro$year <- suppressWarnings(as.integer(ifelse(!is.na(year_from_8digit), year_from_8digit, ifelse(!is.na(year_from_iso), year_from_iso, year_from_us))))

main5 <- c("Listeria monocytogenes", "Salmonella", "Escherichia coli", "Clostridium botulinum", "Hepatitis A")
df_micro$pathogen_main <- ifelse(df_micro$pathogen_cat %in% main5, df_micro$pathogen_cat, NA)
df_micro$pathogen_main <- factor(df_micro$pathogen_main, levels = main5)
df_main5 <- df_micro[!is.na(df_micro$pathogen_main), , drop = FALSE]

write.csv(df_micro, file.path(out_dir, "food_enforcement_microbial_analysis.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(df_main5, file.path(out_dir, "food_enforcement_microbial_main5.csv"), row.names = FALSE, fileEncoding = "UTF-8")

table_main_pathogen <- aggregate(
  cbind(n = rep(1, nrow(df_main5)), class1_n = df_main5$class1 == 1, nationwide_n = df_main5$nationwide == 1) ~ pathogen_main,
  data = df_main5,
  FUN = sum
)
table_main_pathogen$pct <- round(table_main_pathogen$n / sum(table_main_pathogen$n) * 100, 2)
table_main_pathogen$class1_pct <- round(table_main_pathogen$class1_n / table_main_pathogen$n * 100, 2)
table_main_pathogen$nationwide_pct <- round(table_main_pathogen$nationwide_n / table_main_pathogen$n * 100, 2)
table_main_pathogen <- table_main_pathogen[, c("pathogen_main", "n", "pct", "class1_n", "class1_pct", "nationwide_n", "nationwide_pct")]
write.csv(table_main_pathogen, file.path(out_dir, "table_main_pathogen_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 2. Product category mapping --------------------------------------------------

desc_lc <- tolower(norm_text(df_main5$product_description))
reason_lc <- tolower(norm_text(df_main5$reason_for_recall))
all_lc <- paste(desc_lc, reason_lc)
has_desc <- function(pattern) grepl(pattern, desc_lc, perl = TRUE)
has_any <- function(pattern) grepl(pattern, all_lc, perl = TRUE)

product_cat <- rep("Other foods", nrow(df_main5))
dietary <- "dietary supplement|supplement|vitamin|multivitamin|capsule|tablet|softgel|soft gel|pill|protein powder|whey protein|kratom|botanical|herbal supplement|moringa powder|spirulina|chlorella"
meat <- "\\bbeef\\b|\\bpork\\b|\\bham\\b|\\bbacon\\b|sausage|salami|pepperoni|prosciutto|mortadella|bologna|hot dog|frankfurter|jerky|\\bmeat\\b|meatball|meatloaf|steak|veal|lamb|chorizo|chicken|turkey|poultry|\\bduck\\b|\\bgoose\\b|roast beef|deli meat|lunch meat|charcuterie"
seafood <- "seafood|\\bfish\\b|salmon|tuna|cod|haddock|whitefish|trout|tilapia|catfish|sardine|anchov|herring|mackerel|shrimp|prawn|crab|lobster|crawfish|crayfish|oyster|clam|mussel|scallop|squid|octopus|surimi"
dairy <- "\\bmilk\\b|cheese|cheddar|mozzarella|parmesan|pecorino|romano|gouda|brie|queso|yogurt|yoghurt|ice cream|gelato|sorbet|butter|\\bcream\\b|creamer|sour cream|cottage cheese|ricotta|whey|casein|dairy"
produce <- "fresh produce|lettuce|romaine|spinach|kale|arugula|cabbage|collard|greens|spring mix|salad mix|iceberg|celery|carrot|broccoli|cauliflower|cucumber|zucchini|squash|pepper|tomato|onion|mushroom|sprout|alfalfa|cilantro|parsley|basil|fruit|melon|cantaloupe|watermelon|honeydew|apple|pear|peach|plum|nectarine|grape|berry|berries|strawberry|blueberry|raspberry|blackberry|mango|papaya|pineapple|avocado|orange|lemon|lime|bean sprout|vegetable|veggie"
rte <- "ready[ -]?to[ -]?eat|\\brte\\b|prepared food|prepared meal|meal kit|entree|dinner|lunch|sandwich|wrap|submarine|\\bsub\\b|hoagie|panini|burrito|taco|quesadilla|pizza|salad|potato salad|pasta salad|chicken salad|egg salad|tuna salad|coleslaw|slaw|salsa|\\bdip\\b|spread|hummus|pate|sushi|meal|kit|tray|platter|side dish|stuffed|soup|chili|casserole|macaroni|lasagna|noodle|pasta"
produce_to_rte <- "salad kit|salad bowl|salad cup|potato salad|pasta salad|chicken salad|egg salad|tuna salad|veggie tray|vegetable tray|fruit tray|party tray|platter|with dip|with ranch|salsa|coleslaw|slaw|sandwich|wrap|meal kit"
bakery <- "bread|bagel|roll|bun|bakery|cake|cupcake|muffin|cookie|brownie|pie|pastry|doughnut|donut|croissant|tortilla|pita|waffle|pancake|batter mix|cake mix|brownie mix|muffin mix|bread mix|\\bflour\\b"
beverage <- "beverage|\\bdrink\\b|juice|smoothie|tea|coffee|kombucha|lemonade|cider|water beverage|bottled water|soda|shake|drink mix|non-alcoholic drink|cocktail mix"
nuts <- "peanut|almond|cashew|pistachio|walnut|pecan|hazelnut|macadamia|brazil nut|pine nut|mixed nut|\\bnuts\\b|sunflower seed|pumpkin seed|sesame|chia|flax|hemp seed|\\bseed\\b|\\bseeds\\b|nut butter|peanut butter|almond butter|tahini"
frozen <- "frozen|freezer|keep frozen|individually quick frozen|\\biqf\\b"
condiment <- "sauce|salsa|dressing|condiment|marinade|seasoning|spice|spices|herb|herbs|pepper powder|chili powder|curry|paprika|cumin|fennel|garlic powder|onion powder|relish|chutney|mustard|ketchup|mayonnaise|\\bmayo\\b|vinegar|gravy|taco seasoning|dip mix"

product_cat[has_any(dietary)] <- "Dietary supplements"
product_cat[product_cat == "Other foods" & has_desc(meat)] <- "Meat and poultry"
product_cat[product_cat == "Other foods" & has_desc(seafood)] <- "Seafood"
product_cat[product_cat == "Other foods" & has_desc(dairy)] <- "Dairy"
product_cat[product_cat == "Other foods" & has_desc(produce)] <- "Fresh produce"
product_cat[product_cat == "Other foods" & has_desc(rte)] <- "Ready-to-eat / prepared foods"
product_cat[product_cat == "Fresh produce" & has_desc(produce_to_rte)] <- "Ready-to-eat / prepared foods"
product_cat[product_cat == "Other foods" & has_desc(bakery)] <- "Bakery"
product_cat[product_cat == "Other foods" & has_desc(beverage)] <- "Beverages"
product_cat[product_cat == "Other foods" & has_desc(nuts)] <- "Nuts and seeds"
product_cat[product_cat == "Other foods" & has_desc(frozen)] <- "Frozen foods"
product_cat[product_cat == "Other foods" & has_desc(condiment)] <- "Condiments / sauces"
product_cat[product_cat == "Other foods" & has_any("raw milk|pasteurized milk|cheese|ice cream|dairy")] <- "Dairy"
product_cat[product_cat == "Other foods" & has_any("seafood|fish|shrimp|oyster|clam|salmon|tuna")] <- "Seafood"
product_cat[product_cat == "Other foods" & has_any("peanut|almond|pistachio|sesame|seed")] <- "Nuts and seeds"

product_levels <- c("Dairy", "Meat and poultry", "Seafood", "Fresh produce", "Ready-to-eat / prepared foods", "Bakery", "Beverages", "Nuts and seeds", "Frozen foods", "Condiments / sauces", "Dietary supplements", "Other foods")
df_main5$product_cat <- factor(product_cat, levels = product_levels)
write.csv(df_main5, file.path(out_dir, "food_enforcement_microbial_main5_productcat.csv"), row.names = FALSE, fileEncoding = "UTF-8")

table_product_cat_summary <- as.data.frame(sort(table(df_main5$product_cat), decreasing = TRUE))
names(table_product_cat_summary) <- c("product_cat", "n")
table_product_cat_summary$pct <- round(table_product_cat_summary$n / nrow(df_main5) * 100, 2)
write.csv(table_product_cat_summary, file.path(out_dir, "table_product_cat_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pathogen_product_crosstab <- as.data.frame.matrix(table(df_main5$pathogen_main, df_main5$product_cat))
pathogen_product_crosstab <- data.frame(pathogen_main = row.names(pathogen_product_crosstab), pathogen_product_crosstab, row.names = NULL, check.names = FALSE)
write.csv(pathogen_product_crosstab, file.path(out_dir, "table_pathogen_product_crosstab.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 3. Final 2012-2024 analysis --------------------------------------------------

df_sub <- df_main5[!is.na(df_main5$year) & df_main5$year >= 2012 & df_main5$year <= 2024, , drop = FALSE]
df_sub$pathogen_main <- factor(as.character(df_sub$pathogen_main), levels = main5)
df_sub$product_cat <- factor(as.character(df_sub$product_cat), levels = product_levels)
write.csv(df_sub, file.path(out_dir, "food_enforcement_microbial_main5_productcat_2012_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")

make_grid <- function() expand.grid(pathogen_main = main5, product_cat = product_levels, stringsAsFactors = FALSE)
fmt_pct <- function(x) ifelse(is.na(x), "", sprintf("%.1f%%", x))

theme_heatmap <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_blank(), legend.position = "right",
          plot.title = element_text(face = "bold", size = 13))
}

fig2_counts <- as.data.frame(table(pathogen_main = df_sub$pathogen_main, product_cat = df_sub$product_cat), stringsAsFactors = FALSE)
names(fig2_counts) <- c("pathogen_main", "product_cat", "n")
fig2_data <- merge(make_grid(), fig2_counts, by = c("pathogen_main", "product_cat"), all.x = TRUE)
fig2_data$n[is.na(fig2_data$n)] <- 0
row_totals <- aggregate(n ~ pathogen_main, data = fig2_data, sum)
names(row_totals)[2] <- "pathogen_total"
fig2_data <- merge(fig2_data, row_totals, by = "pathogen_main", all.x = TRUE)
fig2_data$row_pct <- round(fig2_data$n / fig2_data$pathogen_total * 100, 2)
fig2_data$label <- fmt_pct(fig2_data$row_pct)
fig2_data$pathogen_main <- factor(fig2_data$pathogen_main, levels = main5)
fig2_data$product_cat <- factor(fig2_data$product_cat, levels = product_levels)
fig2_data <- fig2_data[order(fig2_data$pathogen_main, fig2_data$product_cat), ]

combo_n <- aggregate(list(n = rep(1, nrow(df_sub))), by = list(pathogen_main = df_sub$pathogen_main, product_cat = df_sub$product_cat), FUN = sum)
combo_class1 <- aggregate(list(class1_n = df_sub$class1 == 1), by = list(pathogen_main = df_sub$pathogen_main, product_cat = df_sub$product_cat), FUN = function(x) sum(x, na.rm = TRUE))
fig3_data <- merge(combo_n, combo_class1, by = c("pathogen_main", "product_cat"), all.x = TRUE)
fig3_data$class1_pct <- round(fig3_data$class1_n / fig3_data$n * 100, 2)
fig3_data <- fig3_data[fig3_data$n >= 20, , drop = FALSE]
fig3_data$label <- fmt_pct(fig3_data$class1_pct)

combo_nat <- aggregate(list(nationwide_n = df_sub$nationwide == 1), by = list(pathogen_main = df_sub$pathogen_main, product_cat = df_sub$product_cat), FUN = function(x) sum(x, na.rm = TRUE))
fig4_data <- merge(combo_n, combo_nat, by = c("pathogen_main", "product_cat"), all.x = TRUE)
fig4_data$nationwide_pct <- round(fig4_data$nationwide_n / fig4_data$n * 100, 2)
fig4_data <- fig4_data[fig4_data$n >= 20, , drop = FALSE]
fig4_data$label <- fmt_pct(fig4_data$nationwide_pct)

write.csv(fig2_data, file.path(out_dir, "figure2_distribution_data_2012_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fig3_data, file.path(out_dir, "figure3_class1_data_2012_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fig4_data, file.path(out_dir, "figure4_nationwide_data_2012_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p2 <- ggplot(fig2_data, aes(product_cat, pathogen_main, fill = row_pct)) + geom_tile(color = "white") + geom_text(aes(label = label), size = 3) + scale_fill_gradient(low = "#f7fbff", high = "#2166ac", name = "Row %") + labs(title = "Figure 2. Product distribution by pathogen, 2012-2024") + theme_heatmap()
p3 <- ggplot(fig3_data, aes(product_cat, pathogen_main, fill = class1_pct)) + geom_tile(color = "white") + geom_text(aes(label = label), size = 3) + scale_fill_gradient(low = "#fff7ec", high = "#b30000", name = "Class I %") + labs(title = "Figure 3. Class I percentage, 2012-2024") + theme_heatmap()
p4 <- ggplot(fig4_data, aes(product_cat, pathogen_main, fill = nationwide_pct)) + geom_tile(color = "white") + geom_text(aes(label = label), size = 3) + scale_fill_gradient(low = "#f7fcf5", high = "#238b45", name = "Nationwide %") + labs(title = "Figure 4. Nationwide percentage, 2012-2024") + theme_heatmap()
ggsave(file.path(out_dir, "Figure2_pathogen_product_distribution_heatmap_2012_2024.png"), p2, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "Figure3_class1_heatmap_2012_2024.png"), p3, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "Figure4_nationwide_heatmap_2012_2024.png"), p4, width = 12, height = 5.5, dpi = 300)

# Logistic models
model_df <- df_sub[complete.cases(df_sub[, c("class1", "nationwide", "pathogen_main", "product_cat", "year")]), ]
model_class1 <- glm(class1 ~ pathogen_main + product_cat + year, data = model_df, family = binomial())
model_nationwide <- glm(nationwide ~ pathogen_main + product_cat + year, data = model_df, family = binomial())

clean_term <- function(term) {
  variable <- ifelse(grepl("^pathogen_main", term), "pathogen_main",
                     ifelse(grepl("^product_cat", term), "product_cat",
                            ifelse(term == "year", "year", "Intercept")))
  level <- term
  level <- sub("^pathogen_main", "", level)
  level <- sub("^product_cat", "", level)
  level[level == "year"] <- "per 1-year increase"
  data.frame(variable = variable, level = level, stringsAsFactors = FALSE)
}

model_table <- function(model, outcome) {
  sm <- summary(model)$coefficients
  out <- data.frame(outcome = outcome, term = row.names(sm), estimate = sm[, 1], std_error = sm[, 2], z_value = sm[, 3], p_value = sm[, 4], row.names = NULL)
  out$OR <- exp(out$estimate)
  out$CI_low <- exp(out$estimate - 1.96 * out$std_error)
  out$CI_high <- exp(out$estimate + 1.96 * out$std_error)
  out <- cbind(out[, c("outcome", "term")], clean_term(out$term), out[, c("estimate", "std_error", "z_value", "OR", "CI_low", "CI_high", "p_value")])
  out <- out[out$term != "(Intercept)", , drop = FALSE]
  out$OR_round <- round(out$OR, 3)
  out$CI_low_round <- round(out$CI_low, 3)
  out$CI_high_round <- round(out$CI_high, 3)
  out$p_value_round <- signif(out$p_value, 3)
  out
}

table_model_class1 <- model_table(model_class1, "class1")
table_model_nationwide <- model_table(model_nationwide, "nationwide")
write.csv(table_model_class1, file.path(out_dir, "table_model_class1_2012_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table_model_nationwide, file.path(out_dir, "table_model_nationwide_2012_2024.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 4. Supplementary analysis ----------------------------------------------------

summarise_rate <- function(dat, group_vars, outcome_var) {
  dat$y <- dat[[outcome_var]]
  dat$n_unit <- 1L
  n_tab <- aggregate(as.formula(paste("n_unit ~", paste(group_vars, collapse = " + "))), data = dat, FUN = sum)
  y_tab <- aggregate(as.formula(paste("y ~", paste(group_vars, collapse = " + "))), data = dat, FUN = function(x) sum(x == 1, na.rm = TRUE))
  names(n_tab)[names(n_tab) == "n_unit"] <- "n"
  names(y_tab)[names(y_tab) == "y"] <- paste0(outcome_var, "_n")
  out <- merge(n_tab, y_tab, by = group_vars, all.x = TRUE)
  out[[paste0(outcome_var, "_pct")]] <- round(out[[paste0(outcome_var, "_n")]] / out$n * 100, 2)
  out
}

trend_cp <- summarise_rate(df_sub, c("year", "pathogen_main"), "class1")
trend_np <- summarise_rate(df_sub, c("year", "pathogen_main"), "nationwide")
trend_cprod <- summarise_rate(df_sub, c("year", "product_cat"), "class1")
trend_nprod <- summarise_rate(df_sub, c("year", "product_cat"), "nationwide")
write.csv(trend_cp, file.path(out_dir, "Supplementary_TimeTrend_ClassI_by_pathogen_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(trend_np, file.path(out_dir, "Supplementary_TimeTrend_Nationwide_by_pathogen_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(trend_cprod, file.path(out_dir, "Supplementary_TimeTrend_ClassI_by_product_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(trend_nprod, file.path(out_dir, "Supplementary_TimeTrend_Nationwide_by_product_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

theme_line <- function() theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "right", plot.title = element_text(face = "bold"))
g1 <- ggplot(trend_cp, aes(year, class1_pct, color = pathogen_main, group = pathogen_main)) + geom_line() + geom_point() + scale_x_continuous(breaks = 2012:2024) + scale_y_continuous(limits = c(0, 100)) + labs(title = "Class I percentage by pathogen", y = "Class I (%)") + theme_line()
g2 <- ggplot(trend_np, aes(year, nationwide_pct, color = pathogen_main, group = pathogen_main)) + geom_line() + geom_point() + scale_x_continuous(breaks = 2012:2024) + scale_y_continuous(limits = c(0, 100)) + labs(title = "Nationwide percentage by pathogen", y = "Nationwide (%)") + theme_line()
g3 <- ggplot(trend_cprod, aes(year, class1_pct, color = product_cat, group = product_cat)) + geom_line() + geom_point() + scale_x_continuous(breaks = 2012:2024) + scale_y_continuous(limits = c(0, 100)) + labs(title = "Class I percentage by product", y = "Class I (%)") + theme_line()
g4 <- ggplot(trend_nprod, aes(year, nationwide_pct, color = product_cat, group = product_cat)) + geom_line() + geom_point() + scale_x_continuous(breaks = 2012:2024) + scale_y_continuous(limits = c(0, 100)) + labs(title = "Nationwide percentage by product", y = "Nationwide (%)") + theme_line()
ggsave(file.path(out_dir, "TimeTrend_ClassI_by_pathogen.png"), g1, width = 10, height = 6, dpi = 300)
ggsave(file.path(out_dir, "TimeTrend_Nationwide_by_pathogen.png"), g2, width = 10, height = 6, dpi = 300)
ggsave(file.path(out_dir, "TimeTrend_ClassI_by_product.png"), g3, width = 12, height = 7, dpi = 300)
ggsave(file.path(out_dir, "TimeTrend_Nationwide_by_product.png"), g4, width = 12, height = 7, dpi = 300)

# Combination OR and interaction heatmaps
mod_c_int <- glm(class1 ~ pathogen_main * product_cat + year, data = model_df, family = binomial())
mod_n_int <- glm(nationwide ~ pathogen_main * product_cat + year, data = model_df, family = binomial())
combo_counts <- as.data.frame(table(pathogen_main = model_df$pathogen_main, product_cat = model_df$product_cat), stringsAsFactors = FALSE)
names(combo_counts) <- c("pathogen_main", "product_cat", "n")
combo_events <- aggregate(cbind(class1_n = class1, nationwide_n = nationwide) ~ pathogen_main + product_cat, data = model_df, FUN = function(x) sum(x == 1, na.rm = TRUE))
combo_counts <- merge(combo_counts, combo_events, by = c("pathogen_main", "product_cat"), all.x = TRUE)

combo_or <- function(model, outcome) {
  ref <- data.frame(pathogen_main = factor(main5[1], levels = main5), product_cat = factor(product_levels[1], levels = product_levels), year = median(model_df$year))
  grid <- make_grid()
  grid$year <- median(model_df$year)
  grid$pathogen_main <- factor(grid$pathogen_main, levels = main5)
  grid$product_cat <- factor(grid$product_cat, levels = product_levels)
  xg <- model.matrix(delete.response(terms(model)), grid)
  xr <- model.matrix(delete.response(terms(model)), ref)
  beta <- coef(model)
  keep <- !is.na(beta)
  beta <- beta[keep]
  vc <- vcov(model)[keep, keep, drop = FALSE]
  xg <- xg[, keep, drop = FALSE]
  xr <- xr[, keep, drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    contrast <- xg[i, ] - xr[1, ]
    log_or <- as.numeric(sum(contrast * beta))
    se <- as.numeric(sqrt(t(contrast) %*% vc %*% contrast))
    p <- ifelse(se == 0, NA, 2 * pnorm(abs(log_or / se), lower.tail = FALSE))
    data.frame(outcome = outcome, pathogen_main = as.character(grid$pathogen_main[i]), product_cat = as.character(grid$product_cat[i]), log_OR = log_or, SE = se, OR = exp(log_or), CI_low = exp(log_or - 1.96 * se), CI_high = exp(log_or + 1.96 * se), p_value = p)
  }))
  out <- merge(out, combo_counts, by = c("pathogen_main", "product_cat"), all.x = TRUE)
  out$event_n <- if (outcome == "class1") out$class1_n else out$nationwide_n
  out$OR_round <- round(out$OR, 3)
  out$CI_low_round <- round(out$CI_low, 3)
  out$CI_high_round <- round(out$CI_high, 3)
  out$p_value_text <- ifelse(is.na(out$p_value), NA, ifelse(out$p_value < 0.001, "<0.001", sprintf("%.3f", out$p_value)))
  out[order(factor(out$pathogen_main, levels = main5), factor(out$product_cat, levels = product_levels)), ]
}

sup_or_class1 <- combo_or(mod_c_int, "class1")
sup_or_nat <- combo_or(mod_n_int, "nationwide")
write.csv(sup_or_class1, file.path(out_dir, "SupplementaryTable_OR_ClassI.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sup_or_nat, file.path(out_dir, "SupplementaryTable_OR_Nationwide.csv"), row.names = FALSE, fileEncoding = "UTF-8")

interaction_data <- function(model, outcome) {
  cn <- names(coef(model))
  vc <- vcov(model)
  out <- make_grid()
  out$interaction_OR <- 1
  out$p_value <- NA_real_
  for (i in seq_len(nrow(out))) {
    if (out$pathogen_main[i] == main5[1] || out$product_cat[i] == product_levels[1]) next
    t1 <- paste0("pathogen_main", out$pathogen_main[i], ":product_cat", out$product_cat[i])
    t2 <- paste0("product_cat", out$product_cat[i], ":pathogen_main", out$pathogen_main[i])
    term <- if (t1 %in% cn) t1 else if (t2 %in% cn) t2 else NA
    if (!is.na(term) && !is.na(coef(model)[term])) {
      b <- coef(model)[term]
      se <- sqrt(vc[term, term])
      out$interaction_OR[i] <- exp(b)
      out$p_value[i] <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
    }
  }
  out <- merge(out, combo_counts[, c("pathogen_main", "product_cat", "n")], by = c("pathogen_main", "product_cat"), all.x = TRUE)
  out$outcome <- outcome
  out$label <- sprintf("%.2f", out$interaction_OR)
  out$pathogen_main <- factor(out$pathogen_main, levels = main5)
  out$product_cat <- factor(out$product_cat, levels = product_levels)
  out[order(out$pathogen_main, out$product_cat), ]
}

int_class1 <- interaction_data(mod_c_int, "class1")
int_nat <- interaction_data(mod_n_int, "nationwide")
write.csv(int_class1, file.path(out_dir, "InteractionHeatmap_ClassI_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(int_nat, file.path(out_dir, "InteractionHeatmap_Nationwide_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
ih1 <- ggplot(int_class1, aes(product_cat, pathogen_main, fill = interaction_OR)) + geom_tile(color = "white") + geom_text(aes(label = label), size = 3) + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 1, name = "Interaction OR") + labs(title = "Interaction OR heatmap for Class I") + theme_heatmap()
ih2 <- ggplot(int_nat, aes(product_cat, pathogen_main, fill = interaction_OR)) + geom_tile(color = "white") + geom_text(aes(label = label), size = 3) + scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 1, name = "Interaction OR") + labs(title = "Interaction OR heatmap for nationwide recalls") + theme_heatmap()
ggsave(file.path(out_dir, "InteractionHeatmap_ClassI.png"), ih1, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "InteractionHeatmap_Nationwide.png"), ih2, width = 12, height = 5.5, dpi = 300)

# 5. Publication-ready dictionaries -------------------------------------------

pathogen_pub <- data.frame(
  Mapping_type = "Pathogen classification",
  Source_fields = "reason_for_recall and product_description",
  Priority_order = seq_len(14),
  Mapped_category = c(main5, "Norovirus", "Campylobacter", "Staphylococcus aureus", "Cyclospora", "Cronobacter", "Bacillus cereus", "Vibrio", "Other specified pathogens", "Unspecified microbial contamination"),
  Main_analysis_category = c(main5, rep("Not included in main 5-pathogen analysis", 9)),
  Keywords_used_for_mapping = c(
    "Listeria monocytogenes; L. monocytogenes",
    "Salmonella",
    "Escherichia coli; E. coli; STEC; Shiga toxin; O157; O26; O45; O103; O111; O121; O145",
    "Clostridium botulinum; botulinum; botulism",
    "Hepatitis A; HAV",
    "Norovirus; noro virus",
    "Campylobacter",
    "Staphylococcus aureus; Staph. aureus",
    "Cyclospora",
    "Cronobacter; Cronobacter sakazakii; sakazakii",
    "Bacillus cereus",
    "Vibrio; Vibrio vulnificus; Vibrio parahaemolyticus; Vibrio cholerae",
    "Yersinia; Shigella; Brucella; parasites; Cryptosporidium; Giardia; Toxoplasma; Mycobacterium; Pseudomonas; Enterobacter; Klebsiella; Proteus; scombrotoxin; histamine",
    "Microbial contamination, pathogen growth, spoilage, mold/yeast/fungal contamination, temperature abuse, refrigeration/cold-chain problems, commercial sterility, underprocessing, or other microbial-risk descriptions without a specific pathogen"
  ),
  Included_in_main_analysis = c(rep("Yes", 5), rep("No", 9))
)
write.csv(pathogen_pub, file.path(out_dir, "dictionary_pathogen_mapping_publication_ready.csv"), row.names = FALSE, fileEncoding = "UTF-8")

product_pub <- data.frame(
  Mapping_type = "Product category classification",
  Source_fields = "product_description primarily; reason_for_recall used as supplemental information when needed",
  Priority_order = seq_len(13),
  Product_category = c("Dietary supplements", "Meat and poultry", "Seafood", "Dairy", "Fresh produce", "Ready-to-eat / prepared foods", "Ready-to-eat / prepared foods override", "Bakery", "Beverages", "Nuts and seeds", "Frozen foods", "Condiments / sauces", "Other foods"),
  Keywords_used_for_mapping = c(
    "Dietary supplement; supplement; vitamin; multivitamin; capsule; tablet; softgel; pill; protein powder; whey protein; kratom; botanical; herbal supplement; moringa powder; spirulina; chlorella",
    "Beef; pork; ham; bacon; sausage; salami; pepperoni; prosciutto; mortadella; bologna; hot dog; frankfurter; jerky; meat; meatball; meatloaf; steak; veal; lamb; chorizo; chicken; turkey; poultry; duck; goose; roast beef; deli meat; lunch meat; charcuterie",
    "Seafood; fish; salmon; tuna; cod; haddock; whitefish; trout; tilapia; catfish; sardine; anchovy; herring; mackerel; shrimp; prawn; crab; lobster; crawfish; crayfish; oyster; clam; mussel; scallop; squid; octopus; surimi",
    "Milk; cheese; cheddar; mozzarella; parmesan; pecorino; romano; gouda; brie; queso; yogurt; yoghurt; ice cream; gelato; sorbet; butter; cream; creamer; sour cream; cottage cheese; ricotta; whey; casein; dairy",
    "Fresh produce; lettuce; romaine; spinach; kale; arugula; cabbage; collard greens; spring mix; salad mix; iceberg lettuce; celery; carrot; broccoli; cauliflower; cucumber; zucchini; squash; pepper; tomato; onion; mushroom; sprouts; alfalfa; cilantro; parsley; basil; fruit; melon; cantaloupe; watermelon; honeydew; apple; pear; peach; plum; nectarine; grape; berries; strawberry; blueberry; raspberry; blackberry; mango; papaya; pineapple; avocado; orange; lemon; lime; bean sprouts; vegetables; veggie",
    "Ready-to-eat; RTE; prepared food; prepared meal; meal kit; entree; dinner; lunch; sandwich; wrap; submarine sandwich; hoagie; panini; burrito; taco; quesadilla; pizza; salad; potato salad; pasta salad; chicken salad; egg salad; tuna salad; coleslaw; slaw; salsa; dip; spread; hummus; pate; sushi; meal; kit; tray; platter; side dish; stuffed products; soup; chili; casserole; macaroni; lasagna; noodle; pasta",
    "Salad kit; salad bowl; salad cup; potato salad; pasta salad; chicken salad; egg salad; tuna salad; veggie tray; vegetable tray; fruit tray; party tray; platter; with dip; with ranch; salsa; coleslaw; slaw; sandwich; wrap; meal kit",
    "Bread; bagel; roll; bun; bakery; cake; cupcake; muffin; cookie; brownie; pie; pastry; doughnut; donut; croissant; tortilla; pita; waffle; pancake; batter mix; cake mix; brownie mix; muffin mix; bread mix; flour",
    "Beverage; drink; juice; smoothie; tea; coffee; kombucha; lemonade; cider; bottled water; soda; shake; drink mix; non-alcoholic drink; cocktail mix",
    "Peanut; almond; cashew; pistachio; walnut; pecan; hazelnut; macadamia; Brazil nut; pine nut; mixed nuts; sunflower seed; pumpkin seed; sesame; chia; flax; hemp seed; seeds; nut butter; peanut butter; almond butter; tahini",
    "Frozen; freezer; keep frozen; individually quick frozen; IQF",
    "Sauce; salsa; dressing; condiment; marinade; seasoning; spice; herbs; pepper powder; chili powder; curry; paprika; cumin; fennel; garlic powder; onion powder; relish; chutney; mustard; ketchup; mayonnaise; mayo; vinegar; gravy; taco seasoning; dip mix",
    "No product category keyword identified"
  ),
  Final_product_cat = c("Dietary supplements", "Meat and poultry", "Seafood", "Dairy", "Fresh produce", "Ready-to-eat / prepared foods", "Ready-to-eat / prepared foods", "Bakery", "Beverages", "Nuts and seeds", "Frozen foods", "Condiments / sauces", "Other foods")
)
write.csv(product_pub, file.path(out_dir, "dictionary_product_category_mapping_publication_ready.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 6. Manifest and validation counts -------------------------------------------

counts <- data.frame(
  step = c("Raw JSON records", "Food subset", "Microbial recall screening", "Main 5 pathogen dataset", "Product categorized dataset", "Final 2012-2024 dataset"),
  n = c(n_raw, n_food, n_micro, nrow(df_main5), nrow(df_main5), nrow(df_sub))
)
write.csv(counts, file.path(out_dir, "validation_counts_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

files <- list.files(out_dir, full.names = FALSE)
manifest <- data.frame(file_name = files, file_size_bytes = file.info(file.path(out_dir, files))$size, stringsAsFactors = FALSE)
write.csv(manifest, file.path(out_dir, "validation_file_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\nValidation export complete.\n")
cat("Output directory:", normalizePath(out_dir, winslash = "/"), "\n\n")
print(counts, row.names = FALSE)
cat("\nNumber of output files:", nrow(manifest), "\n")
