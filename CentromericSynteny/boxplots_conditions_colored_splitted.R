#!/usr/bin/env Rscript
# plot_compare_genomes_colored.R

# ------------------------------------------------------------
# 1. Setup e Dependências
# ------------------------------------------------------------
if (!suppressWarnings(require("ggplot2", quietly = TRUE))) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript boxplots_conditions_colored.R <controle.tsv> <tratamento.tsv>")
}

file_ctrl <- args[1]
file_treat <- args[2]

message("--- Processando ---")
message("Controle:   ", file_ctrl)
message("Tratamento: ", file_treat)

# ------------------------------------------------------------
# 2. Funções Auxiliares
# ------------------------------------------------------------
read_input <- function(path) {
  try_csv <- try(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (!inherits(try_csv, "try-error") && ncol(try_csv) > 1) return(try_csv)
  try_tsv <- try(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  if (!inherits(try_tsv, "try-error")) return(try_tsv)
  stop("Erro ao ler arquivo: ", path)
}

clean_genome_name <- function(x) {
  if (is.na(x)) return(NA)
  x <- basename(as.character(x))
  x <- sub("\\.gz$", "", x)
  x <- sub("\\.(fasta|fa|fna|txt|csv|tsv)$", "", x, ignore.case = TRUE)
  # Padronização extra: remove underscores e espaços para facilitar o match
  return(x)
}

# ------------------------------------------------------------
# 3. Carregar e Unir Dados
# ------------------------------------------------------------
df_c <- read_input(file_ctrl)
df_t <- read_input(file_treat)

df_c$Condition <- "Random Regions"
df_t$Condition <- "Centromeric Regions"

df_all <- rbind(df_c, df_t)
cols_num <- c("alen", "match", "qcov")
df_all[cols_num] <- lapply(df_all[cols_num], function(x) as.numeric(as.character(x)))

df_all$identity_pct <- (df_all$match / df_all$alen) * 100
df_all$genome_clean <- sapply(df_all$other_genome, clean_genome_name)

plot_df <- df_all[!is.na(df_all$identity_pct) & !is.na(df_all$qcov) & df_all$alen > 0, ]

# ------------------------------------------------------------
# 4. Ordenação Customizada (Ajuste de nomes aqui)
# ------------------------------------------------------------
# DICA: Verifique se os nomes abaixo batem com os nomes dos seus arquivos.
# Tentei prever variações comuns (como underscores em vez de espaços)
all_found_genomes <- unique(plot_df$genome_clean)
message("Genomas encontrados no arquivo: ", paste(all_found_genomes, collapse=", "))

# Esta lista deve conter os nomes EXATOS que apareceram na mensagem acima
custom_order <- c(
  "Dm25", "BrazilA4", "SilvioX10", # Grupo 1
  "Berenice", "YC6",                            # Grupo 2
  "ClBrenerEsmeraldoLike", "TCC"               # Grupo 3
)

# Filtra apenas os que existem no dado para não criar colunas vazias
existing_order <- custom_order[custom_order %in% all_found_genomes]

if(length(existing_order) == 0) {
  message("AVISO: Nenhum nome da lista custom_order foi encontrado. Usando ordem alfabética.")
  plot_df$genome_clean <- factor(plot_df$genome_clean)
} else {
  plot_df <- plot_df[plot_df$genome_clean %in% existing_order, ]
  plot_df$genome_clean <- factor(plot_df$genome_clean, levels = existing_order)
}

plot_df$Condition <- factor(plot_df$Condition, levels = c("Random Regions", "Centromeric Regions"))

# ------------------------------------------------------------
# 5. Estatística
# ------------------------------------------------------------
#calc_stats <- function(data, variable) {
#  genomes <- levels(data$genome_clean)
#  stats_res <- data.frame()
#
#  for (g in genomes) {
#    sub_d <- data[data$genome_clean == g, ]
#    if (nrow(sub_d) > 0 && length(unique(sub_d$Condition)) == 2) {
#      test <- wilcox.test(as.formula(paste(variable, "~ Condition")), data = sub_d)
#      p_val <- test$p.value
#      sig <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else "ns"
#
#      y_pos <- max(sub_d[[variable]], na.rm = TRUE) * 1.05
#      stats_res <- rbind(stats_res, data.frame(
#        genome_clean = g, label = sig, y_pos = y_pos, Condition = "Control"
#      ))
#    }
#  }
#  return(stats_res)
#}
#
#stats_ident <- calc_stats(plot_df, "identity_pct")
#stats_cov   <- calc_stats(plot_df, "qcov")

# ------------------------------------------------------------
# 6. Cores e Plotagem
# ------------------------------------------------------------
n_genomes <- length(unique(plot_df$genome_clean))
my_colors <- grDevices::rainbow(n_genomes, s = 0.6, v = 0.9)

plot_box <- function(data, stats, y_var, title, y_label) {
  ggplot(data, aes(x = genome_clean, y = .data[[y_var]])) +
    facet_grid(Condition ~ ., scales = "fixed") +
    geom_boxplot(aes(fill = genome_clean), outlier.shape = NA, color = "black", alpha = 0.8) +
    geom_point(position = position_jitter(width = 0.15), size = 0.5, alpha = 0.2) +
    #geom_text(data = stats, aes(x = genome_clean, y = y_pos, label = label),
    #          inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 0) +
    scale_fill_manual(values = my_colors) +
    theme_bw(base_size = 14) +
    labs(title = title, x = "Genome", y = y_label) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()  # CORRIGIDO AQUI
    )
}

p_id <- plot_box(plot_df, stats_ident, "identity_pct", "Identity (%): Random Regions vs Centromeric Regions", "Identity (%)")
p_cov <- plot_box(plot_df, stats_cov, "qcov", "Coverage (%): Random Regions vs Centromeric Regions", "Coverage (%)")

ggsave("genome_grouped_identity.pdf", p_id, width = 11, height = 9, dpi = 600)
ggsave("genome_grouped_coverage.pdf", p_cov, width = 11, height = 9, dpi = 600)

message("Concluído.")
