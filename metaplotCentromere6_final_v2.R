###################################################################################################
# R Script: Centromere PCA Profile Analysis (Bootstrap HDI Edition)
###################################################################################################

# =============================
# Step 0: Packages
# =============================
# Adicionado HDInterval aqui para instalação automática
required_packages <- c("data.table", "ggplot2", "cowplot", "ggsci", "dplyr", "tools", "janitor", "HDInterval")
bioc_packages <- c("GenomicRanges", "rtracklayer")

check_install_packages <- function(pkg_list, is_bioc = FALSE) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (is_bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

check_install_packages(required_packages)
check_install_packages(bioc_packages, is_bioc = TRUE)

# =============================
# Step 1: Inputs & parameters
# =============================
setwd("/home/vinicius_souza/Documents/prj_Jcunha/PrjNatalia/paper_natalia/centromere_analyses/centromeMetaplot/metaplotCentromeres/")

centromere_bed <- "./input/centromeros.bed"

bw_files <- c(
  "PCA21" = "./input/25Kb/EPI_R1_PCA2_25K.bigwig",
  "PCA22" = "./input/25Kb/EPI_R2_PCA2_25K.bigwig"
)

flank_size <- 50000
nbins <- 100
window_width <- flank_size * 2

plot_title <- "Average Centromeric PCA Profile (Oriented)"
x_label <- "Distance from Center (bp)"
y_label <- "Relative PCA Strength (Region - Flanks)"

# =============================
# Step 2: Core function
# =============================
process_pca_profile <- function(bw_path, pca_name, centromere_bed, flank_size, nbins) {
  
  message("Processing: ", pca_name, " -> ", bw_path)
  
  if (!file.exists(bw_path)) {
    message("  BigWig not found.")
    return(NULL)
  }
  
  bw <- BigWigFile(bw_path)
  seqlens <- seqlengths(bw)
  
  bed <- fread(centromere_bed, header = FALSE)[, 1:3]
  setnames(bed, c("chr", "start", "end"))
  bed[, `:=`(start = as.integer(start), end = as.integer(end))]
  
  gr_cent <- GRanges(seqnames = bed$chr, ranges = IRanges(start = bed$start + 1L, end = bed$end))
  gr_cent <- resize(gr_cent, width = 1, fix = "center")
  gr_win <- resize(gr_cent, width = window_width, fix = "center")
  
  chr_names <- as.character(seqnames(gr_win))
  chr_lens <- seqlens[chr_names]
  valid_chr <- !is.na(chr_lens) & chr_lens >= window_width
  if (sum(valid_chr) == 0) return(NULL)
  
  gr_win <- gr_win[valid_chr]
  chr_lens <- chr_lens[valid_chr]
  starts <- start(gr_win); ends <- end(gr_win)
  shift_right <- pmax(1L - starts, 0L); starts <- starts + shift_right; ends <- ends + shift_right
  shift_left <- pmax(ends - chr_lens, 0L); starts <- starts - shift_left; ends <- ends - shift_left
  ranges(gr_win) <- IRanges(start = starts, end = ends)
  
  mat_vals <- NULL
  try({
    m <- summary(bw, gr_win, size = nbins, type = "mean", as = "matrix")
    mat_vals <- as.matrix(m)
    if (ncol(mat_vals) != nbins) mat_vals <- t(mat_vals)
  }, silent = TRUE)
  
  if (is.null(mat_vals)) {
    raw <- summary(bw, gr_win, size = nbins, type = "mean")
    mat_vals <- do.call(rbind, lapply(raw, function(x) {
      v <- as.numeric(unlist(x)); length(v) <- nbins; v
    }))
  }
  
  mat_vals <- as.matrix(mat_vals)
  mode(mat_vals) <- "numeric"
  mat_vals <- mat_vals[rowSums(is.na(mat_vals)) < nbins * 0.5, , drop = FALSE]
  
  n_flank <- ceiling(nbins * 0.15)
  flank_cols <- c(seq_len(n_flank), seq(nbins - n_flank + 1, nbins))
  flank_med <- apply(mat_vals[, flank_cols], 1, median, na.rm = TRUE)
  flip <- flank_med < 0
  if (any(flip)) mat_vals[flip, ] <- -mat_vals[flip, ]
  
  baseline <- apply(mat_vals[, flank_cols], 1, median, na.rm = TRUE)
  mat_vals <- sweep(mat_vals, 1, baseline, "-")
  
  # --- Agregação com Bootstrap HDI ---
  message("  Calculating Bootstrap HDI (95%)...")
  n_boot <- 1000
  boot_stats <- apply(mat_vals, 2, function(x) {
    vals <- x[!is.na(x)]
    if(length(vals) < 3) return(c(mean=mean(vals), low=NA, high=NA))
    b_means <- replicate(n_boot, mean(sample(vals, replace=TRUE)))
    hdi_r <- HDInterval::hdi(b_means, credMass=0.95)
    return(c(mean=mean(b_means), low=as.numeric(hdi_r[1]), high=as.numeric(hdi_r[2])))
  })
  boot_stats <- t(boot_stats)
  
  dt <- data.table(
    bin_id = seq_len(nbins),
    mean_signal = boot_stats[, "mean"],
    ci_lower = boot_stats[, "low"],
    ci_upper = boot_stats[, "high"],
    pca = pca_name
  )
  return(dt)
}

# =============================
# Step 3: Run
# =============================
res_list <- lapply(names(bw_files), function(nm) {
  tryCatch(process_pca_profile(bw_files[[nm]], nm, centromere_bed, flank_size, nbins),
           error = function(e) { message("Erro: ", e$message); NULL })
})
df_final <- rbindlist(Filter(Negate(is.null), res_list))
df_final <- janitor::clean_names(df_final)

bin_width <- (2 * flank_size) / nbins
positions <- seq(-flank_size + bin_width/2, flank_size - bin_width/2, length.out = nbins)
df_final[, position := positions[bin_id]]

# =============================
# Step 5: Metaplot (CORRIGIDO)
# =============================
# Verificação atualizada: removemos 'se' e adicionamos 'ci_lower'/'ci_upper'
expected_cols <- c("mean_signal", "ci_lower", "ci_upper", "pca", "position")
if (!all(expected_cols %in% names(df_final))) {
  stop("Colunas faltando. Encontradas: ", paste(names(df_final), collapse=", "))
}

nature_theme <- theme_classic(base_size = 12) + theme(plot.title = element_text(face="bold", hjust=0.5))

p_metaplot <- ggplot(df_final, aes(x = position, y = mean_signal, color = pca, fill = pca)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  scale_x_continuous(breaks = c(-flank_size, 0, flank_size),
                     labels = c(paste0("-", flank_size/1000, "kb"), "Center", paste0("+", flank_size/1000, "kb"))) +
  scale_color_npg() + scale_fill_npg() +
  labs(title = plot_title, subtitle = "Bootstrap 95% HDI Confidence Interval", x = x_label, y = y_label) +
  nature_theme

# =============================
# Step 6: Boxplot de verificação + testes de Wilcoxon MANUAL
# =============================
# Definir regiões (20% left, 20% right, centro 20% central)
idx_L_end <- floor(nbins * 0.2)
idx_C_start <- floor(nbins * 0.4)
idx_C_end <- floor(nbins * 0.6)
idx_R_start <- floor(nbins * 0.8)

df_final[, region_type := data.table::fcase(
  bin_id <= idx_L_end, "Left Flank",
  bin_id >= idx_C_start & bin_id <= idx_C_end, "Center (Centromere)",
  bin_id >= idx_R_start, "Right Flank",
  default = NA_character_
)]

df_stats <- df_final[!is.na(region_type)]
# fator com ordem desejada
df_stats[, region_type := factor(region_type, levels = c("Left Flank", "Center (Centromere)", "Right Flank"))]

# Convert to data.frame for ggplot text handling
df_stats_df <- as.data.frame(df_stats)

# Encontrar nome da coluna pca (depois do clean_names provavelmente 'pca')
pca_col <- if ("pca" %in% names(df_stats_df)) "pca" else grep("^pca", names(df_stats_df), value = TRUE)[1]
if (is.na(pca_col) || length(pca_col) == 0) stop("Coluna 'pca' não encontrada após limpeza de nomes.")

# Calcular Wilcoxon (Left vs Center; Right vs Center) por pca
library(stats)
dt_tests <- data.table()
for (p in unique(df_stats_df[[pca_col]])) {
  sub <- df_stats_df[df_stats_df[[pca_col]] == p, , drop = FALSE]
  left_vals <- sub$mean_signal[sub$region_type == "Left Flank"]
  center_vals <- sub$mean_signal[sub$region_type == "Center (Centromere)"]
  right_vals <- sub$mean_signal[sub$region_type == "Right Flank"]
  
  p_left <- if (length(left_vals) > 0 && length(center_vals) > 0) wilcox.test(left_vals, center_vals)$p.value else NA_real_
  p_right <- if (length(right_vals) > 0 && length(center_vals) > 0) wilcox.test(right_vals, center_vals)$p.value else NA_real_
  
  dt_tests <- rbind(dt_tests, data.table(pca = p, region_type = "Left Flank", comp = "Left_vs_Center", p_value = p_left))
  dt_tests <- rbind(dt_tests, data.table(pca = p, region_type = "Right Flank", comp = "Right_vs_Center", p_value = p_right))
}

# Função p-value -> significância
p_value_to_sig <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}
dt_tests[, p_sig := sapply(p_value, p_value_to_sig)]

# Definir posições y para os labels (um pouco acima do máximo por pca)
y_max_by_pca <- df_stats[, .(ymax = max(mean_signal, na.rm = TRUE), ymin = min(mean_signal, na.rm = TRUE)), by = "pca"]
dt_tests <- merge(dt_tests, y_max_by_pca, by = "pca", all.x = TRUE)
# offset relativo
dt_tests[, y_pos := ymax + (ymax - ymin) * 0.08]

# Plot boxplot com anotações (geom_text)
p_boxplot <- ggplot(df_stats_df, aes(x = region_type, y = mean_signal, fill = .data[[pca_col]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_npg() +
  labs(title = "Statistical Check", y = "Relative Signal", x = NULL) +
  nature_theme +
  geom_text(
    data = dt_tests,
    aes(x = region_type, y = y_pos, label = p_sig, group = pca),
    position = position_dodge(width = 0.8),
    size = 3,
    vjust = 0
  )

# =============================
# Step 7: Salvar
# =============================
final_panel <- plot_grid(p_metaplot, p_boxplot, ncol = 1, rel_heights = c(1.5, 1))
ggsave("./output/Centromere_PCA_HDI_Bootstrap.pdf", final_panel, width = 18, height = 20, units = "cm")
fwrite(df_final, "./output/Centromere_PCA_HDI_Table.tsv", sep = "\t")
message("Concluído com sucesso!")