library(tidyverse)
library(ggpubr)
library(ggsci)

ffff <- readRDS("~/weihaowang/project2/simufinal/fdr_sqrt3_depdep.rds")
pppp <- readRDS("~/weihaowang/project2/simufinal/power_sqrt3_depdep.rds")

aaaa <- cbind(ffff[, c(1, 2, 4, 7, 5)], pppp[, c(1, 2, 4, 7, 5)])
aaaa <- rbind(aaaa, matrix(0, nrow = 2, ncol = 10))
mat <- aaaa[1:12, ]

df <- as.data.frame(mat)
colnames(df) <- paste0("V", 1:10)

methods <- c("Oracle", "UCDA", "LinDA", "ANCOM-BC2", "CATE")
n_labels <- c("50", "100", "200")

# ★ 关键改动：levels 直接写成表达式字符串
nu_labels <- c("nu == 0.05", "nu == 0.1", "nu == 0.15", "nu == 0.2")

dat_long <- df %>%
  mutate(
    row_id = row_number(),
    nu_idx = (row_id - 1) %/% 3 + 1,
    nu = factor(nu_labels[nu_idx], levels = nu_labels),
    n_idx = (row_id - 1) %% 3 + 1,
    n = factor(n_labels[n_idx], levels = n_labels)
  ) %>%
  pivot_longer(V1:V10, names_to = "Var", values_to = "Value") %>%
  mutate(
    VarIdx = as.integer(sub("V", "", Var)),
    Metric = ifelse(VarIdx <= 5, "FDR", "Power"),
    MethodIdx = ifelse(VarIdx <= 5, VarIdx, VarIdx - 5),
    Method = factor(methods[MethodIdx], levels = methods),
    Metric = factor(Metric, levels = c("FDR", "Power"))
  )

p12 <- ggplot(dat_long, aes(x = n, y = Value, fill = Method)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(
    data = subset(dat_long, Metric == "FDR"),
    aes(yintercept = 0.05),
    linetype = "dashed",
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  # ★ 关键改动：用 label_parsed
  facet_grid(Metric ~ nu, labeller = label_parsed) +
  scale_fill_nejm() +
  labs(x = "n", y = NULL, fill = "Method") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90")
  )

p12

# ggsave("weihaowang/project2/simupng/sqrt3_dep.png", plot = p12,
#        width = 10, height = 8, bg = "white")