# Bonferroni correction

p = c(0.7676, 0.7940, 0.2535, 0.9053, 0.2988, 0.7333, 0.0088)

p_bonferroni <- p.adjust(p, method = "bonferroni", n = length(p))

p_holm <- p.adjust(p, method = "holm", n = length(p))

