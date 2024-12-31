# Forest plot MR

install.packages("ggplot2")

library(ggplot2)

setwd("/Users/guillermocomesanacimadevila/Desktop/Mendelian Randomisation (PA & ID)")

data_iron <- data.frame(
  SNP = c("rs1799945", "rs1800562", "rs57659670", "rs855971"),
  Beta = c(0.17, 0.27, -0.042, -0.17),
  Lower_CI = c(0.159, 0.0076, -0.0042, -0.17),
  Upper_CI = c(0.181, 0.27, 0.0074, 0.0046),
  P_value = c("1.26e-187", "1e-200", "1.08e-8", "1e-200")
)

ggplot(data_iron, aes(x = Beta, y = SNP)) +
  geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI), color = "blue", size = 0.8) +
  geom_point(size = 4, shape = 21, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16)
  ) +
  labs(
    title = "Enhanced Forest Plot",
    x = "Effect Size (Beta)",
    caption = "Beta values with 95% Confidence Intervals"
  ) +
  geom_text(aes(label = P_value), nudge_x = 0.03, hjust = 0, size = 4, color = "darkgreen")