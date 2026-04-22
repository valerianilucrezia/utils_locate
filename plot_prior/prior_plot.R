# Load required libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# Function to calculate beta distribution density
dbeta_mix <- function(x, weights, alphas, betas) {
  result <- 0
  for (i in 1:length(weights)) {
    result <- result + weights[i] * dbeta(x, alphas[i], betas[i])
  }
  return(result)
}

# Create data frame for plotting
x <- seq(0, 1, length.out = 1000)

# Cell line parameters (single beta)
cl_alpha <- 50
cl_beta <- 2

# Clinical sample parameters (mixture of 3 betas)
cs_weights <- c(0.4, 0.4, 0.2)  # High, medium, low purity
cs_alphas <- c(15, 8, 2)
cs_betas <- c(3, 8, 8)

# Generate density values
cell_line_density <- dbeta(x, cl_alpha, cl_beta)
clinical_density <- sapply(x, dbeta_mix, weights=cs_weights, alphas=cs_alphas, betas=cs_betas)

# Also calculate individual components of the mixture for visualization
clinical_comp1 <- cs_weights[1] * dbeta(x, cs_alphas[1], cs_betas[1])
clinical_comp2 <- cs_weights[2] * dbeta(x, cs_alphas[2], cs_betas[2])
clinical_comp3 <- cs_weights[3] * dbeta(x, cs_alphas[3], cs_betas[3])

# Create combined data frame
df <- data.frame(
  x = rep(x, 5),
  density = c(cell_line_density, clinical_density, clinical_comp1, clinical_comp2, clinical_comp3),
  Distribution = c(
    rep("Cell Line", length(x)),
    rep("Clinical", length(x)),
    rep("Clinical - High Purity Component", length(x)),
    rep("Clinical - Medium Purity Component", length(x)),
    rep("Clinical - Low Purity Component", length(x))
  ),
  Type = c(
    rep("Single Beta", length(x)),
    rep("Mixture", length(x)),
    rep("Component", length(x)),
    rep("Component", length(x)),
    rep("Component", length(x))
  )
)

# Create plot
ggplot(df, aes(x = x, y = density, color = Distribution, linetype = Type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c(
    "Cell Line" = "#E41A1C",
    "Clinical (Combined)" = "#377EB8",
    "Clinical - High Purity Component" = "#4DAF4A",
    "Clinical - Medium Purity Component" = "#984EA3",
    "Clinical - Low Purity Component" = "#FF7F00"
  )) +
  scale_linetype_manual(values = c("Single Beta" = "solid", "Mixture" = "solid", "Component" = "dashed")) +
  labs(
    title = "Tumor Purity Prior Distributions",
    subtitle = "Cell Line (Single Beta) vs Clinical (Mixture of Betas)",
    x = "Purity",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Create a separate plot showing just the two main distributions for clarity
main_df <- df %>% filter(Distribution %in% c("Cell Line", "Clinical"))

purity <- ggplot(main_df, aes(x = x, y = density, color = Distribution)) +
  geom_line(size = 1) +
  scale_color_manual('', values = c("Cell Line" = "sienna2", "Clinical" = "plum4")) +
  labs(
    x = "Purity",
    y = "Density"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )


