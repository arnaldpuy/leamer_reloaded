#rm(list = ls())

## ----setup, include=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "tikz", cache = TRUE)


## ----preliminary, warning=FALSE, message=FALSE-----------------------------------------

# LOAD PACKAGES ################################################################

library(sensobol)
library(dplyr)
library(ggplot2)
library(grid)
hgutils::load_packages(c("haven", "tidyverse", "data.table", "margins", "miceadds",
                         "doParallel", "foreach", "cowplot", "benchmarkme", "sandwich", 
                         "countrycode"))

# THEME FOR PLOTTING ###########################################################

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA), 
          strip.background = element_rect(fill = "white"), 
          legend.margin = margin(0.5, 0.1, 0.1, 0.1),
          legend.box.margin = margin(0.2,-4,-7,-7), 
          plot.margin = margin(3, 4, 0, 4), 
          legend.text = element_text(size = 8), 
          axis.title = element_text(size = 10),
          legend.key.width = unit(0.4, "cm"), 
          legend.key.height = unit(0.4, "cm"), 
          legend.title = element_text(size = 9)) 
}


## ----data------------------------------------------------------------------------------

# READ DATA AND STANDARDIZE ####################################################

# these data are taken from the study of Puy et al
#https://github.com/arnaldpuy/universe_of_uncertainty
data <- readRDS("df.rds")
setDT(data)

# create Western Europe and Rich Western Democracies variables

europe <- c(40, 56, 756, 276, 250, 826, 442, 528, 438, 724, 620, 380, 470, 208, 246, 352, 578, 752, 492, 20, 674, 336)

west_rich_dem <- c(40, 56, 756, 276, 250, 826, 442, 528, 438, 724, 620, 380, 470, 208, 246, 352, 578, 752, 492, 20, 674, 336, 840, 124, 392, 372, 36, 554)


data <- data %>%
  mutate(cont_europe = ifelse(iso_country %in% europe, 1, 0),
         rich_dems = ifelse(iso_country %in% west_rich_dem, 1, 0))

# Columns with migration data  -------------------------------------------------

cols <- c("foreignpct", "migstock_wb", "migstock_un", "migstock_oecd")

# Standardize migration data ---------------------------------------------------

data[, (cols) := lapply(.SD, function(x) as.numeric(scale(as.numeric(x)))), .SDcols = cols]

# TURN DOTS AND EMPTY CELLS INTO NA ############################################

data[data == "."] <- NA
data[data == ""] <- NA

# IDENTIFY NUMERIC COLUMNS AND TRANSFORM #######################################

numeric_columns <- sapply(data, is.numeric)
numeric_columns <- names(numeric_columns[numeric_columns])
data[, (numeric_columns):= lapply(.SD, as.numeric), .SDcols = (numeric_columns)]
data[, iso_country:= as.numeric(iso_country)]

to.numeric <- c("gdp_wb", "gdp_twn", "gni_wb", "gini_wid", "ginim_dolt", "top10_wid", 
                "pop_wb", "al_ethnic", "dpi_tf", "wdi_empprilo", "mcp", "mignet_un")

data[, (to.numeric):= lapply(.SD, as.numeric), .SDcols = (to.numeric)]

# IDENTIFY THE PROPORTION OF MISSING VALUES PER COLUMN AND YEAR ################

data[, lapply(.SD, function(x) sum(is.na(x)) / nrow(.SD)), year]

# COUNTRY VARIABLE #############################################################

country <- unique(data$country)
setDT(data)
data[, country_cl:= countrycode(country, "country.name", "continent")]
data[, country_cl:= ifelse(country_cl == "Oceania" | country_cl == "Asia", "Asia_Oceania", 
                           ifelse(country_cl == "Americas", "America", country_cl))]


# DEFINE VARIABLES #############################################################

# Dependent variable (logit and linear) ----------------------------------------

Y.vector <- c("jobs", "jobs_c") 

# Policy variables -------------------------------------------------------------

X <- c("foreignpct", "migstock_un", "netmigpct", "migstock_oecd")         

# Controls variables -----------------------------------------------------------

controls.ind = c("female", "education", "income", 
                 "pop_wb", "al_ethnic","socx_oecd", "country")     

# Controls country -------------------------------------------------------------

controls.country <- c("America - Asia_Oceania", "Asia_Oceania - Europe", 
                      "America - Europe", "America - Asia_Oceania - Europe", "Europe")

# Controls age  ----------------------------------------------------------------

controls.age <- c("age", "age_sq")  

# Controls gini  ---------------------------------------------------------------

controls.gini <- c("ginid_solt", "ginim_dolt")           

# Robust standard errors -------------------------------------------------------

se.vec <- c("no", "country")

# Yes / no vector --------------------------------------------------------------
yes.no.vec <- c("yes", "no")   


## ----sample_matrix, dependson="data"---------------------------------------------------

# DEFINE SAMPLE MATRIX #########################################################

# Settings ---------------------------------------------------------------------

matrices <- c("A", "B", "AB")
N <- 2^10
params <- c("Y", "X1", controls.ind, "age", "SE", "bootstrap", "GINI", "cont_europe", "rich_dems")

# Sample matrix ----------------------------------------------------------------

mat <- data.table(sobol_matrices(matrices = matrices, params = params, N = N))

# To remove -------------------------------------------------------------------

to.remove <- "country"
specific.columns <- setdiff(controls.ind, to.remove)

# Transform to appropriate columns ---------------------------------------------

mat[, Y:= floor(Y * 2) + 1]
mat[, Y:= Y.vector[Y]]
mat[, X1:= floor(X1 * length(X)) + 1]
mat[, X1:= X[X1]]
mat[, (specific.columns):= lapply(.SD, function(x) ifelse(x > 0.5, "yes", "no")), .SDcols = (specific.columns)]
mat[, country:= floor(country * length(controls.country)) + 1]
mat[, country:= controls.country[country]]
mat[, age:= floor(age * length(yes.no.vec)) + 1]
mat[, age:= yes.no.vec[age]]
mat[, SE:= floor(SE * length(se.vec)) + 1]
mat[, SE:= se.vec[SE]]
mat[, bootstrap:= floor(bootstrap * length(yes.no.vec)) + 1]
mat[, bootstrap:= yes.no.vec[bootstrap]]
mat[, GINI:= floor(GINI * length(controls.gini)) + 1]
mat[, GINI:= controls.gini[GINI]]
mat[, cont_europe:= floor(cont_europe * length(yes.no.vec)) + 1]
mat[, cont_europe:= yes.no.vec[cont_europe]]
mat[, rich_dems:= floor(rich_dems * length(yes.no.vec)) + 1]
mat[, rich_dems:= yes.no.vec[rich_dems]]

# DEFINE THE A, B AND AB MATRICES FOR THE BOOTSTRAP SEED #######################

set.seed(40)
bootstrap.seed <- sample(1:(N*3), N * 2)
A.bootstrap.seed <- bootstrap.seed[1:N]
B.bootstrap.seed <- bootstrap.seed[(N + 1):length(bootstrap.seed)]
n.reps <- which(params == "bootstrap") -1
bootstrap.seed.mat <- c(A.bootstrap.seed, B.bootstrap.seed, 
                        rep(A.bootstrap.seed, n.reps), 
                        B.bootstrap.seed, 
                        rep(A.bootstrap.seed, length(params) - which(params == "bootstrap")))

################################################################################

# DEFINE THE A, B AND AB MATRICES FOR THE AGE_SQ VARIABLE ######################

set.seed(40)
sample.age.sq <- runif(N*2)
sample.age.sq.trans <- floor(sample.age.sq * length(yes.no.vec)) + 1
sample.age.sq.trans <- yes.no.vec[sample.age.sq.trans]
A.sample.age.sq <- sample.age.sq.trans[1:N]
B.sample.age.sq <- sample.age.sq.trans[(N + 1):length(sample.age.sq.trans)]
n.reps <- which(params == "age") -1
age_sq <- c(A.sample.age.sq, 
             B.sample.age.sq, 
             rep(A.sample.age.sq, n.reps), 
             B.sample.age.sq, 
             rep(A.sample.age.sq, length(params) - which(params == "age")))

# CREATE FINAL SAMPLING MATRIX #################################################

mat <- cbind(mat, bootstrap.seed.mat, age_sq)

# DEFINE MODEL #################################################################

# Function to create the glm formula -------------------------------------------

formula_fun <- function(mat) {
  
  formula = paste0(mat[, Y], " ~ ", mat[,X1], collapse = " ")
  
  # Add GINI -------------------------------------------------------------------
  
  formula <- paste0(formula, " + ", paste0(mat[, GINI], collapse = " + "))
  
  # List of additional variables to include if they are marked as "yes" --------
  
  additional_vars <- c("female", "education", "income",
                       "pop_wb", "al_ethnic", "socx_oecd", "country", "age", "age_sq")
  
  # Search variables marked as "yes" in mat ------------------------------------
 
  included_vars <-  names(mat[, ..additional_vars])[which(mat[, ..additional_vars] == "yes")]
  
  # Append included variables to the formula -----------------------------------
  
  if (length(included_vars) > 0) {
    
    formula <- paste0(formula, " + ", paste(included_vars, collapse = " + "))
  }

  return(formula)

}

# Wrap up function -------------------------------------------------------------

full_fun <- function(data, mat, Y.vector) {
  
  formula <- as.formula(formula_fun(mat = mat))
  

  
  # Country filtering ----------------------------------------------------------
  
  countries_to_filter <- switch(as.character(mat[, country]),
                                "America - Asia_Oceania" = c("America", "Asia_Oceania"),
                                "Asia_Oceania - Europe" = c("Asia_Oceania", "Europe"),
                                "America - Europe" = c("America", "Europe"),
                                "America - Asia_Oceania - Europe" = c("America", "Asia_Oceania", "Europe"),
                                "Europe" = "Europe",
                                NULL)
  
  if (!is.null(countries_to_filter)) {
    
    data <- data[country_cl %in% countries_to_filter]
  }
  
  # cont_europe and rich_dems ----------------------------------------------
  
  if (mat[, cont_europe] == "yes") {
    
    if (length(table(data$cont_europe)) > 1) {
    
      data = subset(data, cont_europe == 1)
    
    }

  }
  
  if (mat[, rich_dems] == "yes") {
    
    if (length(table(data$rich_dems)) > 1) {
    
      data = subset(data, rich_dems == 1)
      
    }
    
  }
  
  # Bootstrap ------------------------------------------------------------------
  
   if (mat[, bootstrap] == "yes") {
    
     set.seed(mat[, bootstrap.seed.mat])
     data <- data[sample(1:nrow(data), nrow(data), replace = TRUE), ]
    
   }
  
  # Determine model family and cluster variable --------------------------------
  
  family <- if (mat[, Y] %in% Y.vector[1]) "binomial" else "gaussian"
  
  cluster_var <- switch(as.character(mat[, SE]),
                        "no" = NULL,
                        "country" = "country",
                        "year")
  
  # Use only complete cases ----------------------------------------------------
  
  cols_formula <- all.vars(formula)
  
  if(!is.null(cluster_var)) {
    
    colnames_selected <- c(cols_formula, cluster_var)
    
  } else {
    
    colnames_selected <- cols_formula
  }
  
  data <- data[, ..colnames_selected][complete.cases(data[, ..colnames_selected])]
  
  # Execute model --------------------------------------------------------------
  
  if (family == "binomial") {
    
    if (is.null(cluster_var)) {
      
      model <- glm(formula, family = binomial, data = data)
      out <- margins(model, variables = mat[, X1])
      
    } else {
      
      model <- glm.cluster(data = data[, ..cols_formula], formula = formula, 
                           cluster = data[[cluster_var]], family = "binomial")
      out <- margins(model$glm_res, variables = mat[, X1], vcov = model$vcov)
      
    }
    
  } else {
    
    if (is.null(cluster_var)) {
      
      model <- lm(formula = formula, data = data)
      out <- margins(model, variables = mat[, X1])
      
    } else {
      
      model <- lm.cluster(data = data[, ..cols_formula], formula = formula, 
                          cluster = data[[cluster_var]])
      out <- margins(model$lm_res, variables = mat[, X1], vcov = model$vcov)
      
    }
  }
  
  out <- as.numeric(summary(out)[c("AME", "lower", "upper")])
  
  return(out)
}


## ----parallel, dependson="sample_matrix"-----------------------------------------------

# PARALLEL COMPUTING ###########################################################

# Define parallel computing ----------------------------------------------------

cl <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(cl)

# Run the simulations ----------------------------------------------------------

out <-  foreach(i = 1:nrow(mat), 
                .packages = c("margins", "miceadds", "data.table"), 
                .combine = "rbind") %dopar% {
                  
                  full_fun(data = data, mat = mat[i, ], Y.vector = Y.vector)
                  
                }

# Stop the parallel backend ----------------------------------------------------

stopCluster(cl)


## ----arrange_data, dependson="parallel"------------------------------------------------

# ARRANGE DATA #################################################################

out.dt <- data.table(out) %>%
  setnames(., colnames(.), c("AME", "low", "high")) 

AME.dt <- out.dt[1:(2 * N), ]

# EXPORT DATA ##################################################################

full.dt <- cbind(mat, out.dt)
fwrite(full.dt, "full.dt.csv")


## ----plot_uncertainty, dependson="arrange_data", dev = "pdf", fig.height=2, fig.width=3.7----

# PLOT UNCERTAINTY #############################################################

tmp <- out.dt[order(AME.dt)] %>%
  .[, Model:= .I] %>%
  .[, color:= ifelse(low < 0 & high < 0, "negative", 
                     ifelse(low < 0 & high > 0, "includes zero", "positive"))] %>%
  .[, color:= factor(color, levels = c("negative", "includes zero", "positive"))]

tmp[, .N, color] %>%
  .[, total.n:= nrow(tmp)] %>%
  .[, fraction:= N / total.n] %>%
  print()

# Calculate 2.5% and 97.5% quantiles -------------------------------------------

quantiles <- quantile(AME.dt$AME, c(0.025, 0.975))
quantiles

var_AME = var(AME.dt$AME)
tmp <- tmp %>%
  na.omit() %>%
  arrange(AME) %>%
  mutate(Model = row_number())
perc <- tmp %>%
  count(color) %>%
  mutate(perc = 100 * n / sum(n))
cols <- c("orange", "grey", "blue")
plot.uncertainty <- ggplot(tmp, aes(x = Model, y = AME, color = color)) +
  geom_errorbar(aes(ymin = low, ymax = high), size = 0.3) +
  geom_point(size = 0.4) +
  geom_hline(yintercept = 0, lty = 2, color = "black") +
  scale_color_manual(values = cols,
                     labels = c("NEGATIVE (stat. sig.)", 
                                "INCLUDES ZERO (n.s.)", 
                                "POSITIVE (stat. sig.)")) +
  scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.5, 0.5, 0.1)) +
  labs(x = "Models Ordered by AME",
       y = "Average Marginal Effect (AME)",
       color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())
posizioni <- data.frame(
  x = c(nrow(tmp)*0.2, nrow(tmp)*0.5, nrow(tmp)*0.8),
  y = rep(0.55, 3),
  color = c("orange", "grey", "blue"),
  perc = round(perc$perc, 1)
)
plot.uncertainty = plot.uncertainty +
  geom_point(data = posizioni,
             aes(x = x, y = y),
             shape = 21, size = 30, stroke = 1.2,
             color = posizioni$color, fill = "white", inherit.aes = FALSE) +
  geom_text(data = posizioni,
            aes(x = x, y = y, label = perc),
            color = posizioni$color, size = 6, fontface = "bold", inherit.aes = FALSE) +
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(20, 20, 40, 20)) #+
  #annotation_custom(
  #  grid::textGrob(
  #    label = paste0("V[AME] = ", round(var_AME, 3)),
  #    x = unit(0.98, "npc"),  # 98% a destra
  #    y = unit(0.02, "npc"),  # 2% dal basso
  #    just = c("right", "bottom"),
  #    gp = gpar(col = "black", fontsize = 10)
  #  )
  #)
plot.uncertainty

## ----plot_scatter, dependson="arrange_data", dev = "pdf"-------------------------------

# PLOT SCATTERPLOT #############################################################

plot_scatter(data = mat, N = N, Y = out, params = params) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(size = 6)) +
  labs(x = "", y = "AMS") + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2))


## ----sa, dependson="parallel", dev = "pdf", fig.height=2.5, fig.width=4.5--------------

# SENSITIVITY ANALYSIS #########################################################

# Define setting ---------------------------------------------------------------

first <- "saltelli"
total <- "jansen"
R <- 10^3

# Run SA -----------------------------------------------------------------------

round(var(AME.dt$AME),3)
params[params == "country"] = "continent"
ind <- sobol_indices(matrices = matrices, params = params, N = N, Y = out.dt$AME, 
                     first = first, total = total)

round(sum(ind$results[16:30,1]),3)

ind$results$original[ind$results$original<0] = 0
ind$results$original[ind$results$original>1] = 1

# Plot Sobol' indices ----------------------------------------------------------

plot(ind)

MD = sum(ind$results[ind$results$sensitivity=="Ti",1])
plot.sobol <- plot(ind) + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm")
  ) +
  labs(
    x = "Input Variables",
    y = "Sobol Total Index",
    color = "Sobol indices"
  ) +
  annotate("text",
           x = 0.5, y = 1,            
           label = paste0("MD = ", round(MD, 3)),
           hjust = 0, vjust = 1.2,    
           size = 4, color = "black")
plot.sobol

## ----merge_plots, dependson=c("sa", "plot_uncertainty"), fig.width=4.5, fig.height=4, dev = "pdf"----

# ADD ORIGINAL BREZNAU ET AL PLOT ##############################################

original.plot<- cowplot::ggdraw() + cowplot::draw_image("fig01.png", scale = 0.9)

original.plot

p2 <- plot_grid(plot.uncertainty, plot.sobol, ncol = 1, labels = c("b", "c"))

p2



## ----all_together, dependson="merge_plots", fig.width=4.5, fig.height=5.5, dev = "pdf"----

# MERGE ALL PLOTS ###############################################################

pdf(file = "Jobs.pdf", width = 14, height = 7)  # più largo del default
plot_grid(plot.uncertainty, plot.sobol, ncol = 2, labels = c("a", "b"))
dev.off()

## ----session_info----------------------------------------------------------------------

# SESSION INFORMATION ##########################################################

sessionInfo()

## Return the machine CPU
cat("Machine:     "); print(get_cpu()$model_name)

## Return number of true cores
cat("Num cores:   "); print(detectCores(logical = FALSE))

## Return number of threads
cat("Num threads: "); print(detectCores(logical = FALSE))

