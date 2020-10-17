# 2020-08-05
# this is to do cluster analysis
# this will generate the data for later analysis

# first set current path
library(rstudioapi)

# set the path
current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))


# load library ------------------------------------------------------------


library(readxl)
library(dplyr)
library(ggplot2)
library(janitor)
library(fpc) # for cluster.stats
library(scales) # for pretty breaks
library(FactoMineR) # for PCA

# remove all the objects
rm(list = ls())

# define functions ------------------------------------------------

count_plot_func <- function(df_data, out_path, width, height){
  
  # this function does the count
  gp <- ggplot(df_data) +
    geom_col(aes(x=reorder(feature, -count), y=count)) +
    # coord_flip() +
    labs(x="",
         y="Count") +
    scale_y_continuous(breaks = seq(0, 100000, by = 5000)) +
    theme_bw() +
    theme(aspect.ratio = 0.4,
          panel.grid.major = element_line(size = 0.1),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size=7, angle = 90, hjust = 1),
          axis.text.y = element_text(size=7),
          axis.title = element_text(size=7),
          axis.ticks = element_line(colour = '#333333', size = 0.1),
          axis.ticks.length = unit(0.1, 'cm'))
  
  print(gp)
  
  ggsave(filename = out_path, width = width, height = height)
  
}


barplot_single_func <- function(df, x, xlab, ylab, color_par, color_sequence, name_legend, show_legend, out_path, width, height){
  
  # this function counts the garnet types
  gp <- ggplot(df, aes(x = !!sym(x), fill = !!sym(color_par))) +
    geom_bar(color = "white", width=0.5, show.legend = show_legend) +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size=2.5) +
    labs(x = xlab, y = ylab) +
    scale_fill_manual(values=color_sequence) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_bw() +
    theme(
      panel.grid.major = element_line(size = 0.1),
      panel.grid.minor = element_line(size = 0.1),
      axis.title = element_text(size=8),
      axis.text = element_text(size=7),
      axis.text.x = element_text(angle = 0)
    )
  print(gp)
  
  ggsave(filename = out_path, width = width, height = height)
  
}

wss_plot_func <- function(df, out_path, width, height){
  
  # plot the elbow from cluster analysis
  gp <- ggplot(df) +
    geom_line(aes(x=cluster, y=wss), size=0.5) +
    geom_point(aes(x=cluster, y=wss), size=1) +
    geom_vline(xintercept = 6, lty="dashed") +
    labs(x="Number of clusters", y="Total within cluster variation") +
    scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=8),
          axis.title = element_text(size=10)
    )
  
  print(gp)
  
  ggsave(filename = out_path, width = width, height = height)
  
}


colplot_spread_percent_func <- function(df, x, y, xlab, ylab, color_par, color_sequence, name_legend, out_path, width, height){
  
  # this functions plots the distribution of cluster analysis
  if (length(unique(df[[x]])) <= 3){
    size = 1
  }else{
    size = 1
  }
  gp <- ggplot(df) +
    geom_col(aes(x = !!sym(x), y = !!sym(y), fill = !!sym(color_par)), position = "dodge", width = 0.7) +
    geom_text(aes(x = !!sym(x), y = !!sym(y), fill = !!sym(color_par), label = !!sym(y)), position = position_dodge(width = 0.7), vjust = -0.5, size=size) +
    labs(x = xlab,
         y = ylab) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(name = name_legend, values = color_sequence) +
    guides(fill=guide_legend(keywidth = unit(5, "mm"),
                             keyheight = unit(5, "mm"))
    ) + 
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid.major = element_line(size = 0.2),
          panel.grid.minor = element_line(size = 0),
          axis.text = element_text(size=8),
          axis.title = element_text(size=10),
          legend.title = element_text(color = "#333333", size = 9),
          legend.text = element_text(color = "#333333", size=8)
    )
  
  
  
  print(gp)
  
  ggsave(filename = out_path, width = width, height = height)
  
}

# define color ------------------------------------------------------------

# based on rcolorbrewer from low to high
color_sequence = c("#2166ac", "#4393c3", "#abd9e9", "#fee090", "#fdae61", "#d73027")

color_sequence_func <- function(color_num){
  if (color_num == 2){
    return(color_sequence[c(1,6)]) # blue and red for two 
  }else if (color_num == 3){
    return(color_sequence[c(1,4,6)]) # blue, yellow and red for three
  }else{
    return(colorRampPalette(color_sequence)(color_num))
  }
}


# Load data ---------------------------------------------------------------

data_path <- file.path("Final Garnet Dataset.xlsx")

# change guess_max to make sure the blanks don't influence the data type
df <- read_excel(data_path, sheet = "Main Database", guess_max = 9e4)

# data exploration ---------------------------------------------------------------

# clean names
df <- df %>% 
  clean_names()

# select confidence interval of SiO2 that is A or B
df_filter <- df %>% 
  dplyr::filter(confidence_interval_of_si_o2_wt_percent %in% c("A", "B"))

dim(df_filter)

# drop the empty columns that are empty or have less than 1000 non-na values

count_nan_na <- sapply(df_filter, function(x) sum(!is.na(x)))

df_select <- df_filter[, count_nan_na >= 1000]

dim(df_select)


# select major garnet types for preliminary cluster analysis --------

df_select_mineral <- df_select %>% 
  dplyr::filter(mineral %in% c("Uvarovite", "Andradite", "Grossular",
                               "Spessartine", "Pyrope", "Almandine"))

# use 6 types
num_type <- 6

# change from tibble to data.frame format
df_select_mineral <- as.data.frame(df_select_mineral)

# delete repeated samples
table(df_select_mineral['repeat'])

df_select_mineral <- df_select_mineral[df_select_mineral["repeat"] == 0, ]


# count again
# get the non-NA count (dplyr-style)
count_mineral_nan_na <- sapply(df_select_mineral, function(x) sum(!is.na(x)))


df_mineral_count <- data.frame(feature=names(df_select_mineral),
                               count = count_mineral_nan_na)


count_plot_func(df_data = df_mineral_count,
                out_path = file.path("output", "garnet_mineral_col_nanna_count.pdf"),
                width = 12,
                height = 6)

# rank the names based on their count
df_mineral_count <- df_mineral_count %>% 
  arrange(desc(count))

df_mineral_count


# select the parameters for cluster analysis ------------------------------


df_select_mineral_short <- df_select_mineral %>% 
  dplyr::select(si_o2_wt_percent,
                al2o3_wt_percent,
                ca_o_wt_percent,
                mg_o_wt_percent,
                mn_o_wt_percent,
                fe_o_wt_percent,
                cr2o3_wt_percent,
                fe2o3_wt_percent,
  )

# store the number of parameters we want to use for cluster analysis
col_end <- ncol(df_select_mineral_short)

# add the garnet type
df_select_mineral_short$Type <- df_select_mineral$mineral


df_select_mineral_short$Type <- factor(df_select_mineral_short$Type, levels = c("Almandine", "Spessartine", "Pyrope", "Grossular", "Andradite", "Uvarovite"))

# convert the values to numerical format
# df_select_mineral_short_numeric <- sapply(df_select_mineral_short[1:col_end], as.numeric)

df_select_mineral_short_numeric <- data.frame(lapply(df_select_mineral_short[1:col_end], as.numeric))

# bind the remaining columns
df_select_mineral_short_transform <- cbind(df_select_mineral_short_numeric, df_select_mineral_short[(col_end+1):ncol(df_select_mineral_short)])

names(df_select_mineral_short_transform)

# remove the rows that contain any NA values
df_select_mineral_short_transform <- na.omit(df_select_mineral_short_transform)

barplot_single_func(df = df_select_mineral_short_transform,
                    x = "Type",
                    xlab = "Type",
                    ylab = "Count",
                    color_par = "Type",
                    color_sequence = color_sequence_func(num_type),
                    name_legend = "",
                    show_legend = FALSE,
                    out_path = file.path("output", "Mineral_counts.pdf"),
                    width = 4,
                    height = 4)




# assign a short name to our variable
df_final <- df_select_mineral_short_transform

table(df_final$Type)

nrow(df_final)

# transform ---------------------------------------------------------------

# scale

mat_clean_final <- scale(df_final[1:col_end], center = TRUE, scale = TRUE)

df_clean_final <- cbind(mat_clean_final, df_final[(col_end+1):ncol(df_final)])


mat_clean_final_dist_euclidean <- dist(mat_clean_final, method = "euclidean", diag = TRUE, upper = TRUE)

mat_clean_final_dist_manhattan <- dist(mat_clean_final, method = "manhattan", diag = TRUE, upper = TRUE)


# reduce dimension using PCA for later plotting --------------------------------------------------------

pca_out <- PCA(mat_clean_final, scale.unit = FALSE, ncp = col_end, graph = FALSE) # since already scaled, not need to do it again


# Move to cluster analysis (k-means and hierarchical) ------------------------------------------------

# define cluster number ------------------------------------------------------------------

cluster_min <- 2

cluster_max <- 10

num_samples <- nrow(mat_clean_final)


#  Part 1: k-means --------------------------------------------------------

name_algorithm <- "kmeans"


# get the best number of clusters ----------------------------------


# create matrix to store the original cluster
mat_cluster_origin <- matrix(data = NA, nrow = num_samples, ncol = (cluster_max - cluster_min + 1))

cluster_wss_all <- c()

count <- 1

for (k_centers in cluster_min:cluster_max){
  
  origin_cluster <- kmeans(x = mat_clean_final, centers = k_centers, iter.max = 25, nstart = 50)
  
  cluster_origin <- origin_cluster$cluster
  
  cluster_stats <- cluster.stats(d = mat_clean_final_dist_euclidean, 
                                 clustering = cluster_origin)
  
  cluster_wss_all <- c(cluster_wss_all, cluster_stats$within.cluster.ss)
  
  mat_cluster_origin[, count] <- cluster_origin
  
  count <- count + 1
  
  cat(count, "finished...\n")
  
  
}


# create cluster data.frame -------------------------------------------------

df_cluster_origin <- data.frame(mat_cluster_origin)

names(df_cluster_origin) <- paste0(name_algorithm, "_k_", cluster_min:cluster_max)


# plot the elbow ----------------------------------------------------------

df_cluster_wss <- data.frame(cluster=cluster_min:cluster_max, wss=cluster_wss_all) 

wss_plot_func(df = df_cluster_wss,
              out_path = file.path("output", name_algorithm, paste0(name_algorithm, "_wss", ".pdf")),
              width = 4,
              height = 3)


# best number of clusters -----------------------------------------------

k_centers <- 6

cluster_origin <- df_cluster_origin[[paste0(name_algorithm, "_k_", k_centers)]]

df_origin_cluster_pca <- data.frame(x=pca_out$ind$coord[,1], 
                                    y=pca_out$ind$coord[,2],
                                    Cluster=as.factor(cluster_origin))



# compare cluster results with real garnet types

origin_cluster_type <- table(df_origin_cluster_pca$Cluster, df_clean_final$Type)

origin_cluster_type

df_origin_cluster_type <- data.frame(origin_cluster_type) # this will change deposit sytle to factor

names(df_origin_cluster_type)

df_origin_cluster_type <- df_origin_cluster_type %>% 
  dplyr::rename(Cluster=Var1,
                Type=Var2,
                Count=Freq)

# spread percentage of each type
df_origin_cluster_percent_type <- df_origin_cluster_type

df_origin_cluster_percent_type <- df_origin_cluster_percent_type %>% 
  group_by(Type) %>% 
  mutate(count_cluster = sum(Count)) %>% 
  ungroup() %>% 
  mutate(percent_type = round(Count / count_cluster, 2))

df_origin_cluster_percent_type


colplot_spread_percent_func(df = df_origin_cluster_percent_type,
                            x = "Cluster",
                            y = "percent_type",
                            xlab = "Cluster",
                            ylab = "Proportion of garnet type",
                            color_par = "Type",
                            color_sequence = color_sequence_func(6),
                            name_legend = "Garnet type",
                            out_path = file.path("output", name_algorithm, paste0(name_algorithm, "_colplot_origin_cluster_type_spread_count_percent_", k_centers, ".pdf")),
                            width = 6,
                            height = 4)




#  Part 2: hierarchical --------------------------------------------------------

name_algorithm <- "hc"


# get the best number of clusters ----------------------------------

# create matrix to store the original cluster
mat_cluster_origin <- matrix(data = NA, nrow = num_samples, ncol = (cluster_max - cluster_min + 1))

cluster_wss_all <- c()

count <- 1

for (k_centers in cluster_min:cluster_max){
  
  hc_out <- hclust(mat_clean_final_dist_euclidean, method = "ward.D2")
  
  cluster_origin <- cutree(hc_out, k = k_centers)
  
  cluster_stats <- cluster.stats(d = mat_clean_final_dist_euclidean, 
                                 clustering = cluster_origin)
  
  cluster_wss_all <- c(cluster_wss_all, cluster_stats$within.cluster.ss)
  
  mat_cluster_origin[, count] <- cluster_origin
  
  count <- count + 1
  
  cat(count, "finished...\n")
  
}


# create cluster data.frame -------------------------------------------------

df_cluster_origin <- data.frame(mat_cluster_origin)

names(df_cluster_origin) <- paste0(name_algorithm, "_k_", cluster_min:cluster_max)


# plot the elbow ----------------------------------------------------------

wss_plot_func <- function(df, out_path, width, height){
  
  gp <- ggplot(df) +
    geom_line(aes(x=cluster, y=wss), size=0.5) +
    geom_point(aes(x=cluster, y=wss), size=1) +
    geom_vline(xintercept = 6, lty="dashed") +
    labs(x="Number of clusters", y="Total within cluster variation") +
    scale_x_continuous(breaks = pretty_breaks(6), expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size=8),
          axis.title = element_text(size=10)
    )
  
  print(gp)
  
  ggsave(filename = out_path, width = width, height = height)
  
}


df_cluster_wss <- data.frame(cluster=cluster_min:cluster_max, wss=cluster_wss_all) 

wss_plot_func(df = df_cluster_wss,
              out_path = file.path("output", name_algorithm, paste0(name_algorithm, "_wss", ".pdf")),
              width = 4,
              height = 3)


# best number of clusters -----------------------------------------------

k_centers <- 6

cluster_origin <- df_cluster_origin[[paste0(name_algorithm, "_k_", k_centers)]]

df_origin_cluster_pca <- data.frame(x=pca_out$ind$coord[,1], 
                                    y=pca_out$ind$coord[,2],
                                    Cluster=as.factor(cluster_origin))



# compare cluster results with real garnet types

origin_cluster_type <- table(df_origin_cluster_pca$Cluster, df_clean_final$Type)

origin_cluster_type

df_origin_cluster_type <- data.frame(origin_cluster_type) # this will change deposit sytle to factor

names(df_origin_cluster_type)

df_origin_cluster_type <- df_origin_cluster_type %>% 
  dplyr::rename(Cluster=Var1,
                Type=Var2,
                Count=Freq)

# spread percentage of each type
df_origin_cluster_percent_type <- df_origin_cluster_type

df_origin_cluster_percent_type <- df_origin_cluster_percent_type %>% 
  group_by(Type) %>% 
  mutate(count_cluster = sum(Count)) %>% 
  ungroup() %>% 
  mutate(percent_type = round(Count / count_cluster, 2))

df_origin_cluster_percent_type


colplot_spread_percent_func(df = df_origin_cluster_percent_type,
                            x = "Cluster",
                            y = "percent_type",
                            xlab = "Cluster",
                            ylab = "Proportion of garnet type",
                            color_par = "Type",
                            color_sequence = color_sequence_func(6),
                            name_legend = "Garnet type",
                            out_path = file.path("output", name_algorithm, paste0(name_algorithm, "_colplot_origin_cluster_type_spread_count_percent_", k_centers, ".pdf")),
                            width = 6,
                            height = 4)




# now select 30 samples for each type ---------------------------------------------------

# set.seed(1)

# df_final <- df_select_mineral_short_transform %>% 
#   group_by(Type) %>% 
#   sample_n(30) %>% 
#   ungroup()
# 
# table(df_final$Type)