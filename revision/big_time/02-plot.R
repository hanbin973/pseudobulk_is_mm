suppressPackageStartupMessages({
    library(Matrix)
    library(tidyverse)
    library(data.table)
    library(matrixStats)
    library(gridExtra)
    library(svglite)
})

# get parameters
ns <- as.integer(snakemake@params[['ns']])
out_path <- snakemake@output[[1]]
in_path <- snakemake@input

# read file
dfs <- lapply(
              in_path,
              function(x){
                  read.csv(x, row.names=1)
              }
              )
for (i in 1:length(ns)){
    dfs[[i]][,'ns'] <- ns[[i]]
}
df <- dplyr::bind_rows(dfs)
df$ns <- as.character(df$ns)
df$nc <- as.character(df$nc)
df_plot <- df %>% 
    group_by(method,ns,nc) %>%
    summarise(t=mean(time), sd=sqrt(var(time)))

# plot
p <- ggplot(df_plot, aes(fill=method, x=ns, y=t, color=nc)) +
    geom_bar(position='dodge', stat='identity') + 
    theme(text = element_text(size = 15)) +
    labs(
        title = "",
        x ="Number of subjects", 
        y = "Time (minutes)"
    ) +
    geom_errorbar(aes(ymin=t-2*sd, ymax=t+2*sd), width=.6,position=position_dodge(.9)) +
    geom_text(aes(x=ns, y=t+2*sd+9, label=str_c('n=',nc)), size=4, position=position_dodge(width = .9), angle=60) +
    geom_text(aes(x=ns, y=t+2*sd+4, label=str_c('(',round(t,1),')')), size=3, position=position_dodge(width = .9), angle=60) +
    scale_color_manual(values = rep('black',5)) +
    guides(color="none") +
    scale_fill_discrete(name = "Method", labels = c("NB GLMM", "Offset Pb")) +
    ylim(c(0, 65))

ggsave(out_path, height=5, width=10)

