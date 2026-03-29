library(ggplot2)
library(ggthemes)
library(stringr)
library(scales)

theme_Publication <- function(base_size=16, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family = "")
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

setwd("~/EinsteinMed Dropbox/Ronald Cutler/Vijg-lab/Collaborations/Jackson Rogow/transcriptomics/analysis/with-virus-and-stem-loop/batch-2-arc-rep1-removed/network-analysis")

df <- read.csv("clusters 1-6 top GO enrichment biological process.csv",
               header = TRUE)
#df <- df[order(df$false.discovery.rate, decreasing = FALSE),]

df$term.description <- factor(df$term.description, levels = rev(df$term.description))

# Plot with reversed x-axis
pdf("clusters 1-6 top GO enrichment biological process.pdf",
    width = 7,
    height = 5)
ggplot(df, aes(x=term.description, y=false.discovery.rate, color=term.description, size=1.5, color="black")) +
  geom_segment(aes(x=term.description, xend=term.description, y=max(false.discovery.rate), yend=false.discovery.rate),
               size = 1) +
  geom_point(size=7) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) + # Wrap text for x-axis
  scale_y_reverse(labels = label_scientific()) + # Use scientific notation for y-axis
  ggtitle("GO enrichment of clusters") +
  ylab("False Discovery Rate") +
  xlab("") + 
  coord_flip() +
  theme_Publication() +
  scale_color_manual(name = "Cluster",
                     values = rev(c("#ffa9a4", "#ffebc4", "#c7e1a9", "#c1ffc5", "#99e0ce", "#c1e5ff")),
                     labels = rev(c(1,2,3, 4, 5, 6))) +
  guides(color = guide_legend(reverse = TRUE)) +      # Reverse the order of the legend
  theme(
    plot.title.position = "plot",           # Align title with plot area
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),  # Increase size of x-axis text
    legend.position = "right",              # Move legend to the right
    legend.direction = "vertical",           # Set legend direction to vertical
  )
dev.off()

png("clusters 1-6 top GO enrichment biological process.png",
    width = 7,
    height = 5,
    units = "in",
    res = 300)
ggplot(df, aes(x=term.description, y=false.discovery.rate, color=term.description, size=1.5, color="black")) +
  geom_segment(aes(x=term.description, xend=term.description, y=max(false.discovery.rate), yend=false.discovery.rate),
               size = 1) +
  geom_point(size=7) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) + # Wrap text for x-axis
  scale_y_reverse(labels = label_scientific()) + # Use scientific notation for y-axis
  ggtitle("GO enrichment of clusters") +
  ylab("False Discovery Rate") +
  xlab("") + 
  coord_flip() +
  theme_Publication() +
  scale_color_manual(name = "Cluster",
                     values = rev(c("#ffa9a4", "#ffebc4", "#c7e1a9", "#c1ffc5", "#99e0ce", "#c1e5ff")),
                     labels = rev(c(1,2,3, 4, 5, 6))) +
  guides(color = guide_legend(reverse = TRUE)) +      # Reverse the order of the legend
  theme(
    plot.title.position = "plot",           # Align title with plot area
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),  # Increase size of x-axis text
    legend.position = "right",              # Move legend to the right
    legend.direction = "vertical",           # Set legend direction to vertical
  )
dev.off()
