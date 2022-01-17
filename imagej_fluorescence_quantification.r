# R script for visualization of differences in cellular fluorescence 
# Written by Thomas E. Smith in Dec 2020
# Last updated by TES Jan 2022

#### NOTE
# This script assumes the User is in the same directory as 
# TableS9_fluorescence_microscopy.csv, also available at 
# https://zenodo.org/record/5484415#.YTob355KgWp

# setting the stage
library("plyr")
library("tidyverse")
library("stats")
library("FSA")
library("rcompanion")
library("ggplot2")



### Import ImageJ data
fl.dat <- read.csv("TableS9_fluorescence_microscopy.csv", 
                   header = T, stringsAsFactors = F)



### Calculate the average fluorescence intensity per cell area
## This should account for differences in cell size
fl.dat <- fl.dat %>%
  mutate( AFI.area = Mean/Area )



### Set up some generic plotting parameters
condition.order <- unique(fl.dat$Condition)
ch.order <- c("AF488", "WGA640", "DAPI")
ch.colors <- c("green3", "magenta3", "cyan3")
outlier.colors <- c("C7E7E8", "BDE6BD", "E3B6E3")

# Order the data
# Convert Condition and Ch columns into numbers to facilitate dodging at each
# position on the x-axis. Dodging resets the column order alphabetically, so 
# converting to numbers first will maintain the current order.
fl.dat <- fl.dat %>%
  mutate( Condition = as.numeric(factor(Condition, levels = condition.order)),
          Ch = as.numeric(factor(Ch, levels = ch.order)) ) %>%
  # Arrange the dataset first by compound and then by sample type
  arrange( Ch, Condition )



### Perform statistical tests 

## Look at the distribution of data
# Right-skewed data
hist(log(fl.dat$AFI.area), labels = T, breaks = 20)
# Is it normally distributed?
ntest <- shapiro.test(log(fl.dat$AFI.area))
if (ntest$p.value < 0.05) { print("Data is not normally distributed") } 

## Perform Kruskal-Wallis test for each Ch, then extract pairwise comparisons
fl.test <- fl.dat %>%
  select( Ch, Condition, AFI.area ) %>%
  group_by( Ch ) %>% 
  dplyr::summarize( Comparison = dunnTest( AFI.area ~ as.factor(Condition), 
                                           method = "bh" )$res$Comparison,
                    p.adj = dunnTest( AFI.area ~ as.factor(Condition), 
                                      method = "bh" )$res$P.adj )



### Construct box plot of data distributions for each channel

# Loop through channels
for (i in 1:length(ch.order)) {
  # Filter data to include only values for this channel
  df <- fl.dat %>%
    filter( Ch == i )
  
  # Establish the max y-axis value divisible by 10
  y.max <- ceiling( max(df$AFI.area) )
  
  # Construct a compact letter display for the statistical tests
  # of this channel
  df.test <- fl.test %>%
    # Filter by channel
    filter( Ch == i )
  let.disp <- cldList(comparison = df.test$Comparison,
                      p.value = df.test$p.adj,
                      threshold = 0.01 ) %>% 
    as_tibble() %>%
    # Change the column Group to Condition
    mutate( Condition = as.numeric(Group) ) %>%
    select( -c(MonoLetter, Group) )
  
  
  # Calculate sample sizes and upper whisker cutoff value for each Condtion to 
  # determine where to place sample size labels and letter displays
  df.labels <- df %>%
    group_by( Condition ) %>%
    dplyr::summarize( N = paste0( "n=", n() ), 
                      Upper.whisker = max( AFI.area[ AFI.area < ( 1.5*IQR(AFI.area) + quantile(AFI.area, 0.75) ) ] ) ) %>%
    left_join( let.disp, by = "Condition" ) %>%
    group_by( Condition ) %>%
    # Determine y-value for text label by adding to Max
    mutate( n.Height = Upper.whisker + 0.1*y.max,
            let.Height = Upper.whisker + 0.2*y.max )
  
  # Construct plot
  q <- ggplot() +
    geom_boxplot( data = df, 
                  aes(y = AFI.area, x = as.factor(Condition)),
                  color = ch.colors[i], size = 1, 
                  outlier.alpha = 0.2, outlier.size = 3, 
                  show.legend = F ) +
    # Add custom x-axis labels
    scale_x_discrete(labels = condition.order) +
    # Add sample sizes above boxplots
    #geom_text( data = df.labels,
    #           aes( x = Condition, y = n.Height, label = N ),
    #           show.legend = F, size = 6 ) +
    # Add compact letter display above sample sizes
    geom_text( data = df.labels,
               aes( x = Condition, y = let.Height, label = Letter ),
               show.legend = F, size = 6 ) +
    theme(plot.title = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.grid.major =  element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_line(size = 1, colour = "black"),
          axis.ticks.length = unit(7, "pt"),
          axis.line.y = element_line(color = "black", size = 1),
          axis.line.x = element_line(color = "black", size = 1),
          axis.text = element_text(color= "black", size = 18),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 14, angle = 45, 
                                     hjust = 1, vjust = 1),
          axis.title = element_blank() )
  
  # Establish plot width
  plot.width <- 3
  
  # Export as pdf 
  pdf(paste0("Comparison_", ch.order[i], "_plot.pdf"), 
      width = plot.width, height = 3)
  plot(q)
  dev.off()
  
}



