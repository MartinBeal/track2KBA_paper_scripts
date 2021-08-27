## patchworking together track2KBA plots ## 

pacman::p_load(patchwork)

fs_map <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\fur_seals\\analysis\\FS_KBA_ggmap.rds")
gt_map <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\green_turtles\\analysis\\GT_KBA_ggmap_NINDs.rds")
ws_map <- readRDS("C:\\Users\\Martim Bill\\Documents\\mIBA_package\\data\\white_storks\\analysis\\WS_KBA_ggmap.rds")


## example workflow figure ## 

Tripmap / KDEmap + KBAmap



p <- (Tripmap + plot_layout(heights = 1.3)) + KDEmap + repPlot + KBAmap + plot_annotation(tag_levels = 'A') & 
  theme( 
    plot.tag.position = c(0, .95), 
    plot.tag = element_text(size = 20)
    )

ggsave("C:\\Users\\Martim Bill\\Desktop\\test\\AA7.png", p, width = 10, height=6)

gtab <- patchwork:::plot_table(p, 'auto')
overall_width <- grid::convertWidth(sum(gtab$widths) + unit(1, "cm"), unitTo = "cm", valueOnly = TRUE)
overall_height <- grid::convertHeight(sum(gtab$heights) + unit(1, "cm"), unitTo = "cm", valueOnly = TRUE)

ggsave("C:\\Users\\Martim Bill\\Desktop\\test\\AA9.png", p, width = overall_width, height = overall_height+2)


## ggmap maps together ## 


layout <- "
AAABB
AAABB
CCCBB
CCCBB
CCCBB
"

layout <- "
AAAABBB
AAAABBB
AAAABBB
CCCCBBB
CCCCBBB
CCCCBBB
CCCCBBB
CCCCBBB
"

fs_map + gt_map + ws_map + 
  plot_layout(design = layout) #+ plot_annotation(tag_levels = 'A')

ggsave("C:\\Users\\Martim Bill\\Desktop\\test\\A15.png", width = 10, height = 11)


### Green Turtle supplementary figures ###
layout <- "
AAA#####
AAABBBBB
AAABBBBB
AAABBBBB
AAA#####
"

# UDPLOT + repPlot_gt + 
#   plot_layout(design = layout) + 
#   plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30), plot.tag.position = c(.01, .98))

UDPLOT + repPlot_gt +  
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30), plot.tag.position = c(.90, .95))


ggsave("C:\\Users\\Martim Bill\\Desktop\\test\\B1.png", width = 9, height = 7.25)


### Egyptian Vulture supplementary figures ###

# layout <- "
# AAAAAA#######
# AAAAAABBBBBBB
# AAAAAABBBBBBB
# AAAAAABBBBBBB
# AAAAAABBBBBBB
# AAAAAABBBBBBB
# AAAAAA#######
# "

# UDPLOT + repPlot_gt + 
#   plot_layout(design = layout) + 
#   plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30), plot.tag.position = c(.01, .98))

UDPLOT + repPlot_eg +  
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30), plot.tag.position = c(.92, .96))


ggsave("C:\\Users\\Martim Bill\\Desktop\\test\\C1.png", width = 9, height = 7.25)

