library(hexSticker)

imgurl <- "man/figures/hex_graphic.png"

s <- sticker(imgurl, package="pcsstools", 
             p_size=80, p_color = "#b7b7b7ff", 
             s_x=1, s_y=1, s_width=0.8,
             h_color = "#880000ff",
             h_fill = "#0D0630ff",
             filename = "man/figures/logo.png",
             dpi = 1000)
# plot(s)
