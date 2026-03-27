####################################################################################
###################################### MAGICK ######################################
####################################################################################

https://cran.r-project.org/web/packages/magick/vignettes/intro.html#Animation

install.packages("magick")
library(magick)

The principle is : image_read() %>% image_join() %>% image_animate()


####################################################################################
###### Exemple with already existing images


path <- "report/4_hiPAR_further_analysis/Highlights/"

## list file names and read in
imgs <- list.files(path, full.names = TRUE, pattern = ".png$")
# Remmettre dans l'ordre des jours (pour pas que 10 sorte avant 2)
imgs <- imgs[order(as.numeric(sub(".*_D([0-9]+)_.*", "\\1", imgs)))]
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_morph(img_joined, frames = 4) %>% image_animate(fps = 4) 

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "report/4_hiPAR_further_analysis/Highlights/tx-sales.gif")
