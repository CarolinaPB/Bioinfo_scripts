# Save GT tables to word (or other office programs)
Use [gto package](https://github.com/GSK-Biostatistics/gto)
```{r}
#load officer and gt
library(officer)
library(gt)
library(gto)

## create simple gt table
gt_tbl <- your_data %>% gt()

## Create docx and add gt table
doc <- read_docx() # add path to file to which you want to add a table. Leave empty if you want to add table to emprty document
doc <- body_add_gt(doc, value = gt_tbl)

## Save docx
fileout <- tempfile(tmp_dir = "~/Downloads", fileext = ".docx") # A file with a unique ID containing your table will be saved to the downloads folder
print(doc, target = fileout)
```

# Combine several .png into one figure
```{r}
library(cowplot)

plt1 <- ggdraw() +  draw_image(<filename1.png>)
plt2 <- ggdraw() +  draw_image(<filename2.png>)
plt_combined <- plot_grid(plt1, plt2, ncol = 1)
save_plot(glue("combined.pdf"), plt_combined)
```
