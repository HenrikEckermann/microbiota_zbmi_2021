library(tidyverse)
library(Hmisc)
library(corrplot)
library(qvalue)
library(RColorBrewer)



# you must first run import_data.R to create the objects we load here:
load(here::here("rdata/mb_import_H.Rds"))

# create correlation matrices using rcorr
cor_fb_bmi <- select(df_cor, zBMI_1m:zBMI_12y, FB_ratio_1m:FB_ratio_10y) %>% 
  as.matrix() %>% 
  Hmisc::rcorr(type = "spearman")


cor_scfa_bmi <- select(df_cor, zBMI_1m:zBMI_12y, SCFA_prod_1m:SCFA_prod_10y) %>% 
  as.matrix() %>% 
  Hmisc::rcorr(type = "spearman")


# color palette for cor matrix
col <- colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))


# visualize correlation matrix fb ratio
plot1 <- corrplot(cor_fb_bmi$r, method = "color", 
         col = col(200), 
         type = "upper", 
         addCoef.col = "black", 
         p.mat = cor_fb_bmi$P, 
         tl.col = "black", tl.srt=30, 
         sig.level = 0.05, insig = "blank", 
         diag =F)

  # add text of insignificant correlations
text(plot1$corrPos$x, plot1$corrPos$y, round(plot1$corrPos$corr,2), font = 1)


# visualize correlation matrix scfa producers
plot2 <- corrplot::corrplot(cor_scfa_bmi$r, method = "color", 
         col = col(200),
         type = "upper", 
         addCoef.col = "black",
         p.mat = cor_scfa_bmi$P,
         tl.col = "black", tl.srt=30, 
         sig.level = 0.05, insig = "blank", 
         diag =F)

  # add text of insignificant correlations
text(plot2$corrPos$x, plot2$corrPos$y, round(plot2$corrPos$corr,2), font = 1)


# apply Benjamini-Hochberg procedure using qvalue
rownames <- c(row.names(cor_fb_bmi$P), row.names(cor_scfa_bmi$P))[c(9:13,22:26)]

qvalues <- full_join(data.frame(cor_fb_bmi$P[c(9:13),c(1:8)]), data.frame(cor_scfa_bmi$P[c(9:13),c(1:8)])) %>% 
  mutate(var1 = rownames) %>%
  pivot_longer(1:8, names_to = "var2", values_to = "p") %>% 
  select(var1, var2, p) %>% 
  mutate(q = qvalue::qvalue(p, lambda = 0)$qvalues) %>%
  arrange(q)

qvalues
