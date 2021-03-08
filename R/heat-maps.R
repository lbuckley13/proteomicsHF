library(ggplot2)
library(tidyverse)

load('~/Desktop/lassoversions.RData')
load('data-saves/Results.RData')


strain.data <- haven::read_dta('~/Dropbox (Partners HealthCare)/Teams/drshah_shared/ARIC V5 LA strain/Master_ARIC_LA_Analysis.dta')
strain.vars <- names(strain.data)[-1]
strain.data <- strain.data %>% mutate_at(strain.vars, scale, T, T)

data.st <- strain.data %>% right_join(fifth.visit.echo.scaled, by = "id")


echo.data <- fifth.visit.echo.scaled %>% select(all_of(echo.vars))

eigen.data <- bind_cols(eigenData %>% select(-fuptime, -hfdiag), echo.data)
mods <- grep("^ME", names(eigen.data), value = TRUE)
eigen.all.scaled <- lin.arry.aggr(mods, eigen.data, echo.vars, adjust, labels)


heatmapper <- eigen.all.scaled$all.res %>% select(term, outcome, estm)


png(filename = "~/Desktop/modules.png", width=1200, height=600)
ggplot(data = heatmapper, aes(term, outcome, fill = estm))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Linear Regression") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 8, hjust = 1))+
  coord_fixed()
dev.off()

echo.all.sign <- univariate.echo.models.unscaled.alladj$all.res %>%
  filter(pval <= 0.05 / 4877)

png(filename = "~/Desktop/allsign.png", width=12000, height=1600)
ggplot(data = echo.all.sign, aes(term, outcome, fill = estm))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="Linear Regression") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 8, hjust = 1))+
  coord_fixed()
dev.off()
