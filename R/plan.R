source("packages.R")
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["regular"]]

the_plan <-  drake_plan(
  f2 = create_fig2(),
  f3 = create_fig3(),
  f4 = create_fig4(),
  f6 = create_fig6(),
  f7 = create_fig7(),
  geosub = create_geo_data()
)


