message('Installing autocalb R packages...')
install.packages(
  c(
    "nwsrfsr",
    "box",
    "dplyr",
    "data.table",
    "dtplyr",
    "hydroGOF",
    "digest",
    "lubridate",
    "readr",
    "tibble",
    "ggplot2",
    "ggthemes",
    "crayon",
    "argparser",
    "rtop",
    "stringr",
    "tidyr",
    "vctrs",
    "plotly",
    "gridExtra",
    "rlang",
    "rmarkdown"
  ),
  repos = "https://packagemanager.posit.co/cran/latest"
)
