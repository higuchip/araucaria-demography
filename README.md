# Araucaria demography

Data and R scripts for **Adult longevity and a persistent seedling bank in Araucaria angustifolia: evidence from historically logged forests** (Higuchi and co-authors; manuscript under review).

## Run

Using R 4.5.3, open the repository folder as the working directory and run:

```r
install.packages(c("data.table", "popbio", "Rage", "ggplot2",
                   "sf", "ggspatial", "cowplot", "geodata", "rnaturalearth", "rnaturalearthdata"))
source("R/run_all.R")
```

The model uses 12,000 posterior draws and seed 20260912. Figure 1 requires internet access on the first run to download the geographic layers; subsequent runs reuse the local cache.

## Data and results

Input data are in `data/`; variables and sampling areas are described in the [data dictionary](data/README.md). Analysis and figure scripts are in `R/`.

The `out/` folder is created locally when the scripts are run and is not included in the repository. It contains:

- `out/tables/`: Tables 1–2 and Supplementary Tables S1.1, S1.2, and S2.1–S2.4 (CSV).
- `out/figures/`: Figures 1–6 and Supplementary Figure S1 (TIFF, 600 dpi).

## Sources and licenses

Figure 1 uses GADM 4.1 boundaries, Natural Earth country outlines (1:50 million), and the IBGE vegetation map (1:5,000,000). Geographic layers are downloaded from their sources and retain their original terms: [GADM](https://gadm.org/license.html), [Natural Earth](https://www.naturalearthdata.com/about/terms-of-use/), and [IBGE](https://geoftp.ibge.gov.br/informacoes_ambientais/vegetacao/vetores/brasil_5000_mil/).

R scripts: [MIT](LICENSE). The four input CSV datasets: [CC BY 4.0](data/LICENSE). Please cite the associated article when using these materials.
