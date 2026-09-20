# Data dictionary

The CSV datasets are licensed under [CC BY 4.0](LICENSE).

The files contain the focal inventories, regional araucaria records, and remnant
attributes, standardized from the field records. The focal adult inventory retains
all species for the community-level structural analysis and species-abundance
table. The demographic model uses the ***Araucaria angustifolia*** subset.

Run the analyses with `R/run_all.R`.

> Column names and categorical values are kept in Portuguese, as recorded in the
> field. Their meanings are given below.

---

## `focal_arboreo.csv`

Adult component (DBH ≥ 5 cm), all species; 50 plots of 10 × 20 m (≈ 1 ha),
inventoried in 2008, 2012, 2016 and 2021.

| column | type | description |
|---|---|---|
| `transecto` | integer | transect of the plot network |
| `parcela` | integer | plot within the transect |
| `arvore` | text | tree identifier (may carry a stem suffix) |
| `especie` | text | scientific name |
| `ano` | integer | inventory year (2008/2012/2016/2021) |
| `dap_cm` | numeric | diameter at breast height (cm) |
| `ab_cm2` | numeric | basal area of the individual (cm²; summed over stems) |
| `status` | text | `viva` (alive), `morta` (dead) or `recrutamento` (recruit) |

Unique tree identifier = `transecto_parcela_arvore`.

## `focal_regeneracao.csv`

Natural regeneration (size classes), all species, nested subplots, inventoried
annually from **2011 to 2014**. (The 2010 establishment inventory was excluded
because small seedlings were undercounted; it is not present in the data.)

| column | type | description |
|---|---|---|
| `transecto` | integer | transect code (T in the field spreadsheet) |
| `parcela` | integer | plot code within the transect (P in the field spreadsheet) |
| `subparcela` | integer | continuous sampling-unit code (PC, 1–50), shared by the nested size-class sampling areas |
| `individuo` | text | field tag (N in the field spreadsheet; may repeat between species) |
| `especie` | text | scientific name |
| `ano` | integer | inventory year (2011–2014) |
| `classe` | text | state in that year: `1` (C1, 0.20 ≤ height < 1 m), `2` (C2, 1–3 m), `3` (C3, > 3 m and DBH < 5 cm), `M` (dead), `AD` (graduated to the adult class), `NA` (absent) |

Within the araucaria subset, `subparcela` and `individuo` uniquely identify each
plant. The analysis scripts use `transecto_parcela_subparcela_individuo`, which
identifies the same plants. Across all species, include `especie` with
`subparcela` and `individuo` to distinguish repeated tags.
Nested sampling areas (per plot): C1 = 5 m², C2 = 10 m², C3 = 20 m².

## `multisitio_araucaria.csv`

Adult ***Araucaria angustifolia*** (DBH ≥ 5 cm) from the eight remnants of the
regional comparison (Figure 6). Araucaria only. One record per tree × inventory.

| column | type | description |
|---|---|---|
| `sitio` | text | remnant code (PB1 = focal area; PB2, PAI, BJS, Vargem, Epagri, Alu, PARNAMUL) |
| `individuo` | text | unique identifier within the site; at PB1, `PB1_transecto_parcela_arvore` links to the focal inventory |
| `ano` | integer | inventory year |
| `dap_cm` | numeric | DBH (cm). At PB1, the individual DBH is taken directly from `focal_arboreo.csv`, including equivalent DBH for multi-stemmed individuals. At PARNAMUL, where girth was measured per stem, this is the **equivalent DBH** = √Σ(girthᵢ/π)², which preserves basal area |

PARNAMUL has **four inventories** (2011, 2015, 2019, 2023), taken from columns
CAP1–CAP4 of the source spreadsheet (years in ANO1–ANO4). PB1 uses living and newly
recruited individuals from the focal inventory. The other six sites come from the
multi-site source file. Regional basal area is calculated from the published DBH
values; minor rounding differences from the focal inventory's stored basal area
may occur.

## `sites.csv`

**Observational** (field) attributes of each remnant in the regional comparison —
the input for Table S1.1 (site characterization). It is not derived: it is the
source of the coordinates, sampling effort and elevation.

| column | type | description |
|---|---|---|
| `code` | text | site code (matches `sitio` in `multisitio_araucaria.csv`) |
| `nome` | text | display name of the remnant |
| `municipio` | text | municipality (Santa Catarina) |
| `lat`, `lon` | numeric | centroid coordinates (decimal degrees, WGS84) |
| `parcelas` | integer | number of plots sampled |
| `area_parcela_m2` | numeric | area of each plot (m²) |
| `altitude_m_srtm` | integer | elevation (m), SRTM 30 m — cached; recorded value used in Table S1.1 |
| `foco` | logical | `TRUE` only for the focal population (PB1) |

`R/04_regional.R` combines these attributes with the regional inventory to produce Table S1.1. Elevations use the recorded SRTM values in `sites.csv`.
