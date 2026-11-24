# Data sources and terms of use

The input data are public. They are not included in this repository, except the two small files noted below. Check the terms of each source before reuse.

| Data | Content and year | Source | Credit line |
|---|---|---|---|
| Population Census, small-area tabulation and boundaries | 2020; Tables 3 (age), 6-1 and 9 (households), 7-1 (tenure of housing), 8-1 (type of building), 12 (labor force and occupation); boundaries with population | Statistics Bureau of Japan, e-Stat (https://www.e-stat.go.jp). Government of Japan Standard Terms of Use 2.0 | Source: e-Stat, 2020 Population Census small-area tabulation and boundary data, processed |
| Road Traffic Census, section tables and hourly traffic tables | Fiscal year 2021 | Ministry of Land, Infrastructure, Transport and Tourism (MLIT) (https://www.mlit.go.jp/road/census/r3/) | Source: MLIT, Road Traffic Census (FY2021), processed |
| Census section geometry (`traffic_census_2021_converted.parquet`) | Fiscal year 2021 | GitHub shiwaku/mlit-road-traffic-census-converter, release data-v1 (MIT). Original data: MLIT | Source: MLIT, Road Traffic Census (FY2021), processed |
| Continuous roadside noise monitoring | Fiscal years 2019–2023 | Ministry of the Environment, from the Environmental Observatory of the National Institute for Environmental Studies (https://tenbou.nies.go.jp) | Source: Ministry of the Environment, results of continuous monitoring of road traffic noise (via NIES), processed |
| Air quality monitoring | Fiscal year 2022 (used by the NO2 model) | Ministry of the Environment, from the Environmental Observatory of NIES | Source: NIES Environmental Numerical Database, processed |
| Vital Statistics | 2022; deaths by five-year age group (statInfId 000040067346) | Ministry of Health, Labour and Welfare, e-Stat | Source: Vital Statistics 2022, processed |
| Abridged Life Table | 2022 | Ministry of Health, Labour and Welfare | Source: Abridged Life Table 2022, processed |
| National Land Numerical Information, land use subdivision mesh (L03-b) | 2021; 100 m grid | MLIT (https://nlftp.mlit.go.jp/ksj/) | Source: MLIT National Land Numerical Information (land use subdivision mesh), processed |
| Historical weather data | 2018–2023, July–September | Japan Meteorological Agency (https://www.data.jma.go.jp) | Source: Japan Meteorological Agency, processed |
| Receptor-level nighttime noise and NO2 | Traffic of fiscal year 2021; NO2 of fiscal year 2022 | Matsushita, D. Zenodo https://doi.org/10.5281/zenodo.24299015 (CC BY 4.0) | Cite the dataset |
| NOx emission factors | NILIM Technical Note No. 671 (2012), Table 6.4 | `inputs/nox_ef_nilim671.csv` of japan-traffic-no2 | Source: NILIM Technical Note No. 671 |

## Small input files included in `data/`

| File | Content | Note |
|---|---|---|
| `muni_to_area335.csv` | Municipality code and the secondary medical care area (335 areas) that contains it | Read by `burden_v3.py` and `fig2_v8.py` for a secondary calculation on tropical nights. Its results are not reported in the paper |
| `ndb_insufficient_rest_area_2018_2023.csv` | Proportion reporting insufficient rest from sleep by area, from the NDB Open Data of the Ministry of Health, Labour and Welfare (fiscal years 2018–2023) | Same use as above. Source: https://github.com/m138/japan-hotnight-sleep |
| `yll_per_death_30plus_lifetable2022.csv`, `age2_national.csv` | Years of life lost per death at age 30 and over; national crude death rate at age 30 and over | Derived from the 2022 abridged life table, the 2022 Vital Statistics and the 2020 Population Census |

## Methods and coefficients taken from the literature

- WHO exposure–response relationship for high sleep disturbance from road traffic noise: Basner & McGuire, Int. J. Environ. Res. Public Health 15, 519 (2018).
- Japanese exposure–response curve (sensitivity analysis): Morinaga, Yokoshima & Morihara, Int. J. Environ. Res. Public Health 22, 1263 (2025).
- Source power levels and porous asphalt equations: Sakamoto, ASJ RTN-Model 2018, Acoust. Sci. Technol. 41, 529–589 (2020).
- Rolling and propulsion noise coefficients: CNOSSOS-EU, Commission Delegated Directive (EU) 2021/1226.
- Weights of the areal deprivation index: Nakaya, Proc. Inst. Stat. Math. 59, 239–265 (2011).
