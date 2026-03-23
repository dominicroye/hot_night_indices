# Hot-Night Indices for R

R functions to compute **hot-night duration (HNd)** and **hot-night excess (HNe)** from hourly temperature series, using astronomical sunset/sunrise windows as the definition of night. A companion helper computes the location-specific temperature threshold as a fixed or centred moving percentile.

These indices were introduced in [Royé (2017)](https://doi.org/10.1007/s00484-016-1268-6) and applied at a global scale in:

> **Royé D**, Sera F, Tobías A, Hashizume M, Honda Y, Kim H, Vicedo-Cabrera AM, Tong S, Lavigne E, Kyselý J, Pascal M, de'Donato F, Pereira da Silva SN, Madureira J, Huber V, Urban A, Schwartz J, Bell ML, Armstrong B, Iñiguez C, et al. (2025).
> Short-term association between hot nights and mortality: a multicountry analysis in 178 locations considering hourly ambient temperature.
> *Environment International*.
> <https://doi.org/10.1016/j.envint.2025.109526>

---

## Indices

### Hot-night duration — HNd

The proportion of night hours during which air temperature exceeds a threshold *T*thr:

$$\text{HNd}_j = \frac{\sum_{i} \mathbf{1}(T_{ij} > T_{\text{thr}})}{N_j} \times 100$$

where *N*ⱼ is the total number of night hours on night *j* (from sunset to the following sunrise), and *T*ᵢⱼ is the observed hourly temperature. HNd is expressed as a percentage.

### Hot-night excess — HNe

The cumulative thermal excess above the threshold during the night:

$$\text{HNe}_j = \sum_{i} (T_{ij} - T_{\text{thr}}) \cdot \mathbf{1}(T_{ij} > T_{\text{thr}})$$

HNe is expressed in degree-hours (°C·h).

### Threshold

Following Royé et al. (2025), *T*thr is defined as the **moving 95th percentile of daily minimum temperature** computed over a centred 91-day window. This location-specific, seasonally varying threshold captures local climatological conditions and is provided by `night_percentile()`.

---

## Installation

No package installation is required. Source the script directly:

```r
source("nightheat.R")
```

### Dependencies

| Package | Role |
|---|---|
| `lubridate` | Time-zone handling |
| `data.table` | Vectorised joins and rolling operations |
| `suncalc` | Astronomical sunset/sunrise times |
| `future` *(optional)* | Parallel back-end for `night_percentile()` |
| `future.apply` *(optional)* | Parallel `lapply` over station groups |

Install with:

```r
install.packages(c("lubridate", "data.table", "suncalc"))
# optional
install.packages(c("future", "future.apply"))
```

---

## Functions

### `night_duration(x, cols, tz, coord, thres)`

Computes **HNd** (percentage of night hours above threshold).

| Argument | Type | Description |
|---|---|---|
| `x` | data.frame / data.table | Input data with at least two columns |
| `cols` | character(2) | Names of the datetime and temperature columns |
| `tz` | character | Time-zone string (e.g. `"Europe/Madrid"`) |
| `coord` | numeric(2) | `c(longitude, latitude)` in decimal degrees (W = negative) |
| `thres` | numeric | Temperature threshold (°C) |

Returns a `data.table` with columns `date` and `HNd`.

---

### `night_excess(x, cols, tz, coord, thres)`

Computes **HNe** (degree-hours above threshold).

Arguments identical to `night_duration()`. Returns a `data.table` with columns `date` and `HNe`.

---

### `night_percentile(x, cols, p, by, type, window, workers)`

Computes a **fixed or moving percentile** of daily temperature, intended to derive the threshold for `nighthour()` / `nightdegrees()`.

| Argument | Default | Description |
|---|---|---|
| `x` | — | data.frame / data.table with date and temperature columns |
| `cols` | — | character(2): date and temperature column names |
| `p` | `0.95` | Percentile as a proportion |
| `by` | `NULL` | Grouping column(s), e.g. `"city"` |
| `type` | `"fixed"` | `"fixed"` for a single value per group; `"moving"` for a rolling window |
| `window` | `91` | Window size in days (must be odd; even values are rounded up) |
| `workers` | `1L` | Parallel workers (`> 1` requires `{future}` and `{future.apply}`) |

Returns a `data.table` with columns `date`, optional grouping columns, and `pXX` (e.g. `p95`).

---

## Quick start

```r
source("nightheat.R")

# --- Coordinates and time zone ---
coord <- c(-8.54, 42.88)   # Santiago de Compostela (lon, lat)
tz    <- "Europe/Madrid"

# --- Step 1: moving P95 of daily minimum temperature ---
daily_min <- your_hourly_data[
  , .(Ta = min(Ta)), by = .(date = as.Date(Datetime))
]

p95 <- night_percentile(
  x      = daily_min,
  cols   = c("date", "Ta"),
  p      = 0.95,
  type   = "moving",
  window = 91
)

thres <- median(p95$p95, na.rm = TRUE)

# --- Step 2: compute HNd and HNe ---
hnd <- night_duration(your_hourly_data,
                 cols = c("Datetime", "Ta"),
                 tz = tz, coord = coord, thres = thres)

hne <- night_excess(your_hourly_data,
                    cols = c("Datetime", "Ta"),
                    tz = tz, coord = coord, thres = thres)

result <- merge(hnd, hne, by = "date")
head(result)
#>          date      HNd      HNe
#> 1: 2022-01-01   0.0000   0.0000
#> 2: 2022-01-02  12.5000   1.4230
#> ...
```

A self-contained dummy example (synthetic hourly data for Santiago de Compostela, 2022–2023) is included at the bottom of `nightheat.R` inside an `if (FALSE) { ... }` block — copy and run any section interactively.

---

## Input format

The functions expect an hourly time series in wide format:

| Column | Type | Notes |
|---|---|---|
| Datetime | `POSIXct` | Any time zone; internally converted via `with_tz()` |
| Ta | `numeric` | Air temperature in °C |

Column names are passed via `cols` so any names are accepted.

---

## Implementation notes

Night windows are defined astronomically from sunset on day *d* to sunrise on day *d + 1*, using `suncalc::getSunlightTimes()`. Observations are assigned to their window with `data.table::foverlaps()`, an interval-overlap join that replaces the original day-by-day loop and scales to multi-decadal, multi-station series without the O(n · m) cost.

The moving percentile uses `data.table::frollapply(..., align = "center")`, which is equivalent to the lag-matrix approach of the original code but avoids constructing the N × 90 auxiliary matrix.

---

## Citation

If you use these functions, please cite the paper above. BibTeX:

```bibtex
@article{roye2025hotnights,
  author  = {Royé, Dominic and Sera, Francesco and Tobías, Aurelio and
             Hashizume, Masahiro and Honda, Yasushi and Kim, Ho and
             Vicedo-Cabrera, Ana Maria and Tong, Shilu and Lavigne, Eric and
             Kyselý, Jan and Pascal, Mathilde and {de'Donato}, Francesca and
             {Pereira da Silva}, Susana das Neves and Madureira, Joana and
             Huber, Veronika and Urban, Aleš and Schwartz, Joel and
             Bell, Michelle L. and Armstrong, Ben and Iñiguez, Carmen and
             others},
  title   = {Short-term association between hot nights and mortality:
             a multicountry analysis in 178 locations considering hourly
             ambient temperature},
  journal = {Environment International},
  year    = {2025},
  doi     = {10.1016/j.envint.2025.109526}
}
```

---

## License

This work is licensed under a
[Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

[![CC BY 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by/4.0/)
