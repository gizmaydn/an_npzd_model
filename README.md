# NPZD Marine Ecosystem Model in R

A numerical implementation of a zero-dimensional (0D) **NPZD** (Nutrient–Phytoplankton–Zooplankton–Detritus) biogeochemical model using the [`deSolve`](https://cran.r-project.org/package=deSolve) package in R.

The model simulates seasonal plankton bloom dynamics driven by periodic light and temperature forcing over a multi-year cycle, with mass conservation verification.

---

## Model Structure

The system models nitrogen / mass fluxes between four state variables:


[Nutrient (N)] <----------- (Mineralization) <----------- [Detritus (D)]
         |      ^                                                  ^       ^
 Uptake  |      | Respiration                                      |       |
         v      |                                           Mortality      Mortality
 [Phytoplankton (P)] ----(Grazing)----> [Zooplankton (Z)] ---------+-------+
         |                                     |
         +-------------(Respiration)-----------+





### State Variables

| Variable | Description | Initial Value |
| :--- | :--- | :--- |
| `nut` | Dissolved inorganic nutrient concentration | `0.99` |
| `pp` | Phytoplankton biomass | `0.01` |
| `zoo` | Zooplankton biomass | `0.01` |
| `detr` | Non-living particulate organic matter (detritus) | `0.99` |

---

## Governing Equations

$$\frac{dP}{dt} = \text{Uptake} - \text{Mortality}_P - \text{Grazing} - \text{Respiration}_P$$

$$\frac{dZ}{dt} = \text{Grazing} - \text{Mortality}_Z - \text{Respiration}_Z$$

$$\frac{dD}{dt} = \text{Mortality}_P + \text{Mortality}_Z - \text{Mineralization}$$

$$\frac{dN}{dt} = \text{Respiration}_P + \text{Respiration}_Z + \text{Mineralization} - \text{Uptake}$$

### Key Formulations

* **Light Limitation (Michaelis-Menten / Monod):**
  $$\text{light\_lim} = \frac{I}{k_{sLight} + I}, \quad I(t) = 1 + 0.5 \left|\sin\left(\frac{\pi t}{365}\right)\right|$$

* **Temperature Forcing (Eppley curve):**
  $$\text{temp}(t) = 2 + 18 \left|\sin\left(\frac{\pi t}{365}\right)\right|, \quad f(T) = \exp(0.063 \cdot T)$$

* **Nutrient Uptake:**
  $$\text{Uptake} = r_{\max} \cdot \text{light\_lim} \cdot \left(\frac{N}{k_n + N}\right) \cdot f(T) \cdot P$$

* **Zooplankton Grazing (Ivlev formulation):**
  $$\text{Grazing} = g_{\max} \cdot \left(1 - e^{-I_v \cdot P}\right) \cdot Z \cdot f(T)$$

---

## Parameter Values

| Parameter | Value | Unit | Description |
| :--- | :--- | :--- | :--- |
| `rmax` | 1.2 | $\text{day}^{-1}$ | Maximum specific nutrient uptake rate |
| `kn` | 0.1 | conc. | Half-saturation constant for nutrient uptake |
| `ksLight` | 20.0 | dimensionless | Half-saturation constant for light limitation |
| `gmax` | 0.2 | $\text{day}^{-1}$ | Maximum zooplankton grazing rate |
| `Iv` | 1.8 | $\text{conc.}^{-1}$ | Ivlev grazing saturation coefficient |
| `Ipn` | 0.01 | $\text{day}^{-1}$ | Phytoplankton respiration coefficient |
| `Izn` | 0.01 | $\text{day}^{-1}$ | Zooplankton respiration coefficient |
| `Ipd` | 0.02 | $\text{day}^{-1}$ | Phytoplankton natural mortality rate |
| `Izd` | 0.02 | $\text{day}^{-1}$ | Zooplankton natural mortality rate |
| `Idn` | 0.001 | $\text{day}^{-1}$ | Detritus remineralization rate |

---

## Installation & Requirements

Ensure you have R (version ≥ 3.5.0) installed along with `deSolve`:

```R
install.packages("deSolve")


source("npzd_model.R")

output$total_mass <- output$pp + output$zoo + output$detr + output$nut
# Check mass conservation (should remain constant at initial total sum = 2.0)
summary(output$total_mass)

