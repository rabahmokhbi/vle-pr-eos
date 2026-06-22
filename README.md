# VLE Calculator — Peng-Robinson Equation of State

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Author](https://img.shields.io/badge/Author-Elarix%20Dreson-purple)](https://linktr.ee/elarixdreson)

A thermodynamic VLE (Vapor-Liquid Equilibrium) calculator built on the **Peng-Robinson Equation of State (1976)** — one of the most widely used cubic EOS in the petroleum, natural gas, and chemical process industries.

---

## Features

| Feature | Description |
|---|---|
| **Pure component P-v-T** | Compressibility factor Z, molar volume, fugacity coefficient |
| **Saturation pressure** | Fugacity-equality iterative method |
| **Phase envelope** | Vapor pressure curve from Tr ≈ 0.45 to Tc |
| **Binary bubble point** | Successive-substitution algorithm |
| **P-x-y diagram** | Full bubble & dew curve at fixed T |
| **Component database** | 10 pre-loaded components (hydrocarbons, water, CO₂, N₂) |

---

## Theory

The **Peng-Robinson EOS**:

```
P = RT/(V−b) − a(T)/[V(V+b) + b(V−b)]
```

Where:
- `a(T) = 0.45724 · R²Tc²/Pc · α(T)`
- `α(T) = [1 + κ(1 − √Tr)]²`
- `κ = 0.37464 + 1.54226ω − 0.26992ω²`
- `b = 0.07780 · RTc/Pc`

**Binary mixing rules (van der Waals):**
```
a_mix = ΣΣ xi·xj·aij    where aij = √(ai·aj)·(1 − kij)
b_mix = Σ xi·bi
```

---

## Installation

```bash
git clone https://github.com/<your-username>/vle-pr-eos.git
cd vle-pr-eos
pip install numpy scipy matplotlib
```

---

## Usage

### Run the full demo

```bash
python vle_calculator.py
```

This will:
1. Compute pure-component properties for **propane** and **water**
2. Generate phase envelope plots (saved as `.png`)
3. Compute bubble-point pressure for **ethane/propane** at 300 K
4. Plot P-x-y diagrams for two binary systems

### Use in your own script

```python
from vle_calculator import PengRobinson, BinaryVLE, COMPONENTS, plot_pxy
import numpy as np

# Pure component
co2 = PengRobinson("CO2", **COMPONENTS["CO2"])
Psat = co2.saturation_pressure(280.0)
print(f"CO2 Psat at 280 K = {Psat/1e6:.3f} MPa")

# Binary VLE
methane = PengRobinson("methane", **COMPONENTS["methane"])
ethane  = PengRobinson("ethane",  **COMPONENTS["ethane"])
mix = BinaryVLE(methane, ethane, kij=0.0)

z = np.array([0.4, 0.6])
P_bub, y = mix.bubble_point_pressure(200.0, z)
print(f"Bubble P = {P_bub/1e6:.4f} MPa, y_methane = {y[0]:.4f}")

# P-x-y diagram
plot_pxy(mix, T=200.0, save=True)
```

---

## Component Database

| Component | Tc (K) | Pc (MPa) | ω |
|---|---|---|---|
| Methane | 190.6 | 4.60 | 0.011 |
| Ethane | 305.3 | 4.87 | 0.099 |
| Propane | 369.8 | 4.25 | 0.153 |
| n-Butane | 425.1 | 3.80 | 0.200 |
| n-Pentane | 469.7 | 3.37 | 0.251 |
| Benzene | 562.2 | 4.89 | 0.212 |
| Toluene | 591.8 | 4.11 | 0.263 |
| Water | 647.1 | 22.10 | 0.345 |
| CO₂ | 304.2 | 7.38 | 0.225 |
| Nitrogen | 126.2 | 3.39 | 0.037 |

---

## Sample Outputs

### Phase Envelope — Propane
> Vapor pressure curve from ~170 K to Tc = 369.8 K

### P-x-y — Ethane/Propane at 300 K
> Bubble and dew curves showing VLE behavior of the binary system

---

## Validation

| System | Calculated | Literature | Error |
|---|---|---|---|
| Propane Psat (300 K) | ~0.998 MPa | 0.997 MPa | <0.1% |
| Water Psat (373.15 K) | ~1.01 bar | 1.013 bar | <0.3% |

---

## Future Improvements

- [ ] SRK equation of state (Soave-Redlich-Kwong)
- [ ] Ternary and multicomponent flash (Rachford-Rice)
- [ ] T-x-y diagram generation
- [ ] Export results to Excel/CSV
- [ ] Streamlit web interface

---

## Author

**Rabah Mokhbi** — [Elarix Dreson](https://linktr.ee/elarixdreson)  
Master's Student — Pharmaceutical & Process Engineering, Université de Blida 1, Algeria  
[![ORCID](https://img.shields.io/badge/ORCID-0009--0004--8338--1466-brightgreen?logo=orcid)](https://orcid.org/0009-0004-8338-1466)  
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-blue?logo=linkedin)](https://linkedin.com/in/rabah-mokhbi)

---

## License

MIT © 2026 Rabah Mokhbi
```

---

## References

1. Peng, D.-Y., & Robinson, D. B. (1976). A new two-constant equation of state. *Industrial & Engineering Chemistry Fundamentals*, 15(1), 59–64.
2. Smith, J. M., Van Ness, H. C., & Abbott, M. M. (2005). *Introduction to Chemical Engineering Thermodynamics* (7th ed.). McGraw-Hill.
3. Sandler, S. I. (2006). *Chemical, Biochemical, and Engineering Thermodynamics* (4th ed.). Wiley.
