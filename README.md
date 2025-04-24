# PMMA System Generator

> **Accelerated Kinetic Monte Carlo (KMC) workflow for constructing polymethyl methacrylate (PMMA) ensembles with user‑defined chain‑length distributions, dispersities, and bimodal molecular‑weight distributions (MWDs).**

---

## 1  Overview
This repository implements an accelerated KMC engine integrated with control routines to rapidly generate realistic initial configurations of PMMA for subsequent Molecular Dynamics (MD) simulations. It enables you to

* prescribe **target chain lengths** (DP = 10 – 100),
* tune **dispersity** (Ð ≈ 1.30 – 1.90), and
* create **bimodal MWDs**.

All parameters are exposed through descriptive macros so that new polymer architectures can be produced with a two‑line edit and one compile‑run cycle.

---

## 2  Key Features
| Capability | How it works | Typical settings |
|------------|-------------|------------------|
| **Chain‑length selection** | Fixed‐DP growth schedule; each chain stops once its target degree of polymerization is reached. | `DP_TARGETS = {10, 20, 30, 40, 50, 60, 80, 100}` |
| **Dispersity control** | Incremental monomer dosing using `ADD_TIMES` × `ADD_INTERVAL` pairs to broaden or narrow the MWD. | <br>Ð ≈ 1.30 → `(2, 2000)`<br>Ð ≈ 1.50 → `(4, 1000)`<br>Ð ≈ 1.70 → `(8, 500)`<br>Ð ≈ 1.90 → `(20, 200)` |
| **Bimodal MWD** | Alternating feed strategy—equal mass added in four pulses spanning early to late conversion. | `ADD_TIMES = {2, 2, 2, 2}` and `ADD_INTERVAL = {500, 1000, 1500, 2500}` |

> **Note** All macro definitions reside in **`include/add_schedule.h`** and are read at compile time.

---

## 3  Project Structure
```
PMMA-System-Generator/
├─ include/
│  ├─ add_schedule.h   # dosing & growth controls
│  └─ *.h              # core kinetics & utilities
├─ src/
│  ├─ main.cpp         # entry point (KMC loop)
│  └─ ...
├─ output/           # output file
└─ README.md
```

---

## 4  Getting Started
### 4.1  Prerequisites
* **C++17** compiler (GCC ≥ 9, Clang ≥ 10, or MSVC ≥ 19.28)
* **CMake** ≥ 3.20
* **Boost** (ODEInt) for optional stiff‑ODE coupling
* *Optional*: **Python 3.9+** with *NumPy*, *pandas*, *matplotlib* for post‑processing

### 4.2  Build & Run
```bash
# Configure and build (Release recommended)
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Run with default schedule (DP set, Ð ≈ 1.50)
./build/pmma_kmc
```
Edit **`include/add_schedule.h`** and re‑compile to target a different dispersity or to activate the bimodal recipe.

---

## 5  Output Files
| File | Description |
|------|-------------|
| `conversion.out` | Monomer conversion vs. time |
| `molecular_weight.out` | Number‑average (Mn), weight‑average (Mw) vs. time |
| `all_chain.out` | Final chain‑length list (one DP per number) |
| `chain_length_*.out` | Chain‑length output at specific Dpn values |
| `all_chain_conv_*.out` | Chain‑length output at specific conversion values |

These outputs can be fed directly into **Polyply 1.0** to construct CG coordinates for MD.

---

## 6  Customising Dispersity & MWD
Open **`include/add_schedule.h`** and locate the following block:
```cpp
// --- Growth schedule controls ---
#define ADD_TIMES  { 2, 4, 8, 20 }
#define ADD_INTERVAL { 2000, 1000, 500, 200 }
```
Replace the brace‑enclosed lists with any of the presets in **Table 2** or design your own. Re‑build and re‑run to observe the new distribution.

---

## 7  Citing This Work
If you use this code in academic research, please cite:

