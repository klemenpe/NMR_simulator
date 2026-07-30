# 🟢 Educational NMR Spin Simulator (Dense Engine)

This directory contains the single-file, educational implementation of the NMR spin simulator. It is designed to be highly readable, intuitive, and closely aligned with classical textbook quantum mechanics.

## 🌟 Key Features

* **Single-File Simplicity:** The entire program (inputs, matrix construction, solver, and plotting) lives in a single, transparent script (`nmr_simulator_dense.py`).

* **Classical Kronecker Architecture:** Builds full, dense Hamiltonian matrices using direct Kronecker products (`np.kron`), making it trivial to inspect matrices using standard NumPy array printing.

* **Mixed Spin Systems ($I = 1/2$ and $I = 1$):** Full quantum support for spin-1/2 ($^1\text{H}, ^{13}\text{C}, ^{19}\text{F}, ^{31}\text{P}$) and quadrupolar spin-1 ($^2\text{D}$) nuclei.

* **Numerical Stability Scaling:** Implements a temporary frequency-scaling shift to eliminate numerical precision artifacts when Larmor frequencies approach $0\text{ Hz}$ ($0\text{ ppm}$).

* **Interactive Visualization:** Generates a Matplotlib spectrum plot with interactive peak readout and dynamically scaled PPM axes.

---

## 🚀 Run & Customize

### Run the Script

Ensure you have Python and the required scientific libraries installed:

```bash
pip install numpy scipy matplotlib
```

Run the simulator directly:

```bash
python nmr_simulator_dense.py
```

## ⚛️ Default Verification System ($\text{CHD}_2$)

By default, nmr_simulator_dense.py is configured to simulate a residual proton coupled to two equivalent deuterium nuclei ($\text{CHD}_2$ isotopomer, as seen in DMSO-$d_5$).

The central proton ($^1\text{H}$, $I=1/2$) is coupled to two equivalent deuterons ($^2\text{D}$, $I=1$). According to the quantum $2nI + 1$ rule, this coupling produces a symmetric 1:2:3:2:1 quintet:

$$\text{Lines} = 2 \times (2) \times (1) + 1 = 5$$


**System Configuration**

| Nucleus | Spin | Type | Shift (ppm) |
 | ----- | ----- | ----- | ----- |
| **0** | $1/2$ | $^1\text{H}$ | 2.50 |
| **1** | $1/2$ | $^2\text{D}$ | 2.50 |
| **2** | $1$ | $^2\text{D}$ | 2.50 |

**Couplings:**

Couplings:

$J_{0,1}$ (H-D coupling) = $2.0\text{ Hz}$

$J_{0,2}$ (H-D coupling) = $2.0\text{ Hz}$

$J_{1,2}$ (D-D coupling) = $0.0\text{ Hz}$

### ⚙️ Customizing the Simulation

To adapt the simulation for your own molecule, edit the parameters at the top of nmr_simulator_dense.py under the <<< SIMULATION PARAMETERS >>> block:

| Parameter | Purpose | How to Modify |
 | ----- | ----- | ----- |
| `spins` | Spin quantum numbers ($I$) for each nucleus. | Change the array, e.g., `np.array([1/2, 1/2, 1/2])`. |
| `nuclei_types` | Isotope labels matching spins (for $\gamma$ scaling and plotting). | Ensure this matches `spins`, e.g., `['H', 'H', 'H']`. |
| `ppm_positions` | Chemical shift ($\delta$) for each nucleus (must match array length). | Set your desired $\delta$ values, e.g., `[3.5, 1.2, 0.9]`. |
| `J_COUPLING_PAIRS` | Define couplings by index. **The order of indices matters!** | Add or remove tuples: `(Index i, Index j, J_value)`. |
| `PLOT_NUCLEUS` | Defines the reference nucleus for the spectral X-axis. | Use `'H'`, `'D'`, or `'13C'`. |
| `FREQUENCY_SCALING_MODE` | Enables/disables numerical stability fix for Larmor frequencies near 0 Hz. | Set to True or False.|
| `SCALING_SHIFT_PPM` | The temporary uniform PPM shift applied internally if scaling is enabled. | Keep default (e.g., 1), or change if needed.|


## 🔐 Licene

NMR Simulator is licensed under the MIT License.



