# 🟢 Educational NMR Spin Simulator (Dense Engine)

This directory contains the Educational, Single-File implementation of the NMR spin simulator. It is designed to be highly readable, intuitive, and close to classical textbook quantum mechanics.

## 🌟 Key Features

* **Single-File Simplicity:** The entire program (inputs, math, solver, and plotting) is written in one file (nmr_simulator_dense.py).

* **Classical Kronecker Architecture:** Builds the full, dense Hamiltonian matrix, making it easy to inspect using standard NumPy array printing.

* **Spin 1/2 and Spin 1 Support:** Fully supports mixed spin systems, including spin ($^1\text{H}, ^{13}\text{C}, ^{19}\text{F}, ^{31}\text{P}, ^{29}\text{Si}$) and **spin-1 nuclei** ($^2\text{D}$) nuclei.

* **Numerical Stability Fix:** Implements a frequency scaling mode to eliminate numerical artifacts when Larmor frequencies are close to $0\text{ Hz}$ ($0\text{ ppm}$).

* **Interactive Plotting:** Generates a Matplotlib spectrum plot with selectable peaks and dynamically scaled PPM axis.

## 🚀 Run & Customize

### Run the Script

You need Python and the following scientific libraries installed:
```bash
python nmr_simulator_dense.py
```

## ⚛️ The Default Example System

The `nmr_simulator.py` file is pre-configured to simulate an H-D-D system, often found in $\text{CHD}_2$ isotopomers (like the residual solvent signal in DMSO-d5).

The central proton ($\text{H}$, $I=1/2$) is coupled to two equivalent deuterons ($\text{D}$, $I=1$). This coupling is expected to produce a quintet (five equally spaced lines with an intensity ratio of $1:2:3:2:1$) for the proton signal, which is predicted by the $2nI+1$ rule: $2 \times 2 \times 1 + 1 = 5$.

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
To adapt the simulation for your own molecule, you only need to modify the parameters at the very beginning of the nmr_simulator.py file under the <<< S I M U L A T I O N  P A R A M E T E R S >>> block.

| Parameter | Purpose | How to Modify |
 | ----- | ----- | ----- |
| `spins` | Quantum number for each nucleus in the system. | Change the array, e.g., `np.array([1/2, 1/2, 1/2])`. |
| `nuclei_types` | Corresponding type (used for $\gamma$ ratio and plotting scale). | Ensure this matches `spins`, e.g., `['H', 'H', 'H']`. |
| `ppm_positions` | Chemical shift for each nucleus (must match array length). | Set your desired $\delta$ values, e.g., `[3.5, 1.2, 0.9]`. |
| `J_COUPLING_PAIRS` | Define couplings by index. **The order of indices matters!** | Add or remove tuples: `(Index i, Index j, J_value)`. |
| `PLOT_NUCLEUS` | Sets the reference scale for the X-axis (e.g., `'H'` for $^1\text{H}$ ppm). | Use `'H'`, `'D'`, or `'13C'`. |
| `FREQUENCY_SCALING_MODE` | Enables/disables numerical stability fix for Larmor frequencies near 0 Hz. | Set to True or False.|
| `SCALING_SHIFT_PPM` | The temporary uniform PPM shift applied internally if scaling is enabled. | Keep default (e.g., 1), or change if needed.|


## 🔐 Licene

NMR Simulator is licensed under the MIT License.



