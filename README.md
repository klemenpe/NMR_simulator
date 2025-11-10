# ⚛️ NMR Spin Simulator (v0.1.0)

A minimal, quantum-mechanical simulator for calculating high-resolution NMR spectra. This project is focused on educational clarity and providing a straightforward, functional implementation of the spin simulation logic.

## 🌟 Key Features

* **Customizable Systems:** Easily define chemical shifts and J-couplings for any arbitrary spin system.

* **Spin 1/2 and Spin 1 Support:** Includes full support for both spin-1/2 nuclei ($^1\text{H}, ^{13}\text{C}$) and **spin-1 nuclei** ($^2\text{D}$), which is vital for studying isotopic substitution effects.

* **Interactive Plotting:** Generates a Matplotlib spectrum plot with selectable peaks and dynamically scaled PPM axis.

* **Clear Structure:** All core simulation settings are conveniently located in a single parameter block at the top of the Python file.

## 🚀 Setup & Run

### Prerequisites

You need Python and the following scientific libraries installed:
```bash
pip install numpy matplotlib
```

### Execution

1. Download the `nmr_simulator.py` file to your local computer.

2. Run the script directly from your terminal (or via any IDE/editor that supports Python execution):
```bash
python nmr_simulator.py
```
The script will calculate the J-coupling matrix, solve the Hamiltonian, and display the final spectrum plot and the list of calculated transitions in your console.

## ⚛️ The Default Example System

The `nmr_simulator.py` file is pre-configured to simulate an H-D-D system, often found in $\text{CHD}_2$ isotopomers (like the residual solvent signal in DMSO-d5).

The central proton ($\text{H}$, $I=1/2$) is coupled to two equivalent deuterons ($\text{D}$, $I=1$). This coupling is expected to produce a quintet (five equally spaced lines) for the proton signal, which is predicted by the $2nI+1$ rule: $2 \times 2 \times 1 + 1 = 5$.

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

## ⚙️ Customizing the Simulation
To adapt the simulation for your own molecule, you only need to modify the parameters at the very beginning of the nmr_simulator.py file under the <<< S I M U L A T I O N  P A R A M E T E R S >>> block.

| Parameter | Purpose | How to Modify |
 | ----- | ----- | ----- |
| `spins` | Quantum number for each nucleus in the system. | Change the array, e.g., `np.array([1/2, 1/2, 1/2])`. |
| `nuclei_types` | Corresponding type (used for $\gamma$ ratio and plotting scale). | Ensure this matches `spins`, e.g., `['H', 'H', 'H']`. |
| `ppm_positions` | Chemical shift for each nucleus (must match array length). | Set your desired $\delta$ values, e.g., `[3.5, 1.2, 0.9]`. |
| `J_COUPLING_PAIRS` | Define couplings by index. **The order of indices matters!** | Add or remove tuples: `(Index i, Index j, J_value)`. |
| `PLOT_NUCLEUS` | Sets the reference scale for the X-axis (e.g., `'H'` for $^1\text{H}$ ppm). | Use `'H'`, `'D'`, or `'13C'`. |

## 🔐 Licene

NMR Simulator is licensed under the MIT License.



