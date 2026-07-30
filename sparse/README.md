# 🔵 High-Performance NMR Spin Simulator (Sparse Engine) — v1.2.1

This directory houses the high-performance, production-ready simulator engineered to solve massive spin systems containing **20+ active spins** with extreme speed and memory efficiency.

## ⚡ Architectural Optimizations

To overcome the exponential $2^N$ (or $(2I+1)^N$) matrix growth of quantum mechanics, this engine implements two critical optimizations:

### 1. Compressed Sparse Row Storage (SciPy CSR)
Instead of allocating dense matrices filled with zeros, the Hamiltonian is constructed entirely using `scipy.sparse` matrix operators. This keeps the memory footprint lightweight, scaling with non-zero couplings rather than matrix area.

### 2. $M_z$ Symmetry Block-Diagonalization
Because the total magnetic projection operator $F_z$ commutes with the scalar Hamiltonian ($[H, F_z] = 0$), quantum states belonging to different $M_z$ sectors cannot couple.

Instead of diagonalizing one giant matrix of size $D \times D$ the program:
1. **Groups** quantum basis states into independent $M_z$ subspace blocks.
2. **Extracts and solves** small independent block matrices in parallel using targeted eigensolvers.
3. **Applies Selection Rules** only between adjacent $M_z$ blocks ($\Delta M_z = \pm 1$).

This reduces the diagonalization time from an intractable $O(D^3)$ to a series of lightweight $O(d^3)$ solves, turning calculations that would take days into milliseconds!

---

## 🔬 Rigorous Transition Mechanics ($F_x$ Operator)

The sparse engine implements the fully rigorous physical transition operator $F_x = \sum_i I_{x,i}$ for calculating peak line intensities:

* **Physical State Mixing:** When strong J-coupling mixes quantum states, transition intensities are evaluated through explicit operator sandwich products $\langle \psi_b | F_x | \psi_a \rangle$.

* **Heteronuclear Scaling:** Transition amplitudes automatically incorporate physical spin projection weights ($\frac{1}{2}$ for spin-1/2 vs. $\frac{\sqrt{2}}{2}$ for spin-1), yielding exact quantum-mechanical relative line intensities.

---

## 📂 Modular File Architecture

The solver is structured into dedicated, decoupled modules:

* **`nmr_simulator_sparse.py`** — Main entry point. Configure your chemical shifts, coupling constants ($J$), spectrometer frequencies, and output preferences here.
* **`transition_state_solver.py`** — Block-diagonalization solver that sorts basis states by $M_z$, diagonalizes individual subspaces, and calculates allowed transition frequencies and intensities.
* **`spin_database.py`** — Central repository containing gyromagnetic ratios ($\gamma$)
* **`functions.py`** — Utility routines for unit conversions (PPM to Hz), matrix formatting, and plotting.

---


## 🚀 Running the Simulator & Testing

Configure your target system inside `nmr_simulator_sparse.py` and execute:

```bash
python nmr_simulator_sparse.py
```

### Verification & Regression Tests

To run validation tests verifying this solver:


```bash
pytest -v
```

## 🔐 Licene

NMR Simulator is licensed under the MIT License.