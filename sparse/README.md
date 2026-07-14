# 🔵 High-Performance NMR Spin Simulator (Sparse Engine)

This directory houses the production-ready, high-performance simulator designed to solve massive spin systems containing up to $15+$ active spins.

## ⚡ Architectural Optimizations

To overcome the exponential memory scaling of quantum mechanics ($2^N$), this engine implements two critical software optimizations:

### 1. Compressed Sparse Storage (SciPy CSR)

Instead of storing millions of useless zero elements in RAM, the Hamiltonian is constructed using ```scipy.sparse``` matrices. This reduces the memory footprint from gigabytes to kilobytes.

### 2. $M_z$ Block-Diagonalization

Because the total magnetic projection $M_z$ commutes with the Hamiltonian ($[H, F_z] = 0$), states in different $M_z$ sectors cannot couple.

Instead of diagonalizing one massive matrix of size $D \times D$, the solver:

 1. goups your quantum states into independent $M_z$ subspaces.

 2. Extracts and solves tiny block matrices individually.

 3. Calculates transitions only between adjacent blocks ($\Delta M_z = \pm 1.0$).

This reduces the diagonalization complexity from $O(D^3)$ to a series of tiny $O(d^3)$ solves, speeding up calculations from hours to milliseconds!

## ⚛️ Default Verification System: $\text{CHD}_2$

By default, main_experiment.py is configured with the exact same test system as the dense educational version: a residual proton coupled to two equivalent deuterium nuclei ($\text{CHD}_2$ isotopomer of DMSO).

Spins: $[1/2, 1, 1]$

Expected Result: A classic $1:2:3:2:1$ quintet splitting pattern on the Proton ($\text{H}$) channel, governed by the $2nI + 1$ quantum rule:


$$\text{Lines} = 2 \times (2) \times (1) + 1 = 5$$

Using this identical system allows you to easily cross-verify the accuracy of both engines.

## 🔬 Note on Microscopic Numerical Differences
When simulating mixed-spin systems (such as $H-D$ couplings), you will observe differences in the calculated peak intensities in the $3\text{rd}$ decimal place. This is a feature, not a bug! The Dense engine utilizes a simplified topological transition matrix (treating all allowed flips as having a binary weight of $1.0$), while the Sparse engine implements the mathematically rigorous, physical $F_x = \sum I_{x,i}$ operator. Because J-coupling mixes states, the physical spin-projection weights ($\frac{1}{2}$ for spin-$1/2$ vs $\frac{\sqrt{2}}{2}$ for spin-$1$) correctly scale the transition elements, yielding a more physically authentic spectrum.

## 📂 Modular File Architecture

Unlike the dense version, this engine is fully modular:

- nmr_simulator_sparse.py: Your main program file. This is where you configure your chemical shifts, couplings, and run the simulation.

- hamiltonian.py: The sparse operator compiler that constructs your Zeeman and J-coupling terms.

- transition_state_solver.py: The block-diagonal solver that slices the Hamiltonian, diagonalizes the blocks, and applies quantum selection rules to find spectrum peaks.

- spin_database.py: The unified database of physical gyromagnetic ratios ($\gamma$) and quantum projection matrices ($x, y, z,$ and Identity) for spin-$1/2$ and spin-$1$ systems.

- functions.py: Centralized helper tools (PPM to Hz scaling, pairing matrix generation).


## 🚀 Running the Simulator

Configure your molecules (such as Alanine, Proline, or ATP) inside main_experiment.py under the configuration block, then execute:

```bash
python nmr_simulator_sparse.py
```





