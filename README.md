# ⚛️ NMR Spin Simulator (v1.2.0)

Welcome to the High-Resolution NMR Spin Simulator. This repository contains two completely independent implementations of a quantum-mechanical NMR spectrum solver, engineered for both educational clarity and extreme computational performance.

## 📂 Project Overview

The project is split into two independent modules based on performance and design philosophy:

* **`dense/` — The Educational Engine**  
  A single-file, classical Kronecker-product script designed for learning the fundamentals of NMR quantum mechanics and checking theory.
* **`sparse/` — The High-Performance Engine**  
  A modular, heavily optimized solver leveraging block-diagonalization ($M_z$ symmetry sectors) capable of handling **20+ spins** comfortably.

---

## ⚡ Core Comparison

| Feature                 | 🟢 Dense Simulator                   | 🔵 Sparse Simulator                                |
|-------------------------|--------------------------------------|-----------------------------------------------------|
| **Primary Focus**       | Educational clarity & simplicity     | calability & raw performance                        |
| **Architecture**        | Single-file script                   | Clean, modular library                              |
| **Capacity**            | $N \le 6$ spins                      | Comfortably scales to $N \ge 20$ spins              |
| **Hamiltonian Storage** | Dense NumPy Arrays                   | CSR Sparse SciPy Matrices                           |
| **Under the Hood**      | Global dense NumPy matrices          | $M_z$ block-diagonal decomposition                  |
| **Best For**            | Learning NMR quantum math & testing  | Large spin systems & complex spin networks          |

## 🚀 Getting Started

### Prerequisites

Both versions require Python 3 along with standard scientific computation libraries:
```bash
pip install numpy scipy matplotlib pytest
```

### Running the Simulations

1. **To run the Educational (Dense) version:** Navigate into the dense folder and execute the script:

```bash
python dense/nmr_simulator_dense.py
```

2. **To run the High-Performance (Sparse) version:** Navigate into the sparse folder and execute:

```bash
python sparse/nmr_simulator_sparse.py
```

2. **Running Quality Assurance Tests:** The project includes a robust test suite. You can execute all unit and regression tests from the root directory using pytest:

```bash
pytest -v
```

## 📖 Detailed Documentation
For deeper dives into the physics, configuration parameters, and specific usage of each implementation, please check the dedicated README files located within their respective folders


## 🔐 Licene

NMR Simulator is licensed under the MIT License.



