# ⚛️ NMR Spin Simulator (v1.1.0)

Welcome to the High-Resolution NMR Spin Simulator. This repository contains two completely independent implementations of a quantum-mechanical NMR spectrum solver, tailored for different use cases.

## 📂 Repository Structure

The code is organized into two distinct directories based on performance requirements and educational focus:

```
nmr-spin-simulator/
│
├── README.md                  
├── LICENSE.md                    # Project License
│
├── dense/                        # 🟢 THE EDUCATIONAL VERSION
│   ├── README.md                 # Setup and usage for the dense code
│   └── nmr_simulator_dense.py    # Single-file, classical Kronecker-based code
│
└── sparse/                       # 🔵 THE HIGH-PERFORMANCE VERSION
    ├── README.md                 # Detailed documentation on block-diagonalization
    ├── nmr_simulator_sparse.py   # Main runner & input configuration file
    ├── hamiltonian.py            # Sparse matrix construction engine
    ├── transition_state_solver.py  # Block-diagonal eigensolver & transition engine
    ├── spin_database.py          # Central database for spins and operators
    └── functions.py              # Helper/utility functions (PPM to Hz conversion)
```

## ⚡ Comparison: Dense vs. Sparse

| Feature             | 🟢 Dense Simulator                   | 🔵 Sparse Simulator                                |
|---------------------|--------------------------------------|-----------------------------------------------------|
| Primary Focus       | Educational clarity & simplicity     | calability & raw performance                        |
| Architecture        | Single-file script                   | Clean, modular library                              |
| Max Practical Spins | $N \le 6$ spins                      | $N = 15+$ spins                                     |
| Hamiltonian Storage | Dense NumPy Arrays                   | CSR Sparse SciPy Matrices                           |
| Diagonalization     | Solves $2^N \times 2^N$ matrix       | Solves small independent $M_z$ blocks               |
| Best For            | Learning basic NMR quantum mechanics | Large systems, metabolite databases, stress-testing |

## 🚀 Getting Started

### Prerequisites

Both versions require Python 3 and standard scientific libraries:
```bash
pip install numpy scipy matplotlib
```

### Quick Start

1. For an easy, visual introduction to a simple coupled system: Navigate to the dense/ folder and run:

```bash
python nmr_simulator_dense.py
```

2. For simulating massive systems with modular files: Navigate to the sparse/ folder and run:

```bash
python nmr_simulator_sparse.py
```


## 🔐 Licene

NMR Simulator is licensed under the MIT License.



