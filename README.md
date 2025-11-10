Simple NMR Spin Simulator

A minimal, quantum-mechanical simulator for calculating NMR spectra from spin systems, designed for educational purposes and clarity.

This project implements the core logic for calculating transition frequencies and intensities in a coupled spin system, with a focus on simplicity and ease of tracing the underlying quantum mechanics.

🌟 Features

Customizable Spin Systems: Easily define chemical shifts (ppm_positions) and J-couplings (J_COUPLING_PAIRS).

Spin 1/2 and Spin 1 Support: Includes full support for both spin-1/2 nuclei (like $^1\text{H}$) and spin-1 nuclei (like $^2\text{D}$), which is useful for studying isotopic substitution.

Dynamic Plotting: Generates an interactive Matplotlib plot with selectable peaks and dynamically scaled PPM axis.

Clear Parameter Block: All key simulation settings are organized at the top of the file for quick adjustment.

🚀 Getting Started

Prerequisites

You need Python and the following libraries installed:

pip install numpy matplotlib


Running the Simulator

Download the nmr_simulator.py file.

Run the script directly from your terminal:

python nmr_simulator.py


The script will calculate the Hamiltonian, determine allowed transitions, print the J-coupling matrix and peak list to the console, and display the spectrum plot.

⚙️ Configuration

To modify the simulated system, edit the parameters directly within nmr_simulator.py:

Parameter

Description

spectrometer_1H_MHz

The reference frequency of the magnet (e.g., 600 MHz).

spins

List of spin quantum numbers (e.g., [1/2, 1/2, 1]).

nuclei_types

List of corresponding nucleus types (e.g., ['H', 'H', 'D']).

ppm_positions

Chemical shifts for each nucleus in ppm.

J_COUPLING_PAIRS

List of (index i, index j, J_value) to define couplings in Hz.

PLOT_NUCLEUS

Select which nucleus's spectrum scale to display ('H', 'D', etc.).

⚖️ License (CC BY-NC 4.0)

This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International Public License (CC BY-NC 4.0).

In simple terms, this means:

You are free to use, share, and adapt this code.

You must give appropriate credit.

You may NOT use the material for commercial purposes. This code is intended for personal, educational, and research use only.