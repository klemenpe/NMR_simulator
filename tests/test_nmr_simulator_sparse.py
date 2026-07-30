import numpy as np
import pytest
from transition_state_solver import calculate_transition_states_direct_blocks, generate_spin_states_by_mz
import functions as functions


# run this tests in terminal (root project folder) using command: pytest -v 

def test_h_dd_raw_transitions_regression():
    """
    Regression test for the 3-spin system (H-D-D: spins 0.5, 1, 1) 
    checking the raw uncombined 21 transitions directly from the solver.
    """
    spins = np.array([0.5, 1.0, 1.0])
    nuclei_types = ['H', 'D', 'D']
    ppm_positions = [2.50, 2.50, 2.50]
    
    # Spectrometer setup (600 MHz 1H)
    spectrometer_1H_MHz = 600.0
    v, _ = functions.ppm_to_hz(ppm_positions, nuclei_types, spectrometer_1H_MHz)
    
    # J-Coupling matrix matching your main script configuration
    # (0=H1, 1=D1, 2=D2) -> J_H1D1=2.0, J_H1D2=2.0, J_D1D2=0.0
    J_pairs = [(0, 1, 2.0), (0, 2, 2.0), (1, 2, 0.0)]
    J_matrix = functions.convert_pairs_to_j_matrix(len(spins), J_pairs, verbose=False)
    
    # Run solver
    signals = calculate_transition_states_direct_blocks(spins, v, J_matrix, cutoff=0.001)
    
    # Extract frequencies and intensities
    frequencies = signals[0]
    intensities = signals[1]
    
    # Sort by frequency for deterministic comparison
    sorted_idx = np.argsort(frequencies)
    sorted_freqs = frequencies[sorted_idx]
    sorted_ints = intensities[sorted_idx]
    
    # Expected raw transitions (all 21 peaks)
    expected_freqs = np.array([
        229.2488, 229.2504, 229.2504, 229.2520, 229.2520, 229.2535,
        231.2488, 231.2504, 231.2504, 231.2520, 231.2520, 231.2535,
        1496.0032, 1498.0016, 1498.0079, 1500.0000, 1500.0032, 
        1500.0095, 1502.0016, 1502.0079, 1504.0031
    ])
    
    expected_ints = np.array([
        3.9937, 1.9968, 5.9905, 1.9969, 5.9906, 3.9937,
        4.0063, 6.0094, 2.0031, 6.0095, 2.0032, 4.0063,
        1.0063, 1.0032, 1.0031, 1.0000, 1.0000, 
        1.0000, 0.9969, 0.9968, 0.9937
    ])
    
    # Assert dimensions match expected peak count (21 transitions)
    assert len(sorted_freqs) == 21, f"Expected 21 transitions, got {len(sorted_freqs)}"
    
    # Check that frequencies and intensities match within a tight tolerance (e.g., 1e-3)
    np.testing.assert_allclose(sorted_freqs, expected_freqs, atol=1e-3, err_msg="Raw transition frequencies drifted!")
    np.testing.assert_allclose(sorted_ints, expected_ints, atol=1e-3, err_msg="Raw transition intensities drifted!")


def test_basis_state_generation():
    """Verify that spin states are correctly generated and grouped by Mz for a 0.5 and 1 system."""
    spins = np.array([0.5, 1.0])
    mz_basis = generate_spin_states_by_mz(spins)
    
    # Total states should be 2 * 3 = 6
    total_states = sum(len(states) for states in mz_basis.values())
    assert total_states == 6
    
    # Check specific Mz groupings
    assert 1.5 in mz_basis
    assert len(mz_basis[1.5]) == 1  # Only [(0.5, 1.0)]
    assert len(mz_basis[0.5]) == 2  # [(0.5, 0.0), (-0.5, 1.0)]


def test_hd_two_spin_raw_transitions_regression():
    """
    Regression test for the 2-spin system (H-D: spins 0.5, 1.0) 
    checking the raw uncombined transitions directly from the solver.
    """
    spins = np.array([0.5, 1.0])
    nuclei_types = ['H', 'D']
    ppm_positions = [2.50, 2.50]
    
    # Spectrometer setup (600 MHz 1H)
    spectrometer_1H_MHz = 600.0
    v, _ = functions.ppm_to_hz(ppm_positions, nuclei_types, spectrometer_1H_MHz)
    
    # J-Coupling matrix configuration for 2 spins: J_HD = 2.0 Hz
    J_pairs = [(0, 1, 2.0)]
    J_matrix = functions.convert_pairs_to_j_matrix(len(spins), J_pairs, verbose=False)
    
    # Run solver
    signals = calculate_transition_states_direct_blocks(spins, v, J_matrix, cutoff=0.001)
    
    # Extract frequencies and intensities
    frequencies = signals[0]
    intensities = signals[1]
    
    # Sort by frequency for deterministic comparison
    sorted_idx = np.argsort(frequencies)
    sorted_freqs = frequencies[sorted_idx]
    sorted_ints = intensities[sorted_idx]
    
    # Expected raw transitions (all 7 peaks for H-D system)
    expected_freqs = np.array([
        229.2504, 229.2520, 231.2504, 231.2520, 
        1498.0016, 1500.0032, 1502.0016
    ])
    
    expected_ints = np.array([
        1.9968, 1.9969, 2.0031, 2.0032, 
        1.0032, 1.0000, 0.9969
    ])
    
    # Assert dimensions match expected peak count (7 transitions)
    assert len(sorted_freqs) == 7, f"Expected 7 transitions, got {len(sorted_freqs)}"
    
    # Check that frequencies and intensities match within tight tolerance
    np.testing.assert_allclose(sorted_freqs, expected_freqs, atol=1e-3, err_msg="2-spin transition frequencies drifted!")
    np.testing.assert_allclose(sorted_ints, expected_ints, atol=1e-3, err_msg="2-spin transition intensities drifted!")
