import numpy as np
import itertools
from collections import defaultdict


def generate_spin_states_by_mz(spins):
    """
    Generates all possible spin projection combinations and groups them by total Mz.
    Handles arbitrary spin values dynamically.
    """
    all_m_values = []
    for s in spins:
        # Dynamically generate descending mI projections: e.g. 1.0 -> [1.0, 0.0, -1.0]
        m_vals = list(np.arange(s, -s - 0.1, -1.0))
        all_m_values.append(m_vals)
        
    mz_basis = defaultdict(list)
    
    for state in itertools.product(*all_m_values):
        total_mz = round(sum(state), 2)  # Rounded to prevent float precision issues
        mz_basis[total_mz].append(state)
        
    return mz_basis


def calculate_transition_states_direct_blocks(spins, v, J_matrix, cutoff=0.001):
    """
    Constructs and solves Hamiltonian blocks directly in memory without ever
    allocating a global matrix. Handles arbitrary spin-1/2 and spin-1 mixtures.
    """
    print("\n" + "="*60)
    print(">>> STARTING HIGH-PERFORMANCE DIRECT-BLOCK SOLVER <<<")
    print("="*60)
    
    mz_basis = generate_spin_states_by_mz(spins)
    sorted_mzs = sorted(mz_basis.keys(), reverse=True)
    
    block_eigenvalues = {}
    block_eigenvectors = {}
    nspins = len(spins)
    
    print("Building and diagonalizing Hamiltonian blocks directly...")
    for mz in sorted_mzs:
        states_in_block = mz_basis[mz]
        block_size = len(states_in_block)
        
        # Fast state lookup: state tuple -> local index in this dense block
        state_to_local_idx = {state: idx for idx, state in enumerate(states_in_block)}
        
        # Initialize small dense local block matrix
        H_block = np.zeros((block_size, block_size), dtype=np.complex128)
        
        for local_i, state in enumerate(states_in_block):
            # A. Zeeman Term (Diagonal)
            zeeman_energy = sum(v[n] * state[n] for n in range(nspins))
            H_block[local_i, local_i] += zeeman_energy
            
            # B. J-Coupling Terms
            for n in range(nspins):
                for m in range(n + 1, nspins):
                    J = J_matrix[n, m]
                    if J == 0:
                        continue
                    
                    # Sz * Sz term (Diagonal)
                    H_block[local_i, local_i] += J * state[n] * state[m]
                    
                    # Sx*Sx + Sy*Sy term (Off-Diagonal: Flip-Flop via S+S- / S-S+)
                    Sn, Sm = spins[n], spins[m]
                    mn, mm = state[n], state[m]
                    
                    # Case 1: Raise spin n, Lower spin m (S+_n S-_m)
                    if mn < Sn and mm > -Sm:
                        # Construct flipped state
                        flipped = list(state)
                        flipped[n] += 1.0
                        flipped[m] -= 1.0
                        flipped_tuple = tuple(flipped)
                        
                        if flipped_tuple in state_to_local_idx:
                            local_j = state_to_local_idx[flipped_tuple]
                            # Quantum mechanics coefficient: sqrt(S(S+1) - m(m +/- 1))
                            C_raise = np.sqrt(Sn * (Sn + 1.0) - mn * (mn + 1.0))
                            C_lower = np.sqrt(Sm * (Sm + 1.0) - mm * (mm - 1.0))
                            
                            H_block[local_i, local_j] += 0.5 * J * C_raise * C_lower
                    
                    # Case 2: Lower spin n, Raise spin m (S-_n S+_m)
                    if mn > -Sn and mm < Sm:
                        # Construct flipped state
                        flipped = list(state)
                        flipped[n] -= 1.0
                        flipped[m] += 1.0
                        flipped_tuple = tuple(flipped)
                        
                        if flipped_tuple in state_to_local_idx:
                            local_j = state_to_local_idx[flipped_tuple]
                            C_lower = np.sqrt(Sn * (Sn + 1.0) - mn * (mn - 1.0))
                            C_raise = np.sqrt(Sm * (Sm + 1.0) - mm * (mm + 1.0))
                            
                            H_block[local_i, local_j] += 0.5 * J * C_lower * C_raise
        
        # solving eigenvalues
        E_block, V_block = np.linalg.eigh(H_block)
        
        block_eigenvalues[mz] = E_block
        block_eigenvectors[mz] = V_block
        print(f" - Block Mz = {mz:+.1f}: Size {block_size}x{block_size} solved.")

    print("\nCalculating transitions block-by-block...")
    peak_frequencies = []
    peak_intensities = []
    
    # Constructing the Transition Operator
    for i in range(len(sorted_mzs) - 1):
        mz_A = sorted_mzs[i]
        mz_B = sorted_mzs[i+1]
        
        # Verify selection rule (Delta_Mz = +/- 1)
        if np.isclose(abs(mz_A - mz_B), 1.0):
            states_A = mz_basis[mz_A]
            states_B = mz_basis[mz_B]
            
            # Construct a dense local Fx block connecting Block A and Block B
            # Fx_block[r, c] = <A_r | Fx | B_c>
            Fx_block = np.zeros((len(states_A), len(states_B)), dtype=np.complex128)
            
            for r, state_A in enumerate(states_A):
                for c, state_B in enumerate(states_B):
                    # Fx is the sum of Ix_n = 0.5 * (S+_n + S-_n)
                    # For state_A and state_B to be connected by Fx, 
                    # exactly one spin must change by +/- 1.0, and others must be identical.
                    diff_count = 0
                    connected = True
                    active_spin = -1
                    
                    for n in range(nspins):
                        diff = state_A[n] - state_B[n]
                        if diff != 0:
                            diff_count += 1
                            if abs(diff) == 1.0:
                                active_spin = n
                            else:
                                connected = False
                                break
                    
                    if connected and diff_count == 1:
                        # Calculate raising/lowering coefficient
                        S_act = spins[active_spin]
                        m_B = state_B[active_spin]
                        m_A = state_A[active_spin]
                        
                        if m_A > m_B: # Raising from B to A
                            coeff = np.sqrt(S_act * (S_act + 1.0) - m_B * (m_B + 1.0))
                        else: # Lowering from B to A
                            coeff = np.sqrt(S_act * (S_act + 1.0) - m_B * (m_B - 1.0))
                        
                        Fx_block[r, c] = 0.5 * coeff
            
            V_A = block_eigenvectors[mz_A]
            V_B = block_eigenvectors[mz_B]
            T_transformed = V_A.conj().T @ Fx_block @ V_B
            
            intensities = np.square(np.abs(T_transformed)) * 4.0
            E_A = block_eigenvalues[mz_A]
            E_B = block_eigenvalues[mz_B]
            
            for r in range(len(E_A)):
                for c in range(len(E_B)):
                    intensity = intensities[r, c]
                    if intensity >= cutoff:
                        freq = abs(E_A[r] - E_B[c])
                        peak_frequencies.append(freq)
                        peak_intensities.append(intensity)
                        
    if len(peak_frequencies) > 0:
        spectra_peaks = np.vstack([peak_frequencies, peak_intensities])
    else:
        spectra_peaks = np.empty((2, 0))
        
    print(f"Calculation complete. Found {spectra_peaks.shape[1]} allowed transitions.")
    print("="*60 + "\n")
    
    return spectra_peaks
