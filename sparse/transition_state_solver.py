import numpy as np
from scipy import sparse
import itertools
from collections import defaultdict
from spin_database import get_spin_operators, SPIN_CONFIG


def generate_spin_states(spins):
    """
    Generates all possible combinations of mI values for the spin system.
    Dynamically generates the correct spin projection values using the formula: s, s-1, ..., -s.
    Validates that each spin is supported in the operator database.
    """
    all_m_values = []
    for s in spins:
        # Safety Check: ensure the spin is defined in our database
        if s not in SPIN_CONFIG:
            raise ValueError(
                f"Unsupported spin value: {s}. To use this spin, please add its matrices "
                f"(x, y, z, unit) to SPIN_CONFIG in spin_database.py first."
            )
        
        # Dynamically generate the descending list of mI projections (e.g. 1.5 -> [1.5, 0.5, -0.5, -1.5])
        m_vals = list(np.arange(s, -s - 0.1, -1.0))
        all_m_values.append(m_vals)
        
    return list(itertools.product(*all_m_values))


def calculate_transition_states(H_sparse, spins, cutoff=0.001):
    """
    Solves the eigensystem using Mz Block-Diagonalization and calculates 
    allowed transitions efficiently. Designed to scale up to N=15+ spins.
    """
    print("\n" + "="*50)
    print(">>> STARTING BLOCK-DIAGONAL TRANSITION SOLVER <<<")
    print("="*50)
    
    states = generate_spin_states(spins)
    dim = H_sparse.shape[0]
    
    # Group global state indices by their total Mz value
    mz_groups = defaultdict(list)
    for idx, state in enumerate(states):
        mz_groups[sum(state)].append(idx)
        
    sorted_mzs = sorted(mz_groups.keys(), reverse=True)
    
    block_eigenvalues = {}
    block_eigenvectors = {}
    
    print("Diagonalizing Hamiltonian block-by-block...")
    for mz in sorted_mzs:
        indices = mz_groups[mz]
        # Extract the small submatrix block as a dense array for eigh
        H_block = H_sparse[indices, :][:, indices].toarray()
        
        # This is extremely fast because each block is a fraction of the total size
        E_block, V_block = np.linalg.eigh(H_block)
        
        block_eigenvalues[mz] = E_block
        block_eigenvectors[mz] = V_block
        print(f" - Block Mz = {mz:+.1f}: Size {len(indices)}x{len(indices)} solved.")

    print("Constructing sparse transition operator Fx using Kronecker products...")
    Fx = sparse.csr_matrix((dim, dim), dtype=np.complex128)
    for n in range(len(spins)):
        op_full = 1
        for k in range(len(spins)):
            ops = get_spin_operators(spins[k])
            if k == n:
                op_full = sparse.kron(op_full, ops['x']) if not isinstance(op_full, int) else ops['x']
            else:
                op_full = sparse.kron(op_full, ops['unit']) if not isinstance(op_full, int) else ops['unit']
        Fx += op_full

    print("Calculating transitions using selection rules (Delta_Mz = +/- 1)...")
    peak_frequencies = []
    peak_intensities = []
    
    for i in range(len(sorted_mzs) - 1):
        mz_A = sorted_mzs[i]
        mz_B = sorted_mzs[i+1]
        
        # Check if blocks are adjacent (Delta_Mz = 1.0)
        if np.isclose(abs(mz_A - mz_B), 1.0):
            indices_A = mz_groups[mz_A]
            indices_B = mz_groups[mz_B]
            
            # Get the block of Fx that connects Group A and Group B
            Fx_block = Fx[indices_A, :][:, indices_B].toarray()
            
            V_A = block_eigenvectors[mz_A]
            V_B = block_eigenvectors[mz_B]
            
            # Transform transition block to the eigenbasis: V_A^H * Fx_block * V_B
            T_transformed = V_A.conj().T @ Fx_block @ V_B
            
            # Physical intensity is proportional to the square of the transition matrix element
            # We scale by 4.0 to align standard spin-1/2 intensities with your old scale
            intensities = np.square(np.abs(T_transformed)) * 4.0
            
            E_A = block_eigenvalues[mz_A]
            E_B = block_eigenvalues[mz_B]
            
            # Compute differences between all eigenvalues of block A and block B
            for r in range(len(E_A)):
                for c in range(len(E_B)):
                    intensity = intensities[r, c]
                    if intensity >= cutoff:
                        freq = abs(E_A[r] - E_B[c])
                        peak_frequencies.append(freq)
                        peak_intensities.append(intensity)
                        
    if len(peak_frequencies) > 0:
        koncni_signali = np.vstack([peak_frequencies, peak_intensities])
    else:
        koncni_signali = np.empty((2, 0))
        
    print(f"Calculation complete. Found {koncni_signali.shape[1]} allowed transitions.")
    print("="*50 + "\n")
    
    return koncni_signali