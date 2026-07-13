import numpy as np


def convert_pairs_to_j_matrix(nspins, coupling_pairs, verbose=True):
    """
    Converts a list of (i, j, value) tuples into an N x N J_matrix.
    
    Args:
        nspins: Total number of spins in the system.
        coupling_pairs: List of (i, j, J_value) tuples.
        verbose: If True, prints the matrix to console for verification.
    """
    if verbose:
        print("-" * 30)
        print(">>> J-COUPLING MATRIX CONSTRUCTION <<<")

    J = np.zeros((nspins, nspins))
    
    for i, j, j_val in coupling_pairs:
        if i >= nspins or j >= nspins:
            raise IndexError(f"Coupling pair ({i}, {j}) exceeds spin count {nspins}.")
        if i == j:
            print(f"Warning: Ignoring J-coupling for atom to itself at index {i}.")
            continue
            
        # J-matrices are symmetric (J_ij = J_ji)
        J[i, j] = j_val
        J[j, i] = j_val
    
    if verbose:
        print("\nCalculated J Matrix (Hz):\n", np.array_str(J, precision=1, suppress_small=True))
        print("-" * 30)
        
    return J
