import numpy as np
from scipy import sparse
from spin_operators_registry import get_spin_operators


def construct_sparse_hamiltonian(spins, v, J_matrix):
    nspins = len(spins)
    # Calculate the total dimension of the Hilbert space
    dim = int(np.prod([get_spin_operators(s)['dim'] for s in spins]))
    
    # Initialize the Hamiltonian as a sparse matrix
    H = sparse.lil_matrix((dim, dim), dtype=np.complex128)
    
    # 1. ADD ZEEMAN TERM (Diagonal in the basis: v_i * S_z,i)
    for n in range(nspins):
        Sz_full = 1
        for k in range(nspins):
            ops = get_spin_operators(spins[k])
            # If this is the active spin, use its S_z, else use identity
            if k == n:
                Sz_full = sparse.kron(Sz_full, ops['z']) if not isinstance(Sz_full, int) else ops['z']
            else:
                Sz_full = sparse.kron(Sz_full, ops['unit']) if not isinstance(Sz_full, int) else ops['unit']
        
        # Add Zeeman contribution to H
        H += v[n] * Sz_full

    # 2. ADD COUPLING TERM (J_ij * S_i * S_j)
    for i in range(nspins):
        for j in range(i + 1, nspins):
            if J_matrix[i, j] != 0:
                # Logic to kron your sparse matrices for coupling (x, y, z)
                for op_name in ['x', 'y', 'z']:
                    op_i_full = 1
                    op_j_full = 1
                    for k in range(nspins):
                        # Get operator for k-th spin
                        ops_k = get_spin_operators(spins[k])
                        
                        oi = ops_k[op_name] if k == i else ops_k['unit']
                        op_i_full = sparse.kron(op_i_full, oi) if not isinstance(op_i_full, int) else oi
                        
                        oj = ops_k[op_name] if k == j else ops_k['unit']
                        op_j_full = sparse.kron(op_j_full, oj) if not isinstance(op_j_full, int) else oj
                    
                    H += J_matrix[i, j] * (op_i_full @ op_j_full)
                    
    return H.tocsr()