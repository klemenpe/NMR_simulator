import numpy as np

SQRT2_INV = np.sqrt(2) / 2
SQRT3_INV = np.sqrt(3) / 2

# write spin matices here
SPIN_CONFIG = {
    0.5: {
        'dim': 2,
        'x': 0.5 * np.array([[0, 1], [1, 0]]),
        'y': 0.5 * np.array([[0, -1j], [1j, 0]]),
        'z': 0.5 * np.array([[1, 0], [0, -1]]),
        'unit': np.eye(2)
    },
    1.0: {
        'dim': 3,
        'x': SQRT2_INV * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]]),
        'y': SQRT2_INV * np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]]),
        'z': np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]]),
        'unit': np.eye(3)
    }
}


def get_spin_operators(s):
    """
    Returns the operator dictionary for a given spin value s.
    Uses the SPIN_CONFIG dictionary for clean data separation.
    """

    if s not in SPIN_CONFIG:
        raise ValueError(f"Spin value {s} is not defined in the configuration.")
    
    return SPIN_CONFIG[s]
