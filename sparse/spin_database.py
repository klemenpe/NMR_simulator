import numpy as np


# ==============================================================================
#                       <<< PHYSICAL NUCLEI DATABASE >>>
# ==============================================================================
# Gyromagnetic ratio in (MHz/T) for all isotope physical properties.
# If you need to add a new nucleus to your simulator, simply add it here!

_GAMMA_MAP = {
    'H': 42.577478461,  # Proton
    '13C': 10.7083991,  # Carbon-13
    'D': 6.53569888,    # Deuterium
    '19F': 40.078,     # Fluorine-19
    '31P': 17.2349,     # Phosphorus-31
    '29Si': -8.4650,    # Silicon-29
}

# ==============================================================================
#                       <<< QUANTUM SPIN OPERATORS >>>
# ==============================================================================
# Standard spin projection matrices.

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
