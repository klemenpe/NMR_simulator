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


