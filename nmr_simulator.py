# NMR SIMULATOR
# ---------------------------------------------------------------------------------------------------------------------#
# Simple sample for simulating DMSO CHD2 izotopologue  (H, H, D). adjust values to your problem
#
# 
# ---------------------------------------------------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import matplotlib as mpl
import itertools

mpl.rcParams['pdf.fonttype'] = 42 # za pdf matplotlib
plt.rcParams.update({'font.size': 8}) # fontsize za matplotlib


# =====================================================================================================================#
#                         <<< S I M U L A T I O N   P A R A M E T E R S >>>
# =====================================================================================================================#

# 1. SPECTROMETER AND PLOTTING SETUP
spectrometer_1H_MHz = 600  # Spectrometer frequency (in MHz for 1H)

# Plotting settings
PLOT_NUCLEUS = 'H'             # Nucleus whose spectrum is displayed ('H', 'D', or '13C')
PLOT_COMBINED_SIGNALS = True   # If True, closely spaced transitions are grouped into one peak. Else False
tolerance_Hz = 100.0           # Max deviation (in Hz) from Larmor frequency to include a transition in the spectrum

# 2. MOLECULAR SETUP (H-H-D system)
# Define all nuclei in the spin system
spins = np.array([1/2, 1, 1])
nuclei_types = ['H', 'D', 'D']

# Chemical shifts (in ppm) corresponding to the nuclei defined above
# H1 at 2.50 ppm, H2 at 1.0 ppm, D at 2.50 ppm
ppm_positions = [2.50, 2.50, 2.50]

# --- J-COUPLINGS DEFINITION (Hz) ---
# Format: (Index i, Index j, J_value)
# Indices correspond to the order in the 'nuclei_types' list (0=H1, 1=D1, 2=D2)
J_COUPLING_PAIRS = [
    (0, 1, 2.0),  # J_H1D1
    (0, 2, 2.0),  # J_H1D2
    (1, 2, 0.0)   # J_D1D2
] 

# J-Couplings (in Hz) matrix. Order must match the nuclei_types: H1, H2, D
#J_HD = 2.0  # Coupling H-D in Hz
#J_DD = 0.0  # D-D coupling is zero
#
#J = np.array([[0.0,  J_HD, J_HD],
#              [J_HD, 0.0,  J_DD],
#              [J_HD, J_DD, 0.0]])

# 3. TRANSITION FILTERING
cutoff = 0.000    # Intensity cutoff for raw transitions (lower value shows more peaks)
treshold_hz = 0.5 # Threshold (in Hz) for combining peaks if PLOT_COMBINED_SIGNALS is True

# =====================================================================================================================#


# Constants for Gyromagnetic Ratios (MHz/T)
_GAMMA_MAP = {
    'H': 42.577478461,  # Proton
    '13C': 10.7083991,  # Carbon-13
    'D': 6.53569888,    # Deuterium
}

# ---------------------------------------------------------------------------------------------------------------------#
# --- SIMULATION FUNCTIONS --- (Moved below parameters for better user focus)
# ---------------------------------------------------------------------------------------------------------------------#

def get_spectrometer_frequencies(spectrometer_1H_MHz):
    """
    Calculates the spectrometer frequency in MHz for all known nuclei based on the 1H frequency.
    """
    freq_MHz = {}
    gamma_1H = _GAMMA_MAP['H']
    for key, gamma_nucleus in _GAMMA_MAP.items():
        # Larmor frequency (MHz) for nucleus 'key' at the magnetic field B0
        freq_MHz[key] = spectrometer_1H_MHz * (gamma_nucleus / gamma_1H)
    return freq_MHz

def generate_spin_states(spins):
    """
    Generates all possible combinations of mI values for a given set of spins.
    The order of mI values for each spin is assumed to be descending
    (e.g., for spin 1/2: [1/2, -1/2], for spin 1: [1, 0, -1]).
    """
    all_m_values = []
    for s in spins:
        if s == 0.5:
            all_m_values.append([0.5, -0.5])
        elif s == 1:
            all_m_values.append([1.0, 0.0, -1.0])
        else:
            raise ValueError(f"Unsupported spin value: {s}. Only 1/2 and 1 are supported.")

    # Use itertools.product to get all combinations
    states = list(itertools.product(*all_m_values))
    return states

#  transition matrix
def construct_transition_matrix(spins):
    """
    Constructs the allowed transition matrix based on NMR selection rules.
    An allowed transition means only one nucleus's mI value changes by +/- 1.

    Args:
        spins (np.array): Array of spin values for each nucleus (e.g., [0.5, 1]).

    Returns:
        np.array: The transition matrix where 1 indicates an allowed transition, 0 otherwise.
    """
    all_states = generate_spin_states(spins)
    num_states = len(all_states)
    T = np.zeros((num_states, num_states), dtype=int)

    # Optimized loop: only calculate upper triangle and then add the lower.
    # This also avoids checking i == j, as i < j naturally.
    for i in range(num_states):
        for j in range(i + 1, num_states):  # Start j from i + 1
            state1 = all_states[i]  # final state
            state2 = all_states[j]  # initial state

            diff_count = 0
            allowed_change_magnitude = True

            for k in range(len(spins)):  # Iterate through each nucleus
                diff = state1[k] - state2[k]

                if diff != 0:
                    diff_count += 1
                    # Check if the change is exactly +/- 1.0 for this nucleus
                    if not np.isclose(abs(diff), 1.0):  # Use np.isclose for float comparison
                        allowed_change_magnitude = False
                        break  # Not an allowed +/-1 change for this nucleus

            # If exactly one nucleus changed, and that change was +/-1
            if diff_count == 1 and allowed_change_magnitude:
                T[i, j] = 1  # Set upper triangle
    T += T.T  # Add the lower triangle by transposing and adding
    return T


def ppm_to_hz(ppm_positions, nuclei_types, spectrometer_1H_MHz):
    """
    Converts chemical shift positions in ppm to frequencies in Hz.
    
    Args:
        ppm_positions (list/array): Chemical shifts in ppm.
        nuclei_types (list/array): Type of nucleus for each ppm value (e.g., '13C', 'H', 'D').
        spectrometer_1H_MHz (float): Operating frequency of the spectrometer for 1H.
        
    Returns:
        np.array: Frequencies in Hz.
    """
    if len(ppm_positions) != len(nuclei_types):
        raise ValueError("ppm_positions and nuclei_types must have the same length.")

    v_calculated = []
    
    # Calculate the spectrometer frequency (MHz) for each nucleus type relative to 1H
    freq_MHz = {}
    gamma_1H = _GAMMA_MAP['H']
    for key, gamma_nucleus in _GAMMA_MAP.items():
        freq_MHz[key] = spectrometer_1H_MHz * (gamma_nucleus / gamma_1H)
        
    for i, ppm in enumerate(ppm_positions):
        nucleus_type = nuclei_types[i]
        
        if nucleus_type not in freq_MHz:
            raise ValueError(f"Unknown nucleus type: {nucleus_type}. Must be one of {list(_GAMMA_MAP.keys())}")
        
        # Shift in Hz = ppm * Spectrometer Frequency (MHz)
        shift_in_Hz = ppm * freq_MHz[nucleus_type]
        v_calculated.append(shift_in_Hz)

    return np.array(v_calculated)


def construct_j_matrix(spins, J_COUPLING_PAIRS):
    """
    Constructs the symmetric J-coupling matrix (Hz) from a list of coupling pairs.
    Also prints the resulting matrix to the console for user verification.
    
    Args:
        spins (np.array): Array of spin values for each nucleus.
        J_COUPLING_PAIRS (list): List of (index i, index j, J value) tuples.
        
    Returns:
        np.array: The symmetric J matrix.
    """
    print("-" * 30)
    print(">>> J-COUPLING MATRIX CONSTRUCTION <<<")
    nspins = len(spins)
    J = np.zeros((nspins, nspins))

    for i, j, j_val in J_COUPLING_PAIRS:
        # Set symmetric matrix elements (J_ij = J_ji)
        if i == j:
            print(f"Warning: Ignoring J-coupling for atom to itself at index {i}.")
            continue
        if i >= nspins or j >= nspins:
            raise IndexError(f"J-coupling index out of bounds: ({i}, {j}). Max index is {nspins - 1}.")

        J[i, j] = j_val
        J[j, i] = j_val

    # Display J Matrix for verification
    print("\nCalculated J Matrix (Hz):")
    # Use np.array_str for a cleaner, centered display of the matrix
    print(np.array_str(J, precision=1, suppress_small=True))
    print("-" * 30)
    
    return J


def plot_nmr_spectrum(data_to_plot_T, treshold_hz, cutoff, plot_nucleus, plot_combined, spectrometer_freqs_MHz):
    """
    Plots the NMR spectrum for the filtered nucleus with dynamic x-axis limits.
    """
    # --- LOCAL ONPICK FUNCTION DEFINITION---
    # Defined here to access 'data_to_plot_T' (which holds Hz and Intensity data)
    def onpick(event1):
        """
        Display the picked point's coordinates (ppm, Hz, Intensity) in the console.
        """
        thisline = event1.artist
        # xdata is ppm_skala (ppm)
        xdata = thisline.get_xdata() 
        ydata = thisline.get_ydata() # ydata is Intensity
        ind = event1.ind
        
        # Look up the original Hz value using the index 'ind' from the data array.
        # data_to_plot_T[:, 0] holds the Hz values.
        # 'ind' is an array of indices, typically containing one element for a single pick.
        hz_value = data_to_plot_T[ind, 0][0]
        
        ppm_value = xdata[ind][0]
        intensity_value = ydata[ind][0]
        index_value = ind[0]
        
        # Display the required output format: (ppm, Hz, intensity, index)
        print(f'onpick points (ppm, Hz, intensity, index): ({ppm_value:.4f}, {hz_value:.2f}, {intensity_value:.4f}, {index_value})')
           
    # --- END LOCAL ONPICK FUNCTION DEFINITION ---

    if plot_nucleus not in spectrometer_freqs_MHz:
        print(f"Error: Nucleus '{plot_nucleus}' not configured for plotting scale.")
        return

    ref_freq_MHz = spectrometer_freqs_MHz[plot_nucleus]
    
    # Scaling (Hz -> ppm)
    ppm_skala = data_to_plot_T[:, 0] / ref_freq_MHz
    intensities = data_to_plot_T[:, 1]

    cm = 1 / 2.54
    fig, ax = plt.subplots(figsize=(16 * cm, 10 * cm))

    # Use stem plot for discrete spectrum lines
    markerline, stemlines, baseline = ax.stem(ppm_skala, intensities, linefmt='k-', markerfmt='ko', basefmt='k-')

    # Styles
    plt.setp(markerline, picker=5, color='k', markersize=2)
    plt.setp(stemlines, linewidth=1, color='k')
    plt.setp(baseline, linewidth=1, color='k')
    
    baseline.set_xdata([0, 1])
    baseline.set_transform(ax.get_yaxis_transform())

    # --- DYNAMIC X-AXIS LIMITS ---
    if len(ppm_skala) > 0:
        min_ppm = np.min(ppm_skala)
        max_ppm = np.max(ppm_skala)
        
        # Add 1 ppm buffer and invert the axis for NMR display
        x_min_limit = min_ppm - 1.0
        x_max_limit = max_ppm + 1.0
        
        # Set limits as (highest ppm, lowest ppm) to achieve inversion
        ax.set_xlim(x_max_limit, x_min_limit)
        
        print(f"Plotting range set from {x_min_limit:.2f} ppm to {x_max_limit:.2f} ppm.")
    else:
        ax.set_xlim(10, 0) # Default range if no signals found
        
    # Annotations and Interactivity
    title_mode = "Combined" if plot_combined else "Raw Transitions"
    title_params = f"Treshold: {treshold_hz} Hz" if plot_combined else f"Cutoff: {cutoff}"
    
    plt.title(f'Simulated NMR Spectrum ({plot_nucleus} Scale) - {title_mode} ({title_params})')
    plt.xlabel(fr'Chemical Shift $\delta$ ({plot_nucleus} ppm)')
    plt.ylabel('Relative Intensity')

    # Optional: Display the corresponding frequency scale on a second axis
    ax2 = ax.secondary_xaxis('top', functions=(lambda p: p * ref_freq_MHz, lambda h: h / ref_freq_MHz))
    ax2.set_xlabel('Frequency (Hz)')

    fig.canvas.mpl_connect('pick_event', onpick)

    # CURSOR STYLE (Crosshair)
    cursor_ref = Cursor(ax, useblit=True, color='k', linewidth=1)

    plt.tight_layout()
    plt.show(block=False)


# ---------------------------------------------------------------------------------------------------------------------#
# --- MAIN PROGRAM EXECUTION ---
# ---------------------------------------------------------------------------------------------------------------------#
# Convert PPM positions to Larmor Frequencies (v)
v = ppm_to_hz(ppm_positions, nuclei_types, spectrometer_1H_MHz)

# --- J MATRIX CONSTRUCTION ---
J = construct_j_matrix(spins, J_COUPLING_PAIRS)

# Define Pauli matrices (including spin-1)
sigma_x_half = np.array([[0, 1/2], [1/2, 0]])
sigma_y_half = np.array([[0, -1j/2], [1j/2, 0]])
sigma_z_half = np.array([[1/2, 0], [0, -1/2]])
unit_half = np.array([[1, 0], [0, 1]])

# Spin 1 matrices
sigma_x_1 = np.array([[0, np.sqrt(2)/2, 0], [np.sqrt(2)/2, 0, np.sqrt(2)/2], [0, np.sqrt(2)/2, 0]])
sigma_y_1 = np.array([[0, (-1j*np.sqrt(2))/2, 0], [(1j*np.sqrt(2))/2, 0, (-1j*np.sqrt(2))/2], [0, (1j*np.sqrt(2))/2, 0]])
sigma_z_1 = np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]])
unit_1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

# DEFINE OPERATORS
# ----------------------------------------------------------#
nspins = len(v)  # number of spins
nspins_half = np.count_nonzero(spins == 1/2)  # count number of 1/2 spins
nspins_one = np.count_nonzero(spins == 1)  # count number of 1 spins

# construct spin operators

size_states = 2**nspins_half * 3**nspins_one  # size of matrix for all possible states

L = np.empty((3, nspins, size_states, size_states), dtype=np.complex128)  # preaalocate space

for n in range(nspins):
    Lx_current = 1
    Ly_current = 1
    Lz_current = 1
    for k in range(nspins):
        if k == n:
            if spins[n] == 1/2:
                Lx_current = np.kron(Lx_current, sigma_x_half)
                Ly_current = np.kron(Ly_current, sigma_y_half)
                Lz_current = np.kron(Lz_current, sigma_z_half)
            elif spins[n] == 1:
                Lx_current = np.kron(Lx_current, sigma_x_1)
                Ly_current = np.kron(Ly_current, sigma_y_1)
                Lz_current = np.kron(Lz_current, sigma_z_1)
        else:
            if spins[k] == 1/2:
                Lx_current = np.kron(Lx_current, unit_half)
                Ly_current = np.kron(Ly_current, unit_half)
                Lz_current = np.kron(Lz_current, unit_half)
            elif spins[k] == 1:
                Lx_current = np.kron(Lx_current, unit_1)
                Ly_current = np.kron(Ly_current, unit_1)
                Lz_current = np.kron(Lz_current, unit_1)

    L[0, n] = Lx_current
    L[1, n] = Ly_current
    L[2, n] = Lz_current

# constructing operator for spin-spin coupling interactions

L_T = L.transpose(1, 0, 2, 3)
Lproduct = np.tensordot(L_T, L, axes=((1, 3), (0, 2))).swapaxes(1, 2)

# calculate hamiltonian first part
Lz = L[2]  # spin-spin coupling interactions in z axis Lz
H = np.tensordot(v, Lz, axes=1)

# calculate hamiltonian second part and combine with first part
scalars = 0.5 * J
H += np.tensordot(scalars, Lproduct, axes=2)

# --- AUTOMATIC TRANSITION MATRIX GENERATION FROM FUNCTION---
T = construct_transition_matrix(spins)

# Diagonalize and find transition intensities and energies
E, V = np.linalg.eigh(H)
V = V.real
I = np.square(V.T.dot(T.dot(V)))
I_upper = np.triu(I)  # symmetry makes it possible to use only one half of the matrix for faster calculation
E_matrix = np.abs(E[:, np.newaxis] - E)
E_upper = np.triu(E_matrix)
combo = np.stack([E_upper, I_upper])
iv = combo.reshape(2, I.shape[0] ** 2).T

# Filter raw peaks by intensity
peaklist = iv[iv[:, 1] >= cutoff]

# Raw signals array: [Hz, Intensity]
koncni_signali = peaklist.T

print(f"Total transitions found above cutoff {cutoff}: {koncni_signali.shape[1]}")
print("-" * 30)

# ---------------------------------------------------------------------------------------------------------------------#
# --- SIGNAL FILTERING BASED ON PLOT_NUCLEUS ---
# ---------------------------------------------------------------------------------------------------------------------#

spectrometer_freqs_MHz = get_spectrometer_frequencies(spectrometer_1H_MHz)

# Find ALL center frequencies (v) for the PLOT_NUCLEUS
target_nucleus = PLOT_NUCLEUS
# Indices of the target nucleus type (e.g., indices 0 and 1 for 'H')
target_indices = [i for i, n_type in enumerate(nuclei_types) if n_type == target_nucleus]

if not target_indices:
    raise ValueError(f"Nucleus type {target_nucleus} not found in system setup.")

# Get all Larmor frequencies that belong to the target nucleus type (e.g., [1500 Hz, 600 Hz])
target_v_centers = v[target_indices]

# Use a tolerance large enough to catch all coupled signals, but small enough to exclude other nuclei
tolerance_Hz = 200.0 

# Filter raw signals (koncni_signali is [2, N] Hz, Intensity)
raw_hz = koncni_signali[0]
is_target_signal = np.zeros(len(raw_hz), dtype=bool)

# Check proximity to ALL potential center frequencies
for center_freq in target_v_centers:
    # Use OR (|) to combine results: a signal is included if it's close to EITHER center
    is_target_signal = is_target_signal | np.isclose(raw_hz, center_freq, atol=tolerance_Hz, rtol=0.0)

filtered_raw_signals = koncni_signali[:, is_target_signal]

print(f"Plotting for {PLOT_NUCLEUS}. Found {filtered_raw_signals.shape[1]} transitions near centers: {target_v_centers}")


# ---------------------------------------------------------------------------------------------------------------------#
# COMBINES SIGNALS CLOSE TOGETHER (Group the filtered signals)
# ---------------------------------------------------------------------------------------------------------------------#

koncni_signali2_filt = filtered_raw_signals.T
koncni_signali3_filt = koncni_signali2_filt[koncni_signali2_filt[:,0].argsort()].T  # Sort by frequency

differences = np.abs(np.diff(koncni_signali3_filt[0]))
group_separators = np.where(differences > treshold_hz)[0] + 1
group_separators = np.concatenate([np.array([0]), group_separators, np.array([len(koncni_signali3_filt[0])])])

combined_arrays = []
for i in range(len(group_separators) - 1):
    sub_array = koncni_signali3_filt[:, group_separators[i]:group_separators[i+1]]
    
    # Calculate average frequency and sum intensity
    average_row = np.mean(sub_array[0], axis=0)
    sum_row = np.sum(sub_array[1] , axis=0)
    zdruzen = np.vstack([average_row, sum_row])
    combined_arrays.append(zdruzen)

zdruzeni_signali_filt = np.hstack(combined_arrays).T if combined_arrays else np.empty((0, 2))
print(f"Filtered, combined signals for {PLOT_NUCLEUS} (Hz, Intensity):\n", zdruzeni_signali_filt)
print("-" * 30)

# PLOT GRAPH
# ---------------------------------------------------------------------------------------------------------------------#

# --- EXECUTION OF PLOT ---
if PLOT_COMBINED_SIGNALS:
    data_to_plot = zdruzeni_signali_filt
else:
    data_to_plot = filtered_raw_signals.T

plot_nmr_spectrum(data_to_plot, treshold_hz, cutoff, PLOT_NUCLEUS, PLOT_COMBINED_SIGNALS, spectrometer_freqs_MHz)
