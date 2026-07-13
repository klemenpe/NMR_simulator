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
from collections import defaultdict
from functions import convert_pairs_to_j_matrix
from hamiltonian import construct_sparse_hamiltonian
from transition_state_solver import calculate_transition_states

mpl.rcParams['pdf.fonttype'] = 42 # za pdf matplotlib
plt.rcParams.update({'font.size': 8}) # fontsize za matplotlib

__version__ = "0.1.2" 

# =====================================================================================================================#
#                         <<< S I M U L A T I O N   P A R A M E T E R S >>>
# =====================================================================================================================#

# 1. SPECTROMETER AND PLOTTING SETUP
spectrometer_1H_MHz = 600  # Spectrometer frequency (in MHz for 1H)

# Plotting settings
PLOT_NUCLEUS = 'H'             # Nucleus whose spectrum is displayed ('H', 'D', '13C', '19F', '31P', '29Si')
PLOT_COMBINED_SIGNALS = True   # If True, closely spaced transitions are grouped into one peak. Else False
tolerance_Hz = 200.0            # Max deviation (in Hz) from Larmor frequency to include a transition in the spectrum

# 2. NUMERICAL STABILIZATION SETTINGS
# To prevent numerical artifacts when shifts are near 0 ppm, we apply a temporary uniform shift.
# This trick stabilizes the Full Hamiltonian calculation without changing the physics of the coupling.
FREQUENCY_SCALING_MODE = False  # True
SCALING_SHIFT_PPM = 1  # Artificially shift all signals by 1 ppm for calculation stability. The results will be shifted back by -1 ppm before plotting.


# 3. MOLECULAR SETUP (H-H-D system)
# Define all nuclei in the spin system
spins = np.array([1/2, 1])
nuclei_types = ['H', 'D']

# Chemical shifts (in ppm) corresponding to the nuclei defined above
# H1 at 2.50 ppm, H2 at 1.0 ppm, D at 2.50 ppm
ppm_positions = [2.50, 2.50]

# --- J-COUPLINGS DEFINITION (Hz) ---
# Format: (Index i, Index j, J_value)
# Indices correspond to the order in the 'nuclei_types' list (0=H1, 1=D1, 2=D2)
J_COUPLING_PAIRS = [
    (0, 1, 2.0),  # J_H1D1
] 

# 3. TRANSITION FILTERING
cutoff = 0.001    # Intensity cutoff for raw transitions (lower value shows more peaks)
treshold_hz = 0.5 # Threshold (in Hz) for combining peaks if PLOT_COMBINED_SIGNALS is True

# =====================================================================================================================#


# Constants for Gyromagnetic Ratios (MHz/T)
_GAMMA_MAP = {
    'H': 42.577478461,  # Proton
    '13C': 10.7083991,  # Carbon-13
    'D': 6.53569888,    # Deuterium
    '19F': 40.078,     # Fluorine-19
    '31P': 17.2349,     # Phosphorus-31
    '29Si': -8.4650,    # Silicon-29
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


def generate_spin_states_by_mz(spins):
    """
    Generates spin states and groups them by their total Mz value.
    
    Returns:
        dict: A dictionary where keys are total Mz values, 
              and values are lists of state tuples.
    """
    # 1. Prepare mI values for each spin
    all_m_values = []
    for s in spins:
        if s == 0.5:
            all_m_values.append([0.5, -0.5])
        elif s == 1:
            all_m_values.append([1.0, 0.0, -1.0])
        else:
            raise ValueError(f"Unsupported spin value: {s}")

    # 2. Use a defaultdict to automatically group by Mz
    mz_basis = defaultdict(list)
    
    # 3. Generate products and group on the fly
    for state in itertools.product(*all_m_values):
        total_mz = sum(state)
        mz_basis[total_mz].append(state)
        
    return mz_basis

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


def construct_transition_matrix_mz(mz_groups, states, spins):
    """
    Constructs the transition matrix by only checking pairs of states 
    that could possibly be connected (same Mz or adjacent Mz).
    """
    # Note: For NMR, observable transitions (Ix) connect Mz and Mz +/- 1
    # For now, let's keep it simple and just build the connections
    # between the states you already have.
    
    num_states = len(states)
    T = np.zeros((num_states, num_states), dtype=int)

    # We only iterate through the states we know exist
    # This prevents checking millions of irrelevant combinations
    for i in range(num_states):
        for j in range(i + 1, num_states):
            state1 = states[i]
            state2 = states[j]
            
            # The Selection Rule:
            # 1. Only one nucleus changes
            # 2. That change is +/- 1.0
            diff_count = 0
            allowed_change_magnitude = True
            
            for k in range(len(spins)):
                diff = state1[k] - state2[k]
                if diff != 0:
                    diff_count += 1
                    if not np.isclose(abs(diff), 1.0):
                        allowed_change_magnitude = False
                        break
            
            if diff_count == 1 and allowed_change_magnitude:
                T[i, j] = 1
                T[j, i] = 1  # Symmetry
                
    return T


def ppm_to_hz(ppm_positions, nuclei_types, spectrometer_1H_MHz):
    """
    Converts chemical shift positions in ppm to frequencies in Hz.and returns the
    calculated frequencies AND the frequency map for the current spectrometer.
    
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

    return np.array(v_calculated), freq_MHz


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
#v = ppm_to_hz(ppm_positions, nuclei_types, spectrometer_1H_MHz)

# 1. PREPARE PPM INPUTS and TRACK SHIFTS
initial_ppm_positions = np.array(ppm_positions)
shift_ppm_value = SCALING_SHIFT_PPM if FREQUENCY_SCALING_MODE else 0
scaled_ppm_positions = initial_ppm_positions + shift_ppm_value

print(f"PPM Positions Original: {initial_ppm_positions}")
print(f"PPM Positions Shifted for Calculation: {scaled_ppm_positions}")

v_un_shifted, freq_MHz_map = ppm_to_hz(initial_ppm_positions, nuclei_types, spectrometer_1H_MHz)
v, _ = ppm_to_hz(scaled_ppm_positions, nuclei_types, spectrometer_1H_MHz)  # v for Hamiltonian

# Calculate the precise Hz shift applied to each Larmor frequency
v_shift_hz = v - v_un_shifted 
print(f"Un-shifted Larmor Frequencies (v_un_shifted, Hz): {v_un_shifted}")
print(f"Calculated Larmor Frequencies for Hamiltonian (v, Hz): {v}")
print(f"Hz Shift Applied to Each Nucleus (v_shift_hz): {v_shift_hz}")
print("-" * 30)

# --- J MATRIX CONSTRUCTION ---
J_matrix = convert_pairs_to_j_matrix(len(spins), J_COUPLING_PAIRS, verbose=True)

# Construct hamiltonian
# ----------------------------------------------------------#

H = construct_sparse_hamiltonian(spins, v, J_matrix)
print(f"Sparse Hamiltonian constructed. Shape: {H.shape}")

# Convert to dense if you want to compare with original script for verification only
#H_dense = H.toarray()

# --- AUTOMATIC TRANSITION MATRIX GENERATION FROM FUNCTION---

# 1. Generate the grouped states
mz_basis = generate_spin_states_by_mz(spins)

# 2. Flatten the dictionary into a simple list of states
all_states = []
for mz in sorted(mz_basis.keys(), reverse=True):
    all_states.extend(mz_basis[mz])

# new construct_transition_matrix_mz constructed matrix
T = construct_transition_matrix_mz(mz_basis, all_states, spins)

# 3. Solve transitions using the Mz Block Diagonalization engine
koncni_signali = calculate_transition_states(H, spins, cutoff=cutoff)

print("Calculated Signals (Frequency Hz, Intensity):")
# Show the first 20 transitions as a preview
num_peaks = koncni_signali.shape[1]
for idx in range(min(20, num_peaks)):
    freq = koncni_signali[0, idx]
    intensity = koncni_signali[1, idx]
    print(f"Transition {idx+1:02d}: Freq = {freq:8.3f} Hz | Intensity = {intensity:6.4f}")
if num_peaks > 20:
    print(f"... and {num_peaks - 20} more transitions.")

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
#tolerance_Hz = 300.0 

# Filter raw signals (koncni_signali is [2, N] Hz, Intensity)
raw_hz = koncni_signali[0]
is_target_signal = np.zeros(len(raw_hz), dtype=bool)

# Check proximity to ALL potential center frequencies
for center_freq in target_v_centers:
    # Use OR (|) to combine results: a signal is included if it's close to EITHER center
    is_target_signal = is_target_signal | np.isclose(raw_hz, center_freq, atol=tolerance_Hz, rtol=0.0)

filtered_raw_signals = koncni_signali[:, is_target_signal]

print(f"Plotting for {PLOT_NUCLEUS}. Found {filtered_raw_signals.shape[1]} transitions near centers: {target_v_centers}")

# --- NEW: CRITICAL CHECK FOR ZERO SIGNALS ---
if filtered_raw_signals.shape[1] == 0:
    print("\n" + "="*80)
    print(">>> CRITICAL ERROR: NO NMR TRANSITIONS FOUND FOR PLOTTING <<<")
    print("The simulation successfully calculated the transitions, but the plotting filter found zero signals.")
    print(f"Target Nucleus: {PLOT_NUCLEUS}")
    print(f"Calculated Larmor Frequencies (v) for target nucleus: {target_v_centers}")
    print(f"Current Tolerance: {tolerance_Hz:.1f} Hz")
    print("\nThis happens when the distance between the transitions and the Larmor frequency exceeds 'tolerance_Hz'.")
    print("The simplest fixes are:")
    print("1. INCREASE the 'tolerance_Hz' value (currently set to {tolerance_Hz:.1f} Hz).")
    print("2. CHANGE the 'ppm_positions' of the target nucleus to a value closer to 0.")
    print("="*80 + "\n")
    # Raise an error to halt execution cleanly and prevent the subsequent NaN/Inf plot error
    raise ValueError(f"No NMR transitions found within the {tolerance_Hz} Hz tolerance band for {PLOT_NUCLEUS}. See console output for guidance.")


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


# ---------------------------------------------------------------------------------------------------------------------#
# --- POST-COMBINATION FREQUENCY SHIFT (Final correction of grouped signals) ---
# ---------------------------------------------------------------------------------------------------------------------#
if FREQUENCY_SCALING_MODE and zdruzeni_signali_filt.size > 0:
    print("\nApplying final shift correction to combined signals: Reverting the nucleus-specific Hz shift.")
    
    # Iterate through the combined signals (which are still in the shifted Hz scale)
    for i in range(zdruzeni_signali_filt.shape[0]):
        combined_peak_freq_shifted = zdruzeni_signali_filt[i, 0] # Peak in shifted Hz

        # 1. Find the index of the closest Larmor frequency (using SHIFTED values 'v') 
        #    to determine which nucleus the peak belongs to.
        closest_nucleus_index = np.argmin(np.abs(combined_peak_freq_shifted - v))
        
        # 2. Get the specific Hz shift applied to that nucleus (v_shift_hz)
        correction_hz = v_shift_hz[closest_nucleus_index]
        
        # 3. Apply the correction by subtracting the shift. This moves the peak back 
        #    to its correct, un-shifted position relative to the original Larmor center.
        zdruzeni_signali_filt[i, 0] -= correction_hz
        
    print("Combined signals successfully reverted to the correct, un-shifted scale.")

    # Custom print message for when scaling was performed
    final_print_message = f"Filtered, combined and shifted to original scale signals for {PLOT_NUCLEUS} (Hz, Intensity):\n"
else:
    # Standard print message when scaling was NOT performed
    final_print_message = f"Filtered, combined signals for {PLOT_NUCLEUS} (Hz, Intensity):\n"

# Print the final result using the customized message
print(final_print_message)
print(f"{'Frequency (Hz)':<20} {'Intensity':<20}")
for freq, intensity in zdruzeni_signali_filt:
    print(f"{freq:<20.10f} {intensity:<20.10f}")
print("-" * 40)

# testi
tt = generate_spin_states(spins)
tt2 = generate_spin_states_by_mz(spins)

tm = construct_transition_matrix(spins)


# PLOT GRAPH
# ---------------------------------------------------------------------------------------------------------------------#

# --- EXECUTION OF PLOT ---
if PLOT_COMBINED_SIGNALS:
    data_to_plot = zdruzeni_signali_filt
else:
    data_to_plot = filtered_raw_signals.T

plot_nmr_spectrum(data_to_plot, treshold_hz, cutoff, PLOT_NUCLEUS, PLOT_COMBINED_SIGNALS, spectrometer_freqs_MHz)
