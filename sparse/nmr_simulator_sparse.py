# NMR SIMULATOR
# ---------------------------------------------------------------------------------------------------------------------#
# Simple sample for simulating DMSO CHD2 izotopologue  (H, H, D). adjust values to your problem
#
# 
# ---------------------------------------------------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import functions
import matplotlib as mpl
from hamiltonian import construct_sparse_hamiltonian
from transition_state_solver import calculate_transition_states
import time

mpl.rcParams['pdf.fonttype'] = 42 # za pdf matplotlib
plt.rcParams.update({'font.size': 8}) # fontsize za matplotlib

__version__ = "1.1.0" 

# start program execution time
start_time = time.time()

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

# 3. TRANSITION FILTERING
cutoff = 0.001    # Intensity cutoff for raw transitions (lower value shows more peaks)
treshold_hz = 0.5  # Threshold (in Hz) for combining peaks if PLOT_COMBINED_SIGNALS is True

# =====================================================================================================================#

# 1. PREPARE PPM INPUTS and TRACK SHIFTS
# ----------------------------------------------------------#
initial_ppm_positions = np.array(ppm_positions)
shift_ppm_value = SCALING_SHIFT_PPM if FREQUENCY_SCALING_MODE else 0
scaled_ppm_positions = initial_ppm_positions + shift_ppm_value

print(f"PPM Positions Original: {initial_ppm_positions}")
print(f"PPM Positions Shifted for Calculation: {scaled_ppm_positions}")

v_un_shifted, freq_MHz_map = functions.ppm_to_hz(initial_ppm_positions, nuclei_types, spectrometer_1H_MHz)
v, _ = functions.ppm_to_hz(scaled_ppm_positions, nuclei_types, spectrometer_1H_MHz)  # v for Hamiltonian

# Calculate the precise Hz shift applied to each Larmor frequency
v_shift_hz = v - v_un_shifted 
print(f"Un-shifted Larmor Frequencies (v_un_shifted, Hz): {v_un_shifted}")
print(f"Calculated Larmor Frequencies for Hamiltonian (v, Hz): {v}")
print(f"Hz Shift Applied to Each Nucleus (v_shift_hz): {v_shift_hz}")
print("-" * 30)

# ---------------------------------------------------------------------------------------------------------------------#
#                            <<< M A I N  P R O G R A M  E X E C U T I O N >>>
# ---------------------------------------------------------------------------------------------------------------------#

# --- J MATRIX CONSTRUCTION ---
# ----------------------------------------------------------#
J_matrix = functions.convert_pairs_to_j_matrix(len(spins), J_COUPLING_PAIRS, verbose=True)

# --- CONSTRUCT HAMILTONIAN ---
# ----------------------------------------------------------#

H = construct_sparse_hamiltonian(spins, v, J_matrix)
print(f"Sparse Hamiltonian constructed. Shape: {H.shape}")

# Convert to dense if you want to compare with original script for verification only
#H_dense = H.toarray()

# --- AUTOMATIC TRANSITION MATRIX GENERATION ---
# ----------------------------------------------------------#
# 1. Generate the grouped states
mz_basis = functions.generate_spin_states_by_mz(spins)

# 2. Flatten the dictionary into a simple list of states
all_states = []
for mz in sorted(mz_basis.keys(), reverse=True):
    all_states.extend(mz_basis[mz])

# new construct_transition_matrix_mz constructed matrix
T = functions.construct_transition_matrix_mz(mz_basis, all_states, spins)

#--- USE TRANSITION STATE SOLVER FOR FINAL TRANSITION CALCULATION ---
# ----------------------------------------------------------#
# 3. Solve transitions using the Mz Block Diagonalization engine
koncni_signali = calculate_transition_states(H, spins, cutoff=cutoff)

print("Calculated Signals (Frequency Hz, Intensity):")
num_peaks = koncni_signali.shape[1]
for idx in range(num_peaks):
    freq = koncni_signali[0, idx]
    intensity = koncni_signali[1, idx]
    print(f"Transition {idx+1:02d}: Freq = {freq:8.3f} Hz | Intensity = {intensity:6.4f}")


# ---------------------------------------------------------------------------------------------------------------------#
#                            <<< P O S T  P R O C E S S I N G >>>
# ---------------------------------------------------------------------------------------------------------------------#

# --- SIGNAL FILTERING BASED ON PLOT_NUCLEUS ---
# ---------------------------------------------------------------------------------------------------------------------#

spectrometer_freqs_MHz = functions.get_spectrometer_frequencies(spectrometer_1H_MHz)

# Find ALL center frequencies (v) for the PLOT_NUCLEUS
target_nucleus = PLOT_NUCLEUS
# Indices of the target nucleus type (e.g., indices 0 and 1 for 'H')
target_indices = [i for i, n_type in enumerate(nuclei_types) if n_type == target_nucleus]

if not target_indices:
    raise ValueError(f"Nucleus type {target_nucleus} not found in system setup.")

# Get all Larmor frequencies that belong to the target nucleus type (e.g., [1500 Hz, 600 Hz])
target_v_centers = v[target_indices]

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


# end execution time
# ---------------------------------------------------------------------------------------------------------------------#
end_time = time.time()    # 3. Record the end time
elapsed_time = (end_time - start_time)/60

print(f"\nExecution Time: {elapsed_time:.4f} minutes")

# PLOT GRAPH
# ---------------------------------------------------------------------------------------------------------------------#

# --- EXECUTION OF PLOT ---
if PLOT_COMBINED_SIGNALS:
    data_to_plot = zdruzeni_signali_filt
else:
    data_to_plot = filtered_raw_signals.T

functions.plot_nmr_spectrum(data_to_plot, treshold_hz, cutoff, PLOT_NUCLEUS, PLOT_COMBINED_SIGNALS, spectrometer_freqs_MHz)
