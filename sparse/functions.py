import numpy as np
from spin_database import _GAMMA_MAP
import itertools
from collections import defaultdict
from matplotlib.widgets import Cursor
import matplotlib.pyplot as plt


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


# --- PLOTTING---

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