from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq, irfft

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from pp.colors import colors
from scipy.signal import welch


def detect_fundamentals(peak_freqs, tol=0.03):
    """
    Detect fundamental frequencies f0, f1,... such that all other peaks
    can be written as integer combinations n*f0 + m*f1.
    """
    fundamentals = []

    for f in peak_freqs:
        is_combination = False

        # --------------------------------------------------------
        # 1) First peak → automatic fundamental
        # --------------------------------------------------------
        if len(fundamentals) == 0:
            fundamentals.append(f)
            continue

        # --------------------------------------------------------
        # 2) Check pure harmonics: f ≈ n * f0
        # --------------------------------------------------------
        for f0 in fundamentals:
            n = round(f / f0)
            if np.abs(f - n * f0) < tol * f:
                is_combination = True
                break

        if is_combination:
            continue

        # --------------------------------------------------------
        # 3) Check mixed harmonics: f ≈ n*f0 + m*f1
        # --------------------------------------------------------
        if len(fundamentals) >= 2:
            for i in range(len(fundamentals)):
                for j in range(i + 1, len(fundamentals)):
                    f0 = fundamentals[i]
                    f1 = fundamentals[j]

                    # Solve approximately for n, m
                    n = round(f / f0)
                    m = round((f - n * f0) / f1)

                    f_est = n * f0 + m * f1
                    if np.abs(f - f_est) < tol * f:
                        is_combination = True
                        break
                if is_combination:
                    break

        if is_combination:
            continue

        # --------------------------------------------------------
        # 4) If not matched → new fundamental
        # --------------------------------------------------------
        fundamentals.append(f)

    return fundamentals


def refine_peak(omega, amp, k):
    """
    Parabolic interpolation around index k to estimate sub-bin peak location.
    Returns refined_frequency, refined_amplitude.
    """
    if k <= 0 or k >= len(amp) - 1:
        return omega[k], amp[k]  # cannot refine at the edges

    y1, y2, y3 = amp[k - 1], amp[k], amp[k + 1]

    # Parabolic vertex offset relative to k
    denom = y1 - 2 * y2 + y3
    if denom == 0:
        delta = 0
    else:
        delta = 0.5 * (y1 - y3) / denom

    # refine frequency
    domega = omega[k + 1] - omega[k]
    refined_freq = omega[k] + delta * domega

    # refined amplitude at vertex:
    refined_amp = y2 - 0.25 * (y1 - y3) * delta

    return refined_freq, refined_amp


def plot_harmonics(
    omega,
    fft_data,
    width: float,
    depth: float,
    point: np.ndarray,
    variable_label: str,
    prominence: float = 0.15,
) -> np.ndarray:
    """
    Detect fundamental harmonics and plot their superharmonics
    in the FFT spectrum. Each fundamental gets its own color.
    Args:
        omega: array of angular frequencies
        fft_data: FFT complex spectrum
        prominence: relative prominence for peak detection
    """
    amp = np.abs(fft_data)
    if np.max(amp) == 0:
        print(
            colors.WARNING
            + f"No variable signal detected for variable {variable_label} at point x={point[0]:.2f}, y={point[1]:.2f}, z={point[2]:.2f} (d={depth:.2f}, w={width:.2f}). Skipping FFT analysis."
            + colors.ENDC
        )
        return np.array([])
    amp_norm = amp / np.max(amp)

    # -------- 1) Detect all peaks above prominence ----------
    peaks, _ = find_peaks(amp_norm, prominence=prominence)
    peak_freqs = []
    peak_amps = []

    for p in peaks:
        f_ref, a_ref = refine_peak(omega, amp_norm, p)
        peak_freqs.append(f_ref)
        peak_amps.append(a_ref)

    peak_freqs = np.array(peak_freqs)
    peak_amps = np.array(peak_amps)

    # Sort peaks from low to high frequency
    idx = np.argsort(peak_freqs)
    peak_freqs = peak_freqs[idx]
    peak_amps = peak_amps[idx]

    # -------- 2) Remove harmonics: keep only fundamentals ------
    fundamentals = []
    tol = 0.03  # tolerance for identifying integer multiples

    fundamentals = detect_fundamentals(peak_freqs, tol=tol)

    fundamentals = np.array(fundamentals)

    print(
        "Detected fundamental frequencies (ω) for variable "
        f"{variable_label} at point x={point[0]:.2f}, y={point[1]:.2f}, z={point[2]:.2f} (d={depth:.2f}, w={width:.2f}):"
    )
    for omega0 in fundamentals:
        print(f"   ω₀ = {omega0:.6f}")
    try:
        print("Frequency with highest amplitude:")
        print(f"   ω_max = {peak_freqs[np.argmax(peak_amps)]:.6f}")
        print("")
    except ValueError:
        print(
            colors.WARNING
            + "No peaks detected above the prominence threshold."
            + colors.ENDC
        )
        return fundamentals

    # -------- 3) Plot FFT ----------------
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlabel("ω")
    ax.set_ylabel("normalized FFT")
    ax.set_title(
        f"Spectrum of {variable_label} at x = {point[0]:.2f}, y = {point[1]:.2f}, z = {point[2]:.2f} (d{depth:.2f}_w{width:.2f})"
    )
    # -------- 4) Plot superharmonics ----------
    # skip first color (already used in the line plot)

    # colors = plt.cm.tab10(np.linspace(0.1, 1, len(fundamentals)))

    # for color, omega0 in zip(colors, fundamentals):
    #     # Fundamental
    #     ax.axvline(omega0, color=color, linewidth=2.5, label=f"ω₀ = {omega0:.3f}", alpha=0.3)

    #     # Superharmonics
    #     n = 2
    #     while n * omega0 < omega[-1]:
    #         ax.axvline(n*omega0, color=color, linestyle="--", alpha=0.7)
    #         n += 1

    ax.plot(omega, amp_norm, label="FFT amplitude")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend()
    return fundamentals


def fftFreqs(
    time: np.ndarray,
    signal: np.ndarray,
    width: float,
    depth: float,
    point: np.ndarray,
    variable_label: str,
) -> np.ndarray:
    """
    Perform FFT on the given signal and print frequencies with amplitudes above a threshold.
        Args:
            time (np.ndarray): Array of time values.
            signal (np.ndarray): Signal data to perform FFT on.
            variable_label (str): Label for the variable being analyzed (e.g., 'u', 'v').
        Returns:
            np.ndarray: Reconstructed signal from the highest frequencies.
    """
    # substract mean to avoid the zero frequency peak
    # hanning window already reduces the mean, so we can skip this step to avoid losing the mean information
    mean = np.mean(signal)
    signal = signal - mean

    # check if dt is constant
    dt_array = np.diff(time)
    if not np.allclose(dt_array, dt_array[0]):
        print(
            colors.WARNING
            + "Warning: time steps are not constant. Taking only the last data with the last dt."
            + colors.ENDC
        )
        dt = dt_array[-1]
        idx = np.where(dt_array == dt)[0][0] + 1
        print(idx)
        time = time[idx:]
        signal = signal[idx:]
        print(
            colors.OKCYAN
            + f"Using dt = {dt} from t = {time[0]} to t = {time[-1]}"
            + colors.ENDC
        )
    else:
        dt = time[1] - time[0]

    print(colors.OKGREEN + f"Using constant dt = {dt}" + colors.ENDC)

    # Apply a Hann window to reduce leakage
    window = np.hanning(len(signal))
    signal = signal * window

    # use fft to get the frequency
    fft_data = rfft(signal)

    freq = rfftfreq(len(signal), d=dt)
    omega = 2 * np.pi * freq

    # eps = 0.1
    # threshold = eps * np.max(
    #     np.abs(fft_data[1:])
    # )  # we skip the mean value (fft_data[0])
    # # fft_data[np.abs(fft_data) < threshold] = 0

    # idx_to_print = []

    # print(f"Frequencies greater than the {eps*100}% of the highest amplitude (excluding mean mode, thus {threshold:.8f}) for variable {variable_label}")
    # for i in range(len(fft_data)):
    #     if np.abs(fft_data[i]) > threshold:
    #         idx_to_print.append([i, np.abs(fft_data[i])])
    # idx_to_print = np.array(idx_to_print)
    #
    # if len(idx_to_print) == 0:
    #     print("No frequencies found above the threshold.")
    #     return

    # # sort by amplitude
    # idx_to_print = idx_to_print[np.argsort(np.abs(idx_to_print[:, 1]))][::-1]

    # for i, A in idx_to_print:
    #     print(
    #         f"ω (2 π f): {omega[int(i)]:.8f}, Amplitude: {A:.8f}"
    #     )

    # print(f"\nMinimum frequency: ω_min = {np.min(omega[1:]):.8f}")
    # print(f"Maximum frequency: ω_max = {np.max(omega):.8f}\n")

    # fig, ax = plt.subplots()
    # ax.plot(omega, np.abs(fft_data), label=f'FFT Amplitude of {variable_label}')
    # ax.set_xlabel('ω (2 π f)')
    # ax.set_ylabel('Amplitude')
    # ax.set_title(f'FFT Amplitude Spectrum for {variable_label}')

    # set mean to zero to avoid the zero frequency peak

    fundamentals = plot_harmonics(omega, fft_data, width, depth, point, variable_label)
    # # manually add back the mean to the fft data
    # fft_data[0] = mean * len(signal)

    # # reconstruct the signal from the highest frequencies
    # signal_reconstructed = irfft(fft_data, n=len(signal))

    return fundamentals


def make_uniform_time(time: np.ndarray, signal: np.ndarray):
    dt_array = np.diff(time)
    if len(dt_array) == 0:
        raise ValueError("Need at least two time samples to compute a PSD.")

    if np.allclose(dt_array, dt_array[0]):
        return time, signal, dt_array[0]

    dt = np.mean(dt_array)
    time_uniform = np.arange(time[0], time[-1] + 0.5 * dt, dt)
    signal_uniform = np.interp(time_uniform, time, signal)
    print(
        colors.WARNING
        + f"Warning: non-uniform time step detected. Interpolating to dt = {dt:g}."
        + colors.ENDC
    )
    return time_uniform, signal_uniform, dt


def compute_psd(time: np.ndarray, signal: np.ndarray):
    time_uniform, signal_uniform, dt = make_uniform_time(time, signal)
    n_samples = len(time_uniform)
    nperseg = min(4096, n_samples)

    freq, psd_freq = welch(
        signal_uniform,
        fs=1.0 / dt,
        window="hann",
        nperseg=nperseg,
        detrend="constant",
        scaling="density",
    )

    omega = 2.0 * np.pi * freq
    psd_omega = psd_freq / (2.0 * np.pi)
    return omega, psd_omega
