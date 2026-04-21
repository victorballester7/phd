import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, irfft, rfftfreq
from scipy.fft import fft, ifft, fftfreq
from scipy.signal import hilbert


def getTimeFrequency(
    time: np.ndarray,
    signal: np.ndarray,
) -> None:
    # use fft to get the frequency
    fft_data = fft(signal)

    freq = fftfreq(len(signal), d=(time[1] - time[0]))
    omega = 2 * np.pi * freq

    eps = 0.4
    threshold = eps * np.max(
        np.abs(fft_data[1:])
    )  # we skip the mean value (fft_data[0])
    print("Threshold:", threshold)
    # fft_data[np.abs(fft_data) < threshold] = 0

    print("Frequencies with significant amplitude:")
    for i in range(len(fft_data)):
        if np.abs(fft_data[i]) > threshold:
            print(
                f"ω (2 π f): {2 * np.pi * freq[i]:.8f}, Amplitude: {np.abs(fft_data[i]):.8f}"
            )

    # reconstruct the signal from the highest frequencies
    signal_reconstructed = ifft(fft_data, n=len(signal))
    # rescale the reconstructed signal

    plt.plot(time, np.log(signal), label="original signal")
    # plt.plot(time, signal_reconstructed, label="ifft signal", linestyle="--")
    plt.grid()
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("Amplitude")
    plt.show()


def main():
    # create a time vector
    time = np.linspace(0, 1000, 1000)
    w = 0.23
    sigma = 0.004

    x = np.sin(w * time + 2) * np.exp(sigma * time)

    analytic = hilbert(x)
    envelope = np.abs(analytic)
    


    plt.plot(time, np.log(envelope), label="log of envelope")
    sigma, _ = np.polyfit(time, np.log(envelope), 1)

    print(f"Estimated growth rate: {sigma:.6f}")

    getTimeFrequency(time, x)


if __name__ == "__main__":
    main()
