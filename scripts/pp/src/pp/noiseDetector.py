import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks, hilbert
from scipy.stats import variation
import pandas as pd

class NoiseDetector:
    """
    Clase para detectar automáticamente el inicio de oscilaciones en señales CFD
    """
    
    def __init__(self, time, signal_data):
        """
        Inicializar con los datos de tiempo y señal
        
        Parameters:
        -----------
        time : array-like
            Vector de tiempo
        signal_data : array-like
            Datos de la señal (puede ser presión, velocidad, etc.)
        """
        self.time = np.array(time)
        self.signal = np.array(signal_data)
        self.dt = np.mean(np.diff(self.time))
        self.fs = 1.0 / self.dt  # Frecuencia de muestreo
        
    def method_variance_threshold(self, window_size=100, threshold_factor=3.0):
        """
        Método 1: Detección basada en varianza móvil
        
        Parameters:
        -----------
        window_size : int
            Tamaño de la ventana móvil
        threshold_factor : float
            Factor multiplicativo del umbral de varianza
        
        Returns:
        --------
        t_0 : float
            Tiempo de inicio de oscilaciones
        """
        # Calcular varianza móvil
        variance = np.array([np.var(self.signal[max(0, i-window_size):i+1]) 
                           for i in range(len(self.signal))])
        
        # Varianza base (primeros puntos, asumiendo régimen estacionario)
        base_variance = np.mean(variance[:window_size])
        threshold = base_variance * threshold_factor
        
        # Encontrar primer punto donde se supera el umbral
        oscillation_start = np.where(variance > threshold)[0]
        
        if len(oscillation_start) > 0:
            idx_start = oscillation_start[0]
            t_0 = self.time[idx_start]
        else:
            t_0 = self.time[0]
            
        return t_0, variance, threshold
    
    def method_derivative_analysis(self, window_size=50, threshold_factor=2.0):
        """
        Método 2: Análisis de derivadas para detectar cambios abruptos
        
        Parameters:
        -----------
        window_size : int
            Tamaño de ventana para suavizado
        threshold_factor : float
            Factor para el umbral de derivada
        
        Returns:
        --------
        t_0 : float
            Tiempo de inicio de oscilaciones
        """
        # Suavizar la señal
        smoothed = signal.savgol_filter(self.signal, window_size, 3)
        
        # Calcular derivada
        derivative = np.gradient(smoothed, self.dt)
        
        # Derivada móvil RMS
        derivative_rms = np.array([np.sqrt(np.mean(derivative[max(0, i-window_size):i+1]**2)) 
                                 for i in range(len(derivative))])
        
        # Umbral basado en los primeros valores
        base_derivative = np.mean(derivative_rms[:window_size])
        threshold = base_derivative * threshold_factor
        
        # Encontrar inicio de oscilaciones
        oscillation_start = np.where(derivative_rms > threshold)[0]
        
        if len(oscillation_start) > 0:
            idx_start = oscillation_start[0]
            t_0 = self.time[idx_start]
        else:
            t_0 = self.time[0]
            
        return t_0, derivative_rms, threshold
    
    def method_spectral_analysis(self, window_size=200, overlap=0.5, freq_threshold=0.1):
        """
        Método 3: Análisis espectral mediante ventanas deslizantes
        
        Parameters:
        -----------
        window_size : int
            Tamaño de ventana para FFT
        overlap : float
            Solapamiento entre ventanas (0-1)
        freq_threshold : float
            Umbral de energía espectral normalizada
        
        Returns:
        --------
        t_0 : float
            Tiempo de inicio de oscilaciones
        """
        step = int(window_size * (1 - overlap))
        spectral_energy = []
        time_windows = []
        
        for i in range(0, len(self.signal) - window_size, step):
            window_signal = self.signal[i:i + window_size]
            window_time = self.time[i + window_size//2]
            
            # FFT de la ventana
            fft = np.fft.fft(window_signal)
            freqs = np.fft.fftfreq(window_size, self.dt)
            
            # Energía en frecuencias relevantes (excluyendo DC)
            relevant_idx = (freqs > 0) & (freqs < self.fs/4)  # Hasta frecuencia de Nyquist/2
            energy = np.sum(np.abs(fft[relevant_idx])**2)
            
            spectral_energy.append(energy)
            time_windows.append(window_time)
        
        spectral_energy = np.array(spectral_energy)
        time_windows = np.array(time_windows)
        
        # Normalizar energía
        if len(spectral_energy) > 0:
            normalized_energy = spectral_energy / np.max(spectral_energy)
            
            # Encontrar cuando la energía supera el umbral
            oscillation_idx = np.where(normalized_energy > freq_threshold)[0]
            
            if len(oscillation_idx) > 0:
                t_0 = time_windows[oscillation_idx[0]]
            else:
                t_0 = self.time[0]
        else:
            t_0 = self.time[0]
            
        return t_0, time_windows, normalized_energy, freq_threshold
    
    def method_envelope_analysis(self, window_size=100, threshold_factor=2.0):
        """
        Método 4: Análisis de envolvente usando transformada de Hilbert
        
        Parameters:
        -----------
        window_size : int
            Tamaño de ventana móvil
        threshold_factor : float
            Factor para el umbral de envolvente
        
        Returns:
        --------
        t_0 : float
            Tiempo de inicio de oscilaciones
        """
        # Obtener envolvente usando Hilbert
        analytic_signal = hilbert(self.signal)
        envelope = np.abs(analytic_signal)
        
        # Variación de la envolvente
        envelope_variation = np.array([variation(envelope[max(0, i-window_size):i+1]) 
                                     for i in range(len(envelope))])
        
        # Umbral basado en variación inicial
        base_variation = np.mean(envelope_variation[:window_size])
        threshold = base_variation * threshold_factor
        
        # Detectar inicio de oscilaciones
        oscillation_start = np.where(envelope_variation > threshold)[0]
        
        if len(oscillation_start) > 0:
            idx_start = oscillation_start[0]
            t_0 = self.time[idx_start]
        else:
            t_0 = self.time[0]
            
        return t_0, envelope_variation, threshold
    
    def detect_t0_robust(self, methods=['variance', 'derivative', 'spectral', 'envelope']):
        """
        Método robusto que combina múltiples técnicas
        
        Parameters:
        -----------
        methods : list
            Lista de métodos a usar
        
        Returns:
        --------
        t_0_final : float
            Tiempo de inicio robusto
        results : dict
            Resultados de todos los métodos
        """
        results = {}
        t0_values = []
        
        if 'variance' in methods:
            t0_var, var_data, var_thresh = self.method_variance_threshold()
            results['variance'] = {'t0': t0_var, 'data': var_data, 'threshold': var_thresh}
            t0_values.append(t0_var)
        
        if 'derivative' in methods:
            t0_der, der_data, der_thresh = self.method_derivative_analysis()
            results['derivative'] = {'t0': t0_der, 'data': der_data, 'threshold': der_thresh}
            t0_values.append(t0_der)
        
        if 'spectral' in methods:
            t0_spec, time_win, spec_data, spec_thresh = self.method_spectral_analysis()
            results['spectral'] = {'t0': t0_spec, 'time_windows': time_win, 
                                 'data': spec_data, 'threshold': spec_thresh}
            t0_values.append(t0_spec)
        
        if 'envelope' in methods:
            t0_env, env_data, env_thresh = self.method_envelope_analysis()
            results['envelope'] = {'t0': t0_env, 'data': env_data, 'threshold': env_thresh}
            t0_values.append(t0_env)
        
        # Tiempo robusto: mediana de todos los métodos
        t_0_final = np.median(t0_values)
        
        results['final_t0'] = t_0_final
        results['all_t0_values'] = t0_values
        
        return t_0_final, results
    
    def filter_data(self, t_0):
        """
        Filtrar datos desde t_0 hasta el final
        
        Parameters:
        -----------
        t_0 : float
            Tiempo de inicio
        
        Returns:
        --------
        time_filtered : array
            Tiempo filtrado
        signal_filtered : array
            Señal filtrada
        """
        mask = self.time >= t_0
        time_filtered = self.time[mask]
        signal_filtered = self.signal[mask]
        
        return time_filtered, signal_filtered
    
    def plot_analysis(self, results, figsize=(15, 12)):
        """
        Visualizar todos los análisis
        """
        fig, axes = plt.subplots(3, 2, figsize=figsize)
        fig.suptitle('Análisis de Detección de Oscilaciones', fontsize=16)
        
        # Señal original
        axes[0, 0].plot(self.time, self.signal, 'b-', alpha=0.7, label='Señal original')
        if 'final_t0' in results:
            axes[0, 0].axvline(results['final_t0'], color='r', linestyle='--', 
                              label=f't₀ = {results["final_t0"]:.4f}')
        axes[0, 0].set_xlabel('Tiempo')
        axes[0, 0].set_ylabel('Amplitud')
        axes[0, 0].set_title('Señal Original')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Método de varianza
        if 'variance' in results:
            var_res = results['variance']
            axes[0, 1].plot(self.time, var_res['data'], 'g-', label='Varianza móvil')
            axes[0, 1].axhline(var_res['threshold'], color='r', linestyle=':', 
                              label='Umbral')
            axes[0, 1].axvline(var_res['t0'], color='r', linestyle='--', 
                              label=f't₀ = {var_res["t0"]:.4f}')
            axes[0, 1].set_xlabel('Tiempo')
            axes[0, 1].set_ylabel('Varianza')
            axes[0, 1].set_title('Método: Varianza Móvil')
            axes[0, 1].legend()
            axes[0, 1].grid(True, alpha=0.3)
        
        # Método de derivadas
        if 'derivative' in results:
            der_res = results['derivative']
            axes[1, 0].plot(self.time, der_res['data'], 'm-', label='RMS Derivada')
            axes[1, 0].axhline(der_res['threshold'], color='r', linestyle=':', 
                              label='Umbral')
            axes[1, 0].axvline(der_res['t0'], color='r', linestyle='--', 
                              label=f't₀ = {der_res["t0"]:.4f}')
            axes[1, 0].set_xlabel('Tiempo')
            axes[1, 0].set_ylabel('RMS Derivada')
            axes[1, 0].set_title('Método: Análisis de Derivadas')
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        
        # Método espectral
        if 'spectral' in results:
            spec_res = results['spectral']
            axes[1, 1].plot(spec_res['time_windows'], spec_res['data'], 'c-', 
                           label='Energía espectral')
            axes[1, 1].axhline(spec_res['threshold'], color='r', linestyle=':', 
                              label='Umbral')
            axes[1, 1].axvline(spec_res['t0'], color='r', linestyle='--', 
                              label=f't₀ = {spec_res["t0"]:.4f}')
            axes[1, 1].set_xlabel('Tiempo')
            axes[1, 1].set_ylabel('Energía normalizada')
            axes[1, 1].set_title('Método: Análisis Espectral')
            axes[1, 1].legend()
            axes[1, 1].grid(True, alpha=0.3)
        
        # Método de envolvente
        if 'envelope' in results:
            env_res = results['envelope']
            axes[2, 0].plot(self.time, env_res['data'], 'orange', 
                           label='Variación envolvente')
            axes[2, 0].axhline(env_res['threshold'], color='r', linestyle=':', 
                              label='Umbral')
            axes[2, 0].axvline(env_res['t0'], color='r', linestyle='--', 
                              label=f't₀ = {env_res["t0"]:.4f}')
            axes[2, 0].set_xlabel('Tiempo')
            axes[2, 0].set_ylabel('Variación')
            axes[2, 0].set_title('Método: Análisis de Envolvente')
            axes[2, 0].legend()
            axes[2, 0].grid(True, alpha=0.3)
        
        # Comparación de métodos
        if 'all_t0_values' in results:
            method_names = [k for k in results.keys() if k not in ['final_t0', 'all_t0_values']]
            t0_values = results['all_t0_values']
            
            axes[2, 1].bar(method_names, t0_values, alpha=0.7)
            axes[2, 1].axhline(results['final_t0'], color='r', linestyle='--', 
                              label=f'Mediana: {results["final_t0"]:.4f}')
            axes[2, 1].set_ylabel('t₀')
            axes[2, 1].set_title('Comparación de Métodos')
            axes[2, 1].legend()
            axes[2, 1].grid(True, alpha=0.3)
            plt.setp(axes[2, 1].xaxis.get_majorticklabels(), rotation=45)
        
        plt.tight_layout()
        return fig

# Ejemplo de uso
def ejemplo_uso(t: np.ndarray, signal: np.ndarray):
    """
    Ejemplo completo de uso con señal sintética
    """
    # Señal total
    signal_total = signal
    
    print("=== ANÁLISIS AUTOMÁTICO DE DETECCIÓN DE OSCILACIONES ===")
    print()
    
    # Crear detector
    detector = NoiseDetector(t, signal_total)
    
    # Detectar t_0 usando método robusto
    t_0_detected, results = detector.detect_t0_robust()
    
    print(f"Tiempo detectado automáticamente: {t_0_detected:.4f}")
    print()
    
    # Mostrar resultados de cada método
    print("Resultados por método:")
    for method, result in results.items():
        if method not in ['final_t0', 'all_t0_values']:
            print(f"  {method}: t₀ = {result['t0']:.4f}")
    
    # Filtrar datos
    time_filtered, signal_filtered = detector.filter_data(t_0_detected)
    
    print(f"\nDatos originales: {len(t)} puntos")
    print(f"Datos filtrados: {len(time_filtered)} puntos")
    print(f"Porcentaje de datos útiles: {len(time_filtered)/len(t)*100:.1f}%")
    
    # Visualizar
    fig = detector.plot_analysis(results)
    plt.show()
    
    return detector, results, time_filtered, signal_filtered

# Función para aplicar a datos reales
def procesar_datos_cfd(archivo_datos, columna_tiempo='time', columna_senal='pressure'):
    """
    Función para procesar datos CFD reales desde archivo
    
    Parameters:
    -----------
    archivo_datos : str
        Ruta al archivo de datos (CSV, Excel, etc.)
    columna_tiempo : str
        Nombre de la columna de tiempo
    columna_senal : str
        Nombre de la columna de la señal a analizar
    
    Returns:
    --------
    detector : NoiseDetector
        Objeto detector configurado
    t_0 : float
        Tiempo de inicio detectado
    time_filtered : array
        Datos de tiempo filtrados
    signal_filtered : array
        Datos de señal filtrados
    """
    # Leer datos
    if archivo_datos.endswith('.csv'):
        df = pd.read_csv(archivo_datos)
    elif archivo_datos.endswith('.xlsx'):
        df = pd.read_excel(archivo_datos)
    else:
        raise ValueError("Formato de archivo no soportado. Use CSV o Excel.")
    
    # Extraer datos
    time_data = df[columna_tiempo].values
    signal_data = df[columna_senal].values
    
    # Crear detector y procesar
    detector = NoiseDetector(time_data, signal_data)
    t_0, results = detector.detect_t0_robust()
    
    # Filtrar datos
    time_filtered, signal_filtered = detector.filter_data(t_0)
    
    # Mostrar análisis
    detector.plot_analysis(results)
    plt.show()
    
    print(f"Tiempo de inicio detectado: {t_0:.6f}")
    print(f"Datos filtrados: {len(time_filtered)} de {len(time_data)} puntos")
    
    return detector, t_0, time_filtered, signal_filtered
