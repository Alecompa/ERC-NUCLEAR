import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Constants
density_nai = 3.67  # Density of NaI (g/cm³)
H = 7.62  # Detector thickness (cm)
b = 1     # Source and detector separation (cm)
m = 0     # Offset between source and detector center (cm)
R = 3.72  # Radius of detector (cm)
g = 1     # Source Radius (cm)

# Energy and Total attenuation coefficient (without coherent scattering) pairs. From NIST: https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
energy_tau_pairs = [
    (1.000E-03, 7.794E+03), (1.072E-03, 6.736E+03), (1.072E-03, 7.021E+03),
    (1.072E-03, 7.925E+03), (1.072E-03, 7.021E+03), (1.072E-03, 7.924E+03),
    (1.500E-03, 3.801E+03), (2.000E-03, 1.917E+03), (3.000E-03, 7.003E+02),
    (4.000E-03, 3.352E+02), (4.557E-03, 2.387E+02), (4.557E-03, 6.585E+02),
    (4.702E-03, 6.166E+02), (4.852E-03, 5.775E+02), (4.852E-03, 7.719E+02),
    (5.000E-03, 7.277E+02), (5.188E-03, 6.612E+02), (5.188E-03, 7.604E+02),
    (6.000E-03, 5.298E+02), (8.000E-03, 2.489E+02), (1.000E-02, 1.376E+02),
    (1.500E-02, 4.578E+01), (2.000E-02, 2.071E+01), (3.000E-02, 6.714E+00),
    (3.317E-02, 5.081E+00), (3.317E-02, 2.987E+01), (4.000E-02, 1.835E+01),
    (5.000E-02, 1.018E+01), (6.000E-02, 6.228E+00), (8.000E-02, 2.863E+00),
    (1.000E-01, 1.576E+00), (1.500E-01, 5.663E-01), (2.000E-01, 3.020E-01),
    (3.000E-01, 1.534E-01), (4.000E-01, 1.100E-01), (5.000E-01, 9.035E-02),
    (6.000E-01, 7.901E-02), (8.000E-01, 6.571E-02), (1.000E+00, 5.762E-02),
    (1.022E+00, 5.687E-02), (1.250E+00, 5.086E-02), (1.500E+00, 4.644E-02),
    (2.000E+00, 4.119E-02), (2.044E+00, 4.087E-02), (3.000E+00, 3.668E-02),
    (4.000E+00, 3.512E-02), (5.000E+00, 3.472E-02), (6.000E+00, 3.484E-02),
    (7.000E+00, 3.526E-02), (8.000E+00, 3.584E-02), (9.000E+00, 3.650E-02),
    (1.000E+01, 3.722E-02)
]

def mu_1(H, b, s):
    return (H + b) / np.sqrt((H + b)**2 + s**2)

def mu_2(b, s):
    return b / np.sqrt(b**2 + s**2)

def F1(mu, T, H):
    return 1 - np.exp(-T * H / mu)

def F2(mu, phi, T, b, m, R):
    s = -m * np.sin(phi) + np.sqrt(R**2 - m**2 * np.cos(phi)**2)
    return 1 - np.exp(-T * (s / np.sqrt(1 - mu**2) - b / mu))

def integrand(phi, T, H, b, m, R):
    s = -m * np.sin(phi) + np.sqrt(R**2 - m**2 * np.cos(phi)**2)
    mu1 = mu_1(H, b, s)
    mu2 = mu_2(b, s)
    
    inner_integral_1, _ = integrate.quad(F1, mu1, 1, args=(T, H))
    inner_integral_2, _ = integrate.quad(F2, mu2, mu1, args=(phi, T, b, m, R))
    
    return inner_integral_1 + inner_integral_2

def epsilon_AT(T, H, b, m, R):
    outer_integral, _ = integrate.quad(integrand, -np.pi/2, np.pi/2, args=(T, H, b, m, R))
    return outer_integral / (2 * np.pi)

def epsilon_AT_disk(T, H, b, m, R, g):
    def disk_integrand(m):
        return m * epsilon_AT(T, H, b, m, R)
    integral_result, _ = integrate.quad(disk_integrand, 0, g)
    return (2 / g**2) * integral_result


# Convert to linear attenuation coefficients and compute detector efficiency
Egamma = []
linear_attenuation_coeffs = []
efficiencies = []
i = 0
for energy, tau in energy_tau_pairs:
    mu = tau * density_nai  # Convert to 1/cm
    detector_efficiency = epsilon_AT_disk(mu, H, b, m, R, g)
    Egamma.append(energy)
    linear_attenuation_coeffs.append(mu)
    efficiencies.append(detector_efficiency)
    print(f"#: {i} , For Energy: {energy:.4f} MeV, mu: {mu:.3f} 1/cm")
    print(f"   Detector Efficiency: {detector_efficiency*100:.3f} %")
    i = i+1
# Plot efficiency vs energy
plt.figure(figsize=(8, 6))
plt.plot(Egamma, efficiencies, marker='o', linestyle='-', color='b')
plt.xlabel('Energy (MeV)')
plt.ylabel('Detector Total Efficiency')
plt.title('Detector Total Efficiency vs Gamma Energy')
plt.xscale('log')
plt.grid(True)
plt.show()