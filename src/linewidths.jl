const c0 = 3e8  # Speed of light in m/s
const kB = 1.380649e-23  # Boltzmann constant in J/K

@inline sigma_doppler(f0, T; m) = sqrt(kB)/c0 * sqrt(T / m)  * f0
@inline gamma_hard_sphere(p, T; μ, cs) = p/kB/T * cs^2 * sqrt( 8/pi * kB*T/μ )

