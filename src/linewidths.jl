const c0 = 3e8  # Speed of light in m/s
const kB = 1.380649e-23  # Boltzmann constant in J/K

"""
    sigma_doppler(f0, T; m)

Return the Doppler Gaussian standard deviation for a transition at frequency
`f0` and temperature `T` (in kelvin). The molecular mass `m` is supplied as a
keyword argument and should be expressed in kilograms.
"""
@inline sigma_doppler(f0, T; m) = sqrt(kB)/c0 * sqrt(T / m)  * f0

"""
    gamma_hard_sphere(p, T; μ, cs)

Compute the Lorentzian half-width at half-maximum (in cm⁻¹) under the
hard-sphere collisional model. The calculation uses the pressure `p` (Pa),
temperature `T` (K), reduced mass `μ` (kg), and average relative speed `cs`
(m/s).
"""
@inline gamma_hard_sphere(p, T; μ, cs) = p/kB/T * cs^2 * sqrt( 8/pi * kB*T/μ ) / c0/100 # units of cm^-1

