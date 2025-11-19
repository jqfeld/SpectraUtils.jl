module SpectraUtilsMonteCarloMeasurementsExt

using SpectraUtils, MonteCarloMeasurements
import SpectraUtils: gaussian, lorentzian, voigt
import SpecialFunctions: faddeeva

@inline voigt(x, σ::P, γ) where {P <: MonteCarloMeasurements.AbstractParticles} = real(MonteCarloMeasurements.ℂ2ℂ_function(faddeeva, (x + im * γ) / σ / sqrt(2))) / σ / sqrt(2π)
@inline voigt(x, σ, γ::P) where {P <: MonteCarloMeasurements.AbstractParticles} = real(MonteCarloMeasurements.ℂ2ℂ_function(faddeeva, (x + im * γ) / σ / sqrt(2))) / σ / sqrt(2π)
@inline voigt(x, σ::P1, γ::P2) where {P1 <: MonteCarloMeasurements.AbstractParticles, P2 <: MonteCarloMeasurements.AbstractParticles} = 
  real(MonteCarloMeasurements.ℂ2ℂ_function(faddeeva, (x + im * γ) / σ / sqrt(2))) / σ / sqrt(2π)

end
