using SpecialFunctions

"""Abstract supertype for supported line-shape models."""
abstract type LineShape end

"""Evaluate the Lorentzian profile with half-width at half-maximum `γ`."""
@inline lorentzian(x, γ) = γ / π / (γ^2 + x^2)

"""
    Lorentzian(hwhm)

Callable representation of a Lorentzian profile with half-width at
half-maximum `hwhm`. Invoke the object as `shape(x, p)` to evaluate the
profile at offset `x`, optionally resolving `hwhm` from a parameter set `p`.
"""
struct Lorentzian{T} <: LineShape
  hwhm::T
end

@inline (L::Lorentzian{T})(x, p=NullParameters()) where {T} =
  lorentzian(x, calc_param(L.hwhm, p))

"""Evaluate the Gaussian profile with standard deviation `σ`."""
@inline gaussian(x, σ) = 1 / sqrt(2π) / σ * exp(-x^2 / 2 / σ^2)

"""
    Gaussian(sigma)

Callable representation of a Gaussian profile with standard deviation
`sigma`. Call `shape(x, p)` to evaluate the profile at offset `x`, resolving
`sigma` from `p` when it is provided as a function.
"""
struct Gaussian{T} <: LineShape
  sigma::T
end

@inline (L::Gaussian{T})(x, p=NullParameters()) where {T} =
  gaussian(x, calc_param(L.sigma, p))

"""
    VoigtApprx(sigma, gamma)

Pseudo-Voigt approximation that blends Gaussian and Lorentzian profiles with
standard deviation `sigma` and half-width at half-maximum `gamma`. Evaluate
instances as `shape(x, p)` to compute the approximated Voigt profile at
offset `x`.
"""
struct VoigtApprx{S,G} <: LineShape
  sigma::S
  gamma::G
end

function (L::VoigtApprx)(x, p=NullParameters())
  σ = calc_param(L.sigma, p)
  γ = calc_param(L.gamma, p)

  fG = 2 * sqrt(2 * log(2)) * σ
  fL = 2 * γ
  f = (max(0.0, fG^5 + 2.69269 * fG^4 * fL + 2.42843 * fG^3 * fL^2 + 4.47163 * fG^2 * fL^3 + 0.07842 * fG * fL^4 + fL^5))^(1 / 5)
  η = 1.36603 * (fL / f) - 0.47719 * (fL / f)^2 + 0.11116 * (fL / f)^3
  return η * lorentzian(x, f / 2) + (1 - η) * gaussian(x, f / (2 * sqrt(2 * log(2))))
end


struct Voigt{T1,T2} <: LineShape
  sigma::T1
  gamma::T2
end

"""Evaluate the full Voigt profile for standard deviation `σ` and `γ`."""
@inline voigt(x, σ, γ) = real(faddeeva((x + im * γ) / σ / sqrt(2))) / σ / sqrt(2π)

"""
    Voigt(sigma, gamma)

Callable representation of the full Voigt profile with standard deviation
`sigma` and Lorentzian half-width `gamma`. Evaluate using `shape(x, p)` to
compute the line shape at offset `x`, resolving stored parameters from `p`
when they are provided as callables.
"""
@inline (L::Voigt)(x, p=NullParameters()) =
  voigt(x, calc_param(L.sigma, p), calc_param(L.gamma, p))
