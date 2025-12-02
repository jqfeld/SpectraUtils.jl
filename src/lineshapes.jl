using SpecialFunctions

"""Abstract supertype for supported line-shape models."""
abstract type LineShape end


"""
    Lorentzian(hwhm)

Callable representation of a Lorentzian profile with half-width at
half-maximum `hwhm`. Invoke the object as `shape(x)` to evaluate the
profile at offset `x`.
"""
struct Lorentzian{T} <: LineShape
  hwhm::T
end

"""Evaluate the Lorentzian profile with half-width at half-maximum `γ`."""
@inline lorentzian(x, γ) = γ / π / (γ^2 + x^2)
@inline (L::Lorentzian{T})(x) where {T} = lorentzian(x,L.hwhm)


"""
    Gaussian(sigma)

Callable representation of a Gaussian profile with standard deviation
`sigma`. Call `shape(x)` to evaluate the profile at offset `x`.
"""
struct Gaussian{T} <: LineShape
  sigma::T
end

"""Evaluate the Gaussian profile with standard deviation `σ`."""
@inline gaussian(x, σ) = 1 / sqrt(2π) / σ * exp(-x^2 / 2 / σ^2)
@inline (L::Gaussian{T})(x) where {T} = gaussian(x,L.sigma)

"""
    VoigtApprx(sigma, gamma)

Pseudo-Voigt approximation that blends Gaussian and Lorentzian profiles with
standard deviation `sigma` and half-width at half-maximum `gamma`. Evaluate
instances as `shape(x)` to compute the approximated Voigt profile at
offset `x`.
"""
struct VoigtApprx{S,G} <: LineShape
  sigma::S
  gamma::G
end

function (L::VoigtApprx)(x)
  σ = L.sigma
  γ = L.gamma

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
`sigma` and Lorentzian half-width `gamma`. Evaluate using `shape(x)` to
compute the line shape at offset `x`.
"""
@inline (L::Voigt)(x) =
  voigt(x, L.sigma, L.gamma)



"""
    DopplerFree(depth, envelop, dip)

Callable representation of a Doppler-free spectral profile. The resulting
shape multiplies a broad envelope `envelop(x)` by a saturation dip,
returning `envelop(x) * (1 - depth * dip(x))` when evaluated. The `envelop`
and `dip` arguments can be any callable line-shape models (for example,
`Gaussian` or `Lorentzian`).
"""
struct DopplerFree{D,L1,L2} <: LineShape
  depth::D
  envelop::L1
  dip::L2
end

@inline (L::DopplerFree)(x) = 
  L.envelop(x) * (1 - abs(L.depth)*L.dip(x))


"""
    Smith(cross_relaxation, sigma, gamma)

Cross-relaxation–broadened line shape following Smith et al.
([Phys. Rev. Lett. 26, 740](https://doi.org/10.1103/PhysRevLett.26.740)).
The profile smoothly blends a normalized Lorentzian–Gaussian product with a
Gaussian term according to the dimensionless cross-relaxation factor
`cross_relaxation ∈ [0, 1]`:

```
(1 - c) * e^{-(γ/2σ)^2}/erfc(γ/2σ) * (γ/π)/(x^2 + γ^2) * e^{-x^2/4σ^2}
  + c * Gaussian(x, σ)
```

The parameters are the Gaussian standard deviation `sigma`, the Lorentzian
half-width at half-maximum `gamma`, and the mixing factor `cross_relaxation`.
The analytic prefactor ensures the profile integrates to one for any
parameter values.
"""
struct Smith{T1,T2,T3} <: LineShape
  cross_relaxation::T1
  sigma::T2
  gamma::T3
end

@inline (L::Smith)(x) =
  (1 - L.cross_relaxation)*exp(-(L.gamma/2/L.sigma)^2)/erfc(abs(L.gamma/2/L.sigma)) *
  L.gamma/pi /(x^2 + L.gamma^2) * exp(-x^2/4/L.sigma^2) +
  L.cross_relaxation * gaussian(x, L.sigma)
