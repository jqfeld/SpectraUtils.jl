# Line Shapes

SpectraUtils.jl provides several callable types that model common spectral
profiles. Each line shape can be constructed with either numeric values or
functions that resolve parameters at evaluation time.

## Analytic profiles

- [`LineShape`](@ref SpectraUtils.LineShape)
- [`Lorentzian`](@ref SpectraUtils.Lorentzian)
- [`Gaussian`](@ref SpectraUtils.Gaussian)
- [`Voigt`](@ref SpectraUtils.Voigt)
- [`VoigtApprx`](@ref SpectraUtils.VoigtApprx)
- [`DopplerFree`](@ref SpectraUtils.DopplerFree)
- [`Smith`](@ref SpectraUtils.Smith)
- [`lorentzian`](@ref SpectraUtils.lorentzian)
- [`gaussian`](@ref SpectraUtils.gaussian)
- [`voigt`](@ref SpectraUtils.voigt)

## Example

```julia
using SpectraUtils

shape = Voigt(0.01, 0.002)
shape(0.0)              # peak amplitude at line centre
shape(0.03, NullParameters())
```

Passing `NullParameters()` is optional when the stored parameters are constants,
but it highlights how the same shape can be reused with external parameter
bundles.

## Doppler-free composites

`DopplerFree` combines a broad envelope with a narrow saturation dip to model
classic Doppler-free spectroscopy traces. Any callable line shape can serve as
the envelope or dip, allowing flexible compositions such as a Lorentzian
envelope with a Gaussian dip. For example:

```julia
using SpectraUtils

envelop = Lorentzian(0.5)
dip = Gaussian(0.05)

shape = DopplerFree(0.2, envelop, dip)
shape(0.0)         # central dip on top of the envelope
shape(0.1)
```

## Smith line profile

The `Smith` profile models cross-relaxation broadening following the
formulation of Smith et al. (Phys. Rev. Lett. 26, 740). It blends a Lorentzian
envelope with a Doppler-broadened Gaussian core while enforcing unit
normalization so that the total integrated area is one. The `cross_relaxation`
parameter controls the relative weight of the Gaussian component, with
`cross_relaxation = 0` yielding the pure Lorentzian-derived term and
`cross_relaxation = 1` reducing to the Gaussian profile `gaussian(x, sigma)`.
