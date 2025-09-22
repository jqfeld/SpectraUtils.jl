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
