# Linewidth Estimates

The `linewidths` helpers offer quick calculations for Doppler and collisional
broadening effects. They return parameters that can be fed directly into the
line-shape constructors.

## API

```@docs
sigma_doppler
gamma_hard_sphere
```

## Example

```julia
using SpectraUtils

σ = sigma_doppler(5.74e14, 300; m=28 * 1.66054e-27)
γ = gamma_hard_sphere(50, 300; μ=14 * 1.66054e-27, cs=450)
Voigt(σ, γ)
```

The formulas are intentionally lightweight and are best suited for quick
estimates or as building blocks inside more sophisticated models.
