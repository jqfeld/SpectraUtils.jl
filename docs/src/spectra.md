# Spectra

`Spectrum` bundles a collection of `Line` objects with an optional background
function to produce complete spectral traces. Calling the resulting object
returns the combined background and line contributions at a single position or
broadcasts across an array of positions.

## Types and constructors

- [`Spectrum`](@ref SpectraUtils.Spectrum)

A spectrum can be built from any iterable of lines. Pass a callable background
as the second argument to account for offsets or slow drifts, or omit it to
model only the summed line profiles.

## Usage example

```julia
using SpectraUtils

line1 = Line(0.8, -0.05, Gaussian(0.02))
line2 = Line(1.3, 0.04, Lorentzian(0.01))
background = x -> 0.01 + 0.02x

spec = Spectrum([line1, line2], background)
spec(0.0)               # evaluate at a single point

xs = range(-0.2, 0.2; length=200) |> collect
spec(xs)                # vectorized evaluation
```

When no background is provided, `Spectrum(lines)` defaults to summing the line
profiles alone. The `Spectrum` call overloads optimize scalar and vector cases
so that repeated evaluations during fitting or simulation stay efficient.
