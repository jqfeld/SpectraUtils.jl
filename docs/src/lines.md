# Line Models

`Line` objects wrap an amplitude, position, and line shape into a convenient
callable bundle. The stored fields can be constants or functions of an external
parameter set, making it straightforward to evaluate spectra under different
conditions.

## Types and helpers

- [`Line`](@ref SpectraUtils.Line)

## Usage example

```julia
using SpectraUtils

line = Line(1.0, 0.0, Gaussian(0.01))
line(0.0)           # evaluates to the peak amplitude

xs = range(-0.05, 0.05; length=5) |> collect
line(xs)            # returns a Vector with the sampled profile
```

The helper `calc_param` powers the parameter resolution used throughout the
package and can also be employed when building custom line models.
