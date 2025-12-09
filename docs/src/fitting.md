# Fitting utilities

SpectraUtils exposes lightweight wrappers for probabilistic and deterministic
spectral fitting. The core API revolves around defining a callable model that
maps a parameter set to a `Spectrum`, then passing that model to one of the
provided fit types.

## Types and entry points

- [`fit_spectrum`](@ref SpectraUtils.fit_spectrum)
- [`FitResult`](@ref SpectraUtils.FitResult)
- [`TuringFit`](@ref SpectraUtils.TuringFit)
- [`LevenbergMarquardtFit`](@ref SpectraUtils.LevenbergMarquardtFit)

`fit_spectrum` dispatches to the appropriate implementation once optional
dependencies are loaded:

- Load `using Turing` to enable Bayesian inference via `TuringFit` (see
  `ext/SpectraUtilsTuringExt.jl`).
- Load `using NonlinearSolve` to enable nonlinear least-squares fitting via
  `LevenbergMarquardtFit` (see `ext/SpectraUtilsNonlinearSolveExt.jl`).

The returned [`FitResult`](@ref SpectraUtils.FitResult) integrates with
`StatsAPI`, exposing `coef`, `stderror`, `rss`, and related accessors for
analysis and reporting.

## Usage examples

### Probabilistic fit with Turing

```julia
using SpectraUtils, Turing

model = params -> Spectrum([
  Line(params.amplitude, params.position, Voigt(0.01, 0.005)),
])

priors = (
  amplitude = Normal(1.0, 0.3),
  position = Normal(0.0, 0.05),
)

fit = TuringFit(model, priors)
xs = range(-0.1, 0.1; length=100) |> collect
ys = model((; amplitude=1.1, position=0.02)).(xs) .+ 0.01 .* randn(length(xs))

result = fit_spectrum(fit, xs, ys; num_samples=500, progress=false)
coef(result)           # posterior means for each parameter
stderror(result)       # posterior standard deviations
```

### Deterministic Levenberg–Marquardt fit

```julia
using SpectraUtils, NonlinearSolve

model = params -> Spectrum([
  Line(params.amplitude, params.position, Gaussian(0.02)),
])
initial_guess = (
  amplitude = 1.0,
  position = 0.0,
)

lm = LevenbergMarquardtFit(model, initial_guess)
xs = range(-0.2, 0.2; length=80) |> collect
ys = model((; amplitude=0.9, position=0.015)).(xs)

result = fit_spectrum(lm, xs, ys)
coef(result)           # fitted parameter means
rss(result)            # residual sum of squares
```

`FitResult` stores residuals and optional weights; these values propagate to the
`StatsAPI` interface so that downstream tooling (e.g., reporting or diagnostics)
can consume a consistent representation.
