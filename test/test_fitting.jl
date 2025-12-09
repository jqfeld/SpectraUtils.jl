using Test
using Turing, NonlinearSolve

using SpectraUtils
using Distributions
using Random; Random.seed!(42)


line = Line(0., 1., Gaussian(1.))

spec = Spectrum((line,), x -> 1.)


xs = range(-5, 5, length=1000)
data = spec(xs) .+ rand(Normal(0, 0.1), length(xs))



##

spectrum(p) = Spectrum(
  (Line(0, p.ampl, Gaussian(1)),), 
  x -> p.bg
)

lmfit = LevenbergMarquardtFit(spectrum, (ampl = 1, bg=1))
lmres = fit_spectrum(lmfit, xs, data)


turingfit = TuringFit(spectrum, (ampl = Normal(1, 1), bg=Normal(1,1) ))
turingres = fit_spectrum(turingfit, xs, data; num_chains=5)


@test mean(lmres[:ampl]) ≈ mean(turingres[:ampl]) rtol=1e-3
@test mean(lmres[:bg]) ≈ mean(turingres[:bg]) rtol=1e-3
