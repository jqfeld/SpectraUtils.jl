"""
    SpectraUtils

Utilities for constructing and evaluating spectroscopy line models. The
package provides line-shape definitions, convenience helpers for bundling
parameters into reusable `Line` objects, and utilities for estimating
linewidths from physical conditions.
"""
module SpectraUtils



export Lorentzian, Gaussian, Voigt, VoigtApprx, DopplerFree, Smith
include("./lineshapes.jl")

export Line
include("./line.jl")

export sigma_doppler, gamma_hard_sphere
include("./linewidths.jl")

export Spectrum
include("./spectrum.jl")

export fit_spectrum, FitResult
export TuringFit, LevenbergMarquardtFit

using StatsAPI: coef, dof, nobs, rss, stderror, weights, residuals

include("./fitting.jl")

end
