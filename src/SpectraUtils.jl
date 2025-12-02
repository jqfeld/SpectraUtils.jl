"""
    SpectraUtils

Utilities for constructing and evaluating spectroscopy line models. The
package provides line-shape definitions, convenience helpers for bundling
parameters into reusable `Line` objects, and utilities for estimating
linewidths from physical conditions.
"""
module SpectraUtils

"""
    calc_param(value, params)

Evaluate a parameter specification `value` with the supplied `params`.

If `value` is a number it is returned unchanged. When `value` is a function it
is called with `params` to obtain the numerical parameter value. This helper
allows a `Line` to combine fixed and parameterised quantities.
"""
@inline calc_param(x::Real, _) = x
@inline calc_param(x::Function, p) = x(p)

"""
    NullParameters

A sentinel struct used when no external parameter set is required. Methods that
accept optional parameter bundles default to `NullParameters()`.
"""
struct NullParameters end


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
include("./fitting.jl")

end
