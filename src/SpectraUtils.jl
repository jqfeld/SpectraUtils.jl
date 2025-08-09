module SpectraUtils

calc_param(x::Real, _) = x
calc_param(x::Function, p) = x(p)

export Lorentzian, Gaussian, Voigt, VoigtApprx
include("./lineshapes.jl")

export Line
include("./line.jl")

export sigma_doppler, gamma_hard_sphere
include("./linewidths.jl")


end
