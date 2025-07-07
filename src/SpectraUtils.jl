module SpectraUtils


export Lorentzian, Gaussian, Voigt, VoigtApprx, lineparam
include("./lineshapes.jl")

export sigma_doppler, gamma_hard_sphere
include("./linewidths.jl")


end
