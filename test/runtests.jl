using SafeTestsets

@safetestset "Gaussian Line Shape" begin
  include("test_gaussian.jl")
end

@safetestset "Lorentzian Line Shape" begin
  include("test_lorentzian.jl")
end

@safetestset "Voigt Line Shape" begin
  include("test_voigt.jl")
end

@safetestset "VoigtApprx Line Shape" begin
  include("test_voigt_approx.jl")
end

@safetestset "Smith Line Shape" begin
  include("test_smith.jl")
end

@safetestset "DopplerFree Line Shape" begin
  include("test_dopplerfree.jl")
end

@safetestset "Chainrules" begin
  include("test_chainrules.jl")
end

@safetestset "Fitting methods" begin
  include("test_fitting.jl")
end
