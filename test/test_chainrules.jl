using SpectraUtils
using Test
using ChainRulesCore
using ChainRulesTestUtils
using Random

import SpectraUtils: gaussian, lorentzian, voigt

@testset "chainrules" begin
  Random.seed!(1)

  @testset "lineshapes" begin
    test_points = (0., 0.1, 0.5, 1.0, 1.5, 2.5, 10.0,)
    test_widths = (0.01, 0.1, 0.5, 1.0, 1.5, 5., 10.0, )
    for x in test_points, σ in test_widths
      test_frule(gaussian, x, σ)
      test_rrule(gaussian, x, σ)
      test_frule(gaussian, -x, σ)
      test_rrule(gaussian, -x, σ)
      test_frule(lorentzian, x, σ)
      test_rrule(lorentzian, x, σ)
      test_frule(lorentzian, -x, σ)
      test_rrule(lorentzian, -x, σ)
      for γ in test_widths
        test_frule(voigt, x, σ, γ; rtol=1e-3)
        test_rrule(voigt, x, σ, γ; rtol=1e-3)
        test_frule(voigt, -x, σ, γ; rtol=1e-3)
        test_rrule(voigt, -x, σ, γ; rtol=1e-3)
      end
    end
  end

end
