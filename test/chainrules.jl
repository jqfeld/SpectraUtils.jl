using SpectraUtils
using Test
using ChainRulesCore
using ChainRulesTestUtils
using Random

import SpectraUtils: gaussian, lorentzian

@testset "chainrules" begin
  Random.seed!(1)

  @testset "lineshapes" begin
    test_points = (1.5, 2.5, 10.5, 1.6 + 1.6im, 1.6 - 1.6im, 4.6 + 1.6im)
    for x in test_points, y in test_points
      if x isa Real && y isa Real
        # test_scalar(gaussian, x, y)
        test_frule(gaussian, x, y)
        test_rrule(gaussian, x, y)
      end
    end
  end

end
