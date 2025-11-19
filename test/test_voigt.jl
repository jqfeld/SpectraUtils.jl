using SpectraUtils
using Test
using QuadGK
using SpecialFunctions

@testset "Voigt Line Shape" begin
  @testset "Basic functionality" begin
    V = Voigt(1.0, 0.5)

    @test V isa Voigt{Float64,Float64}
    @test V.sigma == 1.0
    @test V.gamma == 0.5

    @test V(0.0) isa Real
    @test V(0.0) > 0
  end

  @testset "Normalization" begin
    V = Voigt(1.0, 0.5)
    area, _ = quadgk(x -> V(x), -Inf, Inf, rtol=1e-8)
    @test area ≈ 1.0 atol = 1e-6

    V2 = Voigt(2.0, 1.0)
    area2, _ = quadgk(x -> V2(x), -Inf, Inf, rtol=1e-8)
    @test area2 ≈ 1.0 atol = 1e-6
  end

  @testset "Symmetry" begin
    V = Voigt(1.5, 0.8)
    test_points = [0.5, 1.0, 2.0, 3.0]

    for x in test_points
      @test V(x) ≈ V(-x) atol = 1e-12
    end
  end


  @testset "Different parameter types" begin
    V_mixed = Voigt(1.0, 1)
    V_float = Voigt(1.0, 1.0)

    @test V_mixed isa Voigt{Float64,Int}
    @test V_float isa Voigt{Float64,Float64}

    @test V_mixed(0.0) ≈ V_float(0.0) atol = 1e-12
  end

  @testset "Limiting cases" begin
    # When gamma -> 0, should approach Gaussian
    V_gaussian_like = Voigt(1.0, 1e-6)
    G = Gaussian(1.0)

    @test V_gaussian_like(0.0) ≈ G(0.0) rtol = 1e-3
    @test V_gaussian_like(1.0) ≈ G(1.0) rtol = 1e-3

    # When sigma -> 0, should approach Lorentzian
    V_lorentzian_like = Voigt(1e-6, 1.0)
    L = Lorentzian(1.0)

    @test V_lorentzian_like(0.0) ≈ L(0.0) rtol = 1e-3
    @test V_lorentzian_like(1.0) ≈ L(1.0) rtol = 1e-3
  end

  @testset "Mathematical properties" begin
    V = Voigt(1.0, 0.5)

    @test V(0.0) > V(1.0)
    @test V(1.0) > V(2.0)

    @test all(V(x) > 0 for x in [-5, -1, 0, 1, 5])

    center_val = V(0.0)
    @test center_val isa Real
    @test center_val > 0
  end

  @testset "Faddeeva function usage" begin
    # Test that the Voigt profile uses the Faddeeva function correctly
    σs = [0.1, 0.5, 1.0, 5.0, 10]
    γs = [0.1, 0.5, 1.0, 5.0, 10,]
    xs = [-1.0, 0.0, 1.0, 10.0, 100.0,]

    for σ in σs, γ in γs
      for x in xs
        V = Voigt(σ, γ)
        expected = real(faddeeva((x + im * γ) / σ / sqrt(2))) / σ / sqrt(2π)

        @test V(x) ≈ expected rtol = 1e-5
      end
    end

  end
end
