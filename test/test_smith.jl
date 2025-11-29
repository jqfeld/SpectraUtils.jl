using SpectraUtils
using Test
using QuadGK

@testset "Smith Line Shape" begin
    @testset "Basic construction" begin
        s = Smith(0.25, 0.5, 0.3)

        @test s isa Smith{Float64,Float64,Float64}
        @test s.cross_relaxation == 0.25
        @test s.sigma == 0.5
        @test s.gamma == 0.3
    end

    @testset "Normalization" begin
        params = [
            (0.0, 0.4, 0.2),
            (0.35, 0.8, 0.5),
            (0.7, 1.1, 0.9),
        ]

        for (cr, σ, γ) in params
            s = Smith(cr, σ, γ)
            area, _ = quadgk(x -> s(x), -Inf, Inf, rtol=1e-8)
            @test area ≈ 1.0 atol = 1e-6
        end
    end

    @testset "Cross-relaxation limits" begin
        σ = 0.6
        γ = 0.25
        xvals = (-1.0, -0.2, 0.0, 0.3, 0.8)

        smith_gaussian = Smith(1.0, σ, γ)
        g = Gaussian(σ)

        for x in xvals
            @test smith_gaussian(x) ≈ g(x) atol = 1e-12
        end

        smith_lorentz = Smith(0.0, σ, γ)
        @test smith_lorentz(0.0) > smith_gaussian(0.0)
        @test smith_lorentz(γ) > smith_gaussian(γ)
    end

    @testset "Symmetry" begin
        s = Smith(0.3, 0.7, 0.4)
        xs = [0.0, 0.2, 0.6, 1.2]

        for x in xs
            @test s(x) ≈ s(-x) atol = 1e-12
        end
    end
end
