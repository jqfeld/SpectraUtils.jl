using SpectraUtils
using Test
using QuadGK

@testset "Gaussian Line Shape" begin
    @testset "Basic functionality" begin
        g = Gaussian(1.0)
        
        @test g isa Gaussian{Float64}
        @test g.sigma == 1.0
        
        @test g(0.0) ≈ 1 / sqrt(2π) atol=1e-10
        
        @test g(1.0) ≈ 1 / sqrt(2π) * exp(-0.5) atol=1e-10
        @test g(-1.0) ≈ 1 / sqrt(2π) * exp(-0.5) atol=1e-10
    end
    
    @testset "Normalization" begin
        g = Gaussian(1.0)
        area, _ = quadgk(x -> g(x), -Inf, Inf, rtol=1e-8)
        @test area ≈ 1.0 atol=1e-6
        
        g2 = Gaussian(2.0)
        area2, _ = quadgk(x -> g2(x), -Inf, Inf, rtol=1e-8)
        @test area2 ≈ 1.0 atol=1e-6
    end
    
    @testset "Symmetry" begin
        g = Gaussian(1.5)
        test_points = [0.5, 1.0, 2.0, 3.0]
        
        for x in test_points
            @test g(x) ≈ g(-x) atol=1e-12
        end
    end
    
    
    @testset "Different parameter types" begin
        g_int = Gaussian(1)
        g_float = Gaussian(1.0)
        g_rational = Gaussian(1//1)
        
        @test g_int isa Gaussian{Int}
        @test g_float isa Gaussian{Float64}
        @test g_rational isa Gaussian{Rational{Int}}
        
        @test g_int(0.0) ≈ g_float(0.0) atol=1e-12
        @test g_float(0.0) ≈ float(g_rational(0.0)) atol=1e-12
    end
    
    @testset "Width scaling behavior" begin
        narrow = Gaussian(0.5)
        wide = Gaussian(2.0)
        
        @test narrow(0.0) > wide(0.0)
        
        @test narrow(0.25) > wide(0.25)
        @test narrow(1.0) < wide(1.0)
    end
    
    @testset "Mathematical properties" begin
        g = Gaussian(1.0)
        
        @test g(0.0) ≈ 1/sqrt(2π) * exp(0) atol=1e-12
        @test g(1.0) ≈ 1/sqrt(2π) * exp(-0.5) atol=1e-12
        @test g(sqrt(2)) ≈ 1/sqrt(2π) * exp(-1.0) atol=1e-12
        
        inflection_point = 1.0  # σ
        @test g(inflection_point) ≈ 1/sqrt(2π) * exp(-0.5) atol=1e-12
    end
end
