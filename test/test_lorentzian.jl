using SpectraUtils
using Test
using QuadGK

@testset "Lorentzian Line Shape" begin
    @testset "Basic functionality" begin
        L = Lorentzian(1.0)
        
        @test L isa Lorentzian{Float64}
        @test L.hwhm == 1.0
        
        @test L(0.0) ≈ 1 / π atol=1e-10
        
        @test L(1.0) ≈ 1 / π / 2 atol=1e-10
        @test L(-1.0) ≈ 1 / π / 2 atol=1e-10
    end
    
    @testset "Normalization" begin
        L = Lorentzian(1.0)
        area, _ = quadgk(x -> L(x), -Inf, Inf, rtol=1e-8)
        @test area ≈ 1.0 atol=1e-6
        
        L2 = Lorentzian(2.0)
        area2, _ = quadgk(x -> L2(x), -Inf, Inf, rtol=1e-8)
        @test area2 ≈ 1.0 atol=1e-6
    end
    
    @testset "Symmetry" begin
        L = Lorentzian(1.5)
        test_points = [0.5, 1.0, 2.0, 3.0]
        
        for x in test_points
            @test L(x) ≈ L(-x) atol=1e-12
        end
    end
    
    @testset "Different parameter types" begin
        L_int = Lorentzian(1)
        L_float = Lorentzian(1.0)
        L_rational = Lorentzian(1//1)
        
        @test L_int isa Lorentzian{Int}
        @test L_float isa Lorentzian{Float64}
        @test L_rational isa Lorentzian{Rational{Int}}
        
        @test L_int(0.0) ≈ L_float(0.0) atol=1e-12
        @test L_float(0.0) ≈ float(L_rational(0.0)) atol=1e-12
    end
    
    @testset "Width scaling behavior" begin
        narrow = Lorentzian(0.5)
        wide = Lorentzian(2.0)
        
        @test narrow(0.0) > wide(0.0)
        
        @test narrow(0.25) > wide(0.25)
        @test narrow(2.0) < wide(2.0)
    end
    
    @testset "Mathematical properties" begin
        L = Lorentzian(1.0)
        
        @test L(0.0) ≈ 1/π atol=1e-12
        @test L(1.0) ≈ 1/π/2 atol=1e-12
        @test L(2.0) ≈ 1/π/5 atol=1e-12
        
        hwhm = 1.0
        @test L(hwhm) ≈ 1/π/2 atol=1e-12
    end
    
    @testset "Heavy tail behavior" begin
        L = Lorentzian(1.0)
        
        @test L(10.0) ≈ 1/π/101 atol=1e-12
        @test L(100.0) ≈ 1/π/10001 atol=1e-12
        
        @test L(10.0) > L(100.0)
    end
end
