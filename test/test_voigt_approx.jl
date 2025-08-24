using SpectraUtils
using Test
using QuadGK

@testset "VoigtApprx Line Shape" begin
    @testset "Basic functionality" begin
        VA = VoigtApprx(1.0, 0.5)
        
        @test VA isa VoigtApprx
        @test VA.sigma == 1.0
        @test VA.gamma == 0.5
        
        @test VA(0.0) isa Real
        @test VA(0.0) > 0
    end
    
    @testset "Normalization" begin
        VA = VoigtApprx(1.0, 0.5)
        area, _ = quadgk(x -> VA(x), -Inf, Inf, rtol=1e-8)
        @test area ≈ 1.0 atol=1e-5
        
        VA2 = VoigtApprx(2.0, 1.0)
        area2, _ = quadgk(x -> VA2(x), -Inf, Inf, rtol=1e-8)
        @test area2 ≈ 1.0 atol=1e-5
    end
    
    @testset "Symmetry" begin
        VA = VoigtApprx(1.5, 0.8)
        test_points = [0.5, 1.0, 2.0, 3.0]
        
        for x in test_points
            @test VA(x) ≈ VA(-x) atol=1e-12
        end
    end
    
    @testset "Parameter dependence" begin
        VA_const = VoigtApprx(2.0, 1.0)
        VA_func = VoigtApprx(p -> p[1], p -> p[2])
        
        @test VA_const(1.0) ≈ VA_func(1.0, [2.0, 1.0]) atol=1e-12
        
        @test VA_func(0.0, [1.0, 0.5]) > 0
        @test VA_func(0.0, [0.5, 0.25]) > VA_func(0.0, [1.0, 0.5])
    end
    
    @testset "Different parameter types" begin
        VA_mixed = VoigtApprx(1.0, 1)
        VA_float = VoigtApprx(1.0, 1.0)
        
        @test VA_mixed(0.0) ≈ VA_float(0.0) atol=1e-12
    end
    
    @testset "Limiting cases" begin
        # When gamma -> 0, should approach Gaussian-like behavior
        VA_gaussian_like = VoigtApprx(1.0, 1e-6)
        G = Gaussian(1.0)
        
        @test VA_gaussian_like(0.0) ≈ G(0.0) rtol=1e-2
        
        # When sigma -> 0, should approach Lorentzian-like behavior
        VA_lorentzian_like = VoigtApprx(1e-6, 1.0)
        L = Lorentzian(1.0)
        
        @test VA_lorentzian_like(0.0) ≈ L(0.0) rtol=1e-2
    end
    
    @testset "Approximation accuracy" begin
        # Test that VoigtApprx is reasonably close to exact Voigt
        σ = 1.0
        γ = 0.5
        
        VA = VoigtApprx(σ, γ)
        V = Voigt(σ, γ)
        
        test_points = [-2.0, -1.0, 0.0, 1.0, 2.0]
        
        for x in test_points
            @test VA(x) ≈ V(x) rtol=2e-2  # Allow 2% relative error
        end
    end
    
    @testset "Mathematical properties" begin
        VA = VoigtApprx(1.0, 0.5)
        
        @test VA(0.0) > VA(1.0)
        @test VA(1.0) > VA(2.0)
        
        @test all(VA(x) > 0 for x in [-5, -1, 0, 1, 5])
        
        center_val = VA(0.0)
        @test center_val isa Real
        @test center_val > 0
    end
    
    @testset "Width parameter conversion" begin
        # Test that the approximation handles FWHM conversion correctly
        σ = 1.0
        γ = 0.5
        
        VA = VoigtApprx(σ, γ)
        
        # The approximation should produce reasonable values
        @test VA(0.0) > 0
        @test VA(1.0) > 0
        @test VA(2.0) > 0
        
        # Check monotonic decrease from center
        @test VA(0.0) > VA(0.5) > VA(1.0)
    end
    
    @testset "Edge cases" begin
        # Test with very small parameters
        VA_small = VoigtApprx(1e-3, 1e-3)
        @test VA_small(0.0) > 0
        
        # Test with large parameters
        VA_large = VoigtApprx(10.0, 5.0)
        @test VA_large(0.0) > 0
        
        # Test with unequal parameters
        VA_unequal = VoigtApprx(0.1, 5.0)
        @test VA_unequal(0.0) > 0
    end
end