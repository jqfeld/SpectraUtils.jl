using SpectraUtils
using SpectraUtils: gaussian
using SpectraUtils: gaussian
using Test

@testset "DopplerFree Line Shape" begin
    @testset "Basic composition" begin
        envelop = x -> 2.0
        dip = x -> 0.5
        shape = DopplerFree(0.1, envelop, dip)

        @test shape isa DopplerFree
        @test shape.depth == 0.1
        @test shape.envelop === envelop
        @test shape.dip === dip

        @test shape(0.0) ≈ 2.0 * (1 - 0.1 * 0.5)
        @test shape(1.0) ≈ shape(0.0)  # constant envelope and dip
    end

    @testset "Symmetry from components" begin
        envelop = Gaussian(0.8)
        dip = Lorentzian(0.2)
        shape = DopplerFree(0.4, envelop, dip)

        for x in (0.0, 0.5, 1.0, 2.0)
            @test shape(x) ≈ shape(-x) atol=1e-12
        end
    end

    @testset "Depth scaling" begin
        envelop = Lorentzian(1.0)
        dip = Gaussian(0.1)

        shallow = DopplerFree(0.0, envelop, dip)
        deep = DopplerFree(0.9, envelop, dip)

        @test shallow(0.0) ≈ envelop(0.0)
        @test deep(0.0) ≈ envelop(0.0) * (1 - 0.9 * dip(0.0))
        @test deep(0.0) < shallow(0.0)
    end

    @testset "Flexible callables" begin
        envelop = VoigtApprx(0.5, 0.1)
        dip = x -> gaussian(x, 0.2)
        shape = DopplerFree(0.3, envelop, dip)

        @test shape(0.0) ≈ envelop(0.0) * (1 - 0.3 * dip(0.0))
        @test shape(0.5) ≈ envelop(0.5) * (1 - 0.3 * dip(0.5))
        @test shape(0.5) ≤ envelop(0.5)
    end
end
