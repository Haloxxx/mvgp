using mvgp
using Test

Ackley(x, y) = -20*exp(-0.2*sqrt(0.5*(x^2+y^2)))-exp(0.5(cos(2*π*x)+cos(2*π*y)))+ℯ+20
Beale(x,y) = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
Matyas(x,y) = 0.26(x^2 + y^2) - 0.48*x*y
Booth(x,y) = (x+2*y-7)^2 + (2*x+y-5)^2
Sphere2(x,y) = x^2+y^2
Sphere3(x,y,z) = x^2+y^2+z^2
Easom(x,y) = -cos(x)*cos(y)*exp(-((x-π)^2 + (y-π)^2))
Levi(x,y) = sin(3*π*x)^2 + (x-1)^2*(1+sin(3*π*y)^2)+(y-1)^2*(1+sin(2*π*y)^2)
McCormick(x,y) = sin(x+y)+(x-y)^2-1.5*x+2.5*y+1
Booth_plus_z(x,y,z) = (x+2*y-7)^2 + (2*x+y-5)^2 + z^2

@testset "mvgp.jl" begin

    @testset "2 variables" begin

        @testset "Sphere2" begin
            for i in [1, 2, 3]
                result = gradientDesc(Sphere2, (i, 2), 2)
                @test isapprox(result[1], 0.0; atol=1e-6)
                @test isapprox(result[2], 0.0; atol=1e-6)
            end
        end

        @testset "Booth" begin
            for i in [1, 2, 3]
                #println("i: ",i)
                result = gradientDesc(Booth, (i, 0.5), 2)
                @test isapprox(result[1], 1.0; atol=1e-6)
                @test isapprox(result[2], 3.0; atol=1e-6)
            end
        end

        @testset "Matyas" begin
            for i in [-0.55, -0.5, -0.45]
                #println("i: ",i)
                result = gradientDesc(Matyas, (i, -1.5), 2)
                @test isapprox(result[1], 0; atol=1e-3)
                @test isapprox(result[2], 0; atol=1e-3)
            end
        end

    end
    
    @testset "3 variables" begin

        @testset "Sphere3" begin
            for i in [-0.55, -0.5, -0.45]
                #println("i: ",i)
                result = gradientDesc(Sphere3, (i, -1.5, 0), 3)
                @test isapprox(result[1], 0; atol=1e-3)
                @test isapprox(result[2], 0; atol=1e-3)
                @test isapprox(result[3], 0; atol=1e-3)
            end
        end

        @testset "Booth_plus_z" begin
            for i in [-0.55, -0.5, -0.45]
                #println("i: ",i)
                result = gradientDesc(Booth_plus_z, (i, -1.5, 0), 3)
                @test isapprox(result[1], 1; atol=1e-3)
                @test isapprox(result[2], 3; atol=1e-3)
                @test isapprox(result[3], 0; atol=1e-3)
            end
        end
    end
    
end