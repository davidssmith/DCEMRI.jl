# units tests.
using Test
using DCEMRI

ccc4, cccnoisy4 = validate(4, makeplots=false)
@testset "QIBA v4 Extended Tofts Phantom" begin
  @test ccc4[1] >= 0.9998
  @test ccc4[2] >= 0.8904
  @test ccc4[3] >= 0.9999
  @test round(cccnoisy4[1], sigdigits=2) >= 0.97
  @test round(cccnoisy4[2], sigdigits=2) >= 0.70
  @test round(cccnoisy4[3], sigdigits=2) >= 0.97
end

ccc6, cccnoisy6 = validate(6, makeplots=false)
@testset "QIBA v6 Standard Tofts Phantom" begin
  @test ccc6[1] >= 0.9999
  @test ccc6[2] >= 0.9999
  @test round(cccnoisy6[1], sigdigits=2) >= 0.84
  @test round(cccnoisy6[2], sigdigits=2) >= 0.86
end
