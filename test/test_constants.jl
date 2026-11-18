@testset "constants" begin
    @test Merzbild.k_B == 1.380649e-23
    @test Merzbild.c_light == 299792458.0
    @test Merzbild.eV == 1.160451812e4
    @test Merzbild.eV_J == 1.602176634e-19
    @test Merzbild.eV_J_inv == 1.0 / Merzbild.eV_J
    @test Merzbild.twopi == 2 * π
    @test Merzbild.sqrt_pi == sqrt(π)
    @test Merzbild.inv_sqrt_2 == 1.0 / sqrt(2.0)
    @test Merzbild.e_mass_div_electron_volt == 5.685630060187e-12
    @test (Merzbild.direction_signs[1] == -1.0 && Merzbild.direction_signs[2] == 1.0) || (Merzbild.direction_signs[1] == 1.0 && Merzbild.direction_signs[2] == -1.0)
    @test Merzbild.q_e == 1.602176634e-19
end
