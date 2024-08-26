@testset "constants" begin
    @test abs(Merzbild.k_B - 1.380649e-23) < eps()
    @test abs(Merzbild.c_light - 299792458.0) < eps()
    @test abs(Merzbild.eV - 1.160451812e4) < eps()
    @test abs(Merzbild.eV_J - 1.602176634e-19) < eps()
    @test abs(Merzbild.eV_J_inv - 1.0 / Merzbild.eV_J) < eps()
    @test abs(Merzbild.twopi - 2 * π) < eps()
    @test abs(Merzbild.e_mass_div_electron_volt - 5.6856301e-12) < eps()
    @test (Merzbild.direction_signs[1] == -1.0 && Merzbild.direction_signs[2] == 1.0) || (Merzbild.direction_signs[1] == 1.0 && Merzbild.direction_signs[2] == -1.0)
    @test abs(Merzbild.q_e - 1.602176634e-19) < eps()
end