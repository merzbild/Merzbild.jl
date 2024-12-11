"""
Boltzmann constant, J/K
"""
const k_B::Float64 = 1.380649e-23  # J / K

"""
Speed of light, m/s
"""
const c_light::Float64 = 299792458.0  # m/s

"""
Electron-Volt, K
"""
const eV::Float64 = 1.160451812e4  # K

"""
Electron-Volt, J
"""
const eV_J::Float64 = 1.602176634e-19 # 1eV = eV_J [J]

"""
1.0/eV[J], 1/J (or equivalent to how much 1 J is equal to expressed in eV)
"""
const eV_J_inv::Float64 = 1.0 / eV_J # 1 J = eV_J_inv [eV]

"""
``2 \\pi``
"""
const twopi::Float64 = 2 * π

"""
Electron mass divided by 1 eV, kg/J
"""
const e_mass_div_electron_volt = 5.6856301e-12  # m_e / eV, kg/J

"""
Positive and negative direction signs
"""
const direction_signs = [-1.0, 1.0]

"""
Elementary charge, C
"""
const q_e = 1.602176634e-19  # elementary charge, C