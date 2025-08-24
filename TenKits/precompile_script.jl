# /path/to/my_project/precompile_script.jl

using Pkg
Pkg.activate(".")

using TenKits
using HDF5


println("Precompiling function...")

N = 10
spin_type = "1"
H_params = [1.0, 1.0, 1.0]

sites = gen_lattice(N, spin_type)

h = H_MPO(sites, H_params)

@time energy, psi = calculate_GS(sites, h)

println("Ground state energy: ", energy)

state_init = ["Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn"]
@time ee = Evolution_TEBD(sites, state_init, 5, H_params; hz=1.)
println(ee)

@time tet4 = measure_local(psi, "Sz")

@time entanglement_entropy(psi, 5)

@time tet = Greenf(psi, "S+", "S-")
