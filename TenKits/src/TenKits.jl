module TenKits
    ############ Import packages
    using ITensors, ITensorMPS
    using CairoMakie
    using CSV, DataFrames

    export gen_lattice, H_MPO, calculate_GS, measure_local,
    entanglement_entropy, Greenf, Evolution_TEBD

    function gen_lattice(
        N::Int64,
        spin_type::String
    )
        sites = siteinds("S=" * spin_type, N)
        return sites
    end

    function H_MPO(
        sites::Vector{Index{Int64}},
        J::Vector{Float64};
        D::Vector{Int64}=[1, 1, 1],
        K::Vector{Int64}=[1, 1, 1],
        A::Vector{Float64}=[0., 0., 0.],
        α::Vector{Float64}=[1., 1., 1.],
        δ::Vector{Float64}=[0., 0., 0.]
    )

        N = length(sites)
        os = OpSum()
        for j in 1:(N-1)
            for d in 1:D[3]
                if j + d > N
                    break
                end
                Jz = J[3] * (1 + A[3] * cos(2π * α[3] * j + δ[3])) / d^(K[3])
                os += Jz, "Sz", j, "Sz", j+d
            end
            for d in 1:D[1]
                if j + d > N
                    break
                end
                Jx = J[1] * (1 + A[1] * cos(2π * α[1] * j + δ[1])) / d^(K[1])
                os += Jx, "Sx", j, "Sx", j+d
            end
            for d in 1:D[2]
                if j + d > N
                    break
                end
                Jy = J[2] * (1 + A[2] * cos(2π * α[2] * j + δ[2])) / d^(K[2])
                os += Jy, "Sy", j, "Sy", j+d
            end
        end
        
        return MPO(os, sites)
    end

    ############ Ground state calculation
    function calculate_GS(
        sites::Vector{Index{Int64}},
        H::MPO;
        linkdim=10,
        nsweeps=5,
        maxdim=[10,20,100,100,200],
        cutoff=1E-10
    )

        psi0 = random_mps(sites; linkdims=linkdim)
        energy, psi = dmrg(
            H, psi0;
            nsweeps=nsweeps,
            maxdim=maxdim,
            cutoff=cutoff
        )
        return energy, psi
    end

    ############ Measure local observables
    function measure_local(psi::MPS, observable::String)
        N = length(psi)
        values = Float64[]
        for j in 1:N
            orthogonalize!(psi, j)
            val = real(expect(psi, observable; sites=j))
            push!(values, val)
        end
        return values
    end

    ############ entanglement entropy
    function entanglement_entropy(psi::MPS, c::Int)
        psi = orthogonalize(psi, c)
        U,S,V = svd(psi[c], (linkinds(psi, c-1)..., siteinds(psi, c)...))
        return sum([-S[n,n]^2 *log(S[n,n]^2) for n in 1:dim(S,1)])
    end

    ############ Green's function
    function Greenf(psi::MPS,A::String,B::String)
        return correlation_matrix(psi,A,B)
    end

    ############ Evolution
    function Evolution_TEBD(
        sites::Vector{Index{Int64}},
        state_init::Vector{String},
        c::Int64,
        J::Vector{Float64};
        A::Vector{Float64}=[0., 0., 0., 0.],
        α::Vector{Float64}=[1., 1., 1., 1.],
        δ::Vector{Float64}=[0., 0., 0., 0.],
        hz::Float64=0.0,
        ttotal=5,
        dt=0.125,
        cutoff=1e-8
    )
        
        N = length(sites)

        gates = ITensor[]

        # Parameters for 4th order Trotter-Suzuki decomposition (Forest-Ruth Scheme)
        s = 1 / (2 - 2^(1/3))
        a1 = s/2
        a2 = (1-s)/2
        b1 = s
        b2 = 1-2s
        TS = [a1, b1, a2, b2, a2, b1, a1]

        for ii in eachindex(TS)
            start_index = 2 - mod(ii, 2)
            final_index = N - 1
            for j in start_index:2:final_index
                s1 = sites[j]
                s2 = sites[j + 1]
                Jz = J[3] * (1 + A[3] * cos(2π * α[3] * j + δ[3]))
                Jx = J[1] * (1 + A[1] * cos(2π * α[1] * j + δ[1]))
                Jy = J[2] * (1 + A[2] * cos(2π * α[2] * j + δ[2]))
                hzj = hz * (1 + A[4] * cos(2π * α[4] * j + δ[4]))
                hj = Jz * op("Sz", s1) * op("Sz", s2) + 
                    Jx * op("Sx", s1) * op("Sx", s2) + 
                    Jy * op("Sy", s1) * op("Sy", s2) +
                    hzj * op("Sz", s1) * op("Id", s2) 
                if j == N - 1
                    hj += hz * op("Id", s1) * op("Sz", s2)
                end
                Gj = exp(-im * TS[ii] * dt * hj)
                push!(gates, Gj)
            end
        end

        psi = MPS(sites, state_init)

        obs = Float64[]
        # Compute <ob> at each time step
        # then apply the gates to go to the next time
        for t in 0.0:dt:ttotal

            ee = entanglement_entropy(psi, c)
            push!(obs, ee)

            t≈ttotal && break

            psi = apply(gates, psi; cutoff)

            normalize!(psi)
        end
        return hcat(0:dt:ttotal, obs)
    end

end # module TenKits
