import os
os.environ['PYTHON_JULIACALL_SYSIMAGE'] = '/Users/xingzy/Documents/UVS/MCP_TN/TenKits.so'  # Or use -X juliacall-sysimage=mysysimage.so when running Python

import uuid
from mcp.server.fastmcp import FastMCP
from mcp.types import TextContent, ImageContent, BlobResourceContents
import logging


# Set up logging (this just prints messages to your terminal for debugging)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(name)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create the MCP server object
mcp = FastMCP()

# Here’s where you define your tools (functions the AI can use)
@mcp.tool()
def create_lattice_system(N: int, spin_type: str) -> TextContent:

    """Create a lattice system with N sites and given spin type.

    Args:
        N: the number of sites
        spin_type: the type of spin, e.g. "1/2", "1", "3/2", "2", ...

    Return:
        The uri and the file path of the data of lattice sites."""

    from juliacall import Main as jl
    jl.seval("using ITensors, ITensorMPS")
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    jl.N = N
    jl.spin_type = spin_type
    jl.seval("sites = gen_lattice(N, spin_type)")
    sites_id = str(uuid.uuid4())
    sites_uri = "sites://" + sites_id
    file_path = f"/tmp/{sites_id}.h5"
    jl.fname = file_path
    jl.seval("""
        h5open(fname, "w") do file
            write(file, "sites", sites)
        end
    """)
    response = f"Sites data is saved to {file_path}, sites_uri: {sites_uri}"
    return TextContent(type="text", text=response)

@mcp.tool()
def generate_Heisenberg_MPO(
    sites_id: str, J: list[float],
    D: list[int] = [1, 1, 1],
    K: list[int] = [1, 1, 1],
    A: list[float] = [0, 0, 0],
    alpha: list[float] = [1, 1, 1],
    delta: list[float] = [0, 0, 0],
) -> TextContent:

    """Generate the MPO of a generalized Heisenberg model with Hamiltonian:

    $$
    H = \sum_{j=1}^{N-1} \sum_{d_{x}=1}^{D_{x}} J^{x}_{j}(d_{x})S^{x}_{j}S^{x}_{j+d_{x}} + \sum_{d_{y}=1}^{D_{y}} J^{y}_{j}(d_{y}) S^{y}_{j}S^{y}_{j+d_{y}} + \sum_{d_{z}=1}^{D_{z}} J^{z}_{j}(d_{z})S^{z}_{j}S^{z}_{j+d_{z}}
    $$

    where
    $$
    J^{r}_{j} (d)= \Theta(N - j - d) \frac{J_{r}}{d^{K_{r}}_{r}} (1 + A_{r} \cos(2\pi \alpha_{r}j + \delta_{r})), \quad r = x, y, z
    $$

    $\Theta(x)$ is the step function, i.e. $\Theta(x) = 1$ if $x \geq 0$ and $\Theta(x) = 0$ if $x < 0$.

    $S_x, S_y, S_z$ are spin operators, $N$ is the number of sites.

    The Hamiltonian can have long-range interactions, where
    $D_x, D_y, D_z$ are the maximum interaction distance in the x, y, z direction respectively,
    $K_x, K_y, K_z$ are the power of the distance of decay.

    Moreover, the interaction strength can be modulated by a periodic function,
    where $A_x, A_y, A_z$ are the amplitude of the modulation,
    $\alpha_x, \alpha_y, \alpha_z$ are the frequency of the modulation,
    $\delta_x, \delta_y, \delta_z$ are the phase of the modulation.

    Args:
        sites_id: the UUID of site index data, with the form "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
        J: a list [Jx, Jy, Jz] where
        Jx: the parameter $J_x$, float number
        Jy: the parameter $J_y$, float number
        Jz: the parameter $J_z$, float number
        D: a list [Dx, Dy, Dz] where
        Dx: the parameter $D_x$, positive integer
        Dy: the parameter $D_y$, positive integer
        Dz: the parameter $D_z$, positive integer
        K: a list [Kx, Ky, Kz] where
        Kx: the parameter $K_x$, positive integer
        Ky: the parameter $K_y$, positive integer
        Kz: the parameter $K_z$, positive integer
        A: a list [Ax, Ay, Az] where
        Ax: the parameter $A_x$, float number in [0, 1]
        Ay: the parameter $A_y$, float number in [0, 1]
        Az: the parameter $A_z$, float number in [0, 1]
        alpha: a list [alphax, alphay, alphaz] where
        alphax: the parameter $\alpha_x$, float number
        alphay: the parameter $\alpha_y$, float number
        alphaz: the parameter $\alpha_z$, float number
        delta: a list [deltax, deltay, deltaz] where
        deltax: the parameter $\delta_x$, float number
        deltay: the parameter $\delta_y$, float number
        deltaz: the parameter $\delta_z$, float number

    Return:
        The MPS of the desired Hamiltonian."""

    from juliacall import Main as jl, convert
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    jl.seval("using ITensors, ITensorMPS")
    sites_file_path = f"/tmp/{sites_id}.h5"
    jl.sites_fname = sites_file_path
    jl.seval("""
        f = h5open(sites_fname,"r")
        sites = read(f, "sites", Vector{Index{Int}})
        close(f)
        N = length(sites)
    """)
    jl.J= convert(jl.Vector[jl.Float64], J)
    jl.D= convert(jl.Vector[jl.Int64], D)
    jl.K= convert(jl.Vector[jl.Int64], K)
    jl.A= convert(jl.Vector[jl.Float64], A)
    jl.α= convert(jl.Vector[jl.Float64], alpha)
    jl.phi= convert(jl.Vector[jl.Float64], delta)
    jl.seval("""
             mpo = H_MPO(sites, J;
             D=D, K=K, A=A, α=α, δ=phi)
    """)
    mpo_id = str(uuid.uuid4())
    mpo_uri = "mpo://" + mpo_id
    mpo_file_path = f"/tmp/{mpo_id}.h5"
    jl.mpo_fname = mpo_file_path
    jl.seval("""
        h5open(mpo_fname, "w") do file
            write(file, "mpo", mpo)
        end
    """)
    response = f"MPO data is saved to {mpo_file_path}, MPO_uri: {mpo_uri}"
    return TextContent(type="text", text=response)

@mcp.tool()
def ground_energy_calculation(
    sites_id: str,
    mpo_id: str,
) -> TextContent:

    """Calculate the ground state of a generalized Heisenberg model with Hamiltonian:

    Args:
        sites_id: the UUID of site index data, e.g. "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
        mpo_id: the uri of MPO data, e.g. "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"

    Return:
        The ground energy of the desired Hamiltonian."""

    from juliacall import Main as jl
    jl.seval("using ITensors, ITensorMPS")
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    sites_file_path = f"/tmp/{sites_id}.h5"
    jl.sites_fname = sites_file_path
    jl.seval("""
        f = h5open(sites_fname,"r")
        sites = read(f, "sites", Vector{Index{Int}})
        close(f)
    """)
    mpo_file_path = f"/tmp/{mpo_id}.h5"
    jl.mpo_fname = mpo_file_path
    jl.seval("""
        f = h5open(mpo_fname,"r")
        mpo = read(f, "mpo", MPO)
        close(f)
    """)
    jl.seval("energy, psi0 = calculate_GS(sites, mpo)")

    gs_id = str(uuid.uuid4())
    gs_uri = "gs://" + gs_id
    gs_file_path = f"/tmp/{gs_id}.h5"
    jl.gs_fname = gs_file_path
    jl.seval("""
        h5open(gs_fname, "w") do file
            write(file, "mps", psi0)
        end
    """)
    response = f"The ground state energy is {str(jl.energy)}. Ground state's MPS data is saved to {gs_file_path}, gs_uri: {gs_uri}"
    return TextContent(type="text", text=response)

@mcp.tool()
def local_observable_measurement(
    psi_id: str, observable: str,
) -> TextContent:

    """Calculate the local observable of the state psi:

    Args:
        psi_id: the UUID of MPS data of the state, e.g. "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
        observable: the local observable to be measured, e.g. "Sz", "Sx", "Sy", "S+", "S-"

    Return:
        Each site's local observables of ground state the desired Hamiltonian."""

    from juliacall import Main as jl
    jl.seval("using ITensors, ITensorMPS")
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    psi_file_path = f"/tmp/{psi_id}.h5"
    jl.psi_fname = psi_file_path
    jl.seval("""
        f = h5open(psi_fname,"r")
        psi = read(f, "mps", MPS)
        close(f)
    """)
    jl.observable = observable
    obs = jl.seval("measure_local(psi, observable)")
    return TextContent(type="text", text=str(obs))

@mcp.tool()
def entanglement_entropy_calculation(
    psi_id: str, c: int,
) -> TextContent:

    """Calculate the entanglement entropy of the state psi,
     the two subsystems of is seperated at the given site:

    Args:
        psi_id: the UUID of MPS data of the state, e.g. "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
        c: integer in [2, N-1], subsystem A contains sites [1, c-1], subsystem B contains sites [c, N]

    Return:
        The entanglement entropy of the state."""

    from juliacall import Main as jl
    jl.seval("using ITensors, ITensorMPS")
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    jl.seval("using DataFrames, CSV")
    psi_file_path = f"/tmp/{psi_id}.h5"
    jl.psi_fname = psi_file_path
    jl.seval("""
        f = h5open(psi_fname,"r")
        psi = read(f, "mps", MPS)
        close(f)
    """)
    jl.c = c
    ee = jl.seval("entanglement_entropy(psi, c)")

    return TextContent(type="text", text=str(ee))

@mcp.tool()
def correlation_function_calculation(
    psi_id: str, op1: str, op2: str,
) -> TextContent:

    """Calculate the correlation function Matrix $C$,
     where $C_{ij} = \langle op1_i op2_j \rangle$ with respect to the state psi,:

    Args:
        psi_id: the UUID of MPS data of the state, e.g. "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
        op1: the first operator in the correlation function, e.g. "S+", "S-", "Sz", "Sx", "Sy"
        op2: the second operator in the correlation function, e.g. "S+", "S-", "Sz", "Sx", "Sy"

    Return:
        The correlation Matrix $C$."""

    from juliacall import Main as jl
    jl.seval("using ITensors, ITensorMPS")
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    psi_file_path = f"/tmp/{psi_id}.h5"
    jl.psi_fname = psi_file_path
    jl.seval("""
        f = h5open(psi_fname,"r")
        psi = read(f, "mps", MPS)
        close(f)
    """)
    jl.op1 = op1
    jl.op2 = op2
    jl.seval("gf = Greenf(psi, op1, op2)")

    gf_id = str(uuid.uuid4())
    gf_uri = "gs://" + gf_id
    gf_file_path = f"/tmp/{gf_id}.csv"
    jl.gf_fname = gf_file_path
    jl.seval("""
        CSV.write(gf_fname, DataFrame(abs.(gf), :auto))
    """)
    response = f"The norm of the correlation matrix is saved to {gf_file_path}, correlation_uri: {gf_uri}"
    return TextContent(type="text", text=response)

@mcp.tool()
def time_evolution_simulation_entanglement_entropy(
    sites_id: str, state_init: list[str], c: int, J: list[float],
    A: list[float] = [0., 0., 0., 0.],
    alpha: list[float] = [1., 1., 1., 1.],
    delta: list[float] = [0., 0., 0., 0.],
    hz: float = 0.,
    ttotal: float = 5.,
    dt: float = 0.125,
) -> TextContent:

    """Simulate the time evolution of a periodic Heisenberg model with periodic magnetic field along z direction:

        $$
            H = \sum_{j=1}^{N-1} J^{x}_{j} S^{x}_{j}S^{x}_{j+1} +  J^{y}_{j} S^{y}_{j}S^{y}_{j+1} + J^{z}_{j} S^{z}_{j}S^{z}_{j+1}+ \sum_{j=1}^{N}h_{j} S^{z}_{j}
        $$
        where
        $$
            J^{r}_{j} (d)=  J_{r} (1 + A_{r} \cos(2\pi \alpha_{r}j + \delta_{r})), \quad r = x, y, z
        $$
        and
        $$
            h_{j} (d)=  h_{z} (1 + A_{h} \cos(2\pi \alpha_{h}j + \delta_{h}))
        $$
        $S_x, S_y, S_z$ are spin operators, $N$ is the number of sites.

        The interaction and the field strength can be modulated by a periodic function,
        where $A_x, A_y, A_z, A_{h}$ are the amplitude of the modulation,
        $\alpha_x, \alpha_y, \alpha_z, \alpha_{h}$ are the frequency of the modulation,
        $\delta_x, \delta_y, \delta_z, \delta_{h}$ are the phase of the modulation.


    Args:
        sites_id: the UUID of site index data, with the form "xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
        state_init: a list of length N, each element is the initial state on each site, e.g. ["Up", "Dn", "Up", ...]
        c: integer in [2, N-1], subsystem A contains sites [1, c-1], subsystem B contains sites [c, N]
        J: a list [Jx, Jy, Jz] where
        Jx: the parameter $J_x$, float number
        Jy: the parameter $J_y$, float number
        Jz: the parameter $J_z$, float number
        A: a list [Ax, Ay, Az, Ah] where
        Ax: the parameter $A_x$, float number in [0, 1]
        Ay: the parameter $A_y$, float number in [0, 1]
        Az: the parameter $A_z$, float number in [0, 1]
        Ah: the parameter $A_h$, float number in [0, 1]
        alpha: a list [alphax, alphay, alphaz, alphah] where
        alphax: the parameter $\alpha_x$, float number
        alphay: the parameter $\alpha_y$, float number
        alphaz: the parameter $\alpha_z$, float number
        alphah: the parameter $\alpha_h$, float number
        delta: a list [deltax, deltay, deltaz, deltah] where
        deltax: the parameter $\delta_x$, float number
        deltay: the parameter $\delta_y$, float number
        deltaz: the parameter $\delta_z$, float number
        deltah: the parameter $\delta_h$, float number
        hz: the parameter $h_z$, float number
        ttotal: total time of evolution, positive float number
        dt: time step of evolution, positive float number

    Return:
        The entanglement entropy at each time step."""

    from juliacall import Main as jl, convert
    jl.seval("using ITensors, ITensorMPS")
    jl.seval("using .TenKits")
    jl.seval("using HDF5")
    jl.seval("using DataFrames, CSV")
    sites_file_path = f"/tmp/{sites_id}.h5"
    jl.sites_fname = sites_file_path
    jl.seval("""
        f = h5open(sites_fname,"r")
        sites = read(f, "sites", Vector{Index{Int}})
        close(f)
    """)
    jl.J= convert(jl.Vector[jl.Float64], J)
    jl.A= convert(jl.Vector[jl.Float64], A)
    jl.α= convert(jl.Vector[jl.Float64], alpha)
    jl.phi= convert(jl.Vector[jl.Float64], delta)
    jl.state_init = convert(jl.Vector[jl.String], state_init)
    jl.c = c
    jl.hz = hz
    jl.ttotal = ttotal
    jl.dt = dt
    jl.seval("""
            ee = Evolution_TEBD(
                sites, state_init, c, J;
                A=A,
                α=α,
                δ=phi,
                hz=hz,
                ttotal=ttotal,
                dt=dt
            )
    """)
    ee_id = str(uuid.uuid4())
    ee_uri = "gs://" + ee_id
    ee_file_path = f"/tmp/{ee_id}.csv"
    jl.ee_fname = ee_file_path
    jl.seval("""
        CSV.write(ee_fname, DataFrame(eachcol(ee), [:t, :C]))
    """)
    response = f"Evolution of the entanglement entropy is saved to {ee_file_path}, evolution_uri: {ee_uri}"
    return TextContent(type="text", text=response)

# This is the main entry point for your server
def main():
    logger.info('Starting your-new-server')
    mcp.run('stdio')


if __name__ == "__main__":
    main()
