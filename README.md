# MCP Server for Scientific Calculation

A [**Model Context Protocol (MCP)**](https://modelcontextprotocol.io/) server that leverages **Large Language Models (LLMs)** to perform and assist in scientific calculations.

## 📦 Requirement
We successfully test the code with softwares as the followings:
- **Julia** = 1.11.6
	- `ITensors` = "0.9.9"
	- `ITensorMPS` = "0.3.19"
	- `HDF5` = "0.17.2"
	- `CSV.jl` = "0.10.15"
	- `DataFrames.jl` = "1.7.0"
- **Python** = 3.13.5, initialized by `uv`, with `uuid` and `juliacall` installed.

## 🚀 Quick Start
### Clone the repository and install dependencies:
```bash
git clone https://github.com/xing-phys/MCP_TN.git
```

### Precompile the Julia package to get a system image `TenKits.so`:
```bash
cd MCP_TN/TenKits
julia generate_so.jl
```

### Configure the server script
Edit the `main.py` file:
```bash
vim MCP_TN/main.py
```
Replace the value of `os.environ['PYTHON_JULIACALL_SYSIMAGE']` to the actual local path to `TenKits.so`

### Test the server with inspector
The following operations are in folder `MCP_TN` by default.
Activate the virtual environment at first:
```bash
source .venv/bin/activate
```
Use the inspector to check the tools:
```bash
mcp dev main.py
```

## 🖥️ Integration with LLM  Clients
 Open an LLM Client (Claude Desktop, 5ire, etc) and configure the MCP server:
- **Sever name**: `MCP_TN`
- **Command**: `uv`
- **Arguments**: `run mcp run /path/to/your/MCP_TN/main.py` _(Replace with your actual local path to `main.py`)_
- **Environment**:
    - Variable: `UV_PROJECT`
    - Value: `/path/to/your/MCP_NetKet` _(Replace with your actual local path to the `MCP_TN` folder)_

## 🎯 Interactive Demo: Ferromagnetic XXZ Heisenberg Model

Try this step-by-step demo to see the MCP server in action! This demonstrates a work flow for calculate the ground state magnetization of the XXZ model.
### Step 1: Generate the MPO of the Hamiltonian
```
Create the MPO of this XXZ model:
$$
H = \sum_{j=1}^{N-1} J(S^{x}_{j}S^{x}_{j+1} + S^{y}_{j}S^{y}_{j+1} ) + J^{z} S^{z}_{j}S^{z}_{j+1}
$$
set $J = -1$, $J^{z}$ from -1.5 to -0.5, step 0.1. The system size is set as $N = 50$.
```
**Expected Results:**
The MPO for each $J^{z}$ is generated properly and saved to the desired `/tmp/` folder and labeled by an unique `UUID`.

### Step 2: Ground State Calculation
```
For each Jz, calculate the ground state
```
**Expected Results:**
The ground energies are returned. The ground state MPS for each $J^{z}$ is generated properly and saved to the desired `/tmp/` folder and labeled by an unique `UUID`.
### Step 3: Calculate the Magnetization
```
Calculate the average magnetization $J^{z}$
$$
\sum_{j=1}^{N} \frac{\left<S^{z}_{j}\right>}{L}
$$
```
**Expected Results:**
Get the average magnetization with absolute values about 0.5 for $J^{z} < -1$ and nearly 0 for $J^{z} > -1$.

## 📂 Project Structure

``` text
MCP_TN
│── README.md
│── TenKits
│   │── Manifest.toml
│   │── Project.toml
│   │── generate_so.jl
│   │── precompile_script.jl
│   └── src
│       └── TenKits.jl
│── __pycache__
│   └── main.cpython-313.pyc
│── main.py
│── pyproject.toml
└── uv.lock
```
