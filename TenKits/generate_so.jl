using PackageCompiler

# NOTE: It's best to run this from your project's root directory.
# This command tells PackageCompiler to:
# 1. Include `MyModule` in the new system image.
# 2. Use the current active project environment (important!).
# 3. Set the output path for the system image to `MyModuleSysimage.so`.
# 4. Run `precompile_script.jl` to trace the necessary methods.

create_sysimage(
    [:TenKits], # Include the TenKits module;
    project=".", # Use the project in the current directory
    sysimage_path="TenKits.so",
    precompile_execution_file="precompile_script.jl"
)
