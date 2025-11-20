import numpy as np
import os
import os.path as osp
import shutil
import json
import argparse
import mdtraj as md
from omegaconf import OmegaConf
from panda.builder import build_system
from panda.parser import parse_C6_C12


def main(config_path, overrides=None):
    # Load config using OmegaConf
    cfg = OmegaConf.load(config_path)
    # Merge CLI overrides if provided
    if overrides:
        cfg = OmegaConf.merge(cfg, overrides)

    output_path = os.path.join(cfg.output_dir, cfg.exp_folder)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # cfg.substrate = substr_gro_path
    with open(os.path.join(output_path, "config.json"), "w") as f:
        json.dump(OmegaConf.to_container(cfg, resolve=True), f, indent=4)

    # Start building system
    build_system(os.path.join(output_path, "config.json"))
    # build_system(config_path)

    # Update config after building system
    cfg = OmegaConf.load(osp.join(output_path, "config.json"))

    # Generating system.itp
    with open(os.path.join(output_path, "system.itp"), "w") as f:
        for component_cfg in set(cfg.components):
            f.write(f'#include "{component_cfg.name}.itp"\n')

        f.write(f'#include "{osp.split(cfg.substr_itp_path)[-1]}"\n')

        f.write(f"\n[ system ]\n{cfg.system_name}\n")
        f.write("\n[ molecules ]\n; molecule name\tnr.\n")
        # TODO: Make it more general
        if cfg.substrate.get("silanol_density", None):
            f.write(f"SIL\t{1}\n")
        elif cfg.substrate.get("hybridisation", None) >= 0:
            f.write(f"MUSC\t{1}\n")
        else:
            f.write(f"CAL\t{1}\n")
        for i, component_cfg in enumerate(cfg.components):
            f.write(f"{component_cfg.name}\t{component_cfg.numbers}\n")

    # Generating run.sh
    # TODO: add support for more than 3 mdp files
    assert len(cfg.pipeline) == 3, "Pipeline must contain 3 mdp files"
    with open(os.path.join(output_path, "run.sh"), "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -J gromacs
#SBATCH -w node{cfg.node}
#SBATCH -N 1\t# Number of nodes requested
#SBATCH -n {cfg.n_mpi}\t# Total number of mpi tasks requested
#SBATCH --cpus-per-task 1\t# Total number of omp tasks requested
{f"#SBATCH --dependency=afterany:{cfg.dependency}" if cfg.dependency != -1 else ""}

gmx grompp -f {osp.split(cfg.pipeline[0])[1]} -c {cfg.system_name}_init.gro -p {osp.split(cfg.topology)[1]} -n index.ndx -o {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[0])[1])[0].split("_")[-1]} -maxwarn 10
mpirun -np {cfg.n_mpi} --cpu-set {cfg.init_core}-{cfg.init_core + cfg.n_mpi - 1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[0])[1])[0].split("_")[-1]} -ntomp 1 -nb gpu -gpu_id {cfg.gpu_id} -mn index.ndx
rm ./*pdb

gmx grompp -f {osp.split(cfg.pipeline[1])[1]} -c {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[0])[1])[0].split("_")[-1]}.gro -p {osp.split(cfg.topology)[1]} -n index.ndx -o {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[1])[1])[0].split("_")[-1]} -maxwarn 10
mpirun -np {cfg.n_mpi} --cpu-set {cfg.init_core}-{cfg.init_core + cfg.n_mpi - 1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[1])[1])[0].split("_")[-1]} -ntomp 1 -nb gpu -gpu_id {cfg.gpu_id} -dlb yes -mn index.ndx
rm ./*pdb

gmx grompp -f {osp.split(cfg.pipeline[2])[1]} -c {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[1])[1])[0].split("_")[-1]}.gro -p {osp.split(cfg.topology)[1]} -n index.ndx -o {cfg.system_name} -maxwarn 10
mpirun -np {cfg.n_mpi} --cpu-set {cfg.init_core}-{cfg.init_core + cfg.n_mpi - 1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {cfg.system_name} -ntomp 1 -nb gpu -gpu_id {cfg.gpu_id} -dlb yes -mn index.ndx
rm ./*pdb""")

    # Generating mdp files from templates
    with open(cfg.pipeline[0]) as f:
        steep_final_text = f.read().format(
            **{
                "freeze": "freezegrps          =  SILICA\nfreezedim           =  Y Y Y"
                if cfg.substrate.get("silanol_density", None)
                else ""
            }
        )
    with open(osp.join(output_path, osp.split(cfg.pipeline[0])[1]), "w") as f:
        f.write(steep_final_text)

    with open(cfg.pipeline[1]) as f:
        short_final_text = f.read().format(
            **{
                "freeze": "freezegrps          =  SILICA\nfreezedim           =  Y Y Y"
                if cfg.substrate.get("silanol_density", None)
                else "",
                "temp": cfg.temp,
            }
        )
    with open(osp.join(output_path, osp.split(cfg.pipeline[1])[1]), "w") as f:
        f.write(short_final_text)

    with open(cfg.pipeline[2]) as f:
        run_final_text = f.read().format(
            **{
                "nsteps": cfg.nsteps,
                "freeze": "freezegrps          =  SILICA\nfreezedim           =  Y Y Y"
                if cfg.substrate.get("silanol_density", None)
                else "",
                "temp": cfg.temp,
            }
        )
    with open(osp.join(output_path, osp.split(cfg.pipeline[2])[1]), "w") as f:
        f.write(run_final_text)

    # Generating force field .top file from template
    with open(cfg.topology) as f:
        topology_text = f.read()
    calcite_params = parse_C6_C12(topology_text, ["CA", "OCA", "CCA"])
    decane_params = parse_C6_C12(topology_text, ["CH2", "CH3"])
    pairwise_params = {}
    for i, Ci in enumerate(["C6", "C12"]):
        for dec_name, Ci_dec in decane_params.items():
            for cal_name, Ci_cal in calcite_params.items():
                pairwise_params["_".join([Ci, dec_name, cal_name])] = "{:.2e}".format(
                    cfg.scale * np.sqrt(Ci_dec[i] * Ci_cal[i])
                )
    topology_text = topology_text.format(**pairwise_params)
    with open(os.path.join(output_path, osp.split(cfg.topology)[1]), "w") as f:
        f.write(topology_text)

    # Copy all neccesary files to the output directory
    for component_cfg in cfg.components:
        shutil.copy(
            osp.join(cfg.mol_path.replace("gro", "itp"), component_cfg.name + ".itp"),
            osp.join(output_path, component_cfg.name + ".itp"),
        )
    shutil.copy(
        osp.join(output_path, cfg.system_name + ".gro"),
        osp.join(output_path, cfg.system_name + "_init.gro"),
    )
    shutil.copy(
        cfg.substr_itp_path,
        osp.join(output_path, osp.split(cfg.substr_itp_path)[-1]),
    )

    # Generate index.ndx file
    ndx_output_path = os.path.join(output_path, "index.ndx")

    # If substrate has ndx file, use it as a base
    if hasattr(cfg, "substr_ndx_path"):
        print(f"Using substrate NDX file: {cfg.substr_ndx_path}")
        shutil.copy(
            cfg.substr_ndx_path,
            ndx_output_path,
        )
    else:
        # Create a new ndx file
        print("Creating new NDX file")
        os.makedirs(os.path.dirname(ndx_output_path), exist_ok=True)
        with open(ndx_output_path, "w") as f:
            # If no substrate NDX exists, we'll just create the ALL group later
            pass

    # Add System and Other groups to index.ndx
    system_gro_path = osp.join(output_path, cfg.system_name + ".gro")
    system = md.load(system_gro_path)
    all_atoms = list(range(1, system.n_atoms + 1))  # 1-indexed for GROMACS

    # Identify substrate atoms
    substrate_atoms = set()
    if hasattr(cfg, "substr_ndx_path"):
        # Read the substrate NDX file to get substrate atom indices
        with open(cfg.substr_ndx_path, "r") as f:
            lines = f.readlines()
            reading_substrate = False
            for line in lines:
                line = line.strip()
                if line.startswith("["):
                    # Check if we're in a substrate group section (SILICA or CALCITE)
                    reading_substrate = "SILICA" in line or "CALCITE" in line
                    continue
                if reading_substrate and line:
                    # Add atom indices to the substrate set
                    substrate_atoms.update(map(int, line.split()))

    # Create the "Other" group containing all atoms except substrate atoms
    other_atoms = [atom for atom in all_atoms if atom not in substrate_atoms]

    with open(ndx_output_path, "a") as f:
        f.write("\n[ Other ]\n")
        # Write atoms in groups of 15 per line
        for i in range(0, len(other_atoms), 15):
            f.write(" ".join(map(str, other_atoms[i : i + 15])) + "\n")

        f.write("\n[ System ]\n")
        # Write atoms in groups of 15 per line
        for i in range(0, len(all_atoms), 15):
            f.write(" ".join(map(str, all_atoms[i : i + 15])) + "\n")

    # Copy the entire experiment folder to the server
    shutil.copytree(
        output_path, osp.join(cfg.server_path, cfg.exp_folder), dirs_exist_ok=True
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", type=str, help="Path to the config file")
    # Parse known args and capture any additional --key value pairs as overrides
    args, unknown = parser.parse_known_args()

    # Convert unknown arguments of the form --key value into a dictionary
    overrides_dict = {}
    i = 0
    while i < len(unknown):
        token = unknown[i]
        if token.startswith("--"):
            key = token.lstrip("-")
            # Determine value (supports "--key value" and "--key" as True)
            value = True
            if i + 1 < len(unknown) and not unknown[i + 1].startswith("--"):
                value = unknown[i + 1]
                i += 1

            # Best-effort type casting: bool, int, float, else str
            if isinstance(value, str):
                lower_val = value.lower()
                if lower_val in {"true", "false"}:
                    value = lower_val == "true"
                else:
                    try:
                        if "." in value or "e" in lower_val or "E" in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except ValueError:
                        pass

            # Support dotted keys for nested overrides via dot-expansion
            # e.g., --substrate.Lx 25.0 -> {"substrate": {"Lx": 25.0}}
            target = overrides_dict
            parts = key.split(".")
            for part in parts[:-1]:
                if part not in target or not isinstance(target[part], dict):
                    target[part] = {}
                target = target[part]
            target[parts[-1]] = value
        i += 1

    main(args.config, overrides=overrides_dict if overrides_dict else None)
