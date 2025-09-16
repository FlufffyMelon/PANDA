import numpy as np
import os
import os.path as osp
import shutil
import json
import argparse
from omegaconf import OmegaConf
from panda.builder import build_system
from panda.substr import generate_substrate, generate_calcite_itp
from panda.parser import parse_C6_C12

# OmegaConf.register_new_resolver("eval", eval)


def main(config_path):
    # Load config using OmegaConf
    cfg = OmegaConf.load(config_path)

    # Generating substrate with target dimensions
    substr_gro_path = generate_substrate(
        cfg.substrate_unitcell, cfg.WIDTH_X, cfg.WIDTH_Y, cfg.HEIGHT
    )

    # Generating itp for substrate if needed
    substr_itp_path = generate_calcite_itp(substr_gro_path)

    output_path = os.path.join(cfg.output_dir, cfg.exp_folder)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    cfg.substrate = substr_gro_path
    with open(os.path.join(output_path, "config.json"), "w") as f:
        json.dump(OmegaConf.to_container(cfg, resolve=True), f, indent=4)

    # Start building system
    build_system(os.path.join(output_path, "config.json"))

    # Update config after building system
    cfg = OmegaConf.load(osp.join(output_path, "config.json"))

    # Generating system.itp
    with open(os.path.join(output_path, "system.itp"), "w") as f:
        for component_cfg in cfg.components:
            f.write(f'#include "{component_cfg.name}.itp"\n')

        f.write(f'#include "{osp.split(substr_itp_path)[1]}"\n')

        f.write(f"\n[ system ]\n{cfg.system_name}\n")
        f.write("\n[ molecules ]\n; molecule name\tnr.\n")
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

gmx grompp -f {osp.split(cfg.pipeline[0])[1]} -c {cfg.system_name}_init.gro -p {osp.split(cfg.topology)[1]} -o {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[0])[1])[0].split("_")[-1]} -maxwarn 10
mpirun -np {cfg.n_mpi} --cpu-set {cfg.init_core}-{cfg.init_core + cfg.n_mpi - 1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[0])[1])[0].split("_")[-1]} -ntomp 1 -nb gpu -gpu_id {cfg.gpu_id}
rm ./*pdb

gmx grompp -f {osp.split(cfg.pipeline[1])[1]} -c {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[0])[1])[0].split("_")[-1]}.gro -p {osp.split(cfg.topology)[1]} -o {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[1])[1])[0].split("_")[-1]} -maxwarn 10
mpirun -np {cfg.n_mpi} --cpu-set {cfg.init_core}-{cfg.init_core + cfg.n_mpi - 1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[1])[1])[0].split("_")[-1]} -ntomp 1 -nb gpu -gpu_id {cfg.gpu_id} -dlb yes
rm ./*pdb

gmx grompp -f {osp.split(cfg.pipeline[2])[1]} -c {cfg.system_name}_{osp.splitext(osp.split(cfg.pipeline[1])[1])[0].split("_")[-1]}.gro -p {osp.split(cfg.topology)[1]} -o {cfg.system_name} -maxwarn 10
mpirun -np {cfg.n_mpi} --cpu-set {cfg.init_core}-{cfg.init_core + cfg.n_mpi - 1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {cfg.system_name} -ntomp 1 -nb gpu -gpu_id {cfg.gpu_id} -dlb yes
rm ./*pdb""")

    # Generating mdp files from templates
    with open(cfg.pipeline[0]) as f:
        steep_final_text = f.read().format(**{"freeze": ""})
    with open(osp.join(output_path, osp.split(cfg.pipeline[0])[1]), "w") as f:
        f.write(steep_final_text)

    with open(cfg.pipeline[1]) as f:
        short_final_text = f.read().format(
            **{
                "freeze": "",
                "temp": cfg.temp,
            }
        )
    with open(osp.join(output_path, osp.split(cfg.pipeline[1])[1]), "w") as f:
        f.write(short_final_text)

    with open(cfg.pipeline[2]) as f:
        run_final_text = f.read().format(
            **{
                "nsteps": cfg.nsteps,
                "freeze": "",
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
        substr_itp_path,
        osp.join(output_path, osp.split(substr_itp_path)[1]),
    )

    # Copy the entire experiment folder to the server
    shutil.copytree(output_path, osp.join(cfg.server_path, cfg.exp_folder))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", "-c", type=str, help="Path to the config file")
    args = parser.parse_args()

    main(args.config)
