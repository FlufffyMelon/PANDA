import numpy as np
import os
import json
import mdtraj as md
from omegaconf import OmegaConf
from hydra.utils import instantiate
from panda.assembler.build import build
from panda.assembler.mixer import mixer

OmegaConf.register_new_resolver("eval", eval)


def build_system(config_path):
    # Load config using OmegaConf
    cfg = OmegaConf.load(config_path)
    print(f"[Config] Loaded config from: {config_path}")
    print(f"[Config] Keys: {', '.join(list(cfg.keys()))}")
    print("[Config] Config loaded successfully.\n")

    # Substrate loading/generation
    substr_path = cfg.substrate
    traj = md.load(substr_path)

    # Check if unitcell is orthogonal
    box_lengths = traj.unitcell_lengths[0]
    assert np.all(box_lengths[:3] > 0), "Box lengths must be positive."
    if box_lengths.shape[0] > 3:
        assert np.allclose(box_lengths[3:], 0), (
            "Box should be orthorhombic (last 6 box components should be zero)."
        )

    # Update system size with actual values
    cfg.WIDTH_X, cfg.WIDTH_Y, cfg.HEIGHT = map(float, box_lengths[:3])

    # Adding the pore heigth to the substr unitcell
    system_size = np.array([cfg.WIDTH_X, cfg.WIDTH_Y, cfg.HEIGHT + cfg.H])
    traj.unitcell_lengths[0] = system_size

    # Reading names and densities of the compounds
    names = list(cfg.names)
    density = list(cfg.density)

    # Instantiate all regions from config
    regions = [instantiate(r) for r in cfg.regions]
    if "insertion_regions" in cfg:
        insertion_regions = [instantiate(r) for r in cfg.insertion_regions]
    else:
        insertion_regions = [instantiate(r) for r in cfg.regions]

    assert len(names) == len(density) == len(regions), (
        f"The size of the names, density and regions arrays do not match. names ({len(names)}), density ({len(density)}), regions ({len(regions)})"
    )
    assert len(regions) == len(insertion_regions), (
        f"The number of regions and insertation regions do not match. regions ({len(regions)}), insertion_regions ({len(insertion_regions)})"
    )

    # Counting the number of molecules of each type based on their density
    mol_numbers = list(
        np.round(
            np.array([regions[i].get_volume() * density[i] for i in range(len(names))])
        ).astype(int)
    )

    # Building system
    traj = build(
        cfg.mol_path,
        traj,
        names,
        mol_numbers,
        insertion_regions,
        insertion_limit=cfg.insertion_limit,
        rotation_limit=cfg.rotation_limit,
        package=cfg.package,
        insertion_attempts=cfg.insertion_attempts,
        min_dist2=cfg.min2,
        opt_dist2=cfg.opt2,
    )

    # Create output dir if needed
    output_path = os.path.join(cfg.output_dir, cfg.exp_folder)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Loging config
    with open(os.path.join(output_path, "config.json"), "w") as f:
        json.dump(OmegaConf.to_container(cfg, resolve=True), f, indent=4)

    # Writing structure in gro fileformat
    traj.save(os.path.join(output_path, cfg.system_name + ".gro"))

    # Mixing the system
    print()
    print("[Mixing] Mixing system...")
    traj = mixer(traj, 3 * cfg.min2**0.5 / 2, cfg.min2, cfg.opt2)
    traj.save(os.path.join(output_path, cfg.system_name + ".gro"))

    print()
    print(f"[Summary] System build complete. Output: {output_path}")
