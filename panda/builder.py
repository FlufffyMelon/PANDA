import numpy as np
import os
import json
import mdtraj as md
from omegaconf import OmegaConf
from hydra.utils import instantiate
from panda.assembler.build import build
from panda.assembler.mixer import mixer
from panda.utils import serialize_component_config

OmegaConf.register_new_resolver("eval", eval)


def build_system(config_path):
    # Load config using OmegaConf
    cfg = OmegaConf.load(config_path)
    print(f"[Config] Loaded config from: {config_path}")
    print(f"[Config] Keys: {', '.join(list(cfg.keys()))}")
    print("[Config] Config loaded successfully.\n")

    # Substrate loading/generation
    if cfg.get("substrate", None):
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

        # Adding the pore height to the substr unitcell
        if cfg.get("H", None):
            system_size = np.array([[cfg.WIDTH_X, cfg.WIDTH_Y, cfg.HEIGHT + cfg.H]])
        else:
            raise ValueError("H must be specified if substrate is provided.")

        traj.unitcell_lengths[0] = system_size
    else:
        traj = md.Trajectory(
            xyz=np.empty((0, 3)),
            topology=md.Topology(),
            unitcell_lengths=np.array([cfg.WIDTH_X, cfg.WIDTH_Y, cfg.HEIGHT]),
            unitcell_angles=np.array([90, 90, 90]),
        )

    # Collect components from config (manual YAML loading, merging with main config for interpolation)
    components = []
    config_dir = os.path.dirname(config_path)
    for comp_yaml in cfg.components:
        comp_path = os.path.join(config_dir, comp_yaml)
        comp_cfg = OmegaConf.load(comp_path)
        component = instantiate(comp_cfg, **cfg)
        components.append(component)

    assert len(components) > 0, "No components defined."

    # Building system
    traj = build(
        cfg.mol_path,
        traj,
        [c.name for c in components],
        [c.numbers for c in components],
        [c.insertion_region for c in components],
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

    # Helper to serialize region objects
    def region_to_dict(region):
        if region is None:
            return None
        d = {"type": type(region).__name__}
        if hasattr(region, "center"):
            center = getattr(region, "center")
            # Convert numpy arrays to lists
            if hasattr(center, "tolist"):
                center = center.tolist()
            d["center"] = center
        if hasattr(region, "borders"):
            borders = getattr(region, "borders")
            if hasattr(borders, "tolist"):
                borders = borders.tolist()
            d["borders"] = borders
        return d

    # Prepare a resolved, human-readable version of the config for saving
    config_to_save = OmegaConf.to_container(cfg, resolve=True)
    config_to_save["components"] = serialize_component_config(components)

    with open(os.path.join(output_path, "config.json"), "w") as f:
        json.dump(config_to_save, f, indent=4)

    # Writing structure in gro fileformat
    traj.save(os.path.join(output_path, cfg.system_name + ".gro"))

    # Mixing the system
    print()
    print("[Mixing] Mixing system...")
    traj = mixer(traj, 3 * cfg.min2**0.5 / 2, cfg.min2, cfg.opt2)
    traj.save(os.path.join(output_path, cfg.system_name + ".gro"))

    print()
    print(f"[Summary] System build complete. Output: {output_path}")
