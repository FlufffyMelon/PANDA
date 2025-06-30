![Logo](logo.png)

# PANDA: Predicting Angle from Nanoscale Density Analysis

<p align="left">
  <a href="https://doi.org/10.1016/j.colsurfa.2024.135994"><img alt="Journal" src="https://img.shields.io/badge/journal-Colloids%20%26%20Surfaces%20A-blue"></a>
  <a href="https://doi.org/10.1016/j.colsurfa.2024.135994"><img alt="Article" src="https://img.shields.io/badge/article-Elsevier-green"></a>
  <a href="https://img.shields.io/badge/license-MIT-brightgreen"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-brightgreen"></a>
  <a href="https://img.shields.io/badge/python-%3E%3D3.8-blue"><img alt="Python >=3.8" src="https://img.shields.io/badge/python-%3E%3D3.8-blue"></a>
</p>

**A Python package for calculating contact angles from molecular simulations using one-dimensional density profiles.**

## 🚀 Installation

Install PANDA directly from GitHub:

```bash
pip install git+https://github.com/FlufffyMelon/PANDA.git@main
```

## 🧩 Basic Usage

Import PANDA and use its core functions:

```python
import panda
```

### Main Functions

- **`get_each_density_profile()`**
  Extracts the density profile for each frame of a molecular dynamics trajectory (currently only XTC format is supported).
  **Typical usage:**
  ```python
  axises, denses = panda.get_each_density_profile(
      trajectory_file, topology_file, residue, sl, chunk_length,
      begin_time, time, timestep, units
  )
  ```

- **`block_average_density_profile()`**
  Computes block-averaged density profiles to reduce noise and improve statistical reliability.
  **Typical usage:**
  ```python
  axises_avg, denses_avg = panda.block_average_density_profile(axises, denses, block_length)
  ```

- **`profile_approx()`**
  Fits a density profile and extracts the contact angle and other parameters.
  **Typical usage:**
  ```python
  axis_fit, dens_fit, result = panda.profile_approx_from_array(
      denses_avg, axises_avg, rho_bulk, l, phi, H,
      interface_type=interface_type, display=False
  )
  angle_deg = np.rad2deg(result['theta'])
  ```

## 🛠️ System Builder

PANDA includes a **system builder** for assembling molecular systems from configuration files. The builder uses a main setup file (e.g., `setup.yaml`) that references component files (e.g., `decane.yaml`, `tip4p.yaml`).

**Usage:**
```bash
panda build --config <path-to-setup.yaml>
```
- `<path-to-setup.yaml>`: Path to your main YAML configuration file describing the system setup and listing component files.

## 📒 Examples

- **System Generation:**
  Explore the [`examples/building_systems`](examples/building_systems) folder for scripts and configuration files demonstrating automated system assembly with the builder.

- **Contact Angle Calculation:**
  The Jupyter notebook for the full PANDA workflow (density profile extraction, block averaging, contact angle calculation) is now located in [`examples/contact_angle/calcite_decane_contact_angle.ipynb`](examples/contact_angle/calcite_decane_contact_angle.ipynb).

## 📖 Reference

If you use PANDA in your research, please cite:

[**Semenchuk, A. A., N. D. Kondratyuk, and I. V. Kopanichuk. "PANDA: Predicting angle from nanoscale density analysis." Colloids and Surfaces A: Physicochemical and Engineering Aspects 708 (2025): 135994.**](https://doi.org/10.1016/j.colsurfa.2024.135994)

**BibTeX:**
```bibtex
@article{semenchuk2025panda,
  title={PANDA: Predicting angle from nanoscale density analysis},
  author={Semenchuk, AA and Kondratyuk, ND and Kopanichuk, IV},
  journal={Colloids and Surfaces A: Physicochemical and Engineering Aspects},
  volume={708},
  pages={135994},
  year={2025},
  publisher={Elsevier}
}
```

## 📝 License
MIT

