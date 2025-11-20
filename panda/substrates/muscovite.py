import numpy as np
import os
import mdtraj as md
from mdtraj.core.topology import Topology
from panda.utils import apply_pbc
from .utils import base62_encode, CustomSubstrate, ase2mdtraj, angle
from ase.io import read
import ase


class MuscoviteSubstrate(CustomSubstrate):
    def __init__(
        self,
        unitcell_path: str,
        Lx: float,
        Ly: float,
        Lz: float,
        build: bool = True,
        hybridisation: int = 1,
    ):
        self.unitcell_path = unitcell_path
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.build = build
        self.hybridisation = hybridisation

        # These will be set during the generation process
        self.atoms = None
        self.unit_atom_types = None
        self.Nx = None
        self.Ny = None
        self.Nz = None
        self.atom_type_map = None
        self.filename = None

        self.gro_path, self.itp_path, self.ndx_path = self._generate()

    def hybridise(self, n: int):
        """
        Replace n Si atoms with Al atoms in each z-layer of the unit cell.
        The unit cell is split into 4 layers along the z-axis.
        In each layer, Si atoms are sorted by ascending y-coordinate,
        and the first n atoms are replaced with Al.

        Parameters
        ----------
        n : int
            Number of Si atoms to replace with Al in each layer (must be ≤ 4)
        """
        if n > 4:
            raise ValueError("Cannot replace more than 4 Si atoms per layer")

        if n == 0:
            return  # No substitution needed

        print(
            f"Hybridising muscovite structure: replacing {n} Si atoms with Al in each z-layer"
        )

        # Get positions and chemical symbols
        positions = self.atoms.get_positions()
        symbols = self.atoms.get_chemical_symbols()

        # Find all Si atoms
        si_indices = [i for i, s in enumerate(symbols) if s == "Si"]

        # Get the cell height (z-dimension)
        cell_height = self.atoms.cell[2, 2]  # Z-height of the unit cell

        # Define 4 layers along z-axis
        layer_height = cell_height / 4
        z_layers = [
            (0, layer_height),
            (layer_height, 2 * layer_height),
            (2 * layer_height, 3 * layer_height),
            (3 * layer_height, cell_height),
        ]

        # Process each z-layer
        for layer_idx, (z_min, z_max) in enumerate(z_layers):
            # Find Si atoms in this z-layer
            layer_si_indices = [
                i for i in si_indices if z_min <= positions[i][2] < z_max
            ]

            if not layer_si_indices:
                raise ValueError(
                    f"No Si atoms found in z-layer {layer_idx + 1} ({z_min:.2f}-{z_max:.2f} Å)"
                )

            # Sort by ascending y-coordinate
            layer_si_indices.sort(key=lambda i: positions[i][1])

            # Replace the first n Si atoms with Al
            for i in range(min(n, len(layer_si_indices))):
                idx = layer_si_indices[(2 * i) % len(layer_si_indices)]
                symbols[idx] = "Al"

        # Update the atoms object with the new symbols
        self.atoms.set_chemical_symbols(symbols)

    def _generate(self):
        """Generate the substrate GRO, ITP and NDX files."""
        print("\nGenerating substrate...")

        # Load the unit cell
        self.atoms = read(self.unitcell_path)

        # Apply hybridisation if requested
        if self.hybridisation > 0:
            self.hybridise(self.hybridisation)

        # Classify atoms (after hybridisation)
        self.unit_atom_types = self._classify_muscovite_atoms(self.atoms)

        # Calculate substrate dimensions
        unitcell_box = np.diag(self.atoms.cell) * 1e-1  # Convert from Å to nm
        self.Nx, self.Ny, self.Nz = np.round(
            np.clip(np.array([self.Lx, self.Ly, self.Lz]) / unitcell_box, 1, None)
        ).astype(int)

        # Determine file paths
        substr_folder, unitcell_filename = os.path.split(self.unitcell_path)
        unitcell_name = os.path.splitext(unitcell_filename)[0]
        self.filename = "_".join(unitcell_name.split("_")[:-1])

        # Define paths for all three files
        substr_name = f"{self.filename}_{self.Nx}x{self.Ny}x{self.Nz}.gro"
        substr_path = os.path.join(substr_folder, "gro", substr_name)

        substr_itp_name = f"{self.filename}_{self.Nx}x{self.Ny}x{self.Nz}.itp"
        substr_itp_path = os.path.join(substr_folder, "itp", substr_itp_name)

        substr_ndx_name = f"{self.filename}_{self.Nx}x{self.Ny}x{self.Nz}.ndx"
        substr_ndx_path = os.path.join(substr_folder, "ndx", substr_ndx_name)

        # Check if all three files exist
        all_files_exist = (
            os.path.isfile(substr_path)
            and os.path.isfile(substr_itp_path)
            and os.path.isfile(substr_ndx_path)
        )

        if self.build and all_files_exist:
            print(f"Using ready-made substrate files from {substr_folder}")
        else:
            # If any file is missing or build is False, force regeneration of all files
            if not all_files_exist:
                print(
                    "One or more substrate files missing. Forcibly generating all files..."
                )

            # Generate the substrate
            substr_name = self._generate_muscovite_substrate()
            substr_path = os.path.join(substr_folder, "gro", substr_name)

            # Generate ITP file using the same atom_type_map
            substr_itp_name = self._generate_muscovite_itp(substr_path)
            substr_itp_path = os.path.join(substr_folder, "itp", substr_itp_name)

            # Generate NDX file using the same atom_type_map
            substr_ndx_name = self._generate_muscovite_ndx(substr_path)
            substr_ndx_path = os.path.join(substr_folder, "ndx", substr_ndx_name)

        return substr_path, substr_itp_path, substr_ndx_path

    def _classify_muscovite_atoms(self, atoms: ase.Atoms):
        """
        Classify and rename atoms in muscovite structure.

        Parameters
        ----------
        atoms : ase.Atoms
            The atomic structure to process

        Returns
        -------
        dict
            Dictionary with atom names and types
        """
        # Constants for bonding distances
        O_SI_BOND_DIST = 2.0  # Angstrom
        O_AL_BOND_DIST = 2.1  # Angstrom

        # Extract atomic data
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        cell = atoms.get_cell()

        # Group atoms by element
        oxygen_indices = [i for i, s in enumerate(symbols) if s == "O"]
        si_indices = [i for i, s in enumerate(symbols) if s == "Si"]
        al_indices = [i for i, s in enumerate(symbols) if s == "Al"]
        k_indices = [i for i, s in enumerate(symbols) if s == "K"]
        h_indices = [i for i, s in enumerate(symbols) if s == "H"]

        # Extract positions for vectorized operations
        oxygen_pos = positions[oxygen_indices, :]
        si_pos = positions[si_indices, :]
        al_pos = positions[al_indices, :]

        # Calculate distances with PBC for non-orthorhombic cells
        def get_pbc_distances(pos1, pos2):
            cell_inv = np.linalg.inv(cell)
            frac_pos1 = np.dot(pos1, cell_inv)
            frac_pos2 = np.dot(pos2, cell_inv)
            frac_diff = frac_pos1[:, np.newaxis, :] - frac_pos2[np.newaxis, :, :]
            frac_diff = frac_diff - np.floor(frac_diff + 0.5)
            cart_diff = np.dot(frac_diff, cell)
            return np.linalg.norm(cart_diff, axis=2)

        # Determine Al coordination
        al_coordination = {}
        al_o_distances = get_pbc_distances(al_pos, oxygen_pos)

        for i, al_idx in enumerate(al_indices):
            neighbor_count = np.sum(al_o_distances[i] < O_AL_BOND_DIST)
            if neighbor_count == 4:
                al_coordination[al_idx] = "tetrahedral"
            elif neighbor_count == 6:
                al_coordination[al_idx] = "octahedral"
            else:
                raise ValueError(f"Al atom {al_idx} has {neighbor_count} neighbors")

        # Calculate O-Si and O-Al distances
        o_si_distances = get_pbc_distances(oxygen_pos, si_pos)
        o_al_distances = get_pbc_distances(oxygen_pos, al_pos)

        # Initialize atom types
        atom_types = {}

        # Assign types for non-oxygen atoms
        for si_idx in si_indices:
            atom_types[si_idx] = "st"  # Silicon tetrahedral

        for al_idx in al_indices:
            atom_types[al_idx] = (
                "at" if al_coordination[al_idx] == "tetrahedral" else "ao"
            )

        for k_idx in k_indices:
            atom_types[k_idx] = "K"

        # Assign type for hydrogen atoms
        for h_idx in h_indices:
            atom_types[h_idx] = "ho"  # Hydrogen

        # Classify oxygen atoms based on their bonding environment
        for i, o_idx in enumerate(oxygen_indices):
            # Count neighbors
            tetra_si = np.sum(o_si_distances[i] < O_SI_BOND_DIST)

            # Count Al neighbors by coordination
            al_neighbors = np.where(o_al_distances[i] < O_AL_BOND_DIST)[0]
            tetra_al = sum(
                1
                for j in al_neighbors
                if al_coordination[al_indices[j]] == "tetrahedral"
            )
            octa_al = len(al_neighbors) - tetra_al

            # Classify oxygen based on bonding environment
            # if tetra_si == 1 and tetra_al == 1:
            #     atom_types[o_idx] = "obt"  # Si–O–Al(t)
            # elif (tetra_si + tetra_al) >= 1 and octa_al >= 1:
            #     atom_types[o_idx] = "obo"  # (Si/Al(t))–O–Al(o)
            # elif tetra_si + tetra_al == 0 and octa_al >= 2:
            #     atom_types[o_idx] = "oh"  # O–Al(o)–Al(o), гидроксильный кислород
            # else:
            #     raise ValueError(
            #         f"Oxygen atom {o_idx} ({', '.join(map(str, np.round(positions[o_idx, :], 2)))} Å) "
            #         f"has {tetra_si} tetrahedral Si, {tetra_al} tetrahedral Al, {octa_al} octahedral Al neighbors"
            #     )

            # replace your original block with this
            tetra_total = tetra_si + tetra_al

            if tetra_total == 0:
                # No tetrahedral neighbors -> oxygen in octahedral sheet (hydroxyl O)
                if octa_al >= 2:
                    atom_types[o_idx] = "oh"
                else:
                    raise ValueError(
                        "Unusual: no tetra and fewer than 2 octa Al: treat as generic bridging O"
                    )
                    # atom_types[o_idx] = "ob"

            elif tetra_total == 2:
                # Bridging oxygen between two tetrahedra
                if tetra_si == 2:
                    atom_types[o_idx] = "ob"  # Si - O - Si
                elif tetra_si == 1 and tetra_al == 1:
                    atom_types[o_idx] = (
                        "obt"  # Si - O - Al(tetra)  (tetra substitution)
                    )
                elif tetra_al == 2:
                    atom_types[o_idx] = (
                        "obs"  # Al(tetra) - O - Al(tetra) (double substitution)
                    )
                else:
                    atom_types[o_idx] = "ob"  # fallback

            elif tetra_total == 1:
                # Apical oxygen (connects tetrahedral layer and octahedral sheet).
                # Classify by whether the single tetra neighbour is Al (substitution) or Si (normal).
                if tetra_al == 1:
                    atom_types[o_idx] = "obt"  # tetrahedral substitution present
                else:
                    atom_types[o_idx] = "ob"

            else:
                # Rare / unexpected configurations (e.g., tetra_total > 2)
                raise ValueError(
                    f"Oxygen atom {o_idx} ({', '.join(map(str, np.round(positions[o_idx, :], 2)))} Å) "
                    f"has {tetra_si} tetrahedral Si, {tetra_al} tetrahedral Al, {octa_al} octahedral Al neighbors"
                )

        # Calculate statistics
        atom_counts = {
            "st": 0,
            "at": 0,
            "ao": 0,
            "K": 0,
            "obt": 0,
            "ob": 0,
            "oh": 0,
            "ho": 0,
        }

        for atom_type in atom_types.values():
            atom_counts[atom_type] = atom_counts.get(atom_type, 0) + 1

        # Print concise summary
        tetra_al = sum(1 for v in al_coordination.values() if v == "tetrahedral")
        octa_al = sum(1 for v in al_coordination.values() if v == "octahedral")
        print(
            f"Classified muscovite atoms: st {len(si_indices)}, at {tetra_al}, ao {octa_al}, "
            f"K {atom_counts['K']}, obt {atom_counts['obt']}, ob {atom_counts['ob']}, oh {atom_counts['oh']}, ho {atom_counts['ho']}"
        )

        return atom_types

    def _add_hydrogen(self, traj, atom_type_map, residue):
        """
        Add hydrogen atoms to hydroxyl groups based on their position in the unit cell.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            The trajectory containing the structure.
        atom_type_map : dict
            Dictionary mapping atom indices to their types.
        residue : mdtraj.Topology.Residue
            The residue to add the hydrogen atoms to.

        Returns
        -------
        mdtraj.Trajectory
            The trajectory with added hydrogen atoms.
        """
        # Find hydroxyl oxygen atoms (type 'oh')
        oh_indices = [idx for idx, atype in atom_type_map.items() if atype == "oh"]
        if not oh_indices:
            print("No hydroxyl oxygen atoms found, skipping hydrogen addition")
            return traj

        print(f"Adding hydrogen atoms to {len(oh_indices)} hydroxyl groups")

        # OH bond length (in nm)
        OH_BOND_LENGTH = 0.095

        # Create arrays for new hydrogen atoms
        h_coords = np.zeros((len(oh_indices), 3))

        # Calculate unit cell dimensions
        box = np.diag(traj.unitcell_vectors[0, :, :])
        unit_cell_box = box / np.array([self.Nx, self.Ny, self.Nz])
        unit_cell_height = unit_cell_box[2]

        # Layer boundaries within a unit cell
        layer1 = unit_cell_height / 4
        layer2 = unit_cell_height / 2
        layer3 = 3 * unit_cell_height / 4

        # Add H atoms based on O position within their unit cell
        for i, idx in enumerate(oh_indices):
            o_pos = traj.xyz[0, idx]

            # Calculate which unit cell this oxygen belongs to
            cell_x = int(o_pos[0] / unit_cell_box[0])
            cell_y = int(o_pos[1] / unit_cell_box[1])
            cell_z = int(o_pos[2] / unit_cell_box[2])

            # Calculate position within the unit cell
            rel_pos = o_pos - np.array(
                [
                    cell_x * unit_cell_box[0],
                    cell_y * unit_cell_box[1],
                    cell_z * unit_cell_box[2],
                ]
            )

            # Get relative z-position within the unit cell
            rel_z = rel_pos[2]

            # Determine H direction based on relative position within unit cell
            if (0 <= rel_z < layer1) or (layer2 <= rel_z < layer3):
                # Add H above O
                direction = np.array([0, 0, 1])
            else:
                # Add H below O
                direction = np.array([0, 0, -1])

            # Calculate H position
            h_coords[i] = o_pos + direction * OH_BOND_LENGTH

        # Create a new topology with H atoms
        new_top = traj.topology

        # Add hydrogen atoms to the topology with "ho" type name
        for i in range(len(oh_indices)):
            # Use the type name directly as the atom name
            new_top.add_atom(
                "ho",  # Use the type name
                element="H",
                residue=residue,
                serial=traj.n_atoms + i + 1,
            )

        # Combine coordinates
        new_xyz = np.concatenate(
            [traj.xyz, h_coords.reshape(1, len(oh_indices), 3)], axis=1
        )

        # Create new trajectory
        new_traj = md.Trajectory(
            xyz=new_xyz,
            topology=new_top,
            unitcell_lengths=traj.unitcell_lengths,
            unitcell_angles=traj.unitcell_angles,
        )

        return new_traj

    def _generate_muscovite_substrate(self):
        """
        Generate a substrate from a given unitcell.

        Returns
        -------
        str
            The path to the file of the generated substrate.
        """
        unitcell_box = np.diag(self.atoms.cell) * 1e-1  # Convert from Å to nm

        # Set output paths
        folder, _ = os.path.split(self.unitcell_path)
        output_name = f"{self.filename}_{self.Nx}x{self.Ny}x{self.Nz}.gro"
        output_path = os.path.join(folder, "gro", output_name)

        print(f"Generating GRO file for substrate: {output_name}")
        print(
            f"Nearest dimensions of the substrate: {unitcell_box[0] * self.Nx:.1f}x{unitcell_box[1] * self.Ny:.1f}x{unitcell_box[2] * self.Nz:.1f} nm ({self.Nx}x{self.Ny}x{self.Nz})..."
        )

        # Create supercell
        supercell = self.atoms.repeat((self.Nx, self.Ny, self.Nz))

        # Convert to mdtraj for further processing
        traj = ase2mdtraj(supercell, residue_name="MUSC")

        # First, create atoms with their type names to preserve type information
        residue = traj.topology.add_residue("MUSC", traj.topology.add_chain(), resSeq=1)

        # Store temporary atom types mapped to their atom indices
        temp_atom_types = {}
        unit_cell_size = len(self.atoms)

        # Initial atom setup - assign type names as atom names
        for idx in range(traj.n_atoms):
            atom = traj.topology.atom(idx)
            # Find the corresponding atom in the unitcell
            unit_cell_idx = idx % unit_cell_size
            atom_type = self.unit_atom_types[unit_cell_idx]

            # Store the type for this atom
            temp_atom_types[idx] = atom_type

            # Use atom type as the name temporarily
            atom.name = atom_type
            atom.residue = residue

        # Add hydrogen atoms to hydroxyl groups
        traj = self._add_hydrogen(traj, temp_atom_types, residue)

        # Set box
        box = np.diag(traj.unitcell_vectors[0, :, :])
        shift = np.array([0.0, 0.0, -0.66])  # shift by -0.605 nm along z-axis
        all_xyz = traj.xyz[0, :, :] + shift
        all_xyz = apply_pbc(all_xyz, box)

        # Remove top layer of K atoms and create final trajectory
        offset = 0.2
        new_top = Topology()
        new_top.add_residue("", new_top.add_chain())
        residue = new_top.add_residue("MUSC", new_top.add_chain())

        # Keep track of the atoms we keep and their types
        kept_atoms = []
        kept_types = {}

        for idx in range(traj.n_atoms):
            # Skip K atoms above the cutoff
            if traj.topology.atom(idx).name == "K":
                if all_xyz[0, idx, 2] > box[2] - offset:
                    continue

            # Store the original index and its type
            kept_atoms.append(idx)

            # Get type from atom name (which we temporarily set to the type)
            atom_type = traj.topology.atom(idx).name
            kept_types[len(kept_atoms) - 1] = atom_type

            new_top.add_atom(
                atom_type,  # Keep the type name for now
                element=None,
                residue=residue,
                serial=len(kept_atoms),
            )

        print(f"Removed {traj.n_atoms - len(kept_atoms)} K atoms")
        new_xyz = all_xyz[:, all_xyz[0, :, 2] < box[2] - offset, :]
        new_box = np.array([box[0], box[1], np.max(new_xyz[0, :, 2])])

        # Create trajectory with type names
        traj = md.Trajectory(
            new_xyz.reshape(1, -1, 3),
            new_top,
            unitcell_lengths=new_box.reshape(1, 3),
            unitcell_angles=np.array([[90.0, 90.0, 90.0]]),
        )

        # Now rename atoms with element names + unique IDs
        final_top = Topology()
        final_top.add_residue("", final_top.add_chain())
        residue = final_top.add_residue("MUSC", final_top.add_chain())

        # Create counters for each element
        counter_dict = {"Si": 0, "Al": 0, "K": 0, "O": 0, "H": 0}

        # Create a new atom_type_map with final indices
        self.atom_type_map = {}

        # Create atoms with proper element names
        for i, atom in enumerate(traj.topology.atoms):
            atom_type = atom.name  # The type name from earlier
            self.atom_type_map[i] = atom_type

            # Determine element from atom type
            if atom_type.startswith("st"):
                element = "Si"
            elif atom_type.startswith("at") or atom_type.startswith("ao"):
                element = "Al"
            elif atom_type.startswith("o"):  # obt, obo, oh
                element = "O"
            elif atom_type.startswith("ho"):
                element = "H"
            else:
                element = "K"

            # Update counter and create unique name
            counter_dict[element] += 1
            atom_name = element + base62_encode(
                counter_dict[element],
                length=5 - len(element),
            )

            # Add atom with element name
            final_top.add_atom(
                atom_name,
                element=element,
                residue=residue,
                serial=i + 1,
            )

        # Create final trajectory with proper atom names
        final_traj = md.Trajectory(
            traj.xyz,
            final_top,
            unitcell_lengths=traj.unitcell_lengths,
            unitcell_angles=traj.unitcell_angles,
        )

        # Write .gro file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        final_traj.save_gro(output_path)
        print("Substrate successfully created!")

        return output_name

    def _generate_muscovite_ndx(self, substr_path):
        """
        Generate .ndx file for a muscovite substrate, based on the .gro file.
        Creates groups for each atom type: st, at, ao, K, obt, obo, oh, H.

        Parameters
        ----------
        substr_path : str
            The path to the .gro file of the substrate.

        Returns
        -------
        str
            The name of the generated .ndx file.
        """
        # Setup paths
        folder, substr_filename = os.path.split(substr_path)
        folder_root = os.path.split(folder)[0]
        substr_name = os.path.splitext(substr_filename)[0]
        output_name = substr_name + ".ndx"
        output_path = os.path.join(folder_root, "ndx", output_name)

        print(f"Generating NDX file for substrate: {output_name}")

        # Load the substrate to ensure it's valid
        md.load(substr_path)

        # Create atom type groups
        atom_groups = {
            "st": [],  # Silicon tetrahedral
            "at": [],  # Aluminum tetrahedral
            "ao": [],  # Aluminum octahedral
            "K": [],  # Potassium
            "obt": [],  # Oxygen bridging tetrahedral
            "ob": [],  # Oxygen bridging octahedral
            "oh": [],  # Oxygen hydroxyl
            "ho": [],  # Hydrogen
        }

        # Use the atom_type_map that was created during substrate generation
        # which already accounts for hydrogen atoms and removed K atoms
        for i, atom_type in self.atom_type_map.items():
            atom_idx = i + 1  # 1-indexed for GROMACS
            atom_groups[atom_type].append(atom_idx)

        # Write the index file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            for group_name, atoms in atom_groups.items():
                if not atoms:
                    continue

                f.write(f"[ {group_name} ]\n")
                # Write atoms in groups of 15 per line
                for i in range(0, len(atoms), 15):
                    f.write(" ".join(map(str, atoms[i : i + 15])) + "\n")
                f.write("\n")

        return output_name

    def _generate_muscovite_itp(self, substr_path):
        """
        Generate .itp file for a muscovite substrate, based on the .gro file.

        Parameters
        ----------
        substr_path : str
            The path to the .gro file of the substrate.

        Returns
        -------
        str
            The name of the generated .itp file.
        """
        # --- Setup paths ---
        folder, substr_filename = os.path.split(substr_path)
        folder_root = os.path.split(folder)[0]
        substr_name = os.path.splitext(substr_filename)[0]
        output_name = substr_name + ".itp"
        output_path = os.path.join(folder_root, "itp", output_name)

        print(f"Generating ITP file for substrate: {output_name}")

        # Load the substrate
        substr = md.load(substr_path)

        # --- Atom type definitions ---
        # metadata = {
        #     "st": ["st", 28.086, 2.1000],  # Silicon tetrahedral
        #     "at": ["at", 26.982, 1.5750],  # Aluminum tetrahedral
        #     "ao": ["ao", 26.982, 1.5750],  # Aluminum octahedral
        #     "K": ["K", 39.098, 1.0000],  # Potassium
        #     "obt": ["obt", 15.999, -1.1688],  # Oxygen bridging tetrahedral
        #     "ob": ["ob", 15.999, -1.0500],  # Oxygen bridging octahedral
        #     "oh": ["oh", 15.999, -0.9500],  # Oxygen hydroxyl
        #     "ho": ["ho", 1.008, 0.4250],  # Hydrogen
        # }
        metadata = {
            "st": ["st", 28.086, 2.1000],  # Silicon tetrahedral
            "at": ["at", 26.982, 1.575],  # Aluminum tetrahedral
            "ao": ["ao", 26.982, 1.575],  # Aluminum octahedral
            "K": ["K", 39.098, 1.0000],  # Potassium
            "obt": ["obt", 15.999, -1.11875],  # Oxygen bridging tetrahedral
            "ob": ["ob", 15.999, -1],  # Oxygen bridging octahedral
            "oh": ["oh", 15.999, -0.9500],  # Oxygen hydroxyl
            "ho": ["ho", 1.008, 0.4250],  # Hydrogen
        }

        # Use the atom_type_map for topology generation
        bonds = []
        angles = []

        # Find atoms by type
        oh_indices = [i for i, atype in self.atom_type_map.items() if atype == "oh"]
        h_indices = [i for i, atype in self.atom_type_map.items() if atype == "ho"]
        ao_indices = [i for i, atype in self.atom_type_map.items() if atype == "ao"]

        # Verify we have necessary atom types
        if not oh_indices:
            raise ValueError("No hydroxyl oxygen atoms (oh) found in substrate")
        if not h_indices:
            raise ValueError("No hydrogen atoms (ho) found in substrate")
        if not ao_indices:
            raise ValueError("No octahedral aluminum atoms (ao) found in substrate")

        # Get coordinates
        oh_coords = substr.xyz[0, oh_indices]
        h_coords = substr.xyz[0, h_indices]
        ao_coords = substr.xyz[0, ao_indices]

        # Cutoff distances
        OH_CUTOFF = 0.12  # nm
        AL_O_CUTOFF = (
            0.22  # nm for Al-O bonds (slightly larger than typical bond length)
        )

        # Get cell information for PBC calculations
        box = substr.unitcell_lengths[0]
        cell = substr.unitcell_vectors[0]

        # Define PBC distance calculation for non-orthogonal cells
        def calculate_pbc_distances(pos1, pos_array):
            """Calculate PBC distances from one point to an array of points"""
            cell_inv = np.linalg.inv(cell)
            frac_pos1 = np.dot(pos1, cell_inv)
            frac_pos2 = np.dot(pos_array, cell_inv)
            frac_diff = frac_pos1 - frac_pos2
            # Apply minimum image convention
            frac_diff = frac_diff - np.round(frac_diff)
            # Convert back to Cartesian
            cart_diff = np.dot(frac_diff, cell)
            return np.linalg.norm(cart_diff, axis=1)

        # Find O-H bonds and relevant angles
        oh_h_pairs = []  # Store (oh_idx, h_idx) pairs for angle calculation

        for h_idx, h_pos in enumerate(h_coords):
            # Find nearest hydroxyl oxygen with PBC
            distances = calculate_pbc_distances(h_pos, oh_coords)
            nearest_oh_idx = np.argmin(distances)

            if distances[nearest_oh_idx] < OH_CUTOFF:
                oh_atom_idx = oh_indices[nearest_oh_idx]
                h_atom_idx = h_indices[h_idx]
                bonds.append((oh_atom_idx, h_atom_idx))
                oh_h_pairs.append((nearest_oh_idx, h_idx))
            else:
                raise ValueError(
                    f"Hydrogen atom {h_indices[h_idx]} is not bonded to any hydroxyl oxygen (closest distance: {distances[nearest_oh_idx]:.3f} nm)"
                )

        # Find Al-O-H angles
        for oh_idx_local, h_idx_local in oh_h_pairs:
            oh_idx = oh_indices[oh_idx_local]
            h_idx = h_indices[h_idx_local]
            oh_pos = oh_coords[oh_idx_local]

            # Find nearest aluminum octahedral atoms using PBC
            ao_found = []
            if len(ao_indices) > 0:
                # Calculate distances with PBC
                ao_distances = calculate_pbc_distances(oh_pos, ao_coords)
                ao_nearby_indices = np.where(ao_distances < AL_O_CUTOFF)[0]

                for i in ao_nearby_indices:
                    ao_idx = ao_indices[i]
                    ao_found.append((ao_idx, ao_distances[i]))

            # Process octahedral aluminum atoms (ao-oh-ho)
            if len(ao_found) > 0:
                # Sort by y-coordinate if there are multiple
                if len(ao_found) > 1:
                    # We need to consider y-coordinates with PBC
                    ao_y_coords = []
                    for idx, _ in ao_found:
                        # Calculate the minimum image y-coordinate difference with PBC
                        cell_inv = np.linalg.inv(cell)
                        ao_pos = substr.xyz[0, idx]
                        frac_oh = np.dot(oh_pos, cell_inv)
                        frac_ao = np.dot(ao_pos, cell_inv)
                        frac_diff = frac_ao - frac_oh
                        # Apply minimum image convention in fractional coords
                        frac_diff = frac_diff - np.round(frac_diff)
                        # Convert back to Cartesian to get proper y-coordinate
                        cart_diff = np.dot(frac_diff, cell)
                        # We're interested in y-coordinate
                        ao_y_coords.append(oh_pos[1] + cart_diff[1])

                    # Choose the one with smaller y-coordinate
                    chosen_idx = np.argmin(ao_y_coords)
                    ao_idx = ao_found[chosen_idx][0]
                else:
                    ao_idx = ao_found[0][0]

                # Calculate angle to verify it's > 90 degrees
                angle_value = angle(ao_idx, oh_idx, h_idx, substr.xyz[0], box)
                if angle_value < 90:
                    raise ValueError(
                        f"ao-oh-ho angle {ao_idx}-{oh_idx}-{h_idx} = {angle_value:.1f}° is <= 90°"
                    )
                angles.append((ao_idx, oh_idx, h_idx))
            else:
                raise ValueError(
                    f"No octahedral aluminum atoms (ao) found for oh-ho angle {oh_idx}-{h_idx}"
                )

        print(f"Found {len(bonds)} bonds and {len(angles)} angles")

        # --- Format ITP file sections ---
        # Atoms section
        atoms_text = ""
        total_charge = 0.0

        for i, atom_type in self.atom_type_map.items():
            atom_name = substr.topology.atom(i).name
            type_name, mass, charge = metadata[atom_type]
            total_charge += charge
            atoms_text += "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}\n".format(
                i + 1, type_name, 1, "MUSC", atom_name, 1, charge, mass
            )

        # Print total charge for verification
        print(f"Total substrate charge: {total_charge:.6f}")

        # Bonds section
        bonds_text = ""
        for i, j in bonds:
            bonds_text += "{:>8}{:>8}{:>8}\n".format(i + 1, j + 1, 1)

        # Angles section
        angles_text = ""
        for i, j, k in angles:
            angles_text += "{:>8}{:>8}{:>8}{:>8}\n".format(i + 1, j + 1, k + 1, 1)

        # --- Generate final ITP file ---
        # Check if template exists, if not create a basic template
        template_path = os.path.join(folder_root, self.filename + "_template" + ".itp")
        # Read template
        with open(template_path) as f:
            template_text = f.read()

        # Format the template
        final_text = template_text.format(
            atoms=atoms_text, bonds=bonds_text, angles=angles_text
        )

        # Write output file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            f.write(final_text)

        return output_name
