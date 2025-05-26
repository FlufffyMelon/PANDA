import os
from panda.substr import generate_substrate, generate_calcite_itp

class CustomSubstrate:
    def __init__(self, unitcell, Lx, Ly, Lz, freeze_substr=True, build=True):
        self.unitcell = unitcell
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.freeze_substr = freeze_substr
        self.build = build
        self.gro_path, self.itp_path = self._generate()

    def _generate(self):
        substr_name = generate_substrate(
            self.unitcell,
            self.Lx,
            self.Ly,
            self.Lz,
            freeze_substr=self.freeze_substr,
            build=self.build,
        )
        substr_folder, _ = os.path.split(self.unitcell)
        substr_path = os.path.join(substr_folder, "gro", substr_name)

        substr_itp_name = generate_calcite_itp(substr_path, build=self.build)
        substr_itp_path = os.path.join(substr_folder, "itp", substr_itp_name)

        return substr_path, substr_itp_path
