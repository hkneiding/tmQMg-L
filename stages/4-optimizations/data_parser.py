import re


class DataParser:

    """Class for reading relevant data from Gaussian output files into python dictionaries."""

    def __init__(self, file_path):

        """Constructor

        Args:
            file_path (string): Path to the output file.
        """

        self.file_path = file_path
        with open(file_path, 'r') as f:
            self.lines = f.read().strip().split('\n')
        self.n_atoms = self._get_number_of_atoms()

    def _get_number_of_atoms(self):

        for i in range(len(self.lines)):
            # find line that contain atom number
            if 'NAtoms=' in self.lines[i]:

                # get position in line
                line_split = self.lines[i].split()
                n_atomsIndex = line_split.index('NAtoms=') + 1

                # return
                return int(line_split[n_atomsIndex])

        raise Exception('Could not find number of atoms in file.')

    def parse(self):

        """Iterates through a given Gaussian output file and returns a dict of extracted data.

        Returns:
            dict[]: A dict containing different information extracted from the Gaussian output file.
        """

        if self._has_failed():
            return None

        # output variable
        qm_data = {}
        qm_data['n_atoms'] = self.n_atoms

        # get id token from file name
        qm_data['id'] = '/'.join(self.file_path.split('/')[-2:])
        for i in range(len(self.lines)):

            # search for keywords and if found call appropriate functions with start index
            # the start index addition offset is based on the Gaussian output format

            if 'Standard orientation' in self.lines[i]:
                qm_data['atomic_numbers'] = self._extract_atomic_numbers(i + 5)
                qm_data['geometric_data'] = self._extract_geometric_data(i + 5)

            if 'Summary of Natural Population Analysis' in self.lines[i]:
                qm_data['natural_atomic_charges'], qm_data['natural_electron_population'] = self._extract_natural_atomic_charges(i + 6)

            if 'Charge = ' in self.lines[i]:
                qm_data['charge'] = self._extract_charge(i)

            if 'Stoichiometry' in self.lines[i]:
                qm_data['stoichiometry'] = self._extract_stoichiometry(i)

            if 'SCF Done' in self.lines[i]:
                qm_data['electronic_energy'] = self._extract_electronic_energy(i)

            if 'Grimme-D3(BJ) Dispersion energy=' in self.lines[i]:
                qm_data['dispersion_energy'] = self._extract_dispersion_energy(i)

            if 'Dipole moment (field-independent basis, Debye)' in self.lines[i]:
                qm_data['dipole_moment'] = self._extract_dipole_moment(i + 1)

            if 'Atomic contributions to Alpha molecular orbitals:' in self.lines[i]:
                qm_data['molecular_orbital_data'] = self._extract_molecular_orbital_data(i + 1)

            if 'Alpha virt. eigenvalues' in self.lines[i]:
                pass

            if 'Isotropic polarizability' in self.lines[i]:
                qm_data['isotropic_polarisability'] = self._extract_polarisability(i)

            if 'Frequencies -- ' in self.lines[i]:
                if 'frequencies' not in qm_data:
                    qm_data['frequencies'] = self._extract_frequency(i)
                else:
                    qm_data['frequencies'].extend(self._extract_frequency(i))

            if 'Molar volume' in self.lines[i]:
                qm_data['molar_volume'] = self._extract_molar_volume(i)

            if 'Principal axes' in self.lines[i]:
                qm_data['principal_moments'] = self._extract_principal_moments(i + 2)

        return qm_data

    def _has_failed(self):

        if 'Normal termination' not in self.lines[-1]:
            return True

        return False

    # - - - extraction functions - - - #

    # Some of the following extraction functions are redundant in the sense that for some properties
    # the extraction procedures are identical. The distinction between these functions is kept
    # nonetheless to ensure maintainability (e.g. when the Gaussian output format changes).

    def _extract_frequency(self, start_index: int):

        line_split = self.lines[start_index].split()
        return list(map(float, line_split[2:]))

    def _extract_polarisability(self, start_index: int):

        line_split = self.lines[start_index].split()
        return float(line_split[5])

    def _extract_charge(self, start_index: int):

        line_split = self.lines[start_index].split()
        return int(line_split[2])

    def _extract_stoichiometry(self, start_index: int):

        line_split = self.lines[start_index].split()
        return line_split[1]

    def _extract_molecular_mass(self, start_index: int):

        line_split = self.lines[start_index].split()
        return float(line_split[2])

    def _extract_dispersion_energy(self, start_index: int):

        line_split = self.lines[start_index].split()
        return float(line_split[4])

    def _extract_electronic_energy(self, start_index: int):

        line_split = self.lines[start_index].split()
        return float(line_split[4])

    def _extract_dipole_moment(self, start_index: int):

        line_split = self.lines[start_index].split()
        return float(line_split[7])

    def _extract_molar_volume(self, start_index: int):

        line_split = self.lines[start_index].split()
        return float(line_split[3])

    def _extract_principal_moments(self, start_index: int):
        
        return self.lines[start_index]

        line_split = self.lines[start_index].split()

        if len(line_split) != 5:
            return None

        return [float(line_split[2]), float(line_split[3]), float(line_split[4])]

    def _extract_atomic_numbers(self, start_index: int):

        atomic_numbers = [int(self.lines[i].split()[1]) for i in range(start_index, start_index + self.n_atoms, 1)]

        return atomic_numbers

    def _extract_geometric_data(self, start_index: int):

        geometric_data = []

        for i in range(start_index, start_index + self.n_atoms, 1):
            # split line at any white space
            line_split = self.lines[i].split()
            # read out data (index number based on Gaussian output format)
            xyz = [float(line_split[3]), float(line_split[4]), float(line_split[5])]
            geometric_data.append(xyz)

        # also return index i to jump ahead in the file
        return geometric_data

    def _extract_natural_atomic_charges(self, start_index):

        natural_atomic_charges = [float(self.lines[i].split()[2]) for i in range(start_index, start_index + self.n_atoms, 1)]
        natural_electron_population = [[float(self.lines[i].split()[j]) for j in range(3, 6)] for i in range(start_index, start_index + self.n_atoms, 1)]

        # also return index i to jump ahead in the file
        return natural_atomic_charges, natural_electron_population

    def _extract_molecular_orbital_data(self, start_index):

        i = start_index

        molecular_orbital_data = []
        while 'Orbital energies and kinetic energies (alpha):' not in self.lines[i]:

            # make new entry
            if 'occ' in self.lines[i] or 'vir' in self.lines[i]:

                combined = self.lines[i]

                j = i + 1
                while 'occ' not in self.lines[j] and 'vir' not in self.lines[j] and self.lines[j].strip() != '':
                    combined += self.lines[j][1:]
                    j += 1

                line_split = combined.split()

                # extract id and energy
                mo_spin = line_split[0]
                mo_type = line_split[1]
                mo_id = line_split[2]
                energy = float(line_split[3].split('=')[-1])

                # split
                line_split = line_split[5:]
                molecular_orbital_occupations = {}
                for item in line_split:
                    item_split = item.split('=')
                    molecular_orbital_occupations[item_split[0]] = float(item_split[1])

                molecular_orbital_data.append({
                    'id': mo_id,
                    'spin': mo_spin,
                    'type': mo_type,
                    'energy': energy,
                    'occupations': molecular_orbital_occupations
                })
            i += 1

        return molecular_orbital_data

    def _extract_natural_electron_configuration(self, start_index: int):

        natural_electron_configuration = []

        for i in range(start_index, start_index + self.n_atoms, 1):

            # single atom electron configuration ([s, p, d, f])
            electron_configuration = [0.0, 0.0, 0.0, 0.0]

            # split line at any white space
            line_split = self.lines[i].split()
            # remove first two columns of data, rejoin and remove '[core]'
            line_cleaned = ''.join(line_split[2:]).replace('[core]', '')
            # split at '(' and ')' so that orbital type and config can be extracted
            line_cleanedSplit = re.split(r'\(|\)', line_cleaned)

            for j in range(0, len(line_cleanedSplit), 2):

                # add value to appropriate list element
                if 's' in line_cleanedSplit[j]:
                    electron_configuration[0] += float(line_cleanedSplit[j + 1])
                elif 'p' in line_cleanedSplit[j]:
                    electron_configuration[1] += float(line_cleanedSplit[j + 1])
                elif 'd' in line_cleanedSplit[j]:
                    electron_configuration[2] += float(line_cleanedSplit[j + 1])
                elif 'f' in line_cleanedSplit[j]:
                    electron_configuration[3] += float(line_cleanedSplit[j + 1])
                else:
                    continue

            # append to full list
            natural_electron_configuration.append(electron_configuration)

        # also return index i to jump ahead in the file
        return natural_electron_configuration

