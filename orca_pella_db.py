import concurrent.futures
import os
import re
import logging
import uuid
import numpy as np
from tqdm import tqdm
from typing import List
from multiprocessing import cpu_count, Pool
from rdkit import Chem, RDLogger
from rdkit.Chem import Lipinski, rdDetermineBonds
from cclib.parser.orcaparser import ORCA


RDLogger.DisableLog("rdApp.*")


logging.basicConfig(
    format="%(asctime)s - %(levelname)s: %(message)s", filename="error.log"
)


class OrcaPellaReader:
    """
    An adapter for turning Orca data into standardised SQL and JSON formats,
    developed by Alejandro Santana Bonilla and Professor Chris Lorenz of
    King's College London.

    """
    # the fractional digit to which every value will be rounded to
    ROUND_TO = 3
    # the fractional digit to which every 3d coordinate will be rounded to
    ROUND_ATOM_COORDINATES_TO = 6

    def __init__(self, log_path: str):
        self._filename = log_path.split("/")[-2]
        try:
            self._orca_data = ORCA(log_path).parse()
        except (ValueError, IndexError):
            logging.error(f"Orca Reader couldn't read: {log_path}")
            raise AttributeError()
        self._lines_data = [line.strip() for line in open(log_path, "r")]
        self._dipole_moment_index = self._get_index_of_keyword("DIPOLE MOMENT")
        self._quad_polarization_index = self._get_index_of_keyword("QUADRUPOLE MOMENT (A.U.)")

    def _get_index_of_keyword(self, keyword):
        """
        Function to get the index of the first instance of a keyword
        within log file.

        Parameters
        ----------
        keyword: str
            A string that is searched for in the string

        Returns
        -------
        idx: int
            The first index of the keyword

        """
        for idx, line in enumerate(self._lines_data):
            if keyword == line:
                return idx

    @staticmethod
    def _get_atom_coords(arr):
        """
        Function to get the 3d coordinates of an atom

        Parameters
        ----------
        arr: list[list]
            A list of atoms

        Returns
        -------
        atom_coords: dict
            A dict containing coordinate info about atoms

        """
        atom_coords = {
            "element": arr[0],
            **{
                k: round(v, OrcaPellaReader.ROUND_ATOM_COORDINATES_TO)
                for k, v in zip(["x", "y", "z"], arr[1:])
            },
        }
        return atom_coords

    def _get_3d_atom_coordinates(self):
        """
        Function to retrieve the coordinates and elements from an ORCA output file.

        Returns
        -------
        atom_coords: list
            A list with the elements and coordinates from a cclib object

        """
        try:
            start_indexes = self._find_indexes_of_string(
                "CARTESIAN COORDINATES (ANGSTROEM)"
            )[0]
        except IndexError:
            raise AttributeError()

        end_is_found, counter = False, 2
        valid_lines = []
        while not end_is_found:
            line = self._lines_data[start_indexes + counter]
            valid_lines.append(line.split())
            if len(line) == 0:
                end_is_found = True
            counter += 1
        valid_lines = [
            [
                element if ind == 0 else float(element)
                for ind, element in enumerate(_list)
            ]
            for _list in valid_lines[:-1]
        ]
        return valid_lines

    @staticmethod
    def _make_xyz_file(atom_coords):
        xyz_file_head = [str(len(atom_coords)) + "\n", "\n"]
        xyz_file = [
            " ".join([str(_float) for _float in atom_coord]) + "\n"
            for atom_coord in atom_coords
        ]
        xyz_file[-1] = xyz_file[-1][:-2]
        xyz_file = xyz_file_head + xyz_file
        file_path = "xyz/" + str(uuid.uuid4()) + ".xyz"
        with open(file_path, "w") as file:
            file.writelines(xyz_file)
        return file_path

    def _get_rdkit_smiles(self, file_path):
        raw_mol = Chem.MolFromXYZFile(file_path)
        try:
            mol = Chem.Mol(raw_mol)
        except Exception as e:
            logging.error(f"RDKIT couldn't read molecule with error {e} - {self._filename}")
            raise AttributeError()

        rdDetermineBonds.DetermineConnectivity(mol)

        # try:
        #     rdDetermineBonds.DetermineBondOrders(mol, charge=0)
        # except ValueError:
        #     logging.error(f"RDKIT couldn't determine bonding orders {self._filename}")
        #     raise AttributeError()
        try:
            mol = Chem.RemoveHs(mol)
        except Chem.rdchem.AtomValenceException:
            logging.error(f"RDKIT couldn't remove H from {self._filename}")
            raise AttributeError()
        smi = Chem.MolToSmiles(mol)
        smi = Chem.CanonSmiles(smi)
        return mol, smi

    def _encode_molecule(self, mol):
        """
        Given a SMILES string return a list of molecular encodings

        Parameters
        ----------
        mol: rdkit.Chem.Mol
            An rdkit molecule

        Returns
        -------
        mol_encoding: list
            A list of molecule features

        """
        mol_encoding = {
            "fractioncsp3": Lipinski.FractionCSP3,
            "no_count": Lipinski.NOCount,
            "num_aliphatic_rings": Lipinski.NumAliphaticRings,
            "num_h_donors": Lipinski.NumHDonors,
            "num_h_accetors": Lipinski.NumHAcceptors,
            "num_rotatable_bonds": Lipinski.NumRotatableBonds,
            "num_hetero_atoms": Lipinski.NumHeteroatoms,
            "heavy_atom_count": Lipinski.HeavyAtomCount,
            "ring_count": Lipinski.RingCount,
            "nhoh_count": Lipinski.NHOHCount,
            "num_aliphatic_carbocycles": Lipinski.NumAliphaticCarbocycles,
            "num_aliphatic_hetercycles": Lipinski.NumAliphaticHeterocycles,
            "num_aromatic_carbocycles": Lipinski.NumAromaticCarbocycles,
            "num_aromatic_heterocycles": Lipinski.NumAromaticHeterocycles,
            "num_aromatic_rings": Lipinski.NumAromaticRings,
        }

        for k, v in mol_encoding.items():
            try:
                mol_encoding[k] = v(mol)
            except Exception:
                logging.error(f"RDKIT couldn't calculate {k} for {self._filename}")
                mol_encoding[k] = 0

        return mol_encoding

    def _get_element_string(self):
        """
        Function to retrieve the elements from a cclib coordinates object.

        Returns
        -------
        element_string: str
            A string containing the elements of a cclib coordinate object.

        """
        element_string = "".join(
            atom["element"] for atom in self._get_3d_atom_coordinates()
        )
        return element_string

    def _count_elements(self):
        """
        Function to retrieve the elements' counts from a cclib coordinates object.

        Returns
        -------
        element_counts: dict
            A dict containing the elements of a cclib coordinate object and their count

        """
        element_counts = {
            c: self._get_element_string().count(c) for c in self._get_element_string()
        }
        return element_counts

    def _mo_energies(self):
        """
        Function to retrieve the occupied energy levels up to the
        Highest Occupied Molecular Orbital (HOMO).

        Returns
        -------
        _mo_energies_list: list
            A list with the energy levels for occupied orbitals from
            a cclib object.

        """
        try:
            mo_energies = self._orca_data.moenergies
        except AttributeError:
            logging.error(
                f"Orca Reader did not find moenergies in file id: {self._filename}"
            )
            raise AttributeError()
        _mo_energies_list = [
            round(i, OrcaPellaReader.ROUND_TO)
            for i in mo_energies[0][: self._orca_data.homos[0]]
        ]
        return _mo_energies_list

    def _find_indexes_of_string(self, string):
        """
        Function to search a string in a user provided input file.

        Parameters
        ------------
        string: str
            String value to find the indexes of

        Returns
        -------
        indexes: list
            Returns a list of indexes where the String has been found.

        """
        indexes = [idx for idx, line in enumerate(self._lines_data) if string in line]
        return indexes

    @staticmethod
    def _seek_string(string, array):
        """
        Function to search a string within an array.

        Parameters
        ------------
        string: str
          An user-provided string.
        array: list
          An user-provided list.

        Returns
        ---------
        idx:
          List of values

        """
        for idx, values in enumerate(array):
            if string in values:
                return idx

    @staticmethod
    def _clean_arrays(array, value, filename):
        """
        Function to get the eV values of singlet and triplet excited states

        Parameters
        ----------
        array:
            array of arrays to be cleaned with columns: 'state', 'eV', 'nm',
            'fL', 'fV', 'Rl', 'RV'
        value:
            index at which the data is found

        Returns
        -------
        cleaned_array: np.ndarray
            Array with only the values of 'eV' column

        """
        valid_columns = len(array) - value
        # the columns are 'state', 'eV', 'nm', 'fL', 'fV', 'Rl', 'RV',
        # but only 'eV' is wanted
        try:
            cleaned_array = np.array(
                [float(array[value + i][1]) for i in range(1, valid_columns - 1)]
            )
        except IndexError:
            logging.error(
                f"Orca Reader did not find singlet and triplet energies in file id: {filename}"
            )
            raise AttributeError()
        return cleaned_array

    def _excited_states(self):
        """
        Function to get the singlet and triplet excited states.

        Returns
        -------
        exc_states_dict: dict
            A dict of the singlet and triplet excited states

        """
        # get indexes to query from log
        start_indexes = self._find_indexes_of_string("sTD-DFT procedure...")
        end_indexes = self._find_indexes_of_string("sTD-DFT done")

        if len(start_indexes) != 2 or len(end_indexes) != 2:
            logging.error(
                f"Orca Reader did not find singlet and triplet energies in file id: {self._filename}"
            )
            raise AttributeError()

        triplet_indexes, singlet_indexes = list(zip(start_indexes, end_indexes))
        singlets = [
            re.split("\s+", s.strip())
            for s in self._lines_data[slice(*singlet_indexes)]
        ]
        triplets = [
            re.split("\s+", t.strip())
            for t in self._lines_data[slice(*triplet_indexes)]
        ]

        # transform array elements from str to float
        singlets_array = self._clean_arrays(
            singlets, self._seek_string("state", singlets), self._filename
        )
        triplets_array = self._clean_arrays(
            triplets, self._seek_string("state", triplets), self._filename
        )
        exc_states_dict = {
            "singlet": singlets_array.tolist(),
            "triplet": triplets_array.tolist(),
        }
        return exc_states_dict

    @staticmethod
    def _pol_to_matrix(quad):
        """
        Function to store the polarization strings as a numpy matrix.

        Parameters
        ----------
        quad: List
           A list of the polarizations

        Returns
        -------
        final_matrix: list
            A list with the polarization Tensor retrieved from Orca.

        """
        diag, quad = np.eye(3, k=0), np.array(quad, dtype="float64")
        xx_yy_zz = diag[:3] * quad[:3]
        xy = np.asarray([[0, 1, 0], [0, 0, 0], [0, 0, 0]]) * quad[3]
        xz = np.asarray([[0, 0, 1], [0, 0, 0], [0, 0, 0]]) * quad[4]
        yz = np.asarray([[0, 0, 0], [0, 0, 1], [0, 0, 0]]) * quad[5]
        final_matrix = (xx_yy_zz + xy + xz + yz).tolist()
        return final_matrix

    def _dipole_moments(self):
        """
        Function to search the Dipole moment from an Orca output file.

        Returns
        -------
        dip_dict: dict
            A dict with the values of the Dipole moments.

        """
        if self._dipole_moment_index is None:
            logging.error(
                f"Orca Reader did not find dipole moment in file id: {self._filename}"
            )
            raise AttributeError()

        # directions = self.lines_data[ind + 2]  # seems like this isn't used at all
        elements = ["electron", "nucleus", "total"]
        data_indices = [self._dipole_moment_index + i for i in [3, 4, 6]]
        offsets = [2, 3, 4]
        values = [
            self._lines_data[ind].split()[offset:]
            for ind, offset in zip(data_indices, offsets)
        ]
        dip_dict = {
            elem: [round(float(value), OrcaPellaReader.ROUND_TO) for value in values]
            for elem, values in zip(elements, values)
        }
        return dip_dict

    def _quad_polarizations(self):
        """
        Function to browse the Quadrupole moment from an Orca output file.

        Returns
        -------
        quad_pol_dict: dict
            A dict with the values of the Quadrupole moments.

        """
        if self._quad_polarization_index is None:
            logging.error(
                f"Orca Reader did not find quadrupole moment in file id: {self._filename}"
            )
            raise AttributeError()

        # directions = self.lines_data[ind+3].split()  # seems like this isn't used at all
        elements = ["electron", "nucleus", "total"]
        data_indices = [self._quad_polarization_index + i for i in [5, 4, 6]]
        values = [self._lines_data[i].split()[1:7] for i in data_indices]
        processed_values = [self._pol_to_matrix(val) for val in values]
        quad_pol_dict = {
            elem: [
                [round(float(i), OrcaPellaReader.ROUND_TO) for i in sublist]
                for sublist in matrix
            ]
            for elem, matrix in zip(elements, processed_values)
        }
        return quad_pol_dict

    def _qm_charges_loewdin(self):
        """
        Function to retrieve the lowdin charges from an ORCA output file.

        Returns
        -------
        qm_charges_list: list
            A list storing the lowdin charges from a cclib object.

        """
        try:
            atom_charges = self._orca_data.atomcharges
        except AttributeError:
            logging.error(
                f"Orca Reader did not find atom charges in file id: {self._filename}"
            )
            raise AttributeError()

        qm_charges_list = [
            round(i, OrcaPellaReader.ROUND_TO) for i in atom_charges["lowdin"]
        ]
        return qm_charges_list

    @staticmethod
    def run_with_timeout(code_segment, *args, timeout=1):
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(code_segment, *args)
            try:
                result = future.result(timeout=timeout)
            except concurrent.futures.TimeoutError:
                raise TimeoutError("Code segment took too long to execute.")
        return result

    @staticmethod
    def get_json(log_path: str):
        """
        Function to get JSON schema of a log and smiles file

        Parameters
        -------
        log_path: str
            The path to output.log

        Returns
        -------
        filename, schema:
            name of the file, all the data in a dict

        """
        try:
            op_db = OrcaPellaReader(log_path)
            _3d_atom_coordinates = op_db._get_3d_atom_coordinates()
            file_path = op_db._make_xyz_file(_3d_atom_coordinates)
            mol, smi = op_db._get_rdkit_smiles(file_path)

            schema = {
                "smiles": smi,
                "mol_features": op_db._encode_molecule(mol),
                "3d_coords": op_db._get_3d_atom_coordinates(),
                "dipole_moments": op_db._dipole_moments(),
                "polarizations": op_db._quad_polarizations(),
                "loewdin": op_db._qm_charges_loewdin(),
                "excited_states": op_db._excited_states(),
                "mo_energies": op_db._mo_energies(),
            }
            return log_path, schema

        except AttributeError:
            return

    @staticmethod
    def get_all_jsons_from_folder(
        folder_path: str,
        with_multiprocessing: bool = True,
        track_progress: bool = True,
        save_errors_to_file: bool = False,
    ):
        """
        Function to get the JSONs of all logs in the folder_path.

        Parameters
        ----------
        folder_path: str
            The path to the folder with the logs
        with_multiprocessing: bool
            Toggle whether multiprocessing is used to double performance
        track_progress: bool
            Toggle whether the progress is tracked through the processing
        save_errors_to_file: bool
            Toggle whether the errors are saved to a file

        Returns
        -------
        jsons: list[dict]
            A list containing all JSONs of the logs in the folder_path.

        """
        file_list = [
            (folder_path, file)
            for file in tqdm(os.listdir(folder_path))
            if os.path.isdir(os.path.join(folder_path, file))
        ]

        log_files = [
            os.path.join(parent_folder, sub_folder, "output.log")
            for parent_folder, sub_folder in tqdm(file_list)
        ]

        if with_multiprocessing:
            with Pool(cpu_count()) as p:
                jsons = (
                    list(
                        tqdm(
                            p.imap(OrcaPellaReader.get_json, log_files),
                            total=len(file_list),
                        )
                    )
                    if track_progress
                    else p.map(OrcaPellaReader.get_json, log_files)
                )
        else:
            log_files = tqdm(log_files) if track_progress else log_files
            jsons = [OrcaPellaReader.get_json(file_path) for file_path in log_files]

        if not save_errors_to_file:
            os.remove("error.log")

        return jsons
