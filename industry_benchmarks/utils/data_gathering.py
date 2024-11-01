import click
from collections import defaultdict
import pathlib
import json
import tqdm
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFreeSASA, rdMolDescriptors
from rdkit.Chem.AtomPairs import Pairs, Torsions
import gufe
from gufe import SmallMoleculeComponent, LigandAtomMapping, AtomMapping
from gufe.tokenization import JSON_HANDLER
import openfe
from openfe import LigandNetwork
from kartograf.atom_mapping_scorer import (
    MappingRMSDScorer, MappingShapeOverlapScorer, MappingVolumeRatioScorer,
)
import shutil
from openff.units import unit

# define all the result files we want to collect
RESULT_FILES = [
    # data files
    "structural_analysis_data.npz",
    "energy_replica_state.npz",
    "simulation_real_time_analysis.yaml",
    "info.yaml",
    # png analysis files,
    "forward_reverse_convergence.png",
    "ligand_COM_drift.png",
    "ligand_RMSD.png",
    "mbar_overlap_matrix.png",
    "protein_2D_RMSD.png",
    "replica_exchange_matrix.png",
    "replica_state_timeseries.png"
]
import abc

class AtomMappingScorer(abc.ABC):
    """A generic class for scoring Atom mappings.
    this class can be used for example to build graph algorithm based networks.

    Implementations of this class can require an arbitrary and non-standardised
    number of input arguments to create.

    Implementations of this class provide the :meth:`.get_score` method

    """

    def __call__(self, mapping: AtomMapping) -> float:
        return self.get_score(mapping)

    @abc.abstractmethod
    def get_score(self, mapping: AtomMapping) -> float:
        """ calculate the score for an  :class:`.AtomMapping`
            the scoring function returns a value between 0 and 1.
            a value close to 1.0 indicates a small distance, a score close
            to zero indicates a large cost/error.

        Parameters
        ----------
        mapping: AtomMapping
            the mapping to be scored
        args
        kwargs

        Returns
        -------
        float
            a value between [0,1] where one is a very bad score and 0 a very
            good one.

        """
        pass


class MoleculeChargeScorer(AtomMappingScorer):

    def __init__(self, charge_class_change_punishment_factor: float = 0.1,
                 smoothing_factor: float = 3.5):
        self.charge_class_change_punishment_factor = charge_class_change_punishment_factor
        self.smoothing_factor = smoothing_factor

    @staticmethod
    def charge_metric(chargeA: float, chargeB: float,
                      charge_class_change_punishment_factor: float = 0.1,
                      smoothing_factor: float = 3.5) -> float:
        def get_charge_class(charge: float) -> int:
            if charge > 0:
                return 1
            elif charge == 0:
                return 2
            else:
                return 3

        if get_charge_class(chargeA) == get_charge_class(chargeB):
            dist = min(1, abs((chargeA - chargeB) / smoothing_factor))
        elif abs(get_charge_class(chargeA) - get_charge_class(chargeB)) == 1:
            dist = min(1, charge_class_change_punishment_factor + abs(
                (chargeA - chargeB) / smoothing_factor))
        else:
            dist = min(1, 2 * charge_class_change_punishment_factor + abs(
                (chargeA - chargeB) / smoothing_factor))

        return 1 - dist

    def get_score(self, atom_mapping: AtomMapping) -> float:
        chargeA = Chem.GetFormalCharge(atom_mapping.componentA.to_rdkit())
        chargeB = Chem.GetFormalCharge(atom_mapping.componentB.to_rdkit())
        return self.charge_metric(chargeA=chargeA, chargeB=chargeB,
                                  smoothing_factor=self.smoothing_factor,
                                  charge_class_change_punishment_factor=self.charge_class_change_punishment_factor)


def parse_ligand_network(
    input_ligand_network_path: str,
) -> LigandNetwork:
    """
    Read a LigandNetwork
    """
    with open(input_ligand_network_path) as f:
        graphml = f.read()

    network = gufe.LigandNetwork.from_graphml(graphml)
    return network


def get_transformation_network_map(
    input_ligand_network: LigandNetwork,
) -> dict:
    """
    Blinded transformation network: Names of ligands and how they are connected.
    """
    blinded_network = {}
    blinded_network["nodes"] = [node.name for node in input_ligand_network.nodes]
    blinded_network["edge"] = [(edge.componentA.name, edge.componentB.name) for edge in input_ligand_network.edges]

    return blinded_network


def get_number_rotatable_bonds(smc: SmallMoleculeComponent) -> int:
    """
    Calculates the number of rotatable bonds in a molecule.
    """
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    num_rotatable_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(m, strict=True)
    return num_rotatable_bonds

def get_number_ring_systems(smc: SmallMoleculeComponent) -> int:
    """
    Calculates the number of ring systems in a molecule.
    """
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    num_rings = Chem.rdMolDescriptors.CalcNumRings(m)
    return num_rings

def get_number_heavy_atoms(smc: SmallMoleculeComponent) -> int:
    """
    Calculates the number of heavy atoms in a molecule.
    """
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    num_heavy_atoms = Chem.rdMolDescriptors.CalcNumHeavyAtoms(m)
    return num_heavy_atoms

def get_system_element_count(smc: SmallMoleculeComponent) -> int:
    """
    Calculates the number of unique elements in a molecule.
    """
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    atomic_numbers = [atom.GetAtomicNum() for atom in m.GetAtoms()]
    return len(set(atomic_numbers))

def get_solvent_accessible_surface_area(smc: SmallMoleculeComponent) -> float:
    """
    Calculates the solvent accessible surface area of a molecule.
    """
    m = smc.to_rdkit()
    Chem.SanitizeMol(m)
    radii = Chem.rdFreeSASA.classifyAtoms(m)
    sasa = Chem.rdFreeSASA.CalcSASA(m, radii)
    return sasa


def get_lomap_score(mapping: LigandAtomMapping) -> float:
    """
    Calculates the LOMAP score of a LigandAtomMapping.
    """
    score = openfe.setup.lomap_scorers.default_lomap_score(mapping)
    return score


def get_formal_charge(smc: SmallMoleculeComponent) -> int:
    """
    Return the formal charge of a molecule
    """
    return Chem.rdmolops.GetFormalCharge(smc.to_rdkit())


def get_alchemical_charge_difference(mapping) -> int:
    """
    Checks and returns the difference in formal charge between state A and B.

    Parameters
    ----------
    mapping : dict[str, ComponentMapping]
      Dictionary of mappings between transforming components.

    Returns
    -------
    int
      The formal charge difference between states A and B.
      This is defined as sum(charge state A) - sum(charge state B)
    """
    chg_A = get_formal_charge(mapping.componentA)
    chg_B = get_formal_charge(mapping.componentB)

    return chg_A - chg_B


def get_charge_score(mapping: LigandAtomMapping) -> float:
    """
    This function calculates the Charge score of a ligand pair.
    """
    charge_scorer = MoleculeChargeScorer()
    score = charge_scorer(mapping=mapping)
    return score


def get_shape_score(mapping: LigandAtomMapping) -> float:
    """
    This function calculates the ShapeTanimoto distance (as implemented in RDKIt)
    for a given ligand pair.
    """
    shape_scorer = MappingShapeOverlapScorer()
    score = shape_scorer(mapping=mapping)
    return score


def get_volume_score(mapping: LigandAtomMapping) -> float:
    """
    This function calculates a Volume ratio based score
    returns a normalized value between 0 and 1, where 0 is the best
    and 1 ist the worst score.
    """
    volume_scorer = MappingVolumeRatioScorer()
    score = volume_scorer(mapping=mapping)
    return score


def get_mapping_RMSD_score(mapping: LigandAtomMapping) -> float:
    """
    This function calculates a mapping RMSD based score. The RMSD between
    mapped atoms of the edge is calculated.
    The score is a normalized value between 0 and 1, where 1.0 is the best
    and 0.0 is the worst score.
    """
    rmsd_scorer = MappingRMSDScorer()
    score = rmsd_scorer(mapping=mapping)
    return score


def get_number_heavy_dummy_heavy_core_atoms(mapping: LigandAtomMapping) -> int:
    """
    This function calculates the size of the mapped and unmapped regions of
    an edge.

    Returns
    -------
    size_core: int
      The number of heavy atoms in the mapped (core) region of the hybrid topology.
    size_dummy_A: int
      The number of heavy atoms in the unmapped (dummy) region of compound A
    size_dummy_B: int
      The number of heavy atoms in the unmapped (dummy) region of compound B
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    mapped_atomIDs = mapping.componentA_to_componentB

    heavy_core_A = []
    heavy_dummy_A = []
    heavy_core_B = []
    heavy_dummy_B = []

    for atom in molA.GetAtoms():
        # For heavy core: Heavy atom indices that are in the mapping:
        if atom.GetIdx() in mapped_atomIDs.keys() and atom.GetAtomicNum() != 1:
            # Check if the atom is mapped to a hydrogen atom
            index_mapping = list(mapped_atomIDs.keys()).index(atom.GetIdx())
            index_atomB = list(mapped_atomIDs.values())[index_mapping]
            atomic_number_mapped_atom_B = molB.GetAtoms()[index_atomB].GetAtomicNum()
            # If the mapped atom in ligand B is a hydrogen, these atoms are unmapped.
            # Therefore we will could the heavy atom as a dummy atom, not core
            if atomic_number_mapped_atom_B == 1:
                heavy_dummy_A.append(atom.GetIdx())
            else:
                heavy_core_A.append(atom.GetIdx())
        # For heavy dummy: Heavy atom indices that are not in the mapping:
        if atom.GetIdx() not in mapped_atomIDs.keys() and atom.GetAtomicNum() != 1:
            heavy_dummy_A.append(atom.GetIdx())

    for atom in molB.GetAtoms():
        # For heavy core: Heavy atom indices that are in the mapping:
        if atom.GetIdx() in mapped_atomIDs.values() and atom.GetAtomicNum() != 1:
            # Check if the atom is mapped to a hydrogen atom
            index_mapping = list(mapped_atomIDs.values()).index(atom.GetIdx())
            index_atomA = list(mapped_atomIDs.keys())[index_mapping]
            atomic_number_mapped_atom_A = molA.GetAtoms()[index_atomA].GetAtomicNum()
            # If the mapped atom in ligand B is a hydrogen, these atoms are unmapped.
            # Therefore we will could the heavy atom as a dummy atom, not core
            if atomic_number_mapped_atom_A == 1:
                heavy_dummy_B.append(atom.GetIdx())
            else:
                heavy_core_B.append(atom.GetIdx())
        # For heavy dummy: Heavy atom indices that are not in the mapping:
        if atom.GetIdx() not in mapped_atomIDs.values() and atom.GetAtomicNum() != 1:
            heavy_dummy_B.append(atom.GetIdx())

    assert len(heavy_core_A) == len(heavy_core_B)

    size_core = len(heavy_core_A)
    size_dummy_A = len(heavy_dummy_A)
    size_dummy_B = len(heavy_dummy_B)

    return size_core, size_dummy_A, size_dummy_B


def get_fingerprint_similarity_score(mapping: LigandAtomMapping) -> float:
    """
    Morgan fingerprint Tanimoto similarity scores.
    """
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    Chem.SanitizeMol(molA)
    Chem.SanitizeMol(molB)

    # Morgan Fingerprint
    fpgen_morgan = AllChem.GetMorganGenerator(radius=3)
    # Get Fingerprint as bit vector
    fpA = fpgen_morgan.GetFingerprint(molA)
    fpB = fpgen_morgan.GetFingerprint(molB)
    # Get Tanimoto similarity
    morgan_fp_score = DataStructs.TanimotoSimilarity(fpA, fpB)

    return morgan_fp_score

def get_atom_pair_similarity(mapping: LigandAtomMapping) -> float:
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    Chem.SanitizeMol(molA)
    Chem.SanitizeMol(molB)

    fp_a = Pairs.GetAtomPairFingerprint(molA)
    fp_b = Pairs.GetAtomPairFingerprint(molB)
    return DataStructs.DiceSimilarity(fp_a, fp_b)

def get_topological_torsion_similarity(mapping: LigandAtomMapping) -> float:
    molA = mapping.componentA.to_rdkit()
    molB = mapping.componentB.to_rdkit()
    Chem.SanitizeMol(molA)
    Chem.SanitizeMol(molB)

    fp_a = Torsions.GetTopologicalTorsionFingerprintAsIntVect(molA)
    fp_b = Torsions.GetTopologicalTorsionFingerprintAsIntVect(molB)
    return DataStructs.DiceSimilarity(fp_a, fp_b)



def get_changing_number_rotatable_bonds(mapping: LigandAtomMapping) -> int:
    """
    Number of rotatable bonds in the non-core region
    """
    num_rot_bonds_A = get_number_rotatable_bonds(mapping.componentA)
    num_rot_bonds_B = get_number_rotatable_bonds(mapping.componentB)

    return abs(num_rot_bonds_A - num_rot_bonds_B)


def get_changing_number_rings(mapping: LigandAtomMapping) -> int:
    """
    Number of ring systems in the non-core region
    """
    num_rings_A = get_number_ring_systems(mapping.componentA)
    num_rings_B = get_number_ring_systems(mapping.componentB)

    return abs(num_rings_A - num_rings_B)


def get_difference_solvent_accessible_surface_area(mapping: LigandAtomMapping) -> float:
    """
    Difference in solvent accessible surface area between two ligands
    """
    sasa_A = get_solvent_accessible_surface_area(mapping.componentA)
    sasa_B = get_solvent_accessible_surface_area(mapping.componentB)

    return abs(sasa_A - sasa_B)


def gather_transformation_scores(
    input_ligand_network: LigandNetwork,
) -> dict[str, dict[str, any]]:
    """
    Gather all scores for the transformations in a dict.
    Transformations information
        Lomap score
        formal charge score/number of charge changes
        Shape score
        Volume score
        RMSD (mapping RMSD score)
        Number of heavy dummy and core atoms
        2/3D fingerprint similarity scores (e.g. Tanimoto)
        changing number of rotatable bonds
        changing number of rings
    Returns
    -------
    dict[str, dict[str, int]]
      Dictionary of the edge name and a dictionary of score name and value.
    """
    transformations_scores = {}
    for edge in input_ligand_network.edges:

        name = f'edge_{edge.componentA.name}_{edge.componentB.name}'
        transformations_scores[name] = {}
        edge_scores = {}
        edge_scores["ligand_A"] = edge.componentA.name
        edge_scores["ligand_B"] = edge.componentB.name
        lomap_score = get_lomap_score(edge)
        edge_scores["lomap_score"] = lomap_score
        alchemical_charge_difference = get_alchemical_charge_difference(edge)
        edge_scores["alchemical_charge_difference"] = alchemical_charge_difference
        charge_score = get_charge_score(edge)
        edge_scores["charge_score"] = charge_score
        shape_score = get_shape_score(edge)
        edge_scores["shape_score"] = shape_score
        volume_score = get_volume_score(edge)
        edge_scores["volume_score"] = volume_score
        mapping_rmsd_score = get_mapping_RMSD_score(edge)
        edge_scores["mapping_rmsd_score"] = mapping_rmsd_score
        size_core, size_dummy_A, size_dummy_B = get_number_heavy_dummy_heavy_core_atoms(edge)
        edge_scores["num_heavy_core"] = size_core
        edge_scores["num_heavy_dummy_A"] = size_dummy_A
        edge_scores["num_heavy_dummy_B"] = size_dummy_B
        diff_rings_AB = get_changing_number_rings(edge)
        diff_rot_bonds_AB = get_changing_number_rotatable_bonds(edge)
        edge_scores["difference_num_rings_AB"] = diff_rings_AB
        edge_scores["difference_num_rot_bonds_AB"] = diff_rot_bonds_AB
        morgan_fp_score = get_fingerprint_similarity_score(edge)
        edge_scores["morgan_tanimoto_similarity"] = morgan_fp_score
        diff_sasa_AB = get_difference_solvent_accessible_surface_area(edge)
        edge_scores["difference_solvent_accessible_surface_area"] = diff_sasa_AB
        diff_atom_pair = get_atom_pair_similarity(edge)
        edge_scores["atom_pair_tanimoto_similarity"] = diff_atom_pair
        diff_topological_torsions = get_topological_torsion_similarity(edge)
        edge_scores["topological_torsion_tanimoto_similarity"] = diff_topological_torsions

        transformations_scores[name] = edge_scores

    return transformations_scores


def gather_ligand_scores(
    input_ligand_network: LigandNetwork,
) -> dict[str, dict[str, int]]:
    """
    Gather all scores for the ligands in a dict.
    Ligand information
        Number of rotatable bonds
        Number of rings
        Number of heavy atoms
        system element counts
        solvent accessible surface area
        formal charge
    Returns
    -------
    dict[str, dict[str, int]]
      Dictionary of the ligand name and a dictionary of score name and value.
    """
    all_ligand_scores = {}
    for node in input_ligand_network.nodes:
        name = f'{node.name}'
        all_ligand_scores[name] = {}
        ligand_scores = {}
        num_rotatable_bonds = get_number_rotatable_bonds(node)
        ligand_scores["num_rotatable_bonds"] = num_rotatable_bonds
        num_rings = get_number_ring_systems(node)
        ligand_scores["num_rings"] = num_rings
        num_heavy_atoms = get_number_heavy_atoms(node)
        ligand_scores["num_heavy_atoms"] = num_heavy_atoms
        num_elements = get_system_element_count(node)
        ligand_scores["num_elements"] = num_elements
        sasa = get_solvent_accessible_surface_area(node)
        ligand_scores["solvent_accessible_surface_area"] = sasa
        formal_charge = get_formal_charge(node)
        ligand_scores["formal_charge"] = formal_charge

        all_ligand_scores[name] = ligand_scores

    return all_ligand_scores

def load_results_file(file_name: pathlib.Path) -> None | dict:
    """Try and load the JSON file as a results file and check that the cleanup script has been used.

    Raises
    ------
        ValueError: If the results clean up script has not been run.
    """
    with open(file_name, "r") as f:
        results = json.load(f)

    # First we check if someone passed in a network_setup.json
    if results.get("__qualname__") == "AlchemicalNetwork":
        print(f"{file_name} is a AlchemicalNetwork json, skipping")
        return None

    #  Now  we check if someome passed in an input json
    if results.get("__qualname__") == "Transformation":
        print(f"{file_name} is an input Transformation json, skipping")
        return None

    # Check to see if we have already cleaned  up this result
    result_key = next(k for k in results["protocol_result"]["data"].keys())
    if (
            "structural_analysis"
            in results["protocol_result"]["data"][result_key][0]["outputs"]
    ):
        raise ValueError(f"{file_name} has not been cleaned up please make sure you run the `results_cleanup.py` script first.")

    # Check to make sure we don't have more than one proto result
    # We might have ProtocolUnitResult-* and ProtocolUnitFailure-*
    # We only handle the case where we have one ProtocolUnitResult
    protocol_unit_result_count = len(
        [
            k
            for k in results["unit_results"].keys()
            if k.startswith("ProtocolUnitResult")
        ]
    )

    # Check to make sure we don't just have failures
    # if all failures, tell user to re-run
    if protocol_unit_result_count == 0:
        print(f"{file_name} failed to run, traceback and exception below, this will not be included in the results. \n")
        proto_failures = [
            k
            for k in results["unit_results"].keys()
            if k.startswith("ProtocolUnitFailure")
        ]
        for proto_failure in proto_failures:
            print("\n")
            print(results["unit_results"][proto_failure]["traceback"])
            print(results["unit_results"][proto_failure]["exception"])
            print("\n")
        return None

    return results

def remove_first_reversed_sequential_duplicate_from_path(path: pathlib.Path) -> pathlib.Path:
    """
    Remove the first duplicated directory from path
    We reverse the path so we remove the first sequential
    duplicated directory starting from the deepest part of the path.

    This is taken from the clean up script.
    """

    # reverse the path parts so we start from the deepest part first
    reversed_path_parts = list(reversed(path.parts))
    max_idx = len(reversed_path_parts)

    # find index of first dupe
    for idx in range(max_idx - 1):
        if reversed_path_parts[idx] == reversed_path_parts[idx + 1]:
            break
    # no dupes
    else:
        print("Path didn't have any dupes")
        return path

    del reversed_path_parts[idx]
    return pathlib.Path(*reversed(reversed_path_parts))


def find_data_folder(result: dict) -> None | pathlib.Path:
    """
    Find the path to the cleaned up results file for this transformation result.
    """
    # get the name of the key which is a gufe token
    # for the only ProtocolUnitResult-* in unit_results
    # this means we can grab the first that matches since there is only
    # one ProtocolUnitResult-*
    proto_key = next(
        k
        for k in result["unit_results"].keys()
        if k.startswith("ProtocolUnitResult")
    )
    results_dir = (
        pathlib.Path(result["unit_results"][proto_key]["outputs"]["nc"]["path"])
        .resolve()
        .parent
    )
    # if the dir doesn't exist, we should try and fix it
    if not results_dir.is_dir():
        print("Fixing path to results dir")
        # Depending on the relative location to the result dir, we might have
        # to fix a duplicate folder, see this post for more details
        # https://github.com/OpenFreeEnergy/IndustryBenchmarks2024/pull/83#discussion_r1689003616
        results_dir = remove_first_reversed_sequential_duplicate_from_path(
            results_dir
        )
        # Now we should check if the dir exists
        if not results_dir.is_dir():
            error_message = f"Can't find results directory: {results_dir}"
            raise FileNotFoundError(error_message)

    # now check that all of the results files can be found in the folder
    # allow skipping of missing png files
    for f_name in RESULT_FILES:
        if not results_dir.joinpath(f_name).exists():
            if ".png" not in f_name:
                error_message = f"Can't find cleaned results file: {f_name} in {results_dir}"
                raise FileNotFoundError(error_message)
            else:
                error_message = f"Can't find cleaned results file: {f_name} in {results_dir} skipping"
                print(error_message)


    return results_dir

def get_transform_name(result: dict, alchemical_network: gufe.AlchemicalNetwork) -> tuple[str, str, str]:
    """
    Get the name of this transformation, taking into account that the inputs might have accidentally been deleted.

    Returns
    -------
    The name of the transformation as a tuple of (phase, ligand_a name, ligand_b name)
    """
    # grab the gufe key of the chemical systems used in the inputs
    unit_result = list(result["unit_results"].values())[0]
    state_a_key = unit_result["inputs"]["stateA"][":gufe-key:"]
    mapping_key = unit_result["inputs"]["ligandmapping"][":gufe-key:"]
    # work out which system this is in the alchemical network
    system_look_up = dict((str(node.key), node) for node in alchemical_network.nodes)
    mapping_look_up = dict((str(edge.mapping.key), edge.mapping) for edge in alchemical_network.edges)
    # build the transform
    if any([isinstance(comp, gufe.ProteinComponent) for comp in system_look_up[state_a_key].components.values()]):
        phase = "complex"
    else:
        phase = "solvent"

    ligmap = mapping_look_up[mapping_key]
    return phase, ligmap.componentA.name, ligmap.componentB.name

def check_network_is_connected(results_data: dict[tuple[str, str, str], list[tuple[unit.Quantity, unit.Quantity, pathlib.Path]]], alchemical_network: gufe.AlchemicalNetwork, name_mapping: None | dict[str, str] = None) -> bool:
    """
    Build a network from the results and check the network is connected.

    This requires that both the complex and solvent phases have 3 repeats.

    """
    from networkx.exception import NetworkXPointlessConcept

    # reverse the name mapping
    if name_mapping is not None:
        reverse_name_mapping = dict(
            (v, k) for k, v in name_mapping.items()
        )
    else:
        reverse_name_mapping = None

    # create a new version of the results data with the correct edge names
    results_by_edge = {}
    edge_to_ligands = {}
    for (phase, ligand_a, ligand_b), results in results_data.items():
        # create the edge string
        if reverse_name_mapping is not None:
            edge_name = f"{phase}_{reverse_name_mapping[ligand_a]}_{reverse_name_mapping[ligand_b]}"

        else:
            edge_name = f"{phase}_{ligand_a}_{ligand_b}"

        results_by_edge[edge_name] = results
        edge_to_ligands[edge_name] = (ligand_a, ligand_b)


    edges = defaultdict(list)
    # group the transforms by ligands
    for transform in alchemical_network.edges:
        if transform.name in results_by_edge and len(results_by_edge[transform.name]) == 3:
            # get the names of the ligands in this edge
            ligand_a, ligand_b = edge_to_ligands[transform.name]
            edges[(ligand_a, ligand_b)].append(transform)

    # extract edges which have both phases completed
    complete_edges = [t for values in edges.values() if len(values) == 2 for t in values]

    # extract the ligand network and check its connected
    result_network = gufe.AlchemicalNetwork(edges=complete_edges)
    ligand_network = extract_ligand_network(result_network)
    try:
        is_connected = ligand_network.is_connected()
    # handle the case where the graph is empty
    except NetworkXPointlessConcept:
        is_connected = False
    return is_connected

def get_estimate(result: dict) -> tuple[unit.Quantity, unit.Quantity]:
    """Extract the DDG and MBAR error estimate from this run"""
    ddg = result["estimate"]["magnitude"] * getattr(unit, result["estimate"]["unit"])
    # get the unit key
    proto_key = next(
        k
        for k in result["unit_results"].keys()
        if k.startswith("ProtocolUnitResult")
    )
    # extract the MBAR uncertainty
    uncertainty = result["unit_results"][proto_key]["outputs"]["unit_estimate_error"]["magnitude"] * getattr(unit, result["unit_results"][proto_key]["outputs"]["unit_estimate_error"]["unit"])
    return ddg, uncertainty


def process_results(results_folders: list[pathlib.Path], output_dir: pathlib.Path, alchemical_network: gufe.AlchemicalNetwork, name_mapping: None | dict[str, str] = None) -> dict[tuple[str, str, str, str], tuple[unit.Quantity, unit.Quantity]]:
    """
    Loop over the results folders extracting the required information and moving it to the output folder.


    Returns
    -------
        results: dict[str, tuple[unit.Quantity, unit.Quantity]]
        The extracted DDG and uncertainty for each transformation repeat
        > {("solvent", "ligand_a", "ligand_b", "repeat_0"): (ddg, uncertainty) ...

    """
    # workout the expected number of results
    # assuming 3 repeats of each solvent and complex transformation
    all_results = defaultdict(list)
    expected_results = len(alchemical_network.edges) * 3

    # map the transformation to the results files
    for results_folder in results_folders:
        for results_file in results_folder.glob("**/*.json"):
            # run checks on the end results
            result = load_results_file(file_name=results_file)
            if result is not None:
                # get the estimate which will be saved
                ddg, uncertainty = get_estimate(result)
                # find the cleaned up results file
                simulation_data_file = find_data_folder(result=result)
                # work out the name of the transform
                # we use the tuple to avoid splitting on _ as ligands might have _ in the name
                try:
                    transformation_name = get_transform_name(result=result, alchemical_network=alchemical_network)
                    if name_mapping is not None:
                        phase, ligand_a, ligand_b = transformation_name
                        transformation_name = (phase, name_mapping[ligand_a], name_mapping[ligand_b])
                    # collect the paths to the results files and the structural data
                    if simulation_data_file is not None:
                        all_results[transformation_name].append((ddg, uncertainty, simulation_data_file))
                except KeyError:
                    # if we can not find the edge name this means the result was not expected
                    raise ValueError(f"Result {results_file} contains a transformation which was not expected for the "
                               f"alchemical network, if you have had to fix this network pass in the extra "
                               f"`alchemical_network.json` using the `--fixed_alchemical_network` flag.")

    # Write stats on the number of transformations found
    found_results = sum([len(v) for v in all_results.values()])
    print(f"Total results found {found_results}/{expected_results} indicating {expected_results - found_results} failed transformations.")

    # check we have a connected network
    if not check_network_is_connected(results_data=all_results, alchemical_network=alchemical_network, name_mapping=name_mapping):
        raise ValueError("The network built from the complete results is disconnected, some simulations may still be"
                         "running, needed restarting or you may have missing repeats. Reproducible edge failures may require extra edges which "
                         "can be generated using the `fix_networks.py` script.")

    estimates = {}
    # move the results to the output folder and collect the estimates
    for transformation_name, results in tqdm.tqdm(all_results.items(), desc="Collecting edges", total=len(all_results), ncols=80):
        for i, (ddg, uncertainty, results_dir) in enumerate(results):
            # store the per repeat estimates
            repeat = f"repeat_{i}"
            estimates[transformation_name + (repeat, )] = (ddg, uncertainty)

            # copy the files for this estimate
            phase, ligand_a, ligand_b = transformation_name
            output_path = output_dir.joinpath(f"{phase}_{ligand_a}_{ligand_b}_{repeat}")
            output_path.mkdir(parents=True, exist_ok=False)
           # copy the analysis files
            for f_name in RESULT_FILES:
                target_file = results_dir.joinpath(f_name)
                # we have already done error handling so just try and move files which are present
                if target_file.exists():
                    shutil.copy(target_file, output_path.joinpath(f_name))


    return estimates

def parse_alchemical_network(file_name: pathlib.Path) -> gufe.AlchemicalNetwork:
    j_dict = json.load(
        open(file_name, "r"), cls=gufe.tokenization.JSON_HANDLER.decoder
    )
    alchem_network = gufe.AlchemicalNetwork.from_dict(j_dict)
    return alchem_network

def extract_ligand_network(alchemical_network: gufe.AlchemicalNetwork) -> LigandNetwork:
    """Extract a ligand network from an alchemical_network"""
    edges = []
    for e in alchemical_network.edges:
        edges.append(e.mapping)

    network = LigandNetwork(edges=set(edges))
    return network

def replace_ligand_names(ligand_network: LigandNetwork) -> tuple[LigandNetwork, dict[str, str]]:
    """
    Replace the names of the ligands in the network with generic names and return a mapping of the current name to
        the generic one of the form ligand1

    Parameters
    ----------
    ligand_network:
        The ligand network with ligands that should have their name replaced.

    Returns
    -------
    A tuple of the new ligand network and a dict mapping the old name to the new ligand name.
    """

    # store a mapping of names
    name_mapping = {}
    # create new small molecules for the network
    node_mapping = {}
    # create new names for each ligand
    for i, node in enumerate(ligand_network.nodes):
        new_name = f"ligand{i}"
        name_mapping[node.name] = new_name
        node_mapping[node.name] = node.copy_with_replacements(name=new_name)
    edges = []
    for edge in ligand_network.edges:
        new_edge = LigandAtomMapping(
            componentA=node_mapping[edge.componentA.name],
            componentB=node_mapping[edge.componentB.name],
            componentA_to_componentB=edge.componentA_to_componentB,
            annotations=edge.annotations
        )
        edges.append(new_edge)

    return LigandNetwork(edges=edges), name_mapping


@click.command
@click.option(
    '--input_alchemical_network',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./alchemicalNetwork/alchemical_network.json"),
    required=True,
    help=("Path to the alchemical_network.json file that was used to run these "
         "simulations."),
)
@click.option(
    '--output_dir',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("./output_data_gathering"),
    required=True,
    help="Path to the output directory that stores all data.",
)
@click.option(
    '--fixed_alchemical_network',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=None,
    required=False,
    help=("Only needed when a broken network was fixed with additional edges. "
          "Path to the alchemical_network.json file that was used to run the "
         "simulations of fixing the network."),
)
@click.option(
    '--hide-ligand-names',
    is_flag=True,
    show_default=True,
    default=False,
    help="If the ligand names should be replaced by generic labels to hide confidential data."
)
@click.option(
    "--results-folder",
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    multiple=True,
    help="The path to the directory which contains transformation results, can be supplied multiple times for each repeat folder."
)
def gather_data(
    input_alchemical_network: pathlib.Path,
    output_dir: pathlib.Path,
    fixed_alchemical_network: None | pathlib.Path,
    results_folder: list[pathlib.Path],
    hide_ligand_names: bool,
):
    """
    Function that gathers all the data.
    """
    # Make folder for outputs
    output_dir.mkdir(exist_ok=False, parents=True)
    # make the subfolder which will become the zip archive
    results_dir = output_dir.joinpath("results_data")
    results_dir.mkdir(exist_ok=False, parents=True)

    # GATHER Input based information
    alchemical_network = parse_alchemical_network(input_alchemical_network)
    ligand_network = extract_ligand_network(alchemical_network)
    # Combine old + new LigandNetwork if a fixed network is provided
    if fixed_alchemical_network is not None:
        fixed_alchemical_network = parse_alchemical_network(fixed_alchemical_network)
        fixed_network = extract_ligand_network(fixed_alchemical_network)
        ligand_network = ligand_network.enlarge_graph(
            edges=fixed_network.edges, nodes=fixed_network.nodes)
        # combine the alchemical networks
        alchemical_network = gufe.AlchemicalNetwork(edges=[*alchemical_network.edges, *fixed_alchemical_network.edges])

    # work out if we should hide the names of the ligands
    if hide_ligand_names:
        ligand_network, name_mapping = replace_ligand_names(ligand_network)
        with open(output_dir / "ligand_name_mapping_PRIVATE.json", "w") as out:
            json.dump(name_mapping, out)
    else:
        # make some noise that we are not hiding the names
        name_mapping = None
        import warnings
        warnings.warn(message="The names of the ligands will be used as-is in the results, if you want to make them anonymous"
                              "run the script again with the '--hide-ligand-names' flag.")


    # get the ligand and edge scores
    transformation_scores = gather_transformation_scores(ligand_network)
    ligand_scores = gather_ligand_scores(ligand_network)

    # process the results one by one
    collected_results = process_results(
        results_folders=results_folder,
        output_dir=results_dir,
        alchemical_network=alchemical_network,
        name_mapping=name_mapping
    )

    # create a copy of the results using a string as the hash to enable saving to json
    formatted_results = dict(
        (f"{phase}_{ligand_a}_{ligand_b}_{repeat}", value)
        for (phase, ligand_a, ligand_b, repeat), value in collected_results.items()
    )

    # workout which edges must have failed by checking for 6 unique transformation results
    edge_counts = defaultdict(list)
    for name, results in collected_results.items():
        _, lig_a, lig_b, _ = name
        # use the name format which matches the gather_transformation_scores function
        edge_name = f"edge_{lig_a}_{lig_b}"
        if results:
            edge_counts[edge_name].append(results)

    # annotate the failed edges only
    for edge, results in edge_counts.items():
        if len(results) != 6:
            transformation_scores[edge]["failed"] = True

    blinded_network = get_transformation_network_map(ligand_network)

    # Create a single dict of all scores
    network_properties = {
        "Ligands": ligand_scores,
        "Edges": transformation_scores,
        "DDG_estimates": formatted_results,
    }

    # Save this to json
    file = pathlib.Path(results_dir / 'all_network_properties.json')
    with open(file, mode='w') as f:
        json.dump(network_properties, f, cls=JSON_HANDLER.encoder, indent=2)

    # finally zip the folder
    shutil.make_archive(results_dir.as_posix(), "zip", results_dir.as_posix())



if __name__ == "__main__":
    gather_data()
