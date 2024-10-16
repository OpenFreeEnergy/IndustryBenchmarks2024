import click
import pathlib
import json
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import gufe
import gufe import SmallMoleculeComponent, LigandAtomMapping
import openfe
from openfe import LigandNetwork
from kartograf.atom_mapping_scorer import (
    MappingRMSDScorer, MappingShapeOverlapScorer, MappingVolumeRatioScorer,
)


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


def get_lomap_score(mapping: SmallMoleculeComponent) -> float:
    """
    Calculates the LOMAP score of a LigandAtomMapping.
    """
    score = openfe.setup.lomap_scorers.default_lomap_score(mapping)
    return score


def get_formal_charge(smc: openfe.SmallMoleculeComponent) -> int:
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
    morgan_fp_score = rdkit.DataStructs.TanimotoSimilarity(fpA, fpB)

    return morgan_fp_score


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
) -> dict[str, dict[str, int]]:
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
        name = f'ligand_{node.name}'
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


@click.command
@click.option(
    '--input_ligand_network',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./alchemicalNetwork/ligand_network.graphml"),
    required=True,
    help=("Path to the ligand_network.graphml file that was used to run these "
         "simulations."),
)
@click.option(
    '--output_dir',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path("./"),
    required=True,
    help="Path to the output directory that stores all data.",
)
@click.option(
    '--fixed_ligand_network',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path("./ligand_network.graphml"),
    required=False,
    help=("Only needed when a broken network was fixed with additional edges. "
          "Path to the ligand_network.graphml file that was used to run the "
         "simulations of fixing the network."),
)
def gather_data(
    input_ligand_network: pathlib.Path,
    output_dir: pathlib.Path,
    fixed_ligand_network: str =None,
):
    """
    Function that gathers all the data.
    """
    # Make folder for outputs
    output_dir.mkdir(exist_ok=False, parents=True)

    # GATHER Input based information
    ligand_network = parse_ligand_network(input_ligand_network)
    # Combine old + new LigandNetwork if a fixed network is provided
    if fixed_ligand_network:
        fixed_network = parse_ligand_network(fixed_ligand_network)
        ligand_network = ligand_network.enlarge_graph(
            edges=fixed_network.edges, nodes=fixed_network.nodes)
    transformation_scores = gather_transformation_scores(ligand_network)
    ligand_scores = gather_ligand_scores(ligand_network)
    # Create a single dict of all scores
    scores = {
        "transformation_scores": transformation_scores,
        "ligand_scores": ligand_scores,
    }
    # Save this to json
    file = pathlib.Path(output_dir / 'all_network_properties.json')
    with open(file, mode='w') as f:
        json.dump(scores, f)

    blinded_network = get_transformation_network_map(ligand_network)
    # Save this to json
    file = pathlib.Path(output_dir / 'network_map.json')
    with open(file, mode='w') as f:
        json.dump(blinded_network, f)


if __name__ == "__main__":
    gather_data()
