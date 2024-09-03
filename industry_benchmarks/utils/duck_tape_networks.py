import json
import click
import pathlib
from glob import glob
import networkx as nx

from gufe.tokenization import JSON_HANDLER
from openfe import ChemicalSystem, LigandAtomMapping, Protocol
from openfe import Transformation, AlchemicalNetwork, LigandNetwork

from openfe.setup import KartografAtomMapper, lomap_scorers
from konnektor.network_planners import MstConcatenator
from konnektor.network_tools import concatenate_networks

try:
    from konnektor.network_planners import MstConcatenator
    from konnektor.network_tools import concatenate_networks
except ModuleNotFoundError:  # This will attempt to install kartograf
    print("Did not find Konnektor! will attempt install!")
    !{sys.executable} -m pip install git+https://github.com/OpenFreeEnergy/konnektor.git
    from konnektor.network_planners import MstConcatenator
    from konnektor.network_tools import concatenate_networks
    

def parse_result_folders(result_files_regex:str )->AlchemicalNetwork:
    nodes = []
    edges = {}
    i=0
    for resp in glob(result_files_regex):
    
        res =  json.load(open(resp, 'rb'), cls=JSON_HANDLER.decoder)

        #TODO: Handle error of dag - continue
        for k in res["protocol_result"]["data"]:
            units =  res["protocol_result"]["data"][k]
            
            name=units[0]["name"].split(" repeat")[0]
            stateA=ChemicalSystem.from_dict(units[0]["inputs"]["stateA"])
            stateB=ChemicalSystem.from_dict(units[0]["inputs"]["stateB"])
            ligmap=LigandAtomMapping.from_dict(units[0]["inputs"]["ligandmapping"])
            try:
                protocol=Protocol.from_dict(units[0]["inputs"]["protocol"])
            except:
                protocol = None
            t = Transformation(stateA=stateA, stateB=stateB, mapping =ligmap, name=name, protocol=protocol)
            
            nodes.extend([stateA, stateB])
            if name in edges:
                edges[name].repeats += 1
            else: 
                t.repeats = 1
                edges[name] = t

    alchemical_network = AlchemicalNetwork(nodes=nodes, edges=edges.values())
    return alchemical_network


def alchemical_network_to_networkx(alchemical_network)->nx.Graph:
    
    edges = alchemical_network.edges
    wedges = []
    for edge in list(edges)[::2]:
        wedges.append([edge.stateA, edge.stateB, edge.mapping["score"]])
    
    g = nx.Graph(nodes=alchemical_network.nodes)
    g.add_weighted_edges_from(ebunch_to_add=wedges)
    return g


def alchemical_network_to_ligand_network(alchemical_network)->LigandNetwork:
    nodes = [n.components["ligand"] for n in alchemical_network.nodes]
    edges = [e.mapping for e in alchemical_network.edges]
    network = LigandNetwork(nodes=nodes, edges=edges)
    return network


def decomposite_disconnected_alchemical_network(alchemical_network:AlchemicalNetwork)->list[AlchemicalNetwork]:
    # build alchemical sub-networks
    g = alchemical_network_to_networkx(alchemical_network)
    all_edges = list(alchemical_network.edges)
    sub_network_nodes = list(nx.connected_components(g))
    
    alchem_sub_networks = []
    for i, snn in enumerate(sub_network_nodes):
        # collect edges
        sub_edges = []
        for e in all_edges:
            if any(n.components["ligand"] in (e.stateA, e.stateB) for n in snn):
                sub_edges.append(e)
                
        alchem_sub_network = AlchemicalNetwork(nodes=snn, edges=sub_edges, name=f"sub_net_{i}")
        alchem_sub_networks.append(alchem_sub_network)
        
    return alchem_sub_networks

def ligand_network_to_networkx(ligand_network)->nx.Graph:
    
    edges = ligand_network.edges
    wedges = []
    for edge in list(edges)[::2]:
        wedges.append([edge.componentA, edge.componentB, edge.annotations["score"]])
    
    g = nx.Graph(nodes=ligand_network.nodes)
    g.add_weighted_edges_from(ebunch_to_add=wedges)
    return g


def decomposite_disconnected_ligand_network(ligand_network:LigandNetwork)->list[LigandNetwork]:
    # build alchemical sub-networks
    g = ligand_network_to_networkx(ligand_network)
    all_edges = list(ligand_network.edges)
    sub_network_nodes = list(nx.connected_components(g))
    
    ligand_sub_networks = []
    for i, snn in enumerate(sub_network_nodes):
        # collect edges
        sub_edges = []
        for e in all_edges:
            if any(n in (e.componentA, e.componentB) for n in snn):
                sub_edges.append(e)
                
        ligand_sub_network = LigandNetwork(nodes=snn, edges=sub_edges)
        ligand_sub_networks.append(ligand_sub_network)
        
    return ligand_sub_networks


def ducktape_networks(ligand_sub_networks:list[LigandNetwork], n_connecting_edges:int = 3):
    
    mapper = KartografAtomMapper()
    scorer = lomap_scorers.default_lomap_score
    concatenator = MstConcatenator(mapper=mapper, scorer=scorer, n_connecting_edges=n_connecting_edges, n_processes=1)
    
    duck_taped_network = concatenate_networks(ligand_sub_networks, concatenator=concatenator)
    return duck_taped_network
    

def do_taping(result_files_regex:str, ouput_alchemical_network_folder: str="./ducktaping_transformations", input_ligand_network:str=None):
    # Parse Alchemical Network
    print("Parsing input files:")
    
    #res_alchemical_network = parse_result_folders(result_files_regex)
    res_ligand_network = alchemical_network_to_ligand_network(res_alchemical_network)
    
    # decomposite disconnected sets
    ligand_sub_networks = decomposite_disconnected_ligand_network(res_ligand_network)
    
    # Add missing nodes to ligand_sub_networks if ligand network provided
    missing_nodes = []
    input_ligand_network = None
    if isinstance(input_ligand_network_path, str) and os.path.isfile(input_ligand_network_path):
        # Open Ligand Network
        with open(input_ligand_network_path, "r") as ligand_network_file:
            graphml_str = "\n".join(ligand_network_file.readlines())
        input_ligand_network = LigandNetwork.from_graphml(graphml_str)
        
        # get network node differences
        alchemical_network_ligands = set([n for n in res_ligand_network.nodes])
        missing_commponents = list(alchemical_network_ligands.difference(input_ligand_network.nodes))
    
        for missing_commponent in missing_commponents:
            ligand_sub_networks.append(LigandNetwokr(nodes=missing_commponent))
    
    # check status
    if len(ligand_sub_networks) == 1:
        raise ValueError("Did not find disconnected components in alchemical network!")
    else:
        print(f"\tresult AlchemicalNetwork nodes: {len(res_ligand_network.nodes)}")
        print(f"\tresult in AlchemicalNetwork edges: {len(res_ligand_network.edges)}")
        print(f"\tresult disconnected networks: {len(ligand_sub_networks)}")
    
        if input_ligand_network is not None:
            print(f"\tinput LigandNetwork  nodes: {len(input_ligand_network.nodes)}")
            print(f"\tinput LigandNetwork  edges: {len(input_ligand_network.edges)}")
            print(f"\tinput missing components: {len(missing_nodes)}")
        else:
            print(f"Did not find input LigandNetwork")
            
    print()
    
    # Tape the networks together
    print("Taping the Networks together:")
    duck_taped_network = ducktape_networks(ligand_sub_networks=ligand_sub_networks, n_connecting_edges = 3)
    print(f"\tGenerated in DuckTapedNetwork nodes: {len(duck_taped_network.nodes)}")
    print(f"\tGenerated in DuckTapedNetwork edges: {len(duck_taped_network.edges)}")
    print()
    
    # write out
    

@click.command
@click.option(
    '--results_regex_path',
    type=click.Path(dir_okay=False, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path('results*/*json'),
    help="regular expression to find all result files of this system. (default: 'results*/*json')",
)
@click.option(
    '--ligand_network_file',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    default=pathlib.Path('./alchemicalNetwork/ligand_network.graphml'),
    help="Path to the used ligand network plan (default: 'alchemicalNetwork/ligand_network.graphml')",
)
@click.option(
    '--output_alchemical_network',
    type=click.Path(dir_okay=True, file_okay=False, path_type=pathlib.Path),
    default=pathlib.Path('./tapesForAlchemicalNetwork'),
    help="Directory name in which to store the additional taping transformation json files. (default: './tapesForAlchemicalNetwork')",
)
def cli_do_taping(results_regex_path:str, ligand_network_file:str, output_alchemical_network:str):
    do_taping(result_files_regex= results_regex_path,
       input_ligand_network= ligand_network_file,
       ouput_alchemical_network_folder= output_alchemical_network,
      )

# main Runs
if __name__ == "__main__":
    cli_do_taping()


    
    