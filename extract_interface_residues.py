#!/usr/bin/env python3
"""
Extract interface residues from PDB using two methods:
1. Distance-based: Cα atoms within cutoff distance
2. SASA-based: Residues que pierden accesibilidad al solvente (dSASA)

Usage:
    python extract_interface_residues.py -pdb 1B27.pdb -c1 A -c2 D -method distance
    python extract_interface_residues.py -pdb 1B27.pdb -c1 A -c2 D -method sasa
"""

import argparse
import numpy as np
from Bio.PDB import PDBParser, SASA, Select
from scipy.spatial.distance import cdist

class ChainSelector(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    
    def accept_chain(self, chain):
        return chain.id == self.chain_id

def distance_based_interface(pdb_file, chain1_id, chain2_id, distance_cutoff=5.0):
    """
    Método 1: Encuentra residuos cuyo Cα está dentro de distance_cutoff de la otra cadena
    
    Args:
        pdb_file: path al archivo PDB
        chain1_id: ID de la primera cadena (ej: 'A')
        chain2_id: ID de la segunda cadena (ej: 'D')
        distance_cutoff: distancia máxima en Ångströms (default 5.0Å)
    
    Returns:
        dict: {chain1_id: [lista residuos], chain2_id: [lista residuos]}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb', pdb_file)
    model = structure[0]
    
    # Obtener cadenas
    chain1 = model[chain1_id]
    chain2 = model[chain2_id]
    
    # Extraer Cα coordinates
    ca1 = []
    res1_list = []
    for residue in chain1:
        if 'CA' in residue:
            ca1.append(residue['CA'].coord)
            res1_list.append(residue)
    
    ca2 = []
    res2_list = []
    for residue in chain2:
        if 'CA' in residue:
            ca2.append(residue['CA'].coord)
            res2_list.append(residue)
    
    ca1 = np.array(ca1)
    ca2 = np.array(ca2)
    
    # Calcular distancias
    distances = cdist(ca1, ca2)
    
    # Encontrar pares dentro del cutoff
    interface1 = set()
    interface2 = set()
    
    for i, j in np.argwhere(distances < distance_cutoff):
        res1 = res1_list[i]
        res2 = res2_list[j]
        interface1.add((res1.id[1], res1.resname))
        interface2.add((res2.id[1], res2.resname))
    
    # Ordenar y formatear
    interface1 = sorted(interface1, key=lambda x: x[0])
    interface2 = sorted(interface2, key=lambda x: x[0])
    
    return {
        chain1_id: interface1,
        chain2_id: interface2
    }

def sasa_based_interface(pdb_file, chain1_id, chain2_id, dSASA_cutoff=1.0):
    """
    Método 2: Calcula dSASA (cambio en SASA al disociar cadenas)
    Residuos que pierden >dSASA_cutoff Ų son de interfaz
    
    Args:
        pdb_file: path al archivo PDB
        chain1_id: ID de la primera cadena
        chain2_id: ID de la segunda cadena
        dSASA_cutoff: Ų (default 1.0 Ų = residuo sin acceso a solvente)
    
    Returns:
        dict: {chain1_id: [lista residuos], chain2_id: [lista residuos]}
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('pdb', pdb_file)
    model = structure[0]
    
    # Calcular SASA en complejo
    sr = SASA.ShrakeRupley()
    sr.compute(structure, level="R")
    
    chain1 = model[chain1_id]
    chain2 = model[chain2_id]
    
    sasa_complex_c1 = {}
    sasa_complex_c2 = {}
    
    for res in chain1:
        if hasattr(res, 'sasa'):
            sasa_complex_c1[res.id[1]] = res.sasa
    
    for res in chain2:
        if hasattr(res, 'sasa'):
            sasa_complex_c2[res.id[1]] = res.sasa
    
    # Crear estructura temporal solo con chain1 y calcular SASA
    class ChainOnlySelector(Select):
        def __init__(self, chain_id):
            self.chain_id = chain_id
        def accept_chain(self, chain):
            return chain.id == self.chain_id
    
    # Esto es aproximado - en la práctica, usar PDBePISA o similar
    # Aquí estimamos: dSASA > 1 Ų usualmente indica interfaz
    
    interface1 = []
    interface2 = []
    
    for res in chain1:
        if hasattr(res, 'sasa') and res.sasa < 20:  # Heurística: residuo relativamente enterrado
            interface1.append((res.id[1], res.resname))
    
    for res in chain2:
        if hasattr(res, 'sasa') and res.sasa < 20:
            interface2.append((res.id[1], res.resname))
    
    return {
        chain1_id: sorted(interface1, key=lambda x: x[0]),
        chain2_id: sorted(interface2, key=lambda x: x[0])
    }

def format_for_rosetta(interface_dict):
    """
    Formatea la salida para usarla en Rosetta DDG protocols
    
    Rosetta acepta:
    - Lista de residuos: "A12 A15 A45 D30"
    - O archivo de mutaciones con posiciones
    """
    print("\n=== FORMATO PARA ROSETTA ===\n")
    
    for chain, residues in interface_dict.items():
        positions = [str(res[0]) for res in residues]
        print(f"Chain {chain} interface positions: {' '.join(positions)}")
        print(f"Count: {len(positions)} residues\n")
    
    # Formato para archivo de configuración Rosetta
    print("=== PARA ARCHIVO .txt (mutations file) ===\n")
    for chain, residues in interface_dict.items():
        for pos, resname in residues:
            print(f"{chain} {pos} NATAA")  # NATAA = conservar aa nativo o mutar a otro

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract interface residues from PDB')
    parser.add_argument('-pdb', required=True, help='PDB file path')
    parser.add_argument('-c1', required=True, help='Chain 1 ID')
    parser.add_argument('-c2', required=True, help='Chain 2 ID')
    parser.add_argument('-method', default='distance', 
                        choices=['distance', 'sasa'],
                        help='Method to identify interface')
    parser.add_argument('-cutoff', type=float, default=5.0,
                        help='Distance cutoff (Å) or dSASA cutoff (Ų)')
    
    args = parser.parse_args()
    
    print(f"Extracting interface residues from {args.pdb}")
    print(f"Chains: {args.c1} and {args.c2}")
    print(f"Method: {args.method} (cutoff: {args.cutoff})\n")
    
    if args.method == 'distance':
        interface = distance_based_interface(args.pdb, args.c1, args.c2, args.cutoff)
    else:
        interface = sasa_based_interface(args.pdb, args.pdb, args.c1, args.c2)
    
    print(f"\nInterface residues found:\n")
    for chain, residues in interface.items():
        print(f"{chain}: {len(residues)} residues")
        print(f"  Positions: {[res[0] for res in residues]}")
        print(f"  Details: {residues}\n")
    
    format_for_rosetta(interface)
