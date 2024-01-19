# -*- coding: utf-8 -*- 

from rdkit import Chem
from rdkit.Chem import AllChem
import os, sys
import time

#usage: python assign_mol2_ids.py input.mol2 output.mol2

#2024-01-11 12:45
#요약: obabel 이용한 pdb->mol2 변환 이슈
#예시: [NH3+][C@H](CC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12)C(=O)[O-]가 D-form(pdb)에서 L-form(mol2)으로 변환됨
#원인 추정: --gen3d
#해결 방안: rdkit이 아닌 gen3d로 pdb 생성
#해결 방안 예시: smi->rdkit mol->pdb->mol2를 smi->obabel pdb->mol2로
# 해결 방안 실패.. rdkit에서 obabel pdb의 bond order를 인식하지 못함(할 수 있는데 내가 못 찾은 것일 수도..)

#2024-01-11 14:07
#요약: rdkit의 obabel 출신 pdb 인식 이슈
#예시: pdb의 bond가 모두 single bond로 표현됨
#원인 추정: rdkdit으로 pdb 파일을 읽고 쓸때 flavor를 설정하지 않음
#해결 방안: pdb 읽고 쓸 때 flavor 설정
#해결 방안 예시:  Chem.MolToPDBFile(mol, f'{tmp_path}/new_aa.pdb', flavor=1&2&4&16&32)
#                mol = Chem.MolFromPDBFile(f'{tmp_path}/aa.pdb', removeHs=False, flavor=1&2&4&16&32)

def get_atom_name(name):
    return ' ' + name + ' '*(3-len(name))

def obabel_smi_mol2(smi, molname, mol2_dest):
    flavor = 1&2&4&16&32
    tmp_path = 'C:/Users/jungm/Jupyter notebook/common/tmp'

    os.system(f'obabel -:"{smi} {molname}" -O "{tmp_path}/aa.pdb" --gen3d -n') #-n 왜 붙였더라?
    mol = Chem.MolFromPDBFile(f'{tmp_path}/aa.pdb', sanitize=False, removeHs=False, flavor=flavor)
    #print(Chem.MolToSmiles(mol))
    try:
        mol = assign_ids_to_pdbmol(mol)
    except:
        print(f'error with {molname}')
        print(f'error smi {smi}')
        return
    Chem.MolToPDBFile(mol, f'{tmp_path}/new_aa.pdb', flavor=flavor)
    os.system(f'obabel -ipdb "{tmp_path}/new_aa.pdb" -omol2 -O {mol2_dest} --title {molname}')

def rdkit_smi_mol2(smi, molname, mol2_dest):
    flavor = 1&2&4&16&32
    tmp_path = 'C:/Users/jungm/Jupyter notebook/common/tmp'

    os.system(f'obabel -:"{smi} {molname}" -O "{tmp_path}/aa.pdb" --gen3d') #-n 왜 붙였더라?
    mol = Chem.MolFromPDBFile(f'{tmp_path}/aa.pdb', sanitize=False, removeHs=False, flavor=flavor)
    #print(Chem.MolToSmiles(mol))
    try:
        mol = assign_ids_to_pdbmol(mol)
    except:
        print(f'error with {molname}')
        print(f'error smi {smi}')
        return
    Chem.MolToPDBFile(mol, f'{tmp_path}/new_aa.pdb', flavor=flavor)
    os.system(f'obabel -ipdb "{tmp_path}/new_aa.pdb" -omol2 -O {mol2_dest} --title {molname}')

def set_residue_info(atom, elem_idxs):
    elem = ''.join([c for c in atom.GetPDBResidueInfo().GetName() if not c.isdigit() and c != ' ']) #get only char
    if not elem_idxs.get(f'{elem}'): elem_idxs[f'{elem}'] = 1
    elem_idx = elem_idxs[f'{elem}']
    name = f'{elem}{elem_idx}'
    atom.GetPDBResidueInfo().SetName(get_atom_name(name))
    elem_idxs[f'{elem}'] += 1
    return atom, elem_idxs

def assign_ids_to_pdbmol(mol):
    smarts = Chem.MolFromSmarts("[OX1-]-[CX3](=[OX1])-[CX4]-[NX4+]") #assume protonated, single OC(=O)(N)C
    anames = ['OXT', 'C', 'O', 'CA', 'N']
    mat = mol.GetSubstructMatches(smarts)[0]

    for i, idx in enumerate(mat): #label amino acid skeleton atoms
        atom = mol.GetAtomWithIdx(idx)
        atom.GetPDBResidueInfo().SetName(get_atom_name(anames[i]))

    elem_idxs = {}
    #last atom is N
    #atom = mat[-1] #alreay last atom
    #amino acid skeleton amine-H
    amino_idxs = [idx for idx in mat]
    aminoh_idxs = [atom.GetIdx() for atom in atom.GetNeighbors() if atom.GetAtomicNum() == 1]

    for idx in aminoh_idxs:
        atom = mol.GetAtomWithIdx(idx)
        atom, elem_idxs = set_residue_info(atom, elem_idxs)
    
    amino_idxs = set(amino_idxs + aminoh_idxs)
    for atom in mol.GetAtoms():
        if atom.GetIdx() in amino_idxs: continue
        atom, elem_idxs = set_residue_info(atom, elem_idxs)
    
    return mol

obabel_smi_mol2(sys.argv[1], sys.argv[2], sys.argv[3])