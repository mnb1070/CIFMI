from openbabel import openbabel as ob
from openbabel import pybel as pb

def set_unique_ids(obmol):
    smarts = pb.Smarts('[#8X1-]-[#6X3](=[#8X1])-[#6X4]-[#7H3X4+]')
    skeleton_idxs = smarts.findall(obmol)
    
    if len(skeleton_idxs) != 1:
        if len(skeleton_idxs) > 1:
            print('warning: more than one skeletons')
        else:
            print('error: no skeleton')
            return
    
    elem_dict = {}
    skeletonh_idxs = ['H1', 'H2', 'H3']
    skeleton_idxs = skeleton_idxs[0]
    for obres in ob.OBResidueIter(obmol.OBMol):
        for obatom in ob.OBResidueAtomIter(obres):
            id = obres.GetAtomID(obatom)
            if obatom.GetIdx() in skeleton_idxs or id in skeletonh_idxs: continue
            if not elem_dict.get(id):
                if id == 'H': elem_dict[id] = 4
                else: elem_dict[id] = 1
            obres.SetAtomID(obatom, f'{id}{elem_dict[id]}')
            elem_dict[id] += 1
    
    return obmol

def pb_gen3d(pbmol):
    builder = ob.OBBuilder()
    obmol = pbmol.OBMol
    
    builder.Build(obmol)
    obmol.AddHydrogens(False, False)
    pff = pb._forcefields['mmff94']
    pff.Setup(obmol)
    pff.EnableCutOff(True)
    pff.SetVDWCutOff(10.0)
    pff.SetElectrostaticCutOff(20.0)
    pff.SetUpdateFrequency(10)
    
    pff.ConjugateGradients(100, 1.0e-4)
    pff.FastRotorSearch(False)
    pff.ConjugateGradients(100, 1.0e-6)
    pff.UpdateCoordinates(obmol)

    return pbmol

def routine_gen3d(smi, molname, dest):
    print('current', smi)
    pbmol = pb.readstring("smi", smi)
    pbmol.OBMol.CorrectForPH(7.4)
    pbmol = pb_gen3d(pbmol)
    set_unique_ids(pbmol)
    pbmol.OBMol.SetTitle(molname)
    pbmol.write("mol2", dest, overwrite=True)
    print('done', smi)

    return pbmol