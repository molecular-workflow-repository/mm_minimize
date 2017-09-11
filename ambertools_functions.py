import moldesign as mdt
from moldesign import units as u

__DOCKERIMAGE__ = 'docker.io/autodesk/mst:moldesign_ambertools-0.7.4b3'


def prep_ligand(mol, ligand_atom_ids, ligandname):
    """
    Create force field parameters for the chosen ligand
    """
    import moldesign as mdt

    print 'Parameterizing ligand "%s"' % ligandname

    ligand = mdt.Molecule([mol.atoms[idx] for idx in ligand_atom_ids])
    ligh = mdt.set_hybridization_and_ph(ligand, 7.4)
    params = mdt.interfaces.ambertools.parameterize(ligh, charges='am1-bcc')
    return {'ligand_parameters': params,
            'ligand': ligh}


def prep_forcefield(mol, ligand_atom_ids, ligand_params):
    """
    Assign forcefield to the protein/ligand complex
    """

    for residue in mol.residues:
        if residue.resname == 'HIS':
            print 'Guessing histidine protonation based on hydrogens present in file:'
            mdt.guess_histidine_states(mol)
            break

    stripmol = mdt.Molecule([atom for atom in mol.atoms
                             if atom.residue.type in ('dna', 'protein')
                             or atom.index in ligand_atom_ids])
    withff = mdt.interfaces.ambertools.assign_forcefield(stripmol, parameters=ligand_params)

    return {'molecule': withff,
            'prmtop': withff.ff.amber_params.prmtop,
            'inpcrd': withff.ff.amber_params.inpcrd}
