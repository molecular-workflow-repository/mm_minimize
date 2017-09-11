import moldesign as mdt
from moldesign import units as u


__DOCKERIMAGE__ = 'docker.io/autodesk/mst:moldesign_subprocess-0.7.4b3'


def get_ligands(mol):
    """ Return a dict of possible ligands in the molecule.
    dict is of the form {ligand_name: [atom_idx1, atom_idx2, ...], ...}
    where the atom_idxN are the 0-based indices of the atoms comprising each potential ligand

    Ligands can span, at most, 2 different residues
    """
    ligands = [residue for residue in mol.residues if residue.type == 'unknown']

    found_ligands = {}
    for ligand in ligands:
        bound_residues = list(ligand.bonded_residues)

        if len(bound_residues) == 0:
            found_ligands[ligand.name] = [atom.index for atom in ligand.atoms]
        elif len(bound_residues) == 1:
            nbr = bound_residues[0]
            if len(list(nbr.bonded_residues)) != 1:
                continue
            else:
                ligname = '%s - %s' % (ligand.name, nbr.name)
                found_ligands[ligname] = [atom.idx for atom in ligand.atoms+nbr.atoms]

    return found_ligands


def validate(mol, ligands):
    missing = missing_internal_residues(mol)
    all_errors = []
    success = True

    if len(ligands) == 0:
        success = False
        all_errors.append('No ligands found in this structure.')

    for chain_name, reslist in missing.iteritems():
        success = False
        all_errors.append('The following residues are not present in the PDB file: %s'
                          % ','.join(name+num for num, name in reslist.iteritems()))

    return {'success': success,
            'errors': ' '.join(all_errors),
            'pdbstring': mol.write('pdb'),
            'ligands': ligands}


def mm_minimization(mol):
    mol.set_energy_model(mdt.models.OpenMMPotential, implicit_solvent='obc')

    # start with a native gradient descent -
    # mostly to provide some intermediate states for animation
    traj = mdt.min.gradient_descent(mol, frame_interval=250, nsteps=1000)
    traj2 = mol.minimize()  # native openmm optimization
    traj.new_frame() # Add final OpenMM state to the native steepest descent trajectory

    return traj


def get_outputs(traj):
    rmsd_json = traj.rmsd()[-1].to(u.angstrom).to_json()
    rmsd_json['name'] = 'RMSD'

    stabilized_energy = traj.frames[0].potential_energy.to(u.kcalpermol).to_json()
    stabilized_energy['units'] = 'kcal/mol'
    stabilized_energy['name'] = 'Energy stabilization'

    final_energy_json = mol.potential_energy.to(u.kcalpermol).to_json()
    final_energy_json['units'] = 'kcal/mol'
    final_energy_json['name'] = 'Final energy'

    results['output_values'] = [final_energy_json, stabilized_energy, rmsd_json]

    return traj_pdb, results
