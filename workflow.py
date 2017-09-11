import molflow
import molflow.definitions as mf

workflow = mf.WorkflowDefinition('MM Minimization')

# Function definitions
convert = mf.get_workflow('convert', version='0.1.0')
get_ligands = mf.Function(sourcefile='./functions.py', funcname='get_ligands')
validate = mf.Function(sourcefile='./functions.py', funcname='validate')
minimize = mf.Function(sourcefile='./functions.py', funcname='mm_minimization')
prep_ligand = mf.Function(sourcefile='./ambertools_functions.py', funcname='prep_ligand')
prep_forcefield = mf.Function(sourcefile='./ambertools_functions.py', funcname='prep_ligand')
mm_minimize = mf.Function(sourcefile='./functions.py',  funcname='mm_minimize')
create_outputs = mf.Function(sourcefile='./functions.py', funcname='create_outputs')

# Workflow inputs
molfile = workflow.add_input('molecule_file',
                             'input file')

### DAG definition
mol = convert(molfile, to_fmt='mdt-0.8')
ligands = get_ligands(mol)
validate(mol, ligands)
ligand_parameters = molflow.map(prep_ligand, ligands)
mol_with_ff = prep_forcefield(mol, ligand_parameters)
trajectory = mm_minimize(mol_with_ff)
trajfile, resultfile = create_outputs(trajectory)

workflow.set_output(trajectory, 'traj.mdt0_8.pkl')
workflow.set_output(trajfile, 'traj.pdb')
workflow.set_output(resultfile, 'results.json')
