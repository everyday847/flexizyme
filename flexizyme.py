###
### This script is an example of the simulation strategy for modeling
### flexizyme-monomer interactions. For convenience, we illustrate the
### modeling method using the residue available in current released
### versions of Rosetta 'BZO'. New parameter files for new flexizyme
### substrates may be created via standard procedures for creating 
### noncanonical monomers in Rosetta. This parameter file may itself be
### edited to simulate a tetrahedral intermediate instead of a starting
### material mimic.
###
###
### This script simply enumerates the torsions that control the relative
### orientation of the nonstandard monomer and its flexizyme, then does
### a local energy minimization on those enumerated conformations. For
### finer sampling, consider smaller values of 'step' -- nothing less than
### 30 degrees should be necessary.
###
###
### Following this script, the resulting PDBs may be used in further
### Rosetta simulation to generate score distributions that assess the
### propensity to form stable interactions with the adjacent base pair.
### We used 'relax' in conjunction with (or without) constraint files 
### supplied on the command line.


from pyrosetta import *
from pyrosetta.rosetta import *
init('-nblist_autoupdate true -mute all')

pose = pose_from_pdb("3cun.pdb")
three_pa_pose = Pose(pose)

chm = rosetta.core.chemical.ChemicalManager.get_instance()
rts = chm.residue_type_set('fa_standard').get_self_ptr()
three_p_acid_residue = core.conformation.Residue(rts.name_map("BZO"), True)

relevant_seqpos = three_pa_pose.pdb_info().pdb2pose('C',92)
core.pose.remove_variant_type_from_pose_residue(three_pa_pose, core.chemical.UPPER_TERMINUS_VARIANT, relevant_seqpos)

three_pa_pose.append_residue_by_atoms(three_p_acid_residue, True, "C", relevant_seqpos, "O3'")
three_pa_pose.dump_pdb("three_p_acid_unoptimized.pdb")

select_bonded_adenosine = core.select.residue_selector.ResidueIndexSelector(relevant_seqpos)
adenosine_neighborhood = core.select.residue_selector.NeighborhoodResidueSelector(select_bonded_adenosine, 12.0)

movemap_factory = core.select.movemap.MoveMapFactory()
movemap_factory.all_chi(False)
movemap_factory.all_bb(False)
movemap_factory.all_branches(True)
movemap_factory.add_bb_action(core.select.movemap.mm_enable, adenosine_neighborhood)
movemap_factory.add_chi_action(core.select.movemap.mm_enable, adenosine_neighborhood)

movemap = movemap_factory.create_movemap_from_pose(three_pa_pose)
# Also set the specific branch PHI to move
movemap.set(core.id.DOF_ID(core.id.AtomID(three_pa_pose.residue_type(three_pa_pose.size()).atom_index('C'), three_pa_pose.size()), core.id.PHI), True)
sfxn = core.scoring.ScoreFunctionFactory.create_score_function("stepwise/rna/rna_res_level_energy4.wts")
minmover = protocols.minimization_packing.MinMover(movemap, sfxn, "lbfgs_armijo_nonmonotone", 0.001, True)

def optimize(p, sfxn):
    
    start = -180
    stop = 180
    step = 120

    def increment_torsions(torsions, start, stop, step):
        new_torsions = {k: v for k,v in torsions.items()}
        for k,v in torsions.items():
            if v < stop:
                new_torsions[k] += step
                break
            else:
                new_torsions[k] = start

        return new_torsions

    best_score = None
    best_pose = None
    torsion_names = ['phi', 'psi', 'chi1', 'chi2']
    torsions = dict(zip(torsion_names, [start]*len(torsion_names)))
    final_torsions = dict(zip(torsion_names, [stop]*len(torsion_names)))

    while torsions != final_torsions:
        q = Pose(p)
        
        #q.set_torsion(core.id.TorsionID(relevant_seqpos, core.id.BRANCH, 1), torsions['phi'])#q.size(), phi)
        #q.set_torsion(core.id.TorsionID(relevant_seqpos, core.id.BB, 6), torsions['phi'])#q.size(), phi)
        q.set_dof(core.id.DOF_ID(core.id.AtomID(q.residue_type(q.size()).atom_index('C'), q.size()), core.id.PHI), torsions['phi']) #core.id.TorsionID(q.size(), core.id.BB, 1), torsions['phi'])
        q.set_torsion(core.id.TorsionID(q.size(), core.id.BB, 1), torsions['phi'])
        q.set_torsion(core.id.TorsionID(q.size(), core.id.BB, 2), torsions['psi'])
        q.set_torsion(core.id.TorsionID(q.size(), core.id.CHI, 1), torsions['chi1'])
        q.set_torsion(core.id.TorsionID(q.size(), core.id.CHI, 2), torsions['chi2'])
        
        minmover.apply(q)
        score = sfxn(q)
        
        print(torsions)
        #print(q.torsion(core.id.TorsionID(275, core.id.BRANCH, 1)), q.phi(q.size()), q.psi(q.size()),q.chi(1,q.size()), q.chi(2, q.size()), score)
        #print(q.torsion(core.id.TorsionID(relevant_seqpos, core.id.BRANCH, 1)), q.torsion(q.size()), q.psi(q.size()),q.chi(1,q.size()), q.chi(2, q.size()), score)
        #print(q.torsion(core.id.TorsionID(relevant_seqpos, core.id.BRANCH, 1)), 
        #print(q.torsion(core.id.TorsionID(relevant_seqpos, core.id.BB, 6)), 
        #print(q.torsion(core.id.TorsionID(q.size(), core.id.BB, 1)),
        print(q.dof(core.id.DOF_ID(core.id.AtomID(q.residue_type(q.size()).atom_index('C'), q.size()), core.id.PHI)),
                q.torsion(core.id.TorsionID(q.size(), core.id.BB, 2)),
                q.torsion(core.id.TorsionID(q.size(), core.id.CHI, 1)),
                q.torsion(core.id.TorsionID(q.size(), core.id.CHI, 2)), score)
        
        if best_pose is None or score < best_score:
            best_score = score
            best_pose = Pose(q)
        
        torsions = increment_torsions(torsions, start, stop, step)

    return best_pose
        
    
best_three_pa_pose = optimize(three_pa_pose, sfxn)
best_three_pa_pose.dump_pdb("three_p_acid_optimized.pdb")

