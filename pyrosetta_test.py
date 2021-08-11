import os
from sys import argv
from random import randint

import pyrosetta.rosetta
from pyrosetta import *
from pyrosetta.rosetta import *


init(extra_options="-extrachi_cutoff 12 -ex1 -ex2 -ex3 -corrections::restore_talaris_behavior")


#os.environ["MKL_NUM_THREADS"] = "1"
#os.environ["NUMEXPR_NUM_THREADS"] = "1"
#os.environ["OMP_NUM_THREADS"] = "1"


from pyrosetta.toolbox import mutate_residue


#_, template, peptide, n = argv
#template, peptide, n = 'hla_a_0201.pdb',  'SPAPPQEKL', 1
template, peptide, n = 'hla_a_0101.pdb',  'IVTSVLLLY', 1

template_model = pose_from_pdb(template)
scorefxn = create_score_function('talaris2014')
positionlist = range(template_model.pdb_info().pdb2pose('C', 1), template_model.pdb_info().pdb2pose('c', 9) + 1)

output_dir = '_models_{0}'.format(peptide)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
os.chdir(output_dir)

jd = PyJobDistributor(peptide, int(n), scorefxn)
jd.native_pose = template_model

while not jd.job_complete:
    mutant = Pose()
    mutant.assign(template_model)

    for i, res in enumerate(peptide):
        mutate_residue(mutant, positionlist[i], res, 0.0, scorefxn)
    remodel_target = Pose()
    remodel_target.assign(mutant)

    peptide_ft = rosetta.protocols.loops.Loop(positionlist[1],
                                              positionlist[7],
                                              positionlist[randint(2, 6)])
    peptide_loops = rosetta.protocols.Loops.Loops()
    peptide_loops.add_loop(peptide_ft)

    rosetta.protocols.loops.set_single_loop_fold_tree(remodel_target, peptide_ft)

    task_pack = rosetta.core.pack.task.TaskFactory.create_packer_task(remodel_target)
    task_pack.restrict_to_repacking()
    task_pack.or_include_current(True)
    pack = rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task_pack)

    pack.apply(remodel_target)

    loops_refine_CCD = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_CCD(peptide, scorefxn)
    loops_refine_CCD.max_inner_cycles(10)
    loops_refine_CCD.apply(remodel_target)

    pack.apply(remodel_target)
    jd.output_decoy(remodel_target)
