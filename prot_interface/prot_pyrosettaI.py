#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30/01/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group


py-rosetta functions used for protein mutation
"""


import pyrosetta
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.core.scoring import * 
from pyrosetta.rosetta.protocols.docking import *
from prot_interface.logging_config import setup_logging
import logging

# Initialize logging before anything else
setup_logging()
logger = logging.getLogger(__name__)


class prot_mut_py_rosetta:
    
    """ 
    Protein Mutation Class based on Point Mutation Scan Jupyter Notebook (pyrosetta)
    https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.08-Point-Mutation-Scan.ipynb
    """
    def __init__(self,pose,partners):
        """Constructor
        
        Parameters
        ----------
        self
        pose: STRING pdb file path and name
        partners: STRING protein main chain components (ex:"A_B") 
        """
        
        self.pose=pyrosetta.io.pose_from_pdb(pose)
        self.partners=partners

   
    "This function returns the wild type amino acids in a given position."
    
    def wildtype(self,aatype = 'AA.aa_gly'):
        AA = ['G','A','L','M','F','W','K','Q','E','S','P'
            ,'V','I','C','Y','H','R','N','D','T','B']

        AA_3 = ['AA.aa_gly','AA.aa_ala','AA.aa_leu','AA.aa_met','AA.aa_phe','AA.aa_trp'
            ,'AA.aa_lys','AA.aa_gln','AA.aa_glu', 'AA.aa_ser','AA.aa_pro','AA.aa_val'
            ,'AA.aa_ile','AA.aa_cys','AA.aa_tyr','AA.aa_his','AA.aa_arg','AA.aa_asn'
            ,'AA.aa_asp','AA.aa_thr','AA.aa_asx']

        for i in range(0, len(AA_3)):
            if(aatype == AA_3[i]):
                return AA[i]

    """ 
    This function utilizes the ResidueSelectorMover() as demonstrated by previous tutorials. 
    Mutated position is allowed to be designed and repacked while the neighborhood residues are limited to repacked only.
    The final mutation will be performed with a PackMover().   
    """ 
    def pack(self,pose,posi,amino, scorefxn):

        # Select Mutate Position
        mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        mut_posi.set_index(posi)

        # Select Neighbor Position
        nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(mut_posi)
        nbr_selector.set_include_focus_in_subset(True)

        # Select No Design Area
        not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)

        # The task factory accepts all the task operations
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

        # These are pretty standard
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

        # Disable Packing
        prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
        tf.push_back(prevent_subset_repacking)

        # Disable design
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

        # Enable design
        aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        aa_to_design.aas_to_keep(amino)
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

        # Create Packer
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
        packer.task_factory(tf)

        #Perform The Move
        packer.apply(pose)

    
    """
    This function is for binding energy analysis. 
    To compute a binding energy, we need to score the total energy of a bound state structure, 
    separate (unbind) the antigen and antibody and then score the unbound state total energy. 
    The binding energy is given by bound energy - unbound energy. 
    """      

    def unbind(self, pose):
        #STEP_SIZE = 100
        STEP_SIZE = 100
        JUMP = 1 ##for two components A_B
        docking.setup_foldtree(pose, self.partners, pyrosetta.Vector1([-1,-1,-1]))
        trans_mover = rigid.RigidBodyTransMover(pose,JUMP)
        trans_mover.step_size(STEP_SIZE)
        trans_mover.apply(pose)



    def bound_score(self,posefile):
        scorefxn = get_score_function()
        pose=pyrosetta.io.pose_from_pdb(posefile)
        return scorefxn(pose)


    def unbound_score(self,posefile):
        pose=pyrosetta.io.pose_from_pdb(posefile)
        self.unbind(pose)
        scorefxn = get_score_function()
        return scorefxn(pose) 


    def mutate(self, posi, amino,gendir,filename):
        #main function for mutation

        CSV_PREFIX = filename
        PDB_PREFIX = filename


        #Initiate test pose
        testPose = pyrosetta.rosetta.core.pose.Pose()
        testPose.assign(self.pose)

        #Initiate energy function
        scorefxn = pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()
        self.unbind(testPose)
        native_ub = scorefxn(testPose)
        testPose.assign(self.pose)

        #Variables initiation
        content = ''
        name = CSV_PREFIX  + '.csv'
        pdbname = PDB_PREFIX + '.pdb'


        wt = self.wildtype(str(self.pose.aa(posi)))

        self.pack(testPose,posi, amino, scorefxn)
        testPose.dump_pdb(gendir+"/"+pdbname)
        bound = scorefxn(testPose)
        self.unbind(testPose)
        unbound = scorefxn(testPose)
        binding = unbound - bound
        testPose.assign(self.pose)

        if (wt == amino):
            wt_energy = binding
        else:
            self.pack(testPose,posi, wt, scorefxn)
            wtbound = scorefxn(testPose)
            self.unbind(testPose)
            wtunbound = scorefxn(testPose)
            wt_energy = wtunbound - wtbound
            testPose.assign(self.pose)

        if (wt_energy == 0.0):
            bwt = 0.0
        else:
            bwt = binding/wt_energy   


        content=(content+str(self.pose.pdb_info().pose2pdb(posi))+','+str(amino)+','
                  +str(native_ub)+','+str(bound)+','+str(unbound)+','+str(binding)+','
                  +str(wt_energy)+','+str(wt)+','+str(bwt)+'\n')


        f = open(gendir+"/"+name,'w+')
        f.write(content)
        f.close()

    def res_mut_list(input_str):
        res_list = []
        substrings = input_str.split(':')
        for substring in substrings:
            number, letter = substring.split('-')
            res_list.append([int(number), letter])
        return res_list

    def multi_mutate(self, mutate_res, gendir, filename):
        #main function for mutation

        CSV_PREFIX = filename
        PDB_PREFIX = filename


        #Initiate test pose
        testPose = pyrosetta.rosetta.core.pose.Pose()
        testPose.assign(self.pose)

        #Initiate energy function
        scorefxn = pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()
        self.unbind(testPose)
        native_ub = scorefxn(testPose)
        testPose.assign(self.pose)

        #Variables initiation
        content = ''
        name = CSV_PREFIX  + '.csv'
        pdbname = PDB_PREFIX + '.pdb'

        logging.info("Mutate Residues: %s", mutate_res)

        f = open(gendir+"/"+name,'w+')
        
        for res in mutate_res:
            
            posi = res[0]
            amino= res[1]

            #Mutate
            self.pack(testPose,posi, amino, scorefxn)
        

        #Save mutated pose
        testPose.dump_pdb(gendir+"/"+pdbname)
        
        # Compute energy
        bound = scorefxn(testPose)
        self.unbind(testPose)
        unbound = scorefxn(testPose)
        dG = bound - unbound 
        
        
        
        content=(content+str(self.pose.pdb_info().pose2pdb(posi))+','+str(amino)+','
                    +str(bound)+','+str(unbound)+','+str(dG)+'\n')
            
        f.write(content)
        f.close()

        # Assign original pose
        testPose.assign(self.pose)


    def sequence(self):
        return self.pose.sequence() 
