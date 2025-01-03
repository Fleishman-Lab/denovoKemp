#!/bin/bash
cd <work_dir>
<rosetta_exe> @flags -s <home>/working_pdbs/1i4n.pdb -out:path:pdb <out_pdb_dir> -out:path:score <out_scores_dir> 
