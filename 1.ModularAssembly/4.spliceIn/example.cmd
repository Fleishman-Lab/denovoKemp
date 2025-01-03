#!/bin/bash
cd <work_dir>
<Rosetta_exe> @flags -out:prefix <pref> -out:path:pdb <pdb_out_dir> -out:path:score <scores_out_dir> -parser:script_vars entry_blade_tail=1i4n_tail blade_tail=1i4n_tail entry_blade_1_2=1i4n_blade1_2 blade_1_2=1i4n_blade1_2 entry_blade_3_4=1i4n_blade3_4 blade_3_4=1i4n_blade3_4 entry_blade_5_6=1i4n_blade5_6 blade_5_6=1i4n_blade5_6 entry_blade_7_8=1i4n_blade7_8 blade_7_8=1i4n_blade7_8
