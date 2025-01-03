#!/bin/bash
cd <work_dir>
<Rosetta_exe> @flags -out:prefix 1i4n_tail_ -out:path:pdb <pdb_out_dir> -out:path:score <scores_out_dir> -parser:script_vars source=../pairfit/tail/1i4n_tail.pdb db=<db_out_path> before=1-0 after=43-250
