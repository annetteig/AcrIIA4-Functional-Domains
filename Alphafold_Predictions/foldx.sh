#!/usr/bin/env bash
#SBATCH --ntasks 6
#SBATCH --nodes 1
#SBATCH --mem 100G

foldx() {
    # take first (and only) argument as pdbfile
    local pdbfile=${1}
    echo working on $pdbfile

    # Run foldx on file
    path/to/foldx_20241231 --command=Stability --pdb=${pdbfile}
}

# export function, so `parallel` can 'see' it
export -f foldx

# pdb_files.txt should be a text file you make, listing all pdb files you want to iterate over
# with one file name per line. 
pdbfiles=$(cat pdb_files.txt)

# Run the `foldx` function we defined earlier across all files in pdb_files.txt, 6 at a time
parallel -j 6 foldx {} ::: ${pdbfiles}