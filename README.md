# Oropharyngeal Cancer Agent-Based Model Sensitivity Analysis
This repo contains all the scripts required to conduct a sensitivity analysis for an oropharyngeal cancer agent-based model hosted at https://github.com/anniedang1234/Oropharyngeal_Cancer_ABM.

## Generating cell position file
This script is used to select a portion of the spatial transcriptomics data to run sensitivity analysis on, to reduce computational cost. This script finds the subsection whose cell type and cell gene expression most closely matches that of the full data.

## Generating parameters
Run sample_generate.py to generate a csv file containing all the sets of parameters.
Run generate.py to parse those sets of parameters into individual csv files.
```bash
py sample_generate.py
py generate.py
```

## Running the sensitivity analysis
To run the analysis on DRAC, use the following in your slurm script:
```bash
module purge
module load StdEnv/2023 gcc/14.3
module load cuda/12.9
module load python/3.11
module load vtk/9.6.0

virtualenv $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r ~/compucell3d.reqs

export PYTHONPATH=$PYTHONPATH:$HOME/CompuCell3D/lib/site-packages

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

all_tumour_counts="/home/annied/Results_SA_${SLURM_JOB_ID}/all_tumour_counts.csv"

for i in $(seq 0 199); do 
    cp /home/annied/OP_Cancer_2D_SA/samples/sample_${i}.csv /home/annied/OP_Cancer_2D_SA/parameters.csv
    python -X faulthandler -m cc3d.run_script -i /home/annied/OP_Cancer_2D_SA/OP_Cancer_2D_SA.cc3d -f 10 -o /home/annied/Results_SA_${SLURM_JOB_ID}/run_${i}
    
    cd "/home/annied/Results_SA_${SLURM_JOB_ID}"

    if [ -f "run_${i}/tumour_count.txt" ]; then
        read -r tumour_count < "run_${i}/tumour_count.txt" # Read final tumour count CSV file
        echo $tumour_count >> "$all_tumour_counts"
    else
        echo "${i} FAILED" >> $all_tumour_counts
    fi
    rm -rf /home/annied/Results_SA_${SLURM_JOB_ID}/run_${i}
done
```
