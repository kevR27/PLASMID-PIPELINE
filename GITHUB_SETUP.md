# GitHub repository and server deployment

This pipeline is maintained in the independent repository:

https://github.com/kevR27/PLASMID-PIPELINE

The repository is **public**. Review reference/seed FASTAs and code before publication; never push raw FASTQs, patient data, credentials, analysis results or intermediate BAMs.

## Clone to a new location

```bash
git clone https://github.com/kevR27/PLASMID-PIPELINE.git
cd PLASMID-PIPELINE
# Edit config/config.yaml and set reads_folder for your server
conda env update -n plasmid_pipeline -f environment.yaml --prune
python3 -m unittest discover -s tests -v
```

## Update the existing server checkout

The user's latest working copy is at `/DATA/casadei7/tools/plasmid_pipeline`. It may include edits not in this snapshot. Compare before overwriting:

```bash
cd /DATA/casadei7/tools/plasmid_pipeline
git status --short
git remote -v
```

Back up local configuration and uncommitted scripts, then pull or merge the repository only if safe. Avoid `git reset --hard` and avoid deleting `results/` or `work/` to sync code.

The default workflow writes `results/` and `work/` inside its project directory. Selecting a custom result root is not yet implemented through a supported `output_dir` configuration option.
