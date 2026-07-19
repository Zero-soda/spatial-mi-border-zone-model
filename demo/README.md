# SSTBA Real-Data Demonstration

`mouse_master_spot_scores.tsv` is a nine-column extract of
`source_data/Source_Data_Master_Spot_Level_Table.tsv`. It contains all 14,147
tissue spots from the six GSE214611 discovery sections and preserves the
original row order used by the released six-nearest-neighbour analysis.

The extract contains public-data-derived identifiers, array coordinates,
author-provided spatial-domain labels and the three frozen state scores. It
does not contain patient identifiers, hospital data or newly generated
experimental data.

From the code-release root, run:

```bash
python3 -m sstba.cli validate --config config/demo_boundary.json
python3 -m sstba.cli run \
  --config config/demo_boundary.json \
  --output results/demo
```

The demonstration constructs a symmetrized six-nearest-neighbour graph within
each section, identifies only domain 3-domain 4 edges, calculates per-domain
boundary fractions and reports both whole-domain and edge-local score
gradients. `run_manifest.json` records parameters, versions and file hashes.

This demonstration reproduces graph-based spatial summaries from processed
scores. Recomputing scores from raw count matrices requires the public inputs
and ordered scripts documented in the main repository README.
