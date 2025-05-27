### A GitHub repository for the AlphaFold3 informatics in "High-resolution structure of a novel contractile injection system in Salmonella enterica subsp. salamae"

**Authors**: Rooshanie N. Ejaz^1^, Kristin Funke^2^, Claudia S. Kielkopf^1^, Freddie J. O. Martin^1^, Marta Siborová^1^, Ivo A. Hendriks^1^, Nicholas Sofos^1,3^, Tillman Pape^1,3^, Eva M. Steiner-Rebrova^1,4^, Michael L. Nielsen^1,5^, Marc Erhardt^2,6^, and Nicholas M. I. Taylor^1^

1: Novo Nordisk Foundation Centre for Protein Research, Faculty of Health and Medical Sciences, University of Copenhagen, Denmark
2: Institute of Biology–Department of Molecular Microbiology, Humboldt-Universität zu Berlin, Berlin 10115, Max Planck Unit for the Science of Pathogens, Berlin 10117, Germany
3: Core Facility for Integrated Microscopy, Faculty of Health and Medical Sciences, University of Copenhagen, Denmark
4: Biomedical Centre, Department of Experimental Medical Science, Faculty of Medicine, University of Lund, Sweden
5: Evosep Biosystems, Odense, Denmark
6: Max Planck Unit for the Science of Pathogens, Berlin 10117, Germany

*To whom correspondence should be addressed*: <Nicholas M I Taylor <nicholas.taylor@cpr.ku.dk>>

Check out the preprint [here](https://www.biorxiv.org/content/)

## About
The code and work in this repo is inspired by the [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown) package by the Kosinski Lab.
We also use code to calculate the pDockQ score of the models, described [here](https://www.nature.com/articles/s41467-022-28865-w), and code adapted from [here](https://github.com/fteufel/alphafold-peptide-receptors/blob/main/pdockq.py)

Note this is not a full implementation of what AF3pulldown could be.

## Using the repository
- `SalCIS_ORF.xlsx` contains the sequences and names for all for the cassette we were investigating.
- `AF3_pulldown_allVall_prep.ipynb` takes this input and makes all the .json files for AlphaFold3.
- `AF3_batch_folder.sh` was used to run AlphaFold3 on a HPC cluster using slurm.
- `AF3_pulldown_allVall_prep.ipynb` takes the predicted results and analyses them, outputting `SALCIS_AF3_all_v_all_results.xlsx`.
- `AF3_pulldown_allVall_plotting.ipynb` plots the heatmap and other figures.

