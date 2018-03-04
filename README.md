*Data and code associated with the integrated analysis (Figure 4) in paper:*

***Developmental diversification of cortical inhibitory interneurons***

by Christian Mayer#, Christoph Hafemeister#, Rachel C. Bandler#, Robert Machold, Renata Batista Brito, Xavier Jaglin, Kathryn Allaway, Andrew Butler, Gord Fishell\* and Rahul Satija\*

\# Equal contribution
\* Corresponding authors

Published in [Nature (2018)](http://dx.doi.org/10.1038/nature25999), DOI 10.1038/nature25999

Preprint on [bioRxiv (2017)](https://www.biorxiv.org/content/early/2017/09/13/105312), DOI 10.1101/105312

The examples subfolder contains an commented R Markdown file (and compiled HTML) for the integration of E18 and P56 datasets.

The R subfolder contains code used to generate Figure 4 and components of Extended Data Figures 5-10.

### How to run

There are individual scripts for the different parts of the analysis. Run them in this order:

1. Rscript R/1_integrated_analysis.R
2. Rscript R/2_cell-assignment.R
3. Rscript R/3_conserved_marker_genes.R
4. Rscript R/4_subtype_analysis.R
5. Rscript R/5_tsne_mappings.R

