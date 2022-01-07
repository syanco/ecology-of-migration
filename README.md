# Ecology-of Migration
Data and analysis workflow for associated with:  
[Yanco, S. W., B. D. Linkhart, P. P. Marra, M. Mika, M. Ciaglo, A. Carver, and M. B. Wunder. 2021. Niche dynamics suggest ecological factors influencing migration in an insectivorous owl. Ecology:e3617.](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecy.3617?af=R)

## Workflow Description

This workflow is basically built as two RMarkdown documents that can be run end-to-end:  
1. load_and_process_data.Rmd  
2. analysis_and_plots.Rmd

The code anticipates the follwoing minimal directory structure:
```
project root
  |
  +-src       # This repo
  |
  +-data      # All raw data stored her (/data directory from OSF)
  |
  +-out       # Code-generated outputs
  |
  +-figures   # Figures generated from code
  |
  ...
```

Additional details on what the code does and how can be foud within those documents.

The data associated with this project can be found [here](https://osf.io/3n2bg/).

Ther are two additonal scripts that allow some of the model fitting to be run in parallel in a high-performance computing (HPC) environment. This parallel code was not written with generality in mind; thus, some of the code may be specific to the HPC on which the analysis was completed.

## Contact
I welcome questions, corections, or comments.  Please feel free to log issues on this repository or contact me directly by email.

e: [scott.yanco@yale.edu](mailto:scott.yanco@yale.edu?subject=[GitHub]%20ecology%20of%20migration)
