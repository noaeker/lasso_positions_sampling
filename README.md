# lasso_positions_sampling


## Basic usage:

### The most basic usage of the code does the following:
1.The code calls raxml-ng to generate one parsimony tree for the MSA . Then,  raxml-ng is used to estimate the value of the Gamma shape parameter, alpha, on this parsimony tree.
2. The code uses raxml-ng to generate many random trees (default = 800). For each random tree, we evaluate per-site log likelihood on the MSA, assuming the same value of alpha we previously estimated for the parsimony tree. The obtained data is used as training data. The code also generates a test set of random trees with optimized branch lengths (default size = 30)
3. The code then applies Lasso model on the training data and evaluates the entire Lasso path. For each threshold (default =  1%, 2.5%, 5%, 10%), the code generates a corresponding reduced MSA and a weights file, and accuracy is evaluated based on the test-set above.

### In order to run the code, please follow the instructions below:
1. Make sure the "pandas", "sklearn", "Biopython" modules are installed in your python environment and that you have a valid installation of RAxML-NG.
1. Go to the file "main_code/config.py" and change the variable RAXML_NG_EXE  to the path of your raxml-ng executable.
2. By default, the input MSAs on which the procedure is applied are those under the folder "example_files". To change the input MSAs folder, set the variable INPUT_FILES_FOLDER under "main_code/config.py" to a folder of your choice.
3. The default assumed evolutionary model is WAG. To change the default model, set the variable EVO_MODEL to an alternative evolutionary model (e.g., JTT, LG)
4. Run the file ps.py. 

### Output files
The default output file folder is example_results. To change the output file folder, set the variable RESULTS_FOLDER under "main_code/config.py" to a folder of your cohice.
1. Summary of the obtained results (for each MSA, training-size, sample percentage...) can be found in example_results/example_run/lasso_summmary.tsv.
2. The reduced MSAs and corresponding weights file for each sampling threshold and MSA configuration (number of sequences, number of positions, training size) can be found under /example_results/example_run/job_0/example_msa.phy/n_seq_xxx/n_loci_xxx/Lasso_folder/exponential/training_xxx_random_tree_eval/trimmed_xxx/threshold_xxx_outputs.

 
 ## Advanced options:
 Advanced options can be managed based on the main_code/input_parsers.py file.
 
 ### Some of the main advanced options are:
 * The Lasso procedure can be compared to naive approaches (random/high rate sampling)
 * The Lasso procedure can be applied in a partitioned analysis.
 * The Lasso procedure can be used during a SPR search, and be compared to a standard SPR search.
