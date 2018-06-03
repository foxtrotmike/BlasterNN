# BlasterNN
BLAST nearest neighbor

Implementation of homology based search. Given a training and test FASTA file, this function will return a dictionary of the lowest e-values of a test protein against the training dataset. "idx" is used to assign an id to the output file. 
Note: it will creat a blast database file based on the training file so the location of the trainfile must be writeable. It will also generate the output file in the current directory. 
Joblib is used to perform parallel testing.
