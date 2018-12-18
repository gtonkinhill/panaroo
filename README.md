# squeaky
An updated pipeline for pan-genome investigation


### Pipeline

1. Run prodigal on full dataset (so as not to run into issue with its model inference step)
2. Order inferred gene sequences from largest to smallest
3. Cluster from largest to smallest using jackhmmer
4. Classify resulting models using profile hmm vs profile hmm tool that I need to find again.
