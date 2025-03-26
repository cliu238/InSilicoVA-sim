# Replication study of simulation using InSilicoVA

The codes folder contains several scripts to produce the results used in the paper: **link to paper**

+ *fit.R* and *eval-results.R*: fit and generate tables and comparison figures using the PHMRC data.
+ *fit-reduced.R* and *eval-results-reduced.R*: The same analysis but using only the subset of symptoms in the shortened PHMRC data. This analysis is just for comparison and not included in the paper, since the study in Flaxman et al. 2018 uses the full PHMRC dataset rather than the shorterned version.

Notice that the model fit script takes a very long time to run all 1000 experiments sequentially on a single machine. We recommend instead of running the for-loop over 1000 random draws, either reduce the number of simulations, or run the experiment on a high performance computing clusters in parallel. In either case, the codes require only minimal adaptation (basically, change the setup of the for-loop at the very begining).

** **docker run -it -v "$(pwd):/InSilicoVA-sim" insilicova-dev
