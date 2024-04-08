# MarkovTest
Code associated with the paper "General Tests of the Markov Assumption in Multi-state Models", A.C. Titman and H.Putter, Biostatistics (DOI: 10.1093/biostatistics/kxaa030).

The contents of the repository are as follows:

sleep_simulated_analysis.R : A worked example that applying the Markov testing method to a simulated dataset of similar structure to the sleep behaviour dataset in the main paper.
sleep_simulated_analysis.Rmd : R Markdown version of the worked example.
sleep_simulated_analysis.pdf : Pdf document generated from the R Markdown worked example.
simulation_code.R : Routines to simulate data of a similar form to the sleep behaviour dataset.
cox_markov_test.R : General routines to implement the test of the Markov property for Markov and Cox-Markov multi-state models
plotMarkovTest.R : Function to produce log-rank trace plots based on objects produced using the cox_markov_test function.

Note that the functions for implementing the test have been incorporated into the mstate R package via the MarkovTest function.
