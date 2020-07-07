# MarkovTest
Code associated with the paper "General Tests of the Markov Assumption in Multi-state Models", A.C. Titman and H.Putter, Biostatistics (To appear).

The contents of the repository are as follows:

sleep_simulated_analysis.R : A worked example that applying the Markov testing method to a simulated dataset of similar structure to the sleep behaviour dataset in the main paper.
simulation_code.R : Routines to simulate data of a similar form to the sleep behaviour dataset.
cox_markov_test.R : General routines to implement the test of the Markov property for Markov and Cox-Markov multi-state models
plotMarkovTest.R : Function to produce log-rank trace plots based on objects produced using the cox_markov_test function.