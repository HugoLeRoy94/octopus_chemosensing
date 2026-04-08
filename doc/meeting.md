# Meeting 2/13/26
## 1. Plans:
### 1.1 Diversity of Genes expression

A 0-parameters fit of the observed diversity as a function of the number of cells using the empirical expression probability of individual genes confirm that there are no correlation between the expression of different genes.

### 1.2 Fitting of dose response curves

It is hard to really infer 2 parameters / curves, as a result we try to fit $EC_{50}$ the half max activation concentration. It is not working perfectly, but I don't think it just don't work.

-> What about the amplitude

### 1.3 Mutual Information at fixed concentration

This is different from previous approaches, where the agent tries to measure the both the ligand identity AND concentration.

-> Computing the empirical mutual information shows that on the pool of ligand proposed, the array of the octopus is overall pretty bad....
-> The question of the heteromerization, in this case, is really not relevant. Comparing the homo-pentamers arrays, the difference is already massive.
-> We see that we need more than $10^4$ ligands with 5 homo-pentamers to feel the difference.
-> Adding the concentration dimension, would justify more heteromerization.

## 2. What's next ?

-> It seems that some of these components are either not part of the octopus environment or he just doesn't care.
    -> Need to model the environment
    -> Classify the chemical in familly ?
-> Looking for optimal array to specific optimization strategy might reveal why so many homomers do not respond.

# Meeting post 03/30/2026
Pablo:
-> Want to meet next week.
-> Like Fig3 and 5

Agnese:
-> Influence of the model of the environment

-> "Define complex enough" -> When independant receptors managed to behave like a perfect array.
-> "compare the distribution of activities" -> That's what I do, using the entropy. This is too high dimensional, you cannot compare them using, mean or std deviation.
-> "Discuss with LCSL" -> why not, but I would like to have a dialogue with one guy, not an endless meeting to talk about general ideas.
-> "drop the sparse chemical space" -> what does sparse means ?

Nick:
-> want to make one final figure
===> Do they have genes expression of nicotinic receptors in a vertebrate ?
|-> Ok, So comparing the genes expression of nicotinic neurotransmitter, to the sensor of octopus could help us to support the question of whether these sensors have been reused by the octopus for sensory purpose ?