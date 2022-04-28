#read me
omp_GenAlg.cpp is a program that implements a genetic algorithm
to solve simple equations. 

#complie
to complie the program use the command: g++ -fopenmp -o genAlg omp_GenAlg.cpp

#equation formats
    equtions should only use addition and should be formatted as such
    coefficient, variable name, space, +, space, =, space, integer  

    examples:
    "1a + 2b + 3c + 4d = 30"
    "20f + 2x + 50g = 1000"

usage: ./genAlg <num threads> <equation> <output file>

E.g. ./genAlg 4 "1a + 2b + 3c + 4d = 30" resultsFile.csv

Main functions
---------------

evaluateFitness: Evaluates the fitness value of a chromosome by comparing the solution stored inside of it to the target solution
evaluateExpression: Takes the expression that needs to be solved and returns the number of genes that should be in each chromosome
getCoefficients: Grabs the coefficients from the user-input equation and parses them into the program for solving
createPop: Creates a list that contains all the chromosomes needed with placeholder genes inside them
crossoverGenes: Generates random numbers between 1 and the length of the chromosome if a chromosome is chosen, then swaps the genes with other selected chromosomes from a cut point.
generateGenes: Replaces the placeholder genes inside all the chromosomes with actual data related to the problem.
selection: Take the population of chromosomes, compute the fitness value of each chromosome (closer to 0 is better), then calculate the probability of each chromosome moving on to the next evolution stage before getting the cumulative probability
mutate: Selects random genes from the entire population's genepool and replaces those selected genes with new randomized genes. The number of mutations is based on the mutation rate variable.
checkForSolution: Checks the chromosome population for a solution that would fit the equation's target value.