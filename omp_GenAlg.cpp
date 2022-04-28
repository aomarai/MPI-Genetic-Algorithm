/*
* IT 388
* Genetic algorithm by John Skluzacek, Ash Omaraie, and Joe Katalinich
*/
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
#include <random>
#include <ctype.h>
#include <omp.h>
#include <fstream>
using namespace std;
class Chromosome
{
public:
    vector<int> genes;
    double fitness;
    double probability;

    Chromosome(int numGenes)
    {
        //#pragma omp parallel for
        for (int i = 0; i < numGenes; i++)
        {
            //#pragma omp critical
            genes.push_back(0);
        }
    }
    void copy(Chromosome &chrom)
    {
        for (int i = 0; i < chrom.genes.size(); i++)
        {
            this->genes[i] = chrom.genes[i];
        }
        this->fitness = chrom.fitness;
        this->probability = chrom.probability;
    }
};

// Evaluate the fitness of a chromosome by comparing its solution to the target solution
void evaluateFitness(vector<Chromosome> &population);

// evaluates the expression to be solved and returns the number of genes //done
int evaluateExpression(string expression);

// Grabs the coefficients from the user provided expression
void getCoefficients(string expression, vector<int> &coef); // done

// Create a list containing the number of Chromosomes with populated genes
void createPop(vector<Chromosome> &population, int populationSize, int numGenes); // done

// Generates random numbers between 1 - (length of Chromosome - 1). Then cuts the parent chromosome at a cut point and swaps genes with the child chromosome.
void crossoverGenes(vector<Chromosome> &population);

// Generate the genes for the population
void generateGenes(vector<Chromosome> &population); // done

// Take the population, compute the fitness of each chromosome, then calculate the probability of each chromosome before getting the cumulative probability
void selection(vector<Chromosome> &population, int numGenes);

// Helper function for selection
int selectionHelper(vector<double> cumulativeProb, double randomNum);

// Select random genes from the genepool and replace those genes with new randomized genes. Based on mutation rate.
void mutate(vector<Chromosome> &population);

// Iterates through the population and prints the genes for every chromosome
void printChromosomes(vector<Chromosome> &population);

// Checks if a viable solution has been found in the population
int checkForSolution(vector<Chromosome> &population);

int g_numChromosomes;
double g_mutationRate;
double g_crossoverRate;
int g_chromosomeLength;
int targetValue; // Most fit functions will produce this value
vector<int> g_coefficients;
int g_limit;
int g_target;
double g_tolerance;
int g_iterationLimit;

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <num threads> <equation> <output file>" << endl;
        return 1;
    }

    double startMain = omp_get_wtime();
    string filename = "results.csv";

    if (argv[3] != NULL)
        filename = argv[3];

    bool appendColumnNames = false;
    ifstream fileStream;
    fileStream.open(filename);
    if (fileStream.fail())
        appendColumnNames = true;
    ofstream myFile(filename, myFile.app);
    // numThreads,totalTime,totalTime,selectionTime,mutationTime,crossoverTime
    if (appendColumnNames)
        myFile << "numThreads"
               << ","
               << "totalTime"
               << ","
               << "totalTime"
               << ","
               << "selectionTime"
               << ","
               << "mutationTime"
               << ","
               << "crossoverTime" << endl;

    // Analyze pieces of program separately (in terms of speedup)
    string expression = argv[2];
    int numThreads = atoi(argv[1]);
    getCoefficients(expression, g_coefficients);
    cout << "Read in target value:" << g_limit << "\nRead in coefficients:\n";
    // initialize #chromosomes, mutation rate, crossover rate...
    omp_set_num_threads(numThreads);

    // Arguments: "Expression", numChromosomes
    g_chromosomeLength = g_coefficients.size();
    g_numChromosomes = 10000;
    g_mutationRate = .1; // change later
    g_crossoverRate = .25;
    g_tolerance = 0;
    g_iterationLimit = 100;

    // generate chromosomes
    //  was 1000
    int numGenes = evaluateExpression(expression);
    double totalTimeElapsed = 0.0;
    vector<Chromosome> population;

    // Give each processor its index and working range to modify
    double start = omp_get_wtime();

    createPop(population, g_numChromosomes, numGenes);
    double done = omp_get_wtime() - start;

    // Measure time to generate genes
    start = omp_get_wtime();
    generateGenes(population);
    done = omp_get_wtime() - start;

    // loop:
    int indexOfsol = -1;
    bool complete = false;

    // Evaluate fitness
    start = omp_get_wtime();
    evaluateFitness(population);
    done = omp_get_wtime() - start;

    printf("TIME TO EVAL FITNESS: %f\n", done);
    indexOfsol = checkForSolution(population);
    if (indexOfsol > -1)
    {
        complete = true;
    }
    int iterNum = 0;

    // Begin iterating
    double totalTimeSpentOnFitnessEvaluation;
    double totalTimeSpentOnSelection;
    double totalTimeSpentOnCrossover;
    double totalTimeSpentOnMutate;
    for (int i = 0; i < g_iterationLimit; i++)
    {

        // Select chromosomes
        start = omp_get_wtime();
        selection(population, numGenes);
        done = omp_get_wtime() - start;
        totalTimeSpentOnSelection += done;

        // Genecrossover
        start = omp_get_wtime();
        crossoverGenes(population);
        done = omp_get_wtime() - start;
        totalTimeSpentOnCrossover += done;

        // Mutation
        start = omp_get_wtime();
        mutate(population);
        done = omp_get_wtime() - start;
        totalTimeSpentOnMutate += done;
        start = omp_get_wtime();

        //#pragma omp parallel for
        // for (int i = 0; i < population.size(); i++)
        //{
        evaluateFitness(population);
        //}
        done = omp_get_wtime() - start;
        totalTimeSpentOnFitnessEvaluation += done;
        indexOfsol = checkForSolution(population);
        if (indexOfsol != -1)
        {
            complete = true;
        }
        iterNum++;
    }
    printf("Iteration Number: %d\n", iterNum);
    cout << "index of sol: " << indexOfsol << endl;
    if (indexOfsol != -1)
    {
        printf("Coefficients: ");
        for (int i = 0; i < g_coefficients.size(); i++)
            cout << "" << g_coefficients[i] << ", ";
        printf("\n");
        printf("Solution Chromosome %d's Genes: ", indexOfsol);
        for (int j = 0; j < population[indexOfsol].genes.size(); j++)
        {
            printf("%d ", population[indexOfsol].genes[j]);
        }
    }
    else
        printf("No solution found.\n");

    double endMain = omp_get_wtime() - startMain;
    printf("Total time to run main: %f\nTotal time to eval fitness: %f\nTotal Time spent on selection: %f\nTotal Time spent on mutation: %f\nTotal time spent on crossover:%f\n", endMain, totalTimeSpentOnFitnessEvaluation, totalTimeSpentOnSelection, totalTimeSpentOnMutate, totalTimeSpentOnCrossover);
    myFile << numThreads << "," << endMain << "," << totalTimeSpentOnFitnessEvaluation << "," << totalTimeSpentOnSelection << "," << totalTimeSpentOnMutate << "," << totalTimeSpentOnCrossover << endl;
    myFile.close();
    return 0;
}

int checkForSolution(vector<Chromosome> &population)
{
    int solution = -1;
#pragma omp parallel for // reduction (min:solution)
    for (int i = 0; i < population.size(); i++)
    {
        if (population[i].fitness <= (g_tolerance * g_limit))
        {
            solution = i;
        }
    }
    return solution;
}
//  This is responsible for taking each chromosome and calculating its fitness value.
// it assumes that chromosomeLength and coefficients are global (for now)
void evaluateFitness(vector<Chromosome> &population)
{
#pragma omp parallel for
    for (int j = 0; j < population.size(); j++)
    {
        int tempSum = 0;
        for (int i = 0; i < g_chromosomeLength; i++)
        {                                                            // chromsomeLength is the same as the number of coefficients in the function
            tempSum += (population[j].genes[i] * g_coefficients[i]); // this will not work without coefficients
        }
        // the absolute value of the difference between tempSum and the target = fitness.
        population[j].fitness = abs(tempSum - g_limit);
        // A value close to 0 means it is more fit. Farther from 0 means it is less fit.
    }
}

// gets the number of letter variables from the expression
int evaluateExpression(string expression)
{
    int numGenes = 0;
    for (int i = 0; i < expression.length(); i++)
    {
        if (isalpha(expression[i]))
        {
            numGenes++;
        }
    }
    return numGenes;
}

// get coefficients from initially given function
void getCoefficients(string expression, vector<int> &coef)
{
    // i.e if it was 2x + 5y, then we would get {2,5}
    // string expression = "22x + 555y + 66z = 50";
    int startIndex = -1;
    int endIndex = -1;
    for (int i = 0; i < expression.length(); i++)
    {
        if (startIndex == -1 && isdigit(expression[i]))
            startIndex = i;
        else if (endIndex == -1 && isalpha(expression[i]))
            endIndex = i;
        if (startIndex > -1 && endIndex > -1)
        {
            string substring = expression.substr(startIndex, endIndex);
            coef.push_back(stoi(substring));
            startIndex = -1;
            endIndex = -1;
        }
    }
    string afterEquals = expression.substr(expression.find("=") + 1); // get number after = sign. Ie: in 20x + 50y = 10, get " 10"
    g_limit = stoi(afterEquals);                                      // set our fitness target to that value
}

// Switch genes inside the chromosome with a random number based on mutation rate
void mutate(vector<Chromosome> &chromosomeVector)
{
    srand(time(NULL));
    int totalGenes;
    vector<int> chosenGenes;
    // Check if vector is empty
    if (!chromosomeVector.empty())
    {
        // Calculate total number of genes
        totalGenes = g_chromosomeLength * g_numChromosomes;
        int numMutations = g_mutationRate * totalGenes;
        // Keep generating numbers until there is the expected amount of mutations
        while (chosenGenes.size() < numMutations)
        {
            int genes = rand() % totalGenes;
            chosenGenes.push_back(genes);
        }
        // Change the chosen genes inside the chromosomes to a random number between 0 and target (g_limit)
        #pragma omp parallel for
        for (int i = 0; i < chosenGenes.size(); i++)
        {
            int indexOfChromosome = chosenGenes[i] / g_chromosomeLength; // Get the index of the chrosomome.
            int indexOfGene = chosenGenes[i] % g_chromosomeLength;
            chromosomeVector[indexOfChromosome].genes[indexOfGene] = (rand() % g_limit) + 1; // set that targeted gene to a number between 0 and g_limit
        }
    }
    else
        return; // Vector is empty. No chromosomes to mutate
}
// function to create a population of chromosomes
void createPop(vector<Chromosome> &population, int populationSize, int numGenes)
{
    for (int i = 0; i < populationSize; i++)
    {
        Chromosome c(numGenes);
        population.push_back(c);
    }
}

// Takes two chromosomes and crosses their genes at a given crossover index in the gene list
void crossoverGenes(vector<Chromosome> &population)
{
    srand(time(NULL));
    vector<int> parents(population.size(), 0);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<float> distribution(0.0, 1.0);
    /*
    Select the parents by generating a random number between 0 and 1.
    If that number is smaller than the crossover rate, then add the index of the chromosome to the vector of parents.
    */
    int filledParents = 0;
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        float randCutNum = distribution(gen);
        if (randCutNum < g_crossoverRate)
        {
            parents[filledParents] = i;
            filledParents++;
        }
    }
    // Crossross over all the parents with each other from the cut point onwards using the parent indexes
    #pragma omp parallel for
    for (int i = 0; i < filledParents; i++)
    {
        int cutPoint = rand() % (population[0].genes.size() - 1);
        for (int j = i + 1; j < filledParents; j++)
        {
            for (int k = cutPoint; k < population[0].genes.size(); k++)
            {
                int temp = population[parents[i]].genes[k];
                population[parents[i]].genes[k] = population[parents[j]].genes[k];
                population[parents[j]].genes[k] = temp;
            }
        }
    }
}

void generateGenes(vector<Chromosome> &population)
{
    srand(time(NULL));
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        for (int j = 0; j < population[i].genes.size(); j++)
        {
            int randInt = rand() % (g_limit + 1);
            population[i].genes[j] = randInt;
        }
    }
}

void selection(vector<Chromosome> &population, int numGenes)
{
    // Take population, compute fitness of each chromosome,
    // then calculate the probably of each chromosome to be selected
    // Then get the cumulative probability
    srand((unsigned)time(NULL));
    double total = 0.0;
    vector<double> cumulativeProb(population.size(), 0.0);
    vector<Chromosome> pop2 = population;
    // computing the fitness
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        population[i].fitness = 1 / (1 + population[i].fitness);
        total += population[i].fitness;
    }
    // probability of each chromosome
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        population[i].probability = population[i].fitness / total;
    }
    // cumulative probability
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        for (int j = 0; j <= i; j++)
        {
            cumulativeProb[i] += population[j].probability;
        }
    }

    // selection of the chromosomes
    // looping through the population to determine which chromosome we will copy
    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        double randD = ((double)rand() / (RAND_MAX));
        int indexOfSelectedChromosome = selectionHelper(cumulativeProb, randD);
        Chromosome c(population[indexOfSelectedChromosome].genes.size());
        c.copy(population[indexOfSelectedChromosome]);
        pop2[i] = c;
    }

    #pragma omp parallel for
    for (int i = 0; i < population.size(); i++)
    {
        population[i] = pop2[i];
    }
    pop2.clear();
}

int selectionHelper(vector<double> cumulativeProb, double randomNum)
{
    double diff = 99.9;
    int indexOfClosest = 0;
    for (int i = 0; i < cumulativeProb.size(); i++)
    {

        double tempDiff = abs(cumulativeProb[i] - randomNum);
        if (tempDiff < diff)
        {
            diff = tempDiff;
            indexOfClosest = i;
        }
    }
    return indexOfClosest;
}

// Print out the the population's genes by chromosome
void printChromosomes(vector<Chromosome> &population)
{
    for (int i = 0; i < population.size(); i++)
    {
        printf("Chromosome %d's Genes: ", i);
        for (int j = 0; j < population[i].genes.size(); j++)
        {
            printf("%d ", population[i].genes[j]);
        }
        cout << endl;
    }
}
