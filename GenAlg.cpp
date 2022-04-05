
#include <iostream>
#include <mpi.h>
#include <string.h>
#include <cmath>
using namespace std;
struct Chromosome
{
    int genes[];
    Chromosome(int numGenes){
        genes = new int[numGenes];
    }
};

int evaluateFitness(Chromomsome chrom, int target);
//evaluates the expression to be solved and returns the number of genes
int evaluateExpression(string expression);

int* getCoefficients(string expression); 
//will return a list containing the number of Chromosomes with populated genes
Chromosome* createPopulation(int populationSize);

Chromosome crossoverGenes(Chromosome inputChromosome);

//
void mutate(Chromosome inputChromosome);

int g_chromosomeLength;
int g_numChromosomes;
double g_mutationRate;
double g_crossoverRate;
int chromosomeLength;
int targetValue; //Most fit functions will produce this value

int main(int argc, char *argv[]){
    string expression;
    //initialize #chromosomes, mutation rate, crossover rate...

    // HARDCODED! CHANGE LATER!
    g_chromosomeLength = 3;
    g_numChromosomes = 48;
    g_mutationRate = .2;
    g_crossoverRate = .4;
    /*
    int numChromosomes = argv[1];
    int mutationRate = argv[2];
    int crossoverRate = argv[3];
    */
    //generate chromosomes


    //loop:
    /*
    Evaluate fitness
    Select chromosomes 
    Genecrossover
    Mutation
    */
        
    return 0;
}

int evaluateFitness(Chromosome chrom){
    int tempSum = 0; 
    for (int i = 0; i < chromsomeLength; i++){ //chromsomeLength is the same as the number of coefficients in the function
        tempSum = tempSum + (chrom.genes[i] * coefficients[i]); //this will not work without coefficients
    }
    abs(tempSum - target) // the absolute value of the difference between tempSum and the target = fitness.
    // A value close to 0 means it is more fit. Farther from 0 means it is less fit.
    
    return 0;
}

int evaluateExpression(string expression){
    int numGenes=0;
    for(int i=0;i<expression.length();i++){
        if(isalpha(line[i])){
            numGenes++;
        }
    }
    return numGenes;
}

int getCoefficients(String expression){ //get coefficients from initially given function
    
    char* coefficients;
    coefficients = strtok(expression,"+") //strip away the +'s, leaving us with "1a", "44x" etc
}

void mutate(Chromosome inputChromosome)
{
    srand(time(NULL));
    int totalGene = sizeof(inputChromosome.genes) * g_numChromosomes;
    
    //Generate a random number between 1 and totalGene.
    int randNum = rand() % totalGene + 1;
    
    //If the random number is smaller than the gene, then replace it with the random number
    for (int i = 0; i < sizeof(inputChromosome.genes; i++))
    {

    }
}

Chromosome* createPopulation(int populationSize,int numGenes){
    Chromosome *population[populationSize];
    for(int i=0;i<populationSize;i++){
        population[i] = new Chromosome(numGenes);
    }
    return *population;
}
