#include <iostream>
#include <mpi.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
using namespace std;
struct Chromosome
{
    public:
    vector<int> genes;
    Chromosome(int numGenes){
        for(int i = 0; i < numGenes; i++){
            genes.push_back(0);
        }
    }
};

int evaluateFitness(Chromosome chrom, int target);
//evaluates the expression to be solved and returns the number of genes //done
int evaluateExpression(string expression);

void getCoefficients(string expression,vector<int> coef); 
//will return a list containing the number of Chromosomes with populated genes
Chromosome* createPopulation(int populationSize);

Chromosome crossoverGenes(Chromosome inputChromosome);

void generateGenes(vector<chromosome> population, int limit);

//
void mutate(Chromosome inputChromosome);

int g_chromosomeLength;
int g_numChromosomes;
double g_mutationRate;
double g_crossoverRate;
int chromosomeLength;
int targetValue; //Most fit functions will produce this value
vector<int> g_coefficients; 
int main(int argc, char *argv[]){
    string expression;
    //initialize #chromosomes, mutation rate, crossover rate...
    
    
    //TODO: HARDCODED! CHANGE LATER!
    
    // Arguments: "Expression", numChromosomes
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
    int numGenes = evaluateExpression(expression);
    vector<Chromosome> population = createPopulation(10);//just using 10 as a base line
    int limit = 0;//will change this later to be what ever is on the right of =
    generateGenes(population,limit);

    //loop:
    /*
    Evaluate fitness
    Select chromosomes 
    Genecrossover
    Mutation
    */
        
    return 0;
}

//  This is responsible for taking each chromosome and calculating its fitness value. 
// it assumes that chromosomeLength and coefficients are global (for now)
int evaluateFitness(Chromosome chrom){
    int tempSum = 0; 
    for (int i = 0; i < chromsomeLength; i++){ //chromsomeLength is the same as the number of coefficients in the function
        tempSum = tempSum + (chrom.genes[i] * coefficients[i]); //this will not work without coefficients
    }
    return abs(tempSum - target); // the absolute value of the difference between tempSum and the target = fitness.
    // A value close to 0 means it is more fit. Farther from 0 means it is less fit.
}
// gets the number of letter variables from the expression
int evaluateExpression(string expression){
    int numGenes=0;
    for(int i=0;i<expression.length();i++){
        if(isalpha(line[i])){
            numGenes++;
        }
    }
    return numGenes;
}


void getCoefficients(String expression, vector<int> coef){ //get coefficients from initially given function
    // i.e if it was 2x + 5y, then we would get {2,5}
    char* tokens;
    tokens = strtok(expression,"+") //strip away the +'s, leaving us with "1a", "44x" etc
    for( int i =0; i < g_chromosomeLength; i++){
        int tempNum = strtok(tokens[i], "abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUVXYZ")[0];
        //The original expression is separated into terms, for example "5x", "2b", etc...
        //Doing strtok again would give us a number then the letters. Ex: "5x" strtok -> "5", "x"
        coef.push_back(atoi(number)); 
    }

}

//Switch genes inside the chromosome with a random number based on mutation rate
void mutate(vector<Chromosome> chromosomeVector)
{
    srand(time(NULL));
    int totalGenes;
    set<int> chosenGenes;

    //Check if vector is empty
    if (!chromosomeVector.empty())
    {
        //Calculate total number of genes
        totalGenes = chromosomeVector[0].genes.size() * g_numChromosomes;
        int numMutations = g_mutationRate * totalGenes;
        
        //Keep generating numbers until there is the expected amount of mutations
        while (chosenGenes.size() < numMutations)
        {
            int genes = rand() % totalGenes + 1;
            chosenGenes.insert(genes);
        }

        //Change the chosen genes inside the chromosomes to a random number between 0 and 30
        for (int i = 0; i < chosenGenes.size(); i++)
        {
            int chosenChromosome = totalGenes % chromosomeVector[i].genes.size();
        }
    }
    else
    //Vector is empty. No chromosomes to mutate
    {
        return;
    }
}
//function to create a population of chromosomes
vector<Chromsome> createPopulation(int populationSize,int numGenes){
    vector<Chromosome> population;
    for(int i=0;i<populationSize;i++){
        population.push_back = Chromosome c(numGenes);
    }
    return population;
}

void generateGenes(vector<Chromosome> population, int limit){
    srand(time(NULL));
    for(int i = 0; i< population.size(); i++){
        for(int j = 0; j < population[i].genes.size();j++){
            int rand = rand() % limit + 1;
            population[i].genes[j]=rand;
        }
    }
}