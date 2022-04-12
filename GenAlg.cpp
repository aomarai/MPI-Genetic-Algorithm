#include <iostream>
//#include <mpi.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
#include <ctype.h>
using namespace std;
class Chromosome
{
    public:
    vector<int> genes;
    double fitness;
    double probability;

    Chromosome(int numGenes){
        for(int i = 0; i < numGenes; i++){
            genes.push_back(0);
        }
    }
    void copy(Chromosome &chrom){
        for(int i =0; i < chrom.genes.size();i++){
            this->genes[i] = chrom.genes[i];
        }
        this->fitness=chrom.fitness;
        this->probability=chrom.probability;
    }
};

//TODO: Fix errors
//Evaluate the fitness of a chromosome by comparing its solution to the target solution
void evaluateFitness(Chromosome &chrom);

//evaluates the expression to be solved and returns the number of genes //done
int evaluateExpression(string expression);

void getCoefficients(string expression,vector<int> &coef); //done

//will return a list containing the number of Chromosomes with populated genes
void createPop(vector<Chromosome> &population, int populationSize, int numGenes); //done

//Generates random numbers between 1 - (length of Chromsome - 1). Then cuts the parent chromosome at a cut point and swaps genes with the child chromosome.
void crossoverGenes(Chromosome &parentChromosome, Chromosome &childChromosome);

void generateGenes(vector<Chromosome> &population); //done

void selection(vector<Chromosome> &population, int numGenes);
//Select random genes from the genepool and replace those genes with new randomized genes
void mutate(Chromosome inputChromosome);

void printChromosomes(vector<Chromosome> &population);

int g_numChromosomes;
double g_mutationRate;
double g_crossoverRate;
int g_chromosomeLength;
int targetValue; //Most fit functions will produce this value
vector<int> g_coefficients; 
int g_limit;
int g_target;

int main(int argc, char *argv[]){
    string expression = argv[1];
    getCoefficients(expression, g_coefficients);
    cout<<"Read in target value:"<<g_limit<<"\nRead in coefficients:\n";
    for(int i =0; i < g_coefficients.size(); i++) cout<<"Coefficient "<<i<<": "<<g_coefficients[i]<<"\n";
    //initialize #chromosomes, mutation rate, crossover rate...
    
    
    //TODO: HARDCODED! CHANGE LATER!
    
    // Arguments: "Expression", numChromosomes
    g_chromosomeLength = g_coefficients.size();
    g_numChromosomes = 48;
    g_mutationRate = .2;
    g_crossoverRate = .4;

    //generate chromosomes
    int numGenes = evaluateExpression(expression);


    vector<Chromosome> population;
    createPop(population,10,numGenes);//just using 10 as a base line


    int g_limit = 0;//will change this later to be what ever is on the right of =
    generateGenes(population);


    
    printChromosomes(population);
    

    //loop:

    //Evaluate fitness
    for(int i = 0; i < population.size();i++) evaluateFitness(population[i]);
    //Select chromosomes 


        selection(population, numGenes );

    //Genecrossover
    //Mutation
    
        
    return 0;
}

//  This is responsible for taking each chromosome and calculating its fitness value. 
// it assumes that chromosomeLength and coefficients are global (for now)
void evaluateFitness(Chromosome &chrom){
    int tempSum = 0; 
    for (int i = 0; i < g_chromosomeLength; i++){ //chromsomeLength is the same as the number of coefficients in the function
        tempSum = tempSum + (chrom.genes[i] * g_coefficients[i]); //this will not work without coefficients
    }
    chrom.fitness = abs(tempSum - g_limit); // the absolute value of the difference between tempSum and the target = fitness.
    // A value close to 0 means it is more fit. Farther from 0 means it is less fit.
}
// gets the number of letter variables from the expression
int evaluateExpression(string expression){
    int numGenes=0;
    for(int i=0;i<expression.length();i++){
        if(isalpha(expression[i])){
            numGenes++;
        }
    }
    return numGenes;
}

//get coefficients from initially given function
void getCoefficients(string expression, vector<int> &coef){ //works as intended
    // i.e if it was 2x + 5y, then we would get {2,5}

    //string expression = "22x + 555y + 66z = 50";
    int startIndex = -1;
    int endIndex = -1;
    for (int i =0; i < expression.length(); i++){
        if(startIndex ==-1 && isdigit(expression[i])) startIndex = i;
        else if (endIndex == -1 && isalpha(expression[i])) endIndex = i;
        if (startIndex > -1 && endIndex > -1){
            string substring = expression.substr(startIndex,endIndex);
            coef.push_back(stoi(substring)); 
            startIndex = -1;
            endIndex = -1;
        }  
    } 
    string afterEquals = expression.substr(expression.find("=")+1); // get number after = sign. Ie: in 20x + 50y = 10, get " 10"
    g_limit = stoi(afterEquals); // set our fitness target to that value


}

//Switch genes inside the chromosome with a random number based on mutation rate
void mutate(vector<Chromosome> &chromosomeVector)
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
            int specificGene =  totalGenes % chromosomeVector[i].genes.size()-1;
            int chosenChromosome = specificGene % chromosomeVector[i].genes.size();
            chromosomeVector[chosenChromosome].genes[specificGene]=(rand() % g_limit)+1;//TODO: might have to check this later
        }
    }
    else return;     //Vector is empty. No chromosomes to mutate
}
//function to create a population of chromosomes
void createPop(vector<Chromosome> &population, int populationSize,int numGenes){//works as intended
    for(int i=0;i<populationSize;i++){
        Chromosome c(numGenes);
        population.push_back(c);
    }
}

//Takes two chromosomes and crosses their genes at a given crossover index in the gene list
void crossoverGenes(Chromosome &parentChromosome, Chromosome &childChromosome)
{
    srand(time(NULL));
    int crossoverPoint = rand() % (parentChromosome.genes.size() - 1);

    //Replace the parent chromosome's genes with the child chromosome's genes from the crossover point onward
    for (int i = crossoverPoint; i < parentChromosome.genes.size(); i++)
    {
        parentChromosome.genes[i] = childChromosome.genes[i];
    }
}

void generateGenes(vector<Chromosome> &population){//works as intended
    srand(time(NULL));
    for(int i = 0; i< population.size(); i++){
        for(int j = 0; j < population[i].genes.size();j++){
            int randInt = rand() % (g_limit + 1);
            population[i].genes[j]=randInt;
        }
    }
}


void selection(vector<Chromosome> &population,int numGenes){
    srand((unsigned)time(NULL));
    int total;
    vector<double> cumulativeProb;
    vector<Chromosome> pop2;
    for(int i = 0;i<population.size();i++){
        population[i].fitness = 1 / (1+population[i].fitness);
        total = population[i].fitness;
    }

    for(int i=0;i<population.size();i++){
        population[i].probability = population[i].fitness / total;
    }

    for(int i=0;i<population.size();i++){
        cumulativeProb.push_back(0.0);
        for(int j=0; j<=i; j++){
            cumulativeProb[i] += population[j].probability; 
        }
    }


    for(int i = 0;i<population.size();i++){
        double randD = rand() % RAND_MAX;\


        for(int j = 0;j<population.size()-1;j++){


            if(randD <= cumulativeProb[i+1] && randD >= cumulativeProb[i]){
                //might have to change things here with 
                Chromosome c(population[i].genes.size());
                c.copy(population[i]);
            }
        }
    }
}

//Print out the the population's genes by chromosome
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

