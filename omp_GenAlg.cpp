#include <iostream>
//#include <mpi.h>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
#include <random>
#include <ctype.h>
#include <omp.h>
using namespace std;
class Chromosome
{
    public:
    vector<int> genes;
    double fitness;
    double probability;

    Chromosome(int numGenes){
        //#pragma omp parallel for
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
void crossoverGenes(vector<Chromosome> &population);

void generateGenes(vector<Chromosome> &population); //done

void selection(vector<Chromosome> &population, int numGenes);
int selectionHelper(vector<double>  cumulativeProb, double randomNum);
//Select random genes from the genepool and replace those genes with new randomized genes
void mutate(vector<Chromosome> &population);

void printChromosomes(vector<Chromosome> &population);
int checkForSolution(vector<Chromosome> &population);
int g_numChromosomes;
double g_mutationRate;
double g_crossoverRate;
int g_chromosomeLength;
int targetValue; //Most fit functions will produce this value
vector<int> g_coefficients; 
int g_limit;
int g_target;
double g_tolerance;
int g_iterationLimit;

int main(int argc, char *argv[]){
/// Start with openMP
/// Maybe use MPI

// Analyze pieces of program separately (in terms of speedup)



    string expression = argv[1];
    //getCoefficients(expression, g_coefficients);
    cout<<"Read in target value:"<<g_limit<<"\nRead in coefficients:\n";
    //initialize #chromosomes, mutation rate, crossover rate...
    omp_set_num_threads(10);
    // Arguments: "Expression", numChromosomes
    g_chromosomeLength = 10;//g_coefficients.size();
    g_numChromosomes = 100000;
    g_mutationRate = .1;//change later
    g_crossoverRate = .25;
    g_tolerance = 0;
    g_iterationLimit = 100;

    //generate chromosomes
    int numGenes = 1000;//evaluateExpression(expression);
    double totalTimeElapsed = 0.0;

    vector<Chromosome> population;

    // WE SHOULD PARALLELIZE THIS.
    // Give each processor its index and working range to modify
    double start = omp_get_wtime();
    printf("Trying to grab omp time\n");
    createPop(population,1000000,numGenes);
    //cout<<"here\n"<<endl;
    double done = omp_get_wtime() - start;
    //totalTimeElapsed+=done;
    printf("TIME TO GENERATE POP: %f\n\n\n", done);
    start = omp_get_wtime();
    generateGenes(population);
    done = omp_get_wtime() - start;
    //totalTimeElapsed+=done;
    printf("TIME TO GENERATE GENES: %f\n\n\n", done);
    return 0;
    //loop:
    int indexOfsol = -1;
    bool complete = false;
    //Evaluate fitness
    start = omp_get_wtime();
    //#pragma omp parallel for
    for(int i = 0; i < population.size();i++){//this can be parallelized
        evaluateFitness(population[i]);
    }
    done = omp_get_wtime() - start;
    //totalTimeElapsed+=done;
    printf("TIME TO EVAL FITNESS: %f\n\n\n", done);
    indexOfsol = checkForSolution(population);
    if(indexOfsol > -1){
            complete = true;
    }
    int iterNum = 0;
    double startLoop = omp_get_wtime();
    //Begin iterating
    for (int i = 0; i < g_iterationLimit; i++){
        //Select chromosomes 
        selection(population, numGenes );
        //Genecrossover
        crossoverGenes(population);
        //Mutation
        mutate(population);
        for(int i = 0; i < population.size();i++){//this can be parallelized 
            evaluateFitness(population[i]);
        }
        indexOfsol = checkForSolution(population);
        if(indexOfsol!=-1){
            complete = true;
        }
        iterNum++;
    }
    double loopDone = omp_get_wtime() - startLoop;
    printf("TIME TO RUN FOR LOOP: %f\n", loopDone);
    //printChromosomes(population);
    printf("Iteration Number: %d\n", iterNum);
    cout<<"index of sol: "<<indexOfsol<<endl;
    if(indexOfsol !=-1){
        printf("Coefficients: ");
        for(int i =0; i < g_coefficients.size(); i++) cout<<""<<g_coefficients[i]<<", ";
        printf("\n");
        printf("Solution Chromosome %d's Genes: ", indexOfsol);
        for (int j = 0; j < population[indexOfsol].genes.size(); j++)
        {
            printf("%d ", population[indexOfsol].genes[j]);
        }
    } else printf("No solution found.\n");

    return 0;
}


int checkForSolution(vector<Chromosome> &population){//this can be parallelized
    for(int i = 0;i < population.size();i++){
        //printf("Fitness gram pacer test %f\n", population[i].fitness);
        if(population[i].fitness<=(g_tolerance * g_limit)){
            return i;
        }
    }
    return -1;
}
//  This is responsible for taking each chromosome and calculating its fitness value. 
// it assumes that chromosomeLength and coefficients are global (for now)
void evaluateFitness(Chromosome &chrom){ //this can be parallelized
    int tempSum = 0; 
    for (int i = 0; i < g_chromosomeLength; i++){ //chromsomeLength is the same as the number of coefficients in the function
        tempSum += (chrom.genes[i] * g_coefficients[i]); //this will not work without coefficients
    }
    //chrom.fitness = ((double)1.0 / ( 1.0 + abs(tempSum - g_limit))); // the absolute value of the difference between tempSum and the target = fitness.
    chrom.fitness = abs(tempSum - g_limit);
    //printf("My fitness: %f\n", chrom.fitness);
   // double percentInaccurate = abs(tempSum - g_limit) / g_limit;
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
    //srand(0);
    int totalGenes;
    vector<int> chosenGenes;
    //Check if vector is empty
    if (!chromosomeVector.empty())
    {
        //Calculate total number of genes
        totalGenes = g_chromosomeLength * g_numChromosomes;
        int numMutations = g_mutationRate * totalGenes;
        //printf("Genes: %d\n", totalGenes);
        //printf("Num mutations expected: %d\n", numMutations);
        //Keep generating numbers until there is the expected amount of mutations
        while (chosenGenes.size() < numMutations)
        {
            int genes = rand() % totalGenes;
            // int genes = rand() % totalGenes + 1; // old version
            chosenGenes.push_back(genes);
        }

        //Change the chosen genes inside the chromosomes to a random number between 0 and target (g_limit)
        //cout<<"past while loop"<<endl;
        //#pragma omp parallel for
        for (int i = 0; i < chosenGenes.size(); i++)
        {   
            int indexOfChromosome = chosenGenes[i] / g_chromosomeLength; // Get the index of the chrosomome.
            //printf("Index of chromosome: %d\n", indexOfChromosome);
            // IF we have 6 chrosomsomes of length 4 and we want the 17th gene, 17 / 4 = 4
            int indexOfGene = chosenGenes[i] % g_chromosomeLength;
            //printf("Gene index: %d\n", indexOfGene);
            // If we have 6 chromosomes of length 4 and we need to know the index of the exact gene, 17 % 4 = 1
            // 4 * 4 + 1 =17
            chromosomeVector[indexOfChromosome].genes[indexOfGene] = (rand() % g_limit)+1; // set that targeted gene to a number between 0 and g_limit 
            
        }
        //cout<<"past for loop"<<endl;
    }
    else 
        return;     //Vector is empty. No chromosomes to mutate
}
//function to create a population of chromosomes

// GOOD CANDIDATE FOR PARALLELIZATION
void createPop(vector<Chromosome> &population, int populationSize,int numGenes){//works as intended
    #pragma omp parallel for private(population)
        for(int i=0;i<populationSize;i++){
            Chromosome c(numGenes);
            population.push_back(c);
        }
}

//Takes two chromosomes and crosses their genes at a given crossover index in the gene list
void crossoverGenes(vector<Chromosome> &population)
{
    srand(time(NULL));
    //srand(0);
    vector<int> parents;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<float> distribution(0.0, 1.0);

    /*
    Select the parents by generating a random number between 0 and 1.
    If that number is smaller than the crossover rate, then add the index of the chromosome to the vector of parents.
    */
    for (int i = 0; i < population.size(); i++)
    {
        float randCutNum = distribution(gen);
        if (randCutNum < g_crossoverRate)
        {
            parents.push_back(i);
        }
    }

    //Crossross over all the parents with each other from the cut point onwards using the parent indexes
    for (int i = 0; i < parents.size(); i++)
    {
        int cutPoint = rand() % (population[0].genes.size() - 1);
        for (int j = i + 1; j < parents.size(); j++)
        {
            //#pragma omp parallel for
            for (int k = cutPoint; k < population[0].genes.size(); k++)
            {
                int temp = population[parents[i]].genes[k];
                population[parents[i]].genes[k] = population[parents[j]].genes[k];
                population[parents[j]].genes[k] = temp;
            }
        }
    }
}

void generateGenes(vector<Chromosome> &population){//works as intended
    srand(time(NULL));
    //srand(0);
    #pragma omp parallel for
        for(int i = 0; i< population.size(); i++){
            for(int j = 0; j < population[i].genes.size();j++){
                int randInt = rand() % (g_limit + 1);
                population[i].genes[j]=randInt;
            }
        }
}

void selection(vector<Chromosome> &population,int numGenes){
    // Take population, compute fitness of each chromosome,
    //then calculate the probably of each chromosome to be selected
    //Then get the cumulative probability
    srand((unsigned)time(NULL));
    //srand(0);
    double total = 0.0;
    vector<double> cumulativeProb;
    vector<Chromosome> pop2;
    //computing the fitness
    for(int i = 0;i<population.size();i++){
        population[i].fitness = 1 / (1+population[i].fitness);
        total += population[i].fitness;
    }
    //probability of each chromosome
    for(int i=0;i<population.size();i++){
        population[i].probability = population[i].fitness / total;
    }
    //cumulative probability
    for(int i=0;i<population.size();i++){
        cumulativeProb.push_back(0.0);
        for(int j=0; j<=i; j++){
            cumulativeProb[i] += population[j].probability; 
        }
        //cout<<"cumuProb at i = "<<cumulativeProb[i]<<endl;
    }

    //selection of the chromosomes
        //looping through the population to determine which chromosome we will copy
    for(int i = 0;i<population.size();i++){
        double randD = ((double)rand() / (RAND_MAX));
        int indexOfSelectedChromosome = selectionHelper(cumulativeProb, randD);
     //   printf("Selected chromosome: %d\n", indexOfSelectedChromosome);
        Chromosome c(population[indexOfSelectedChromosome].genes.size());
        c.copy(population[indexOfSelectedChromosome]);
        pop2.push_back(c);
    }
    for(int i=0;i<population.size();i++){
        population[i]=pop2[i];
    }
    pop2.clear();
}


int selectionHelper(vector<double> cumulativeProb, double randomNum){
  double diff = 99.9;
  int indexOfClosest = 0;
  for(int i = 0; i < cumulativeProb.size(); i++){

    double tempDiff = abs(cumulativeProb[i]-randomNum);
    if(tempDiff < diff){
      diff = tempDiff;
      indexOfClosest = i;
    }
  }
  return indexOfClosest;
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
