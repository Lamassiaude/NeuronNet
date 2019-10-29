#include "random.h"
#include <algorithm>

// constructeur qui est initialisé a une seed
/*! @name Initializing
 The generator \ref rng is a Mersenne twister *mt19937* engine.
 */

/* A seed *s>0* can be provided, by default it is seeded with a *random_device*.
 */

/*
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
*/


RandomNumbers::RandomNumbers(unsigned long int s)
:seed(s)
{
    if (seed == 0) {
        std::random_device rd;
        seed = rd();
    }
    rng = std::mt19937(seed);

}




void RandomNumbers::uniform_double(std::vector<double>& v, double lower, double upper)
{
    std::uniform_real_distribution<> unif(lower, upper);
    for (auto I = v.begin(); I != v.end(); I++) *I = unif(rng);
   
}

double RandomNumbers::uniform_double(double lower, double upper)
{
    std::uniform_real_distribution<> dis(lower, upper);
    return (dis(rng));
}


void RandomNumbers::normal(std::vector<double>& v, double mean, double sd)
{
    std::normal_distribution<> normal(mean, sd);
    for (auto I = v.begin(); I != v.end(); I++) *I = normal(rng);
}

double RandomNumbers::normal(double mean, double sd)
{
    std::normal_distribution<> normal(mean, sd);
    return normal(rng);
}


void RandomNumbers::poisson(std::vector<int>& v, double mean)
{

    std::poisson_distribution<> poisson(mean);
    for (auto I = v.begin(); I != v.end(); I++) *I = poisson(rng);
}


int RandomNumbers::poisson(double mean)
{
    std::poisson_distribution<> poisson(mean);
    return poisson(rng);
}



