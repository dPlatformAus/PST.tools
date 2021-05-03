#define MAX_FACTORIAL 10000

std::vector <double>  log_factorial;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void initialise_log_factorial(){
    unsigned i;
    log_factorial.push_back(0);
    for (i=1; i<=MAX_FACTORIAL; i++) log_factorial.push_back(log_factorial[i-1]+log(i));
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double factorial(int n) { //this is the slow way! use the datastructure log_factorial defined in factorials.h
    return exp(log_factorial[n]); //if n>170 this will return INF as the result will be too large for a double!
}

