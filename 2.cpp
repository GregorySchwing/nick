int main(void){

        // Sets the number of replicas, based on the number of values the user gives for temperature
        // The only thing different about these parallel simulations are their temperatures
         int num_replicas = sim.replica_temps.size();
        
        // There has to be at least one core on your computer for each thread
        // Number of thread equals number of cores
        if (num_replicas > numThreads){
          std::cout << "Error: Not enough threads!\n";
          std::cout << "Use  a 1:1 ratio of replicas:threads.\n";
          exit(EXIT_FAILURE);
        }
               
        int i;

        // An array of parallel simulations
        Simulation* sim_re[num_replicas];
    

        // This is where the magic happens
        // Each thread will run an iteration of this for loop in parallel, instead of one after the other, aka serially
        // They each share the two shared variables , the array of simulations and number of replicas
        // They have private copies of the i, iteration variable,
        
        #pragma omp parallel for default(none) private(i) shared(num_replicas, sim_re)
            for (i = 0; i < num_replicas; i++) {
                sim_re[i] = new Simulation();
             // overloaded RunSim for RE
                sim_re[i]->RunSimulation();
            }
    }
    // GJS
}

void Simulation::RunSimulation()
{

  // An object to represent a replica, or indepenent simulation
  gmx_repl_ex_t     repl_ex = nullptr;

  // An object to represent the parts of the replica that can be swapped
  Replica_State * state_global = new Replica_State();

  int bDoReplEx, bExchanged;

  // A file to write to
  FILE *fplog = fopen(replica_log.c_str(), "a");  

  // This is how many "steps" or tries we want to run
  for (ulong step = 0; step < totalSteps; step++) {

    // The first MC in MCMC. Monte Carlo
    // It chooses a random move and runs it
    system->ChooseAndRunMove(step);

    // A boolean to calculate whether this is a step to try swapping on
    // The && are "logical ands" everything has to be true for bDoReplEx to be true
    bDoReplEx = (step > 0) && !bLastStep && (step % replExParams->exchangeInterval == 0);

    if (bDoReplEx) {
  // This makes all the paralell simulations wait for everyone to catch up to this point
  #pragma omp barrier          
   
        // Set the current state from the system
        GetSystem(state_global, system);

        // Try to exchange
        bExchanged = replica_exchange(repl_ex, state_global);
        
        // If we exchanged, bExchanged is a boolean, then overwrite the old copy
        // Of the state with the newly swapped one.
        if (bExchanged){
            // Swapping happens in this replica_states array.
            // under the hood during the replica_exchange method
            // If a swap took place the new state is in "my" slot in the array
            // repl_ex->repl is "my" index into an array of states
            state_global = replExParams->replica_states[repl_ex->repl];
            SetSystem(state_global, system);
        }
    }
}

void
replica_exchange(         FILE                          *fplog,
                          struct gmx_repl_ex            *re,
                          Replica_State*                state_global)
{
    /* standard nearest neighbor replica exchange */

    // This will be either 0 or 1
    // Just trust me
    m = (step / re->nst) % 2;
   
    // Ok nick so nrepl : number of replicas aka independent simulations with 
    // Slightly different variables.
    // We are only trying the neighboring ones, hence the ^ standard nearest neighbor comment
    // This for loop starts at index 1, but there is an index 0
    // a is set to the index -1, b is the index of the loop.
    // The first iteration will be with index 0 and 1, second 1 and 2, ect 
    for (i = 1; i < re->nrepl; i++)
    {
        a = re->ind[i-1];
        b = re->ind[i];

        // This is a modulo operator.  It is the remainder of a division.  
        // It lets us check if this is an even or odd call to this method
        // m will be either 0 or 1, so we only try this code on the evens or odds on each call
        
        if (i % 2 == m)
        {
        
        // This method is at the bottom of the file. 
        // This is the "closeness" calculation.
        // When using the card analogy, it would be
        // how close the denominations are
        delta = calc_delta(re, a, b);
            
            if (delta <= 0)
            {
                /* accepted */
                // Probability , ranging 0 to 1
                prob[i] = 1;
                // Boolean, 1 = true, 0 = false
                // bExchanged
                bEx[i]  = 1;
            }
            else
            {
                if (delta > PROBABILITYCUTOFF)
                {
                    prob[i] = 0;
                }
                else
                {
                    prob[i] = exp(-delta);
                }
                // roll a number to determine if accepted. For now it is superfluous to
                // reset, but just in case we ever add more calls in different branches
                // it is safer to always reset the distribution.
                

                // The second MC in MCMC
                // The metropolis criterion
                float randVar = ((float) rand() / (RAND_MAX)) + 1;
                bEx[i] = (int)(randVar < prob[i]);
            }
            
            if (bEx[i])
            {
                /* swap these two */
                tmp       = pind[i-1];
                pind[i-1] = pind[i];
                pind[i]   = tmp;
                re->nexchange[i]++;  /* statistics for back compatibility */
            }
        }
        else
        {
            prob[i] = -1;
            bEx[i]  = 0;
        }
    }
    
    re->nattempt[m]++;
}

static float calc_delta(gmx_repl_ex *re, int a, int b)
{
    float   ediff, delta = 0;
    
    float  *Epot = re->Epot;
  
    float  *beta = re->beta;
   
    // Epot is the potential energy of each replica
    // Think of the potential being the card denomination 
    // Ediff is the difference in denominations
    ediff = Epot[b] - Epot[a];
    
    // The calculation requires multiplying by a factor of the simulation conditions
    // Beta = 1 / temperature
  
    delta = -(beta[b] - beta[a])*ediff;

    return delta;
}

