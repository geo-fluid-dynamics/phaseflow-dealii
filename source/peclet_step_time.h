#ifndef _peclet_step_time_h_
#define _peclet_step_time_h_

template<int dim>
void Peclet<dim>::set_time_step_size(double new_size)
{
    
    if (new_size < this->params.time.min_step_size)
    {
        new_size = this->params.time.min_step_size;
    }
    else if(new_size > this->params.time.max_step_size)
    {
        new_size = this->params.time.max_step_size;
    }
    
    if ((new_size != this->time_step_size))
    {
        std::cout << "Set time step to deltat = " << this->time_step_size << std::endl;
    }
    
    this->time_step_size = new_size;
    
}


/*! Step the simulation from the current time step to the next time step.

This requires iterating through each Newton substep of the timestep, 
assembling and solving a linear system for each substep, until convergence.

*/
template <int dim>
void Peclet<dim>::step_time()
{   
    if (this->time > this->params.time.end)
    {
        this->time_step_size = this->params.time.end - this->time;
        
        this->time = this->params.time.end;
    }
    
    this->old_solution = this->solution;
    
    bool converged;

    do 
    {
        converged = this->solve_nonlinear_problem();
        
        if (!converged)
        {
            this->set_time_step_size(this->time_step_size/TIME_GROWTH_RATE);
        }
        
    } while (!converged);

    this->time += this->time_step_size;
    
    std::cout << "Reached time t = " << this->time + this->time_step_size << std::endl;
    
    if (converged)
    {   
        this->set_time_step_size(TIME_GROWTH_RATE*this->time_step_size);        
    }

}

#endif
