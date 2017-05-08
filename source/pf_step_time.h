#ifndef _pf_step_time_h_
#define _pf_step_time_h_

template<int dim>
void Phaseflow<dim>::set_time_step_size(const double _new_size)
{
    double new_size = _new_size;
    
    if (new_size < this->params.time.min_step_size)
    {
        new_size = this->params.time.min_step_size;
    }
    else if(new_size > this->params.time.max_step_size)
    {
        new_size = this->params.time.max_step_size;
    }
    
    if ((this->time + new_size) > this->params.time.end)
    {
        new_size = this->params.time.end - this->time;
    }
    
    if (numbers::NumberTraits<double>::abs(new_size - this->time_step_size) > EPSILON)
    {
        std::cout << "Set time step to deltat = " << new_size << std::endl;
    }
    
    this->time_step_size = new_size;
    
}


/*! Step the simulation from the current time step to the next time step.

This requires iterating through each Newton substep of the timestep, 
assembling and solving a linear system for each substep, until convergence.

*/
template <int dim>
void Phaseflow<dim>::step_time()
{   
    this->old_solution = this->solution;
    
    bool converged;

    do 
    {
        this->new_time = this->time + this->time_step_size;
        
        if (this->params.physics.prescribe_convection_velocity)
        {
            this->prescribed_convection_velocity_function.set_time(this->new_time);
            
            VectorTools::interpolate(
                this->dof_handler,
                this->prescribed_convection_velocity_function,
                this->solution,
                this->velocity_mask);
        }
        
        converged = this->solve_nonlinear_problem();
        
        if (!converged)
        {
            this->set_time_step_size(this->time_step_size/TIME_GROWTH_RATE);
        }
        
    } while (!converged);

    this->time = this->new_time;
    
    std::cout << "Reached time t = " << this->time << std::endl;
    
    if (this->time >= (this->params.time.end - EPSILON))
    {
        return;
    }
    
    if (converged)
    {   
        this->set_time_step_size(TIME_GROWTH_RATE*this->time_step_size);        
    }

}

#endif
