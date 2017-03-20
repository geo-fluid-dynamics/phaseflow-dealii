#ifndef _peclet_step_time_h_
#define _peclet_step_time_h_


/*!
@brief Solve the linear system.

@author Alexander G. Zimmerman 2016 <zimmerman@aices.rwth-aachen.de>
*/
template<int dim>
SolverStatus Peclet<dim>::solve_linear_system(bool quiet)
{
    double tolerance = this->params.solver.tolerance;
    if (this->params.solver.normalize_tolerance)
    {
        tolerance *= this->system_rhs.l2_norm();
    }
    SolverControl solver_control(
        this->params.solver.max_iterations,
        tolerance);
       
    SolverGMRES<> solver_gmres(solver_control);

    PreconditionSSOR<> preconditioner; // @todo: What is the SSOR preconditioner? Probably not appropriate for Newton linearized Navier-Stokes-Boussinesq.
    
    preconditioner.initialize(this->system_matrix, 1.0);

    std::string solver_name;
    
    if (this->params.solver.method == "GMRES")
    {
        solver_name = "GMRES";
        solver_gmres.solve(
            this->system_matrix,
            this->newton_solution,
            this->system_rhs,
            preconditioner);    
    }
    else
    {
        Assert(false, ExcNotImplemented());
    }

    this->constraints.distribute(this->solution);

    if (!quiet)
    {
        std::cout << "     " << solver_control.last_step()
              << " " << solver_name << " iterations." << std::endl;
    }
    
    SolverStatus status;
    status.last_step = solver_control.last_step();
    
    return status;

}

/*!
@brief Step the simulation from the current time step to the next time step.

@detail

    This requires iterating through each Newton substep of the timestep, 
    assembling and solving a linear system for each substep, until convergence.

@author Alexander G. Zimmerman 2017 <zimmerman@aices.rwth-aachen.de>
*/
template<int dim>
void Peclet<dim>::step_time(bool quiet)
{
    if (!quiet & output_this_step)
    {
        std::cout << "Time step " << this->time_step_counter 
            << " at t=" << this->time << std::endl;    
    }

    bool converged = false;
    double residual = 1.e32; // Initialize to an arbitrarily large number

    unsigned int i;

    this->old_solution = this->solution;

    Vector<double> diff(this->solution.size());

    for (i = 0; i < MAX_NEWTON_ITERATIONS; ++i)
    {
        this->assemble_system();

        this->apply_boundary_values_and_constraints();

        this->old_newton_solution = this->newton_solution;

        this->solve_linear_system();

        diff = this->newton_solution;
        diff -= this->old_newton_solution;

        residual = diff.l2_norm();

        if (residual < NEWTON_TOLERANCE)
        {
            converged = true;
            this->solution = this->newton_solution;
            break;
        }

    }

    assert(converged);

    if (!quiet)
    {
        std::cout << "Newton method converged after " << i + 1 << " iterations." << std::endl;
    }
    
    bool steady = false;
    diff = this->solution;
    diff -= this->old_solution;
    const double l2_norm_diff = diff.l2_norm();
    if ((l2_norm_diff < NEWTON_TOLERANCE) & 
        (l2_norm_diff < this->params.solver.tolerance))
    {
        steady = true;
    }

    if (this->params.time.stop_when_steady & steady)
    {
        if (!quiet)
        {
            std::cout << "Reached steady state at t = " << this->time << std::endl;
        }
        
        this->final_time_step = true;
        this->output_this_step = true;
    }
        
    if (this->output_this_step)
    {
        this->write_solution();
        
        if (this->params.verification.enabled)
        {
            this->append_verification_table();
        }
        
    }

}

#endif
