#ifndef _peclet_step_time_h_
#define _peclet_step_time_h_


/*!
@brief Solve the linear system.

@author Alexander G. Zimmerman 2016 <zimmerman@aices.rwth-aachen.de>
*/
template<int dim>
SolverStatus Peclet<dim>::solve_time_step(bool quiet)
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
            this->solution,
            this->system_rhs,
            preconditioner);    
    }
    else
    {
        Assert(false, NotImplemented());
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
    if (output_this_step)
    {
        std::cout << "Time step " << this->time_step_counter 
            << " at t=" << this->time << std::endl;    
    }

    bool converged = false;
    double residual = 1.e32; // Initialize to an arbitrarily large number

    for (unsigned int i = 0; i < max_newton_iterations; ++i)
    {
        this->assemble_system();

        this->apply_boundary_conditions_and_constraints();

        this->solve_linear_system();

        residual = (this->solution - this->old_solution).l2_norm();

        if (residual < tolerance)
        {
            converged = true;
            break;
        }

    }

    assert(converged);

    std::cout << "Newton method converged after " << i + 1 << " iterations." << std::endl;

    this->solver_status = this->solve_time_step(!output_this_step);
    
    bool steady = false;
    if (solver_status.last_step == 0)
    {
        steady = true;
    }

    if (this->params.time.stop_when_steady & steady)
    {
        std::cout << "Reached steady state at t = " << this->time << std::endl;
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
