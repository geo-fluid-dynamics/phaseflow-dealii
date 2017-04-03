#ifndef _peclet_step_time_h_
#define _peclet_step_time_h_


/*!
@brief Solve the linear system.

@author Alexander G. Zimmerman 2016 <zimmerman@aices.rwth-aachen.de>
*/
template<int dim>
SolverStatus Peclet<dim>::solve_linear_system(bool quiet)
{
    if (write_linear_system)
    {
        Output::write_linear_system(this->system_matrix, this->system_rhs);
    }
    
    double tolerance = this->params.linear_solver.tolerance;
    if (this->params.linear_solver.normalize_tolerance)
    {
        double norm = this->system_rhs.l2_norm(); 
        tolerance *= norm;
    }
    SolverControl solver_control(
        this->params.linear_solver.max_iterations,
        tolerance);
       
    SolverGMRES<> solver_gmres(solver_control);

    SparseDirectUMFPACK A_inv;
    
    PreconditionIdentity preconditioner;

    std::string solver_name;
    
    std::string solver_type;
    
    if (this->params.linear_solver.method == "LU")
    {
        solver_name = "LU";
        solver_type = "direct";
        A_inv.initialize(this->system_matrix);
        A_inv.vmult(this->newton_residual, this->system_rhs);
    }
    else if (this->params.linear_solver.method == "GMRES")
    {
        solver_name = "GMRES";
        solver_type = "iterative";
        solver_gmres.solve(
            this->system_matrix,
            this->newton_residual,
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
        if (solver_type == "iterative")\
        {
            std::cout << "     " << solver_control.last_step()
                << " " << solver_name << " iterations." << std::endl;
        }
        else if (solver_type == "direct")
        {
            std::cout << "Solved linear system" << std::endl;
        }
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

    unsigned int i;

    this->old_solution = this->solution;

    for (i = 0; i < this->params.nonlinear_solver.max_iterations; ++i)
    {
        this->old_newton_solution = this->solution;
        
        this->assemble_system();

        this->apply_boundary_values_and_constraints();

        this->solve_linear_system();

        this->solution -= this->newton_residual;

        if (this->newton_residual.l2_norm() < this->params.nonlinear_solver.tolerance)
        {
            converged = true;
            break;
        }

    }

    assert(converged);

    if (!quiet)
    {
        std::cout << "Newton method converged after " << i + 1 << " iterations." << std::endl;
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
