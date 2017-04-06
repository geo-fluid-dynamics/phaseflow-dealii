#ifndef _peclet_step_time_h_
#define _peclet_step_time_h_


/*!
@brief Solve the linear system.

@author Alexander G. Zimmerman 2016 <zimmerman@aices.rwth-aachen.de>
*/
template<int dim>
SolverStatus Peclet<dim>::solve_linear_system(bool quiet)
{
    if (WRITE_LINEAR_SYSTEM)
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

/*! Setup and solve a Newton iteration */
template<int dim>
void Peclet<dim>::step_newton()
{
    this->old_newton_solution = this->solution;
    
    this->assemble_system();

    this->apply_boundary_values_and_constraints();

    this->solve_linear_system();

    this->solution -= this->newton_residual;
}

/*! Iterate the Newton method to solve the nonlinear problem */
template<int dim>
void Peclet<dim>::solve_nonlinear_problem(bool quiet)
{
    bool converged = false;

    unsigned int i;

    for (i = 0; i < this->params.nonlinear_solver.max_iterations; ++i)
    {
        this->step_newton();
        
        Output::write_solution_to_vtk( // @todo Debugging
            "newton_residual.vtk",
            this->dof_handler,
            this->newton_residual);

        Output::write_solution_to_vtk( // @todo Debugging
            "newton_solution.vtk",
            this->dof_handler,
            this->solution);
            
        double norm_residual = this->newton_residual.l2_norm();
        
        if (!quiet)
        {
            std::cout << "Newton iteration L2 norm residual = " << norm_residual << std::endl;
        }
        
        if (norm_residual < this->params.nonlinear_solver.tolerance)
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
    
}

#endif
