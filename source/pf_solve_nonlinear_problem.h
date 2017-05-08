#ifndef _pf_solve_nonlinear_problem_h_
#define _pf_solve_nonlinear_problem_h_


/*! Setup and solve a Newton iteration */
template<int dim>
void Phaseflow<dim>::step_newton()
{
    this->old_newton_solution = this->newton_solution;
    
    this->assemble_system();

    this->apply_boundary_values_and_constraints();

    this->solve_linear_system();

    this->newton_solution -= this->newton_residual;
}

/*! Iterate the Newton method to solve the nonlinear problem */
template<int dim>
bool Phaseflow<dim>::solve_nonlinear_problem()
{
    this->newton_solution = this->solution;
    
    bool converged = false;

    unsigned int i;
    
    double old_norm_residual = 1.e32;
    
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
            this->newton_solution);
        
        double norm_residual = this->newton_residual.l2_norm()/this->newton_solution.l2_norm();
        
        std::cout << "Newton iteration: L2 norm of relative residual, || w_w || / || w_k || = " << norm_residual << std::endl;
        
        if (norm_residual > old_norm_residual)
        {
            converged = false;
            
            std::cout << "Newton iteration diverged." << std::endl;
            
            Output::write_solution_to_vtk( // @todo Debugging
                "diverged_newton_solution.vtk",
                this->dof_handler,
                this->newton_solution);
                
            if (this->time_step_size == this->params.time.min_step_size)
            {
                assert(converged);
            }
            
            return converged;
        }
        
        old_norm_residual = norm_residual;
        
        if (norm_residual < this->params.nonlinear_solver.tolerance)
        {
            converged = true;
            break;
        }

    }

    if ((this->time_step_size > this->params.time.min_step_size) & !converged)
    {
        return converged;
    }
    
    assert(converged);

    std::cout << "Newton method converged after " << i + 1 << " iterations." << std::endl;
    
    this->solution = this->newton_solution;
    
    return converged;
}

#endif
