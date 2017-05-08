#ifndef _pf_system_h_
#define _pf_system_h_

/*! Setup the linear system objects. */
template<int dim>
void Phaseflow<dim>::setup_system()
{
    
    this->dof_handler.distribute_dofs(this->fe);

    DoFRenumbering::component_wise(this->dof_handler);

    std::cout << std::endl
            << "==========================================="
            << std::endl
            << "Number of active cells: " << this->triangulation.n_active_cells()
            << std::endl
            << "Number of degrees of freedom: " << this->dof_handler.n_dofs()
            << std::endl
            << std::endl;
            
    this->constraints.clear();

    DoFTools::make_hanging_node_constraints(
        this->dof_handler,
        this->constraints);
        
    this->constraints.close();

    DynamicSparsityPattern dsp(this->dof_handler.n_dofs());

    DoFTools::make_sparsity_pattern(
        this->dof_handler,
        dsp,
        this->constraints,
        /*keep_constrained_dofs = */ true);
        
    this->sparsity_pattern.copy_from(dsp);

    this->system_matrix.reinit(this->sparsity_pattern);

    this->solution.reinit(this->dof_handler.n_dofs());

    this->newton_residual.reinit(this->dof_handler.n_dofs());
    
    this->newton_solution.reinit(this->dof_handler.n_dofs());
    
    this->old_solution.reinit(this->dof_handler.n_dofs());
    
    this->old_newton_solution.reinit(this->dof_handler.n_dofs());
    
    this->system_rhs.reinit(this->dof_handler.n_dofs());

}

/*!
 @brief Assemble the system for the Newton linearized Navier-Stokes-Boussinesq equtaions.
 
 @detail
 
    This implements equation (17) from Danaila et al. 2014.

     This is the bouyancy force function from Danaila 2014,

        $f_B(\theta) = theta Ra / (Pr Re^2)$

    in terms of the Reynolds, Prandtl, and Rayleigh numbers, defined as:
        
        $Re = \rho_{ref} V_{ref} L_ref/\mu_l,
            Pr = \nu_l / \alpha_l,
            Ra = g \beta L_{ref}^3 (T_h - T_c) / (\nu_l \alpha_l)$

    Danaila's paper assumes gravity is pointing uniformly downward and treats
    gravity as a scalar. This implementation is generalized to a dim-dimensional
    gravity vector. This requires that we interpret the Rayleigh number also as
    vector-valued.

    This models the bouyancy force as only a function of $\theta$.

    Following the implementation approach in
        
        http://dealii.org/8.4.1/doxygen/deal.II/group__vector__valued.html
    
    and 
    
        http://dealii.org/8.4.1/doxygen/deal.II/step_20.html#Assemblingthelinearsystem
 
 @author Alexander Zimmerman 2016
*/
template<int dim>
void Phaseflow<dim>::assemble_system()
{
    this->system_matrix = 0.;
    
    this->system_rhs = 0.;
    
    /*!
     Local parameters
    */
    const double PENALTY = 1.e-7; // @todo: Expose this to ParameterHandler.

    const double
        Ra = RAYLEIGH_NUMBER,
        Pr = PRANDTL_NUMBER,
        Re = REYNOLDS_NUMBER;
    
    const double K = SOLID_CONDUCTIVITY/LIQUID_CONDUCTIVITY;

    Tensor<1, dim> g; // @todo: Make this const.
    for (unsigned int i = 0; i < dim; ++i)
    {
        g[i] = this->params.physics.gravity[i];
    }

    const double mu_l = this->params.physics.liquid_dynamic_viscosity;

    /*!
     lambda function for classical (linear) Boussinesq bouyancy
    */
    auto f_B = [Ra, Pr, Re, g](const double _theta) 
    {
        return _theta*Ra/(Pr*Re*Re)*g;
    };

    /*!
     Analytical derivative of classical (linear) Boussinesq bouyancy
    */
    const Tensor<1, dim> df_B_over_dtheta(Ra/(Pr*Re*Re)*g);

    /*!
     lambda functions for linear, bilinear, and trilinear operator
    */
    auto a = [](
        const double _mu,
        const Tensor<2, dim> _gradu,
        const Tensor<2, dim> _gradv)
    {
        auto D = [](
            const Tensor<2, dim> _gradw)
        {
            return 0.5*(_gradw + transpose(_gradw));
        };

        return 2.*_mu*scalar_product(D(_gradu), D(_gradv));
    };

    auto b = [](
        const double _divu,
        const double _q)
    {
        return -_divu*_q;
    };

    auto c = [](
        const Tensor<1, dim> _w,
        const Tensor<2, dim> _gradz,
        const Tensor<1, dim> _v)
    {
        return (_v*_gradz)*_w;
    };

    /*!
     Organize data
    */

    QGauss<dim> quadrature_formula(SCALAR_DEGREE + 2);

    FEValues<dim> fe_values(
        this->fe,
        quadrature_formula,
        update_values | update_gradients | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    
    const unsigned int n_quad_points = quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    
    Vector<double> local_rhs(dofs_per_cell);
    
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    
    std::vector<Tensor<1,dim>> old_velocity_values(n_quad_points);
    
    std::vector<double> old_pressure_values(n_quad_points);

    std::vector<double> old_temperature_values(n_quad_points);
    
    std::vector<Tensor<1,dim>> old_newton_velocity_values(n_quad_points);
    
    std::vector<double> old_newton_pressure_values(n_quad_points);
    
    std::vector<double> old_newton_temperature_values(n_quad_points);
    
    std::vector<Tensor<2,dim>> old_newton_velocity_gradients(n_quad_points);
    
    std::vector<Tensor<1,dim>> old_newton_temperature_gradients(n_quad_points);

    std::vector<double> old_newton_velocity_divergences(n_quad_points);
    
    std::vector<Vector<double>> source_values(n_quad_points, Vector<double>(dim + 2));
    
    this->source_function.set_time(this->new_time);
    
    typename DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        endc = this->dof_handler.end();

    
    /*!
        Set local variables to match notation in Danaila 2014
    */
    const double deltat = this->time_step_size;
    
    const double gamma = PENALTY;
    
    for (; cell != endc; ++cell) /*! Assemble element-wise */
    {
        fe_values.reinit(cell);

        fe_values[this->velocity_extractor].get_function_values(
            this->old_solution,
            old_velocity_values);

        fe_values[this->pressure_extractor].get_function_values(
            this->old_solution,
            old_pressure_values);
        
        fe_values[this->temperature_extractor].get_function_values(
            this->old_solution,
            old_temperature_values);
        
        fe_values[this->velocity_extractor].get_function_values(
            this->old_newton_solution,
            old_newton_velocity_values);

        fe_values[this->pressure_extractor].get_function_values(
            this->old_newton_solution,
            old_newton_pressure_values);

        fe_values[this->temperature_extractor].get_function_values(
            this->old_newton_solution,
            old_newton_temperature_values);

        fe_values[this->temperature_extractor].get_function_gradients(
            this->old_newton_solution,
            old_newton_temperature_gradients);

        fe_values[this->velocity_extractor].get_function_gradients(
            this->old_newton_solution,
            old_newton_velocity_gradients);

        fe_values[this->velocity_extractor].get_function_divergences(
            this->old_newton_solution,
            old_newton_velocity_divergences);
                    
        std::vector<Tensor<1, dim>> velocity_fe_values(dofs_per_cell);       
        std::vector<double> pressure_fe_values(dofs_per_cell);
        std::vector<double> temperature_fe_values(dofs_per_cell);
        std::vector<Tensor<1, dim>> grad_temperature_fe_values(dofs_per_cell);
        std::vector<Tensor<2, dim>> grad_velocity_fe_values(dofs_per_cell);
        std::vector<double> div_velocity_fe_values(dofs_per_cell);
        
        local_matrix = 0.;
        
        local_rhs = 0.;
        
        this->source_function.vector_value_list(
            fe_values.get_quadrature_points(),
            source_values);
        
        for (unsigned int quad = 0; quad< n_quad_points; ++quad)
        {
            /* Name local variables to match notation in Danaila 2014 */
            const Tensor<1, dim>  u_n = old_velocity_values[quad];
            const double theta_n = old_temperature_values[quad];
            
            const Tensor<1, dim> u_k = old_newton_velocity_values[quad];
            const double p_k = old_newton_pressure_values[quad];
            const double theta_k = old_newton_temperature_values[quad];
            
            const Tensor<1, dim> gradtheta_k = old_newton_temperature_gradients[quad];
            const Tensor<2, dim> gradu_k = old_newton_velocity_gradients[quad];
            const double divu_k = old_newton_velocity_divergences[quad];
            
            Tensor<1, dim> s_u;
            
            for (unsigned int d = 0; d < dim; ++d)
            {
                s_u[d] = source_values[quad][d];
            }
             
            const double s_p = source_values[quad][dim];
            
            const double s_theta = source_values[quad][dim+1];
            
            for (unsigned int dof = 0; dof < dofs_per_cell; ++dof)
            {
                velocity_fe_values[dof] = fe_values[this->velocity_extractor].value(dof, quad);
                pressure_fe_values[dof] = fe_values[this->pressure_extractor].value(dof, quad);
                temperature_fe_values[dof] = fe_values[this->temperature_extractor].value(dof, quad);
                grad_temperature_fe_values[dof] = fe_values[this->temperature_extractor].gradient(dof, quad);
                grad_velocity_fe_values[dof] = fe_values[this->velocity_extractor].gradient(dof, quad);
                div_velocity_fe_values[dof] = fe_values[this->velocity_extractor].divergence(dof, quad);
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                /* Name local variables to match notation in Danaila 2014 */
                const Tensor<1, dim> v = velocity_fe_values[i];
                const double q = pressure_fe_values[i];
                const double phi = temperature_fe_values[i];
                const Tensor<1, dim> gradphi = grad_temperature_fe_values[i];
                const Tensor<2, dim> gradv = grad_velocity_fe_values[i];
                const double divv = div_velocity_fe_values[i];
                
                /* @todo Here I implemented the form derived in danaila2014newton, where they
                multiplied by test functions from the right. deal.II's tutorials recommend to get
                in the habit of instead multiplying from the left, to avoid a common class of errors.
                If verification fails, then I should try deriving my own form, with the left multiplication, and see if this helps.
                */
                for (unsigned int j = 0; j< dofs_per_cell; ++j)
                {
                    const Tensor<1, dim> u_w = velocity_fe_values[j];
                    const double p_w = pressure_fe_values[j];
                    const double theta_w = temperature_fe_values[j];
                    const Tensor<1, dim> gradtheta_w = grad_temperature_fe_values[j];
                    const Tensor<2, dim> gradu_w = grad_velocity_fe_values[j];
                    const double divu_w = div_velocity_fe_values[j];

                    local_matrix(i,j) += (
                        b(divu_w, q) - gamma*p_w*q // Mass
                        + scalar_product(u_w, v)/deltat + c(u_w, gradu_k, v) + c(u_k, gradu_w, v) + a(mu_l, gradu_w, gradv) + b(divv, p_w) // Momentum: Incompressible Navier-Stokes
                        + scalar_product(theta_w*df_B_over_dtheta, v) // Momentum: Bouyancy (Classical linear Boussinesq approximation)
                        + theta_w*phi/deltat - scalar_product(u_k, gradphi)*theta_w - scalar_product(u_w, gradphi)*theta_k + scalar_product(K/Pr*gradtheta_w, gradphi) // Energy
                        )*fe_values.JxW(quad); /* Map to the reference element */                        

                }
                
                local_rhs(i) += (
                        b(divu_k, q) - gamma*p_k*q // Mass
                        + scalar_product(u_k - u_n, v)/deltat + c(u_k, gradu_k, v) + a(mu_l, gradu_k, gradv) + b(divv, p_k) // Momentum: Incompressible Navier-Stokes
                        + scalar_product(f_B(theta_k), v) // Momentum: Bouyancy (Classical linear Boussinesq approximation)
                        + (theta_k - theta_n)*phi/deltat - scalar_product(u_k, gradphi)*theta_k + scalar_product(K/Pr*gradtheta_k, gradphi) // Energy
                        + s_p*q + scalar_product(s_u, v) + s_theta*phi // Source (MMS)
                        )*fe_values.JxW(quad); 

            }            
            
        }
            
        // Export local contributions to the global system
        cell->get_dof_indices(local_dof_indices);
        
        this->constraints.distribute_local_to_global(
            local_matrix, local_rhs, local_dof_indices,
            this->system_matrix, this->system_rhs);

    }

}

template<int dim>
void Phaseflow<dim>::interpolate_boundary_values(
    Function<dim>* function,
    std::map<types::global_dof_index, double> &boundary_values) const
{    
    for (unsigned int ib = 0; ib < this->params.boundary_conditions.strong_boundaries.size(); ++ib) /* For each boundary */
    {    
        unsigned int b = this->params.boundary_conditions.strong_boundaries[ib];
        
        auto mask = this->params.boundary_conditions.strong_masks[b];
        
        for (auto field_name : FIELD_NAMES) /* For each field variable */
        {
            if (std::find(mask.begin(), mask.end(), field_name) == mask.end()) /* Skip if the field name is not in the mask */
            {
                continue;
            }
            
            /* @todo                
            Is there some way to contain the extractors (or pointers to them) in a single object that can be indexed?
            Neither std::vector<void*> nor tuple (because the tuple could not be indexed with a variable) worked, and I'm out of ideas. */
            if (field_name == "velocity")
            {
                VectorTools::interpolate_boundary_values(
                    this->dof_handler, b, *function, boundary_values,
                    this->fe.component_mask(this->velocity_extractor));
            }
            else if (field_name == "pressure")
            {
                VectorTools::interpolate_boundary_values(
                    this->dof_handler, b, *function, boundary_values,
                    this->fe.component_mask(this->pressure_extractor));
            }
            else if (field_name == "temperature")
            {
                VectorTools::interpolate_boundary_values(
                    this->dof_handler, b, *function, boundary_values,
                    this->fe.component_mask(this->temperature_extractor));
            }
            else
            {
                assert(false);
            }
                        
        }

    }
    
}

/*! Apply the boundary conditions (strong and natural) and apply constraints (including those for hanging nodes */
template<int dim>
void Phaseflow<dim>::apply_boundary_values_and_constraints()
{       
    /* Since we are applying boundary conditions to the Newton linearized system
    to compute a residual, we want to apply the boundary conditions residual, rather
    than the user supplied boundary conditions.
    
    To do this, we evaluate the BC's both at the new time and the current time,
    and we apply the difference. */
    std::map<types::global_dof_index, double> residual_boundary_values, boundary_values, new_boundary_values;
    
    /* @todo Using a FEFieldFunction to interpolate the solution values seems like a terrible idea, 
    since FEFieldFunction is designed to interpolate within the domain. Is there another method when I really just
    need the solution values on boundary nodes?
    
    Another reason this is a terrible approach: I am essentially interpolating all of the finite element functions
    to get velocity, pressure, and temperature values, and then these have to be decomposed onto the finite element functions again with interpolate_boundary_values. There should be an easy way to just use the map<global_dof_index, double> to do this directly. */
    Functions::FEFieldFunction<dim> solution_field_function(this->dof_handler, this->solution);
    
    this->interpolate_boundary_values(&solution_field_function, boundary_values);
    
    this->boundary_function.set_time(this->new_time);
        
    this->interpolate_boundary_values(&this->boundary_function, new_boundary_values);
    
    for (auto m: new_boundary_values)
    {
        residual_boundary_values.insert(m);
        
        residual_boundary_values[m.first] -= boundary_values[m.first];
    }

    MatrixTools::apply_boundary_values(
        residual_boundary_values,
        this->system_matrix,
        this->newton_residual,
        this->system_rhs);
}


/*!
@brief Solve the linear system.

@author Alexander G. Zimmerman 2016 <zimmerman@aices.rwth-aachen.de>
*/
template<int dim>
void Phaseflow<dim>::solve_linear_system()
{
    if (WRITE_LINEAR_SYSTEM)
    {
        Output::write_linear_system(this->system_matrix, this->system_rhs);
    }
    
    SparseDirectUMFPACK A_inv;
    
    A_inv.initialize(this->system_matrix);
    
    A_inv.vmult(this->newton_residual, this->system_rhs);

    this->constraints.distribute(this->newton_residual);

    std::cout << "Solved linear system" << std::endl;

}

#endif
