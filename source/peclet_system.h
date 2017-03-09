#ifndef _peclet_system_h_
#define _peclet_system_h_

/*!
 @brief Setup the linear system objects.
 
 @author Alexander Zimmerman 2016
 */
template<int dim>
void Peclet<dim>::setup_system(bool quiet)
{
    
    this->dof_handler.distribute_dofs(this->e);

    DoFRenumbering::component_wise(this->dof_handler);

    if (!quiet)
    {
        std::cout << std::endl
                << "==========================================="
                << std::endl
                << "Number of active cells: " << this->triangulation.n_active_cells()
                << std::endl
                << "Number of degrees of freedom: " << this->dof_handler.n_dofs()
                << std::endl
                << std::endl;    
    }

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

    this->mass_matrix.reinit(this->sparsity_pattern);

    this->convection_diffusion_matrix.reinit(this->sparsity_pattern);

    this->system_matrix.reinit(this->sparsity_pattern);

    /*!
        @todo
         Why must we reinit the solution and RHS vectors here?
         This is an artifact from the step-26 tutorial.
    */
    this->solution.reinit(this->dof_handler.n_dofs());

    this->old_solution.reinit(this->dof_handler.n_dofs());

    this->newton_solution.reinit(this->dof_handler.n_dofs());

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
void Peclet<dim>::assemble_system(bool quiet)
{
    /*!
     Local parameters
     */
    const double penalty = 1.e-7; // @todo: Expose this to ParameterHandler.

    const Tensor<1, dim> Ra_over_PrRe2 = this->Ra/(this->Pr*this->Re*this->Re));

    /*!
     lambda function for classical (linear) Boussinesq bouyancy
     */
    Tensor<1, dim> f_B = [](const double _theta) 
    {
        return _theta*Ra_over_PrRe2;
    }

    /*!
     Analytical derivative of classical (linear) Boussinesq bouyancy
     */
    const Tensor<1, dim> df_B_over_dtheta = Ra_over_PrRe2;

    /*!
     lambda functions for linear, bilinear, and trilinear operator
     */
    double a = [](
        const double _mu,
        const Tensor<dim, dim> _gradu, const Tensor<dim,dim> _gradv)
    {
        Tensor<dim, dim> D = [](
            const Tensor<dim, dim> _gradu)
        {
            return 0.5*(_gradu + transpose(_gradu));
        }

        return 2.*_mu*scalar_product(D(_gradu), D(_gradv));
    }

    double b = [](
        const double _divu,
        const double _q)
    {
        return -_divu*_q;
    }

    double c = [](
        const Tensor<1, dim> _w,
        const Tensor<dim, dim> _gradz,
        const Tensor<1, dim> _v)
    {
        double sum = 0.;

        for (unsigned int i = 0; i < dim; ++dim)
        {
            for (unsigned int j = 0; j < dim; ++dim)
            {
                /*!
                 @todo: Is this indexing correct? 
                 I don't understand that notation in equation (15) from Danaila et al. 2014.
                 */
                sum += _w[j]*_gradz[i,j]*_v[i];
            }
        }

        return sum;
    }

    /*!
     Organize data
     */

    QGauss<dim>   quadrature_formula(this->scalar_degree + 2);
    QGauss<dim-1> face_quadrature_formula(this->scalar_degree + 2);

    FEValues<dim> fe_values (
        this->fe, quadrature_formula,
        update_values    | update_gradients | update_quadrature_points  | update_JxW_values);

    FEFaceValues<dim> fe_face_values (
        this->fe, face_quadrature_formula,
        update_values    | update_normal_vectors | update_quadrature_points  | update_JxW_values);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_quad_points = quadrature_formula.size();
    const unsigned int n_face_quad_points = face_quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> local_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
    const FEValuesExtractors::Scalar temperature(dim);

    std::vector<Tensor<1,dim>> u_k(n_q_points);
    std::vector<double> p_k(n_q_points);
    std::vector<double> theta_k(n_q_points);

    std::vector<Tensor<1,dim>> old_velocity_values(n_q_points);
    std::vector<double> old_pressure_values(n_q_points);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<Tensor<1,dim>> old_temperature_gradients(n_q_points);

    std::vector<Tensor<1,dim>> newton_velocity_values(n_q_points);
    std::vector<double> newton_pressure_values(n_q_points);
    std::vector<double> newton_temperature_values(n_q_points);
    std::vector<Tensor<1,dim>> newton_temperature_gradients(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        endc = this->dof_handler.end();

    /*!
      Assemble element-wise
     */
    const double deltat = this->time_step_size;

    for (; cell != endc; ++cell)
    {
        fe_values.reinit(cell);
        local_matrix = 0;
        local_rhs = 0;


        fe_values[velocities].get_function_values(
            this->old_solution,
            old_velocity_values);

        fe_values[pressure].get_function_values(
            this->old_solution,
            old_pressure_values);

        fe_values[temperature].get_function_values(
            this->old_solution,
            old_temperature_values);

        fe_values[temperature].get_function_gradients(
            this->old_solution,
            old_temperature_gradients);


        fe_values[velocities].get_function_values(
            this->solution,
            newton_velocity_values);

        fe_values[pressure].get_function_values(
            this->solution,
            newton_pressure_values);

        fe_values[temperature].get_function_values(
            this->solution,
            newton_temperature_values);

        fe_values[temperature].get_function_gradients(
            this->solution,
            newton_temperature_gradients);

        fe_values[velocities].get_function_divergences(
            this->solution,
            newton_velocity_divergences);


        for (unsigned int quad = 0; quad < n_quad_points; ++quad)
        {
            /*!
            Initialize local variables to match notation in Danaila 2014
            */
            const Tensor<1, dim>  u_n = old_velocity_values[q];
            const double p_n = old_pressure_values[q];
            const double theta_n = old_temperature_values[q];
            const Tensor<1, dim> gradtheta_n = old_temperature_gradients[q];

            const Tensor<1, dim> u_k = newton_velocity_values[q];
            const double p_k = newton_pressure_values[q];
            const double theta_k = newton_temperature_values[q];
            const double gradtheta_k = newton_temperature_gradients[q];
            const double divu_k = newton_velocity_divergences[q];

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                const Tensor<1, dim> u_w = fe_values[velocities].value(i, quad);
                const double p_w = fe_values[pressure].value(i, quad);
                const double theta_w = fe_values[temperature].value(i, quad);
                const Tensor<1, dim> gradtheta_w = fe_values[temperature].gradient(i, quad);
                const double divu_w = fe_values[velocity].divergence(i, quad);

                for (unsigned int j = 0; j< dofs_per_cell; ++j)
                {
                    const Tensor<1, dim> q = fe_values[velocities].value(j, quad);
                    const double v = fe_values[pressure].value(j, quad);
                    const double phi = fe_values[temperature].value(j, quad);
                    const Tensor<1, dim> gradphi = fe_values[temperature].gradient(j, quad);
                    const double divv = fe_values[velocity].divergence(j, quad);

                    local_matrix(i,j) += // Mass
                        b(gradu_w, q) - gamma(p_w, q);

                    local_matrix(i,j) += // Momentum: Incompressible Navier-Stokes
                        scalar_product(u_w, v)/deltat
                        + c(u_w, u_k, v) + c(u_k, u_w, v) + b(gradv, p_w);

                    local_matrix(i,j) -= // Momentum: Bouyancy (Classical linear Boussinesq approximation)
                        scalar_product(df_B_over_dtheta*theta_w, v);

                    local_matrix(i,j) += // Energy
                        theta_w*phi/deltat
                        - scalar_product(u_k, gradphi)*theta_w
                        - scalar_product(u_w, gradphi)*theta_k
                        + scalar_product(this->K/this->Pr*gradtheta_w, gradphi);
                        
                    // Map to the reference element
                    local_matrix(i,j) *= fe_values.JxW(quad);
                }

                local_rhs(i) += // Mass
                    b(gradu_k, q) - gamma(p_k, q);

                local_rhs(i) += // Momentum: Incompressible Navier-Stokes
                    scalar_product(u_k - u_n, v) + c(u_k, u_k, v) + a(this->mu_l, u_k, v) 
                    + b(gradv, p_k);
                
                local_rhs(i) -= // Momentum: Bouyancy (Classical linear Boussinesq approximation)
                    scalar_product(f_B(theta_k), v);

                local_rhs(i) += // Energy
                    (theta_k - theta_n)*phi/deltat - scalar_product(u_k, gradphi)*theta_k
                    + this->K/this->Pr*gradtheta_k*gradphi;

                /*! @todo: Add forcing function to RHS, e.g. for method of manufactured solution */

        }
            
        // Export local contributions to the global system
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i=0; i < dofs_per_cell; ++i)
        {
            for (unsigned int j=0; j < dofs_per_cell; ++j)
            {
                this->system_matrix.add(
                    local_dof_indices[i],
                    local_dof_indices[j],
                    local_matrix(i,j));
            }
        }

        for (unsigned int i=0; i < dofs_per_cell; ++i)
        {
            this->system_rhs(local_dof_indices[i]) += local_rhs(i);
        }

    }

}

/*!
 @brief Apply the boundary conditions (strong and natural) and apply constraints (including those for hanging nodes).
 
 @author Alexander Zimmerman 2016
*/
template<int dim>
void Peclet<dim>::apply_boundary_values_and_constraints()
{

    /*! Organize boundary functions to simplify application during the time loop */
    
    unsigned int constant_function_index = 0;
    
    for (unsigned int boundary = 0; boundary < boundary_count; boundary++)        
    {
        std::string boundary_type = this->params.boundary_conditions.implementation_types[boundary];
        std::string function_name = this->params.boundary_conditions.function_names[boundary];

        if (function_name == "constant")
        {
            assert(constant_function_index < constant_functions.size());
            this->boundary_functions.push_back(&constant_functions[constant_function_index]);
            constant_function_index++;
        }
        else if (function_name == "parsed")
        {
            this->boundary_functions.push_back(&parsed_boundary_function);
        }
        
    }

    {
        /*!
        @todo Apply natural boundary conditions
        */    

        /*! Homogeneous Neumann boundary conditions on the adiabatic walls are implied */
    }
    
    {
        /*!
         Apply strong boundary conditions

         @todo How do we apply separate boundary conditions to each variable?
        */

        std::map<types::global_dof_index, double> boundary_values;

        for (unsigned int boundary = 0; boundary < boundary_count; boundary++)
        {
            if (this->params.boundary_conditions.implementation_types[boundary] != "strong") 
            {
                continue;
            }
            
            boundary_functions[boundary]->set_time(this->time);
            
            /*! Apply homogeneous Dirichlet boundary conditions for all components of the velocity. */
            const FEValuesExtractors::Vector velocities(0);
            
            VectorTools::interpolate_boundary_values (
                this->dof_handler,
                boundary,
                *boundary_functions[boundary],
                boundary_values,
                this->fe.component_mask(velocities));

            /*! Apply nonhomogeneous Dirichlet boundary conditions for temperature on the hot and cold walls */
            const FEValuesExtractors::Scalar temperature(dim);

            VectorTools::interpolate_boundary_values (
                this->dof_handler,
                boundary,
                *boundary_functions[boundary],
                boundary_values,
                this->fe.component_mask(temperature));
        }

        MatrixTools::apply_boundary_values(
            boundary_values,
            this->system_matrix,
            this->solution,
            this->system_rhs);
    }

}

#endif
