#ifndef _peclet_system_h_
#define _peclet_system_h_

/*!
    @brief Setup the linear system objects.
 
    @author Alexander Zimmerman 2016
*/
template<int dim>
void Peclet<dim>::setup_system(bool quiet)
{
    
    this->dof_handler.distribute_dofs(this->fe);

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
void Peclet<dim>::assemble_system()
{
    /*!
     Local parameters
    */
    const double penalty = 1.e-7; // @todo: Expose this to ParameterHandler.

    const double
        Ra = RAYLEIGH_NUMBER,
        Pr = PRANDTL_NUMBER,
        Re = REYNOLDS_NUMBER;
    
    const double K = SOLID_CONDUCTIVITY/LIQUID_CONDUCTIVITY;

    Tensor<1, dim> g; // @todo: Make this const.
    for (unsigned int i = 0; i < dim; ++i)
    {
        g[i] = GRAVITY[i];
    }

    const double mu_l = LIQUID_DYNAMIC_VISOCITY;

    const double Ra_over_PrRe2(Ra/(Pr*Re*Re));

    /*!
     lambda function for classical (linear) Boussinesq bouyancy
    */
    auto f_B = [Ra_over_PrRe2, g](const double _theta) 
    {
        return _theta*Ra_over_PrRe2*g;
    };

    /*!
     Analytical derivative of classical (linear) Boussinesq bouyancy
    */
    const Tensor<1, dim> df_B_over_dtheta(Ra_over_PrRe2*g);

    /*!
     lambda functions for linear, bilinear, and trilinear operator
    */
    auto a = [](
        const double _mu,
        const Tensor<2, dim> _gradu, const Tensor<2,dim> _gradv)
    {
        auto D = [_mu](
            const Tensor<2, dim> _gradu)
        {
            return 0.5*(_gradu + transpose(_gradu));
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
        const double _divw,
        const Tensor<1, dim> _z,
        const Tensor<1, dim> _v)
    {
        return scalar_product(_divw*_z, _v);
    };

    /*!
     Organize data
    */

    QGauss<dim>   quadrature_formula(SCALAR_DEGREE + 2);
    QGauss<dim-1> face_quadrature_formula(SCALAR_DEGREE + 2);

    FEValues<dim> fe_values (
        this->fe, quadrature_formula,
        update_values    | update_gradients | update_quadrature_points  | update_JxW_values);

    FEFaceValues<dim> fe_face_values (
        this->fe, face_quadrature_formula,
        update_values    | update_normal_vectors | update_quadrature_points  | update_JxW_values);

    const unsigned int dofs_per_cell = this->fe.dofs_per_cell;
    const unsigned int n_quad_points = quadrature_formula.size();

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> local_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector velocity(0);
    const FEValuesExtractors::Scalar pressure(dim);
    const FEValuesExtractors::Scalar temperature(dim);

    std::vector<Tensor<1,dim>> u_k(n_quad_points);
    std::vector<double> p_k(n_quad_points);
    std::vector<double> theta_k(n_quad_points);

    std::vector<Tensor<1,dim>> old_velocity_values(n_quad_points);
    std::vector<double> old_pressure_values(n_quad_points);
    std::vector<double> old_temperature_values(n_quad_points);
    std::vector<Tensor<1,dim>> old_temperature_gradients(n_quad_points);

    std::vector<Tensor<1,dim>> newton_velocity_values(n_quad_points);
    std::vector<double> newton_pressure_values(n_quad_points);
    std::vector<double> newton_temperature_values(n_quad_points);
    std::vector<Tensor<1,dim>> newton_temperature_gradients(n_quad_points);
    std::vector<Tensor<2,dim>> newton_velocity_gradients(n_quad_points);

    std::vector<double> newton_velocity_divergences(n_quad_points);

    typename DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        endc = this->dof_handler.end();

    
    /*!
        Set local variables to match notation in Danaila 2014
    */
    const double deltat = this->time_step_size;
    const double gamma = penalty;
    
    for (; cell != endc; ++cell) /*! Assemble element-wise */
    {
        fe_values.reinit(cell);
        local_matrix = 0;
        local_rhs = 0;


        fe_values[velocity].get_function_values(
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


        fe_values[velocity].get_function_values(
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

        fe_values[velocity].get_function_gradients(
            this->solution,
            newton_velocity_gradients);

        fe_values[velocity].get_function_divergences(
            this->solution,
            newton_velocity_divergences);


        for (unsigned int quad = 0; quad< n_quad_points; ++quad)
        {
            /*!
            Name local variables to match notation in Danaila 2014
            */
            const Tensor<1, dim>  u_n = old_velocity_values[quad];
            const double theta_n = old_temperature_values[quad];
            const Tensor<1, dim> u_k = newton_velocity_values[quad];
            const double p_k = newton_pressure_values[quad];
            const double theta_k = newton_temperature_values[quad];
            const Tensor<1, dim> gradtheta_k = newton_temperature_gradients[quad];
            const Tensor<2, dim> gradu_k = newton_velocity_gradients[quad];
            const double divu_k = newton_velocity_divergences[quad];

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                const Tensor<1, dim> u_w = fe_values[velocity].value(i, quad);
                const double p_w = fe_values[pressure].value(i, quad);
                const double theta_w = fe_values[temperature].value(i, quad);
                const Tensor<1, dim> gradtheta_w = fe_values[temperature].gradient(i, quad);
                const Tensor<2, dim> gradu_w = fe_values[velocity].gradient(i, quad);
                const double divu_w = fe_values[velocity].divergence(i, quad);

                for (unsigned int j = 0; j< dofs_per_cell; ++j)
                {
                    const Tensor<1, dim> v = fe_values[velocity].value(j, quad);
                    const double q = fe_values[pressure].value(j, quad);
                    const double phi = fe_values[temperature].value(j, quad);
                    const Tensor<1, dim> gradphi = fe_values[temperature].gradient(j, quad);
                    const Tensor<2, dim> gradv = fe_values[velocity].gradient(j, quad);
                    const double divv = fe_values[velocity].divergence(j, quad);

                    local_matrix(i,j) += // Mass
                        b(divu_w, q) - gamma*p_w*q;

                    local_matrix(i,j) += // Momentum: Incompressible Navier-Stokes
                        scalar_product(u_w, v)/deltat
                        + c(divu_w, u_k, v) + c(divu_k, u_w, v) + a(mu_l, gradu_w, gradv) + b(divv, p_w);

                    local_matrix(i,j) -= // Momentum: Bouyancy (Classical linear Boussinesq approximation)
                        scalar_product(df_B_over_dtheta*theta_w, v);

                    local_matrix(i,j) += // Energy
                        theta_w*phi/deltat
                        - scalar_product(u_k, gradphi)*theta_w
                        - scalar_product(u_w, gradphi)*theta_k
                        + scalar_product(K/Pr*gradtheta_w, gradphi);

                    local_matrix(i,j) *= fe_values.JxW(quad);  /*! Map to the reference element */

                    local_rhs(i) += // Mass
                        b(divu_k, q) - gamma*p_k*q;

                    local_rhs(i) += // Momentum: Incompressible Navier-Stokes
                        scalar_product(u_k - u_n, v) + c(divu_k, u_k, v) + a(mu_l, gradu_k, gradv) 
                        + b(divv, p_k);

                    local_rhs(i) -= // Momentum: Bouyancy (Classical linear Boussinesq approximation)
                        scalar_product(f_B(theta_k), v);

                    local_rhs(i) += // Energy
                        (theta_k - theta_n)*phi/deltat - scalar_product(u_k, gradphi)*theta_k
                        + K/Pr*gradtheta_k*gradphi;
                }

                /*! @todo: Add forcing function to RHS, e.g. for method of manufactured solution */

                local_rhs(i) *= fe_values.JxW(quad); /*! Map to the reference element */
            }
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

    std::map<types::global_dof_index, double> boundary_values;

    {
        /*!
        @todo Apply natural boundary conditions
        */    

        /*! Homogeneous Neumann boundary conditions on the adiabatic walls are implied */
    }
    
    {
        /*!
         Apply strong boundary conditions
        */

        std::map<types::global_dof_index, double> boundary_values;

        for (unsigned int b = 0; b < this->boundary_count; ++b) /* For each boundary */
        {                    
            std::vector<std::string> mask = this->params.boundary_conditions.strong_masks[b];

            FEValuesExtractors::Vector velocity_extractor(0);
            FEValuesExtractors::Scalar pressure_extractor(dim);
            FEValuesExtractors::Scalar temperature_extractor(dim + 1);

            std::vector<std::string> field_names({"velocity", "pressure", "temperature"});

            for (unsigned int f = 0; f < field_names.size(); ++f) /* For each field variable */
            {
                std::string field_name = field_names[f];

                if (std::find(mask.begin(), mask.end(), field_name) == mask.end()) /* Skip if the field name is not in the mask */
                {
                    continue;
                }

                /*!
                    @todo: 
                
                    Is there some way to contain the extractors (or pointers to them) in a single object that can be indexed?

                    Neither std::vector<void*> nor tuple (because the tuple could not be indexed with a variable) worked, and I'm out of ideas.
                */

                if (field_name == "velocity")
                {
                    VectorTools::interpolate_boundary_values(
                        this->dof_handler, b, *this->boundary_function_pointers[b], boundary_values,
                        this->fe.component_mask(velocity_extractor));
                }
                else if (field_name == "pressure")
                {
                    VectorTools::interpolate_boundary_values(
                        this->dof_handler, b, *this->boundary_function_pointers[b], boundary_values,
                        this->fe.component_mask(pressure_extractor));
                }
                else if (field_name == "temperature")
                {
                    VectorTools::interpolate_boundary_values(
                        this->dof_handler, b, *this->boundary_function_pointers[b], boundary_values,
                        this->fe.component_mask(temperature_extractor));
                }
                else
                {
                    assert(false);
                }

            }
            
        }

        MatrixTools::apply_boundary_values(
            boundary_values,
            this->system_matrix,
            this->solution,
            this->system_rhs);
    }

}

#endif
