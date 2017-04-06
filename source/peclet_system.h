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

    this->solution.reinit(this->dof_handler.n_dofs());

    this->newton_residual.reinit(this->dof_handler.n_dofs());
    
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
void Peclet<dim>::assemble_system()
{
    /*!
     Local parameters
    */
    const double penalty = 1.e-7; // @todo: Expose this to ParameterHandler.

    const double mu_l = LIQUID_DYNAMIC_VISOCITY;

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
        return _v*_gradz*_w;
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

    std::vector<Tensor<1,dim>> u_k(n_quad_points);
    std::vector<double> p_k(n_quad_points);
    std::vector<double> theta_k(n_quad_points);

    std::vector<Tensor<1,dim>> old_velocity_values(n_quad_points);
    std::vector<double> old_pressure_values(n_quad_points);

    std::vector<Tensor<1,dim>> newton_velocity_values(n_quad_points);
    std::vector<double> newton_pressure_values(n_quad_points);
    std::vector<Tensor<2,dim>> newton_velocity_gradients(n_quad_points);

    std::vector<double> newton_velocity_divergences(n_quad_points);

    typename DoFHandler<dim>::active_cell_iterator
        cell = this->dof_handler.begin_active(),
        endc = this->dof_handler.end();

    
    /*!
        Set local variables to match notation in Danaila 2014
    */
    const double gamma = penalty;
    
    for (; cell != endc; ++cell) /*! Assemble element-wise */
    {
        fe_values.reinit(cell);
        local_matrix = 0;
        local_rhs = 0;

        fe_values[velocity].get_function_values(
            this->old_newton_solution,
            newton_velocity_values);

        fe_values[pressure].get_function_values(
            this->old_newton_solution,
            newton_pressure_values);

        fe_values[velocity].get_function_gradients(
            this->old_newton_solution,
            newton_velocity_gradients);

        fe_values[velocity].get_function_divergences(
            this->old_newton_solution,
            newton_velocity_divergences);

        for (unsigned int quad = 0; quad< n_quad_points; ++quad)
        {
            /*!
            Name local variables to match notation in Danaila 2014
            */
            const Tensor<1, dim> u_k = newton_velocity_values[quad];
            const double p_k = newton_pressure_values[quad];
            const Tensor<2, dim> gradu_k = newton_velocity_gradients[quad];
            const double divu_k = newton_velocity_divergences[quad];

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                const Tensor<1, dim> u_w = fe_values[velocity].value(i, quad);
                const double p_w = fe_values[pressure].value(i, quad);
                const Tensor<2, dim> gradu_w = fe_values[velocity].gradient(i, quad);
                const double divu_w = fe_values[velocity].divergence(i, quad);

                for (unsigned int j = 0; j< dofs_per_cell; ++j)
                {
                    const Tensor<1, dim> v = fe_values[velocity].value(j, quad);
                    const double q = fe_values[pressure].value(j, quad);
                    const Tensor<2, dim> gradv = fe_values[velocity].gradient(j, quad);
                    const double divv = fe_values[velocity].divergence(j, quad);

                    local_matrix(i,j) += // Mass
                        b(divu_w, q) - gamma*p_w*q;

                    local_matrix(i,j) += // Momentum: Incompressible Navier-Stokes
                        + c(u_w, gradu_k, v) + c(u_k, gradu_w, v) + a(mu_l, gradu_w, gradv) + b(divv, p_w);

                    local_matrix(i,j) *= fe_values.JxW(quad);  /*! Map to the reference element */

                    local_rhs(i) += // Mass
                        b(divu_k, q) - gamma*p_k*q;

                    local_rhs(i) += // Momentum: Incompressible Navier-Stokes
                        c(u_k, gradu_k, v) + a(mu_l, gradu_k, gradv) 
                        + b(divv, p_k);

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

         To reproduce results in Danaila 2014, use homogeneous Dirichlet BC's on every boundary
        */

        std::map<types::global_dof_index, double> boundary_values;

        for (unsigned int b = 0; b < this->boundary_count; ++b) /* For each boundary */
        {                    
            VectorTools::interpolate_boundary_values(
                this->dof_handler,
                b,
                ZeroFunction<dim>(dim + 1),
                boundary_values);      
        }

        MatrixTools::apply_boundary_values(
            boundary_values,
            this->system_matrix,
            this->newton_residual,
            this->system_rhs);
    }

}

#endif
