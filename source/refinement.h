#ifndef _refinement_h_
#define _refinement_h_

namespace Refinement
{

    /*!
    @brief Adaptive grid refinement method from deal.II's step-26 tutorial.
    */
    template <int dim>
    void adaptive_refine_mesh(
        Triangulation<dim> &triangulation,
        DoFHandler<dim> &dof_handler,
        Vector<double> &solution,
        SolutionTransfer<dim> &solution_trans,
        const FE_Q<dim> fe,
        const unsigned int min_grid_level,
        const unsigned int max_grid_level,
        const unsigned int max_cells,
        const double refine_fraction,
        const double coarsen_fraction)
    {
        Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
        KellyErrorEstimator<dim>::estimate(
            dof_handler,
            QGauss<dim-1>(fe.degree+1),
            typename FunctionMap<dim>::type(),
            solution,
            estimated_error_per_cell);
        GridRefinement::refine_and_coarsen_fixed_fraction(
            triangulation,
            estimated_error_per_cell,
            refine_fraction,
            coarsen_fraction);
        if (triangulation.n_levels() > max_grid_level)
        {
            for (auto cell = triangulation.begin_active(max_grid_level);
                 cell != triangulation.end(); ++cell)
            {
                cell->clear_refine_flag ();    
            }
        }
        for (auto cell = triangulation.begin_active(min_grid_level);
             cell != triangulation.end_active(min_grid_level); ++cell)
        {
            cell->clear_coarsen_flag ();     
        }
        if ((max_cells > 0) & (triangulation.n_active_cells() > max_cells))
        {
           for (auto cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
            {
                cell->clear_refine_flag ();
            }
        }
        triangulation.prepare_coarsening_and_refinement();
        solution_trans.prepare_for_coarsening_and_refinement(solution);
        triangulation.execute_coarsening_and_refinement();
    }

    /*!
    @brief Adaptive grid refinement for vector-valued PDE
    @todo Implement Danaila's metric
    */
    template <int dim>
    void adaptive_refine_mesh(
        Triangulation<dim> &triangulation,
        DoFHandler<dim> &dof_handler,
        Vector<double> &solution,
        SolutionTransfer<dim> &solution_trans,
        const FESystem<dim,dim> fe,
        const unsigned int min_grid_level,
        const unsigned int max_grid_level,
        const unsigned int max_cells,
        const double refine_fraction,
        const double coarsen_fraction)
    {
        Assert(false, ExcNotImplemented());
    }
    
    template <int dim>
    void refine_mesh_near_boundaries (
        Triangulation<dim> &triangulation,
        const std::vector<unsigned int> boundary_ids,
        const unsigned int refinement_cycles)
    {
        for (unsigned int i = 0; i < refinement_cycles; i++)
        {
            for (auto cell : triangulation.active_cell_iterators())
            {
                if (!(cell->at_boundary()))
                {
                    continue;
                }
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                    if (!cell->face(f)->at_boundary())
                    {
                        continue;
                    }
                    bool flagged = false;
                    for (auto boundary_id : boundary_ids)
                    {
                        if (cell->face(f)->boundary_id() == boundary_id)
                        {
                            cell->set_refine_flag();
                            flagged = true;
                            break;
                        }   
                    }
                    if (flagged)
                    {
                        break;
                    }
                }
            }
            triangulation.execute_coarsening_and_refinement(); 
        }
    }   
    
}

#endif
