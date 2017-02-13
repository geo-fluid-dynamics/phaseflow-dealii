#ifndef my_vector_tools_h
#define my_vector_tools_h

#include <deal.II/numerics/vector_tools.h>

namespace MyVectorTools
{
    
  using namespace dealii;

  using namespace VectorTools;
  
  /*
  This is just copied from deal.II's VectorTools because, for an unknown reason,
  they have asserted that it is impossible to use this when dim == 1.
  */
  template <int dim, int spacedim>
  void
  my_create_boundary_right_hand_side (const Mapping<dim, spacedim>      &mapping,
                                   const DoFHandler<dim,spacedim>   &dof_handler,
                                   const Quadrature<dim-1> &quadrature,
                                   const Function<spacedim>     &rhs_function,
                                   Vector<double>          &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_ids)
  {
    
    const FiniteElement<dim> &fe  = dof_handler.get_fe();
    Assert (fe.n_components() == rhs_function.n_components,
            ExcDimensionMismatch(fe.n_components(), rhs_function.n_components));
    Assert (rhs_vector.size() == dof_handler.n_dofs(),
            ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));

    rhs_vector = 0;

    UpdateFlags update_flags = UpdateFlags(update_values   |
                                           update_quadrature_points |
                                           update_JxW_values);
    FEFaceValues<dim> fe_values (mapping, fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       n_q_points    = fe_values.n_quadrature_points,
                       n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs (dofs_per_cell);
    Vector<double> cell_vector (dofs_per_cell);

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                            endc = dof_handler.end();

    if (n_components==1)
      {
        std::vector<double> rhs_values(n_q_points);

        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary () &&
                (boundary_ids.empty() ||
                 (boundary_ids.find (cell->face(face)->boundary_id())
                  !=
                  boundary_ids.end())))
              {
                fe_values.reinit(cell, face);

                const std::vector<double> &weights   = fe_values.get_JxW_values ();
                rhs_function.value_list (fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;
                for (unsigned int point=0; point<n_q_points; ++point)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    cell_vector(i) += rhs_values[point] *
                                      fe_values.shape_value(i,point) *
                                      weights[point];

                cell->get_dof_indices (dofs);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
    else
      {
        std::vector<Vector<double> > rhs_values(n_q_points, Vector<double>(n_components));

        for (; cell!=endc; ++cell)
          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary () &&
                (boundary_ids.empty() ||
                 (boundary_ids.find (cell->face(face)->boundary_id())
                  !=
                  boundary_ids.end())))
              {
                fe_values.reinit(cell, face);

                const std::vector<double> &weights   = fe_values.get_JxW_values ();
                rhs_function.vector_value_list (fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;

                // Use the faster code if the
                // FiniteElement is primitive
                if (fe.is_primitive ())
                  {
                    for (unsigned int point=0; point<n_q_points; ++point)
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        {
                          const unsigned int component
                            = fe.system_to_component_index(i).first;

                          cell_vector(i) += rhs_values[point](component) *
                                            fe_values.shape_value(i,point) *
                                            weights[point];
                        }
                  }
                else
                  {
                    // And the full featured
                    // code, if vector valued
                    // FEs are used
                    for (unsigned int point=0; point<n_q_points; ++point)
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int comp_i = 0; comp_i < n_components; ++comp_i)
                          if (fe.get_nonzero_components(i)[comp_i])
                            {
                              cell_vector(i)
                              += rhs_values[point](comp_i) *
                                 fe_values.shape_value_component(i,point,comp_i) *
                                 weights[point];
                            }
                  }

                cell->get_dof_indices (dofs);

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
        }
  }
  
  template <int dim, int spacedim>
  void
  my_create_boundary_right_hand_side (const DoFHandler<dim,spacedim>   &dof_handler,
                                   const Quadrature<dim-1> &quadrature,
                                   const Function<spacedim>     &rhs_function,
                                   Vector<double>          &rhs_vector,
                                   const std::set<types::boundary_id> &boundary_ids)
  {
    my_create_boundary_right_hand_side(StaticMappingQ1<dim>::mapping, dof_handler,
                                    quadrature,
                                    rhs_function, rhs_vector,
                                    boundary_ids);
  }

}

#endif
