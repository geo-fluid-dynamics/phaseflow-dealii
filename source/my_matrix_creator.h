#ifndef _my_matrix_creator_h_
#define _my_matrix_creator_h_

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>

#include <deal.II/numerics/matrix_creator.templates.h>


namespace MyMatrixCreator
{
    using namespace dealii;
    
    namespace AssemblerData
    {
      template <int dim,
                typename number>
      struct Scratch
      {
        Scratch (const ::dealii::hp::FECollection<dim> &fe,
                 const UpdateFlags         update_flags,
                 const Function<dim,number> *diffusivity,
                 const Function<dim,number> *convection_velocity,
                 const ::dealii::hp::QCollection<dim> &quadrature,
                 const ::dealii::hp::MappingCollection<dim> &mapping)
          :
          fe_collection (fe),
          quadrature_collection (quadrature),
          mapping_collection (mapping),
          x_fe_values (mapping_collection,
                       fe_collection,
                       quadrature_collection,
                       update_flags),
          diffusivity_values(quadrature_collection.max_n_quadrature_points()),
          convection_velocity_values(quadrature_collection.max_n_quadrature_points(),
                                     Vector<number> (dim)),
          diffusivity (diffusivity),
          convection_velocity (convection_velocity),
          update_flags (update_flags)
        {}

        Scratch (const Scratch &data)
          :
          fe_collection (data.fe_collection),
          quadrature_collection (data.quadrature_collection),
          mapping_collection (data.mapping_collection),
          x_fe_values (mapping_collection,
                       fe_collection,
                       quadrature_collection,
                       data.update_flags),
          diffusivity_values (data.diffusivity_values),
          convection_velocity_values (data.convection_velocity_values),
          diffusivity (data.diffusivity),
          convection_velocity (data.convection_velocity),
          rhs_function (data.rhs_function),
          update_flags (data.update_flags)
        {}

        const ::dealii::hp::FECollection<dim>      &fe_collection;
        const ::dealii::hp::QCollection<dim>       &quadrature_collection;
        const ::dealii::hp::MappingCollection<dim> &mapping_collection;

        ::dealii::hp::FEValues<dim> x_fe_values;

        std::vector<number>                  diffusivity_values;
        std::vector<dealii::Vector<number> > convection_velocity_values;
        std::vector<number>                  rhs_values;
        std::vector<dealii::Vector<number> > rhs_vector_values;

        const Function<dim,number>   *diffusivity;
        const Function<dim,number>   *convection_velocity;
        const Function<dim,number>   *rhs_function;

        const UpdateFlags update_flags;
      };

    }
    
    /*

    @brief Create the convection-diffusion matrix

    @detail

        This is the sum of the convention and diffusions matrices, C + K,
        for the convection-diffusion equation in
        Finite Element Methods for Flow Problems (Donea & Huerta, 2003)
        
        K is equivalent to the Laplace matrix, and so most of this is 
        copied from dealii::MatrixCreator::create_laplace_matrix
        
        Super important note: The convection operator is asymmetric.
        
    @author A. Zimmerman <zimmerman@aices.rwth-aachen.de> 
    
    */
    template <int dim,
              typename CellIterator>
    void convection_diffusion_assembler (
        const CellIterator &cell,
        AssemblerData::Scratch<dim,double> &data,
        MatrixCreator::internal::AssemblerData::CopyData<double> &copy_data)
    {
        data.x_fe_values.reinit (cell);
        const FEValues<dim> &fe_values = data.x_fe_values.get_present_fe_values ();

        const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         n_q_points    = fe_values.n_quadrature_points;
                         
        const FiniteElement<dim>    &fe  = fe_values.get_fe();
             
        Assert(data.diffusivity->n_components==1,
             ::dealii::MatrixCreator::ExcComponentMismatch());
             
        Assert(data.convection_velocity->n_components==dim,
             ::dealii::MatrixCreator::ExcComponentMismatch());

        copy_data.cell_matrix.reinit (dofs_per_cell, dofs_per_cell);

        copy_data.dof_indices.resize (dofs_per_cell);

        cell->get_dof_indices (copy_data.dof_indices);


        data.diffusivity_values.resize (n_q_points);

        data.diffusivity->value_list (fe_values.get_quadrature_points(),
                                    data.diffusivity_values);

        data.convection_velocity_values.resize (n_q_points,
                                              dealii::Vector<double>(dim));

        data.convection_velocity->vector_value_list (fe_values.get_quadrature_points(),
                                            data.convection_velocity_values);
                                    

        const std::vector<double> &JxW = fe_values.get_JxW_values();

        double add_data;

        assert(fe.is_primitive());

        for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
            const Tensor<1,dim> *grad_phi_i = &fe_values.shape_grad(i,0);
            const double *phi_i = &fe_values.shape_value(i,0);
            
            for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
                const Tensor<1,dim> *grad_phi_j = &fe_values.shape_grad(j,0);
                
                add_data = 0;
                
                for (unsigned int point=0; point<n_q_points; ++point)
                {
                    add_data += (grad_phi_i[point]*grad_phi_j[point]) *
                                JxW[point] *
                                data.diffusivity_values[point];  

                    Tensor<1,dim> a;
                    for (unsigned int ia = 0; ia < dim; ia++)
                    { // @todo: What is the proper way to convert this?
                        a[ia] = data.convection_velocity_values[point][ia];
                    }
                    
                    add_data += (phi_i[point] * ( a * grad_phi_j[point]) ) * JxW[point];
                    
                }
                
                copy_data.cell_matrix(i,j) = add_data;
                
            }
            
        }
            
    }
    
    
    template <int dim>
    void create_convection_diffusion_matrix (
        const Mapping<dim> &mapping,
        const DoFHandler<dim> &dof,
        const Quadrature<dim> &q,
        SparseMatrix<double> &matrix,
        const Function<dim> *const diffusivity,
        const Function<dim> *const convection_velocity,
        const ConstraintMatrix & constraints = ConstraintMatrix())
    {
        Assert (matrix.m() == dof.n_dofs(), ExcDimensionMismatch (matrix.m(), dof.n_dofs()));
        Assert (matrix.n() == dof.n_dofs(), ExcDimensionMismatch (matrix.n(), dof.n_dofs()));

        hp::FECollection<dim>      fe_collection (dof.get_fe());
        hp::QCollection<dim>                q_collection (q);
        hp::MappingCollection<dim> mapping_collection (mapping);
        
        AssemblerData::Scratch<dim,double> assembler_data (
                            fe_collection,
                            update_values | update_gradients  |
                                update_JxW_values | update_quadrature_points,
                            diffusivity,
                            convection_velocity,
                            q_collection, mapping_collection);
                        
        MatrixCreator::internal::AssemblerData::CopyData<double> copy_data;
        
        copy_data.cell_matrix.reinit (assembler_data.fe_collection.max_dofs_per_cell(),
                                      assembler_data.fe_collection.max_dofs_per_cell());
        
        copy_data.dof_indices.resize (assembler_data.fe_collection.max_dofs_per_cell());

        copy_data.constraints = &constraints;
        
        WorkStream::run(
            dof.begin_active(),
            static_cast<typename DoFHandler<dim>::active_cell_iterator>(dof.end()),
            &convection_diffusion_assembler<dim, typename DoFHandler<dim>::active_cell_iterator>,
            std_cxx11::bind (&MatrixCreator::internal::
                          copy_local_to_global<double,SparseMatrix<double>, Vector<double> >,
                          std_cxx11::_1,
                          &matrix,
                          (Vector<double> *)NULL),
            assembler_data,
            copy_data);
    }


    template <int dim>
    void create_convection_diffusion_matrix (
        const DoFHandler<dim> &dof,
        const Quadrature<dim> &q,
        SparseMatrix<double> &matrix,
        const Function<dim> *const diffusivity,
        const Function<dim> *const convection_velocity,
        const ConstraintMatrix & constraints = ConstraintMatrix())
    {
        create_convection_diffusion_matrix(StaticMappingQ1<dim>::mapping,
            dof, q, matrix, diffusivity, convection_velocity, constraints);
    }
    
}

#endif
