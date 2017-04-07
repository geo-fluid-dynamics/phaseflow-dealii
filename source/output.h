#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/sparse_matrix.h>

#include <iostream>
#include <fstream>


namespace Output
{
    using namespace dealii;
    
    template<int dim>
    void write_solution_to_vtk(
        const std::string filename,
        DoFHandler<dim> &dof_handler,
        Vector<double> &solution
        )
    {
        /*! Organize the vector-valued solution so that visualization software can understand it*/
        std::vector<std::string> solution_names(dim, "velocity");
        
        solution_names.push_back("pressure");

        if (ENERGY_ENABLED)
        {
            solution_names.push_back("temperature");
        }
        
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
            data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
        
        data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
        
        if (ENERGY_ENABLED)
        {
            data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
        }
        
        DataOut<dim> data_out;
        
        data_out.attach_dof_handler(dof_handler);
        
        data_out.add_data_vector(
            solution,
            solution_names,
            DataOut<dim>::type_dof_data,
            data_component_interpretation);
        
        data_out.build_patches ();

        /*! Write the solution data. */
        std::ofstream output(filename.c_str());
        data_out.write_vtk(output);
    }  
    
    void write_linear_system(SparseMatrix<double> &A, Vector<double> &b)
    {
        
        {
            std::ofstream output("A.txt");
        
            A.print(output);
        }
        
        {
            std::ofstream output("b.txt");
        
            b.print(output);
        }
        
        
    }
    
}
