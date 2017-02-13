#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>

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
        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution, "U");
        data_out.build_patches();
        std::ofstream output(filename.c_str());
        data_out.write_vtk(output);
    }    
    
}
