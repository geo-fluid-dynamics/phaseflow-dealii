#include <iostream>
#include <fstream>
#include <string>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/signaling_nan.h>


using namespace dealii;

const unsigned int dim = 2;


void validate_vertices(Triangulation<dim> &tria)
{   
    Triangulation<dim>::cell_iterator cell;
    
    unsigned int nv = GeometryInfo<dim>::vertices_per_cell;
      
    for (cell = tria.begin(); cell != tria.end(); ++cell)
    {    
        for (unsigned int v = 0; v < nv; ++v)
        {
            for (unsigned int i = 0; i < dim; ++i)
            {
                assert(!numbers::is_nan(cell->vertex(v)[i]));
            }
            
        }
        
    }
    
}


int main(int /*argc*/, char** /*argv*/)
{
    Triangulation<dim> triangulation;
    
    GridGenerator::hyper_ball(triangulation);
   
    GridOut grid_writer;

    for (unsigned int r = 0; r < 2; ++r)
    {
        std::string file_name = "grid_"+std::to_string(r)+".vtk";
        
        std::ofstream file(file_name);
        
        grid_writer.write_vtk(triangulation, file);   
        
        validate_vertices(triangulation);
        
        triangulation.refine_global(1);
    }
    
    triangulation.set_manifold(0);
    
    return 0;
}
