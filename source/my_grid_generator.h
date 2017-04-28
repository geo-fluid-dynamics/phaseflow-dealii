#ifndef my_grid_generator_h
#define my_grid_generator_h

#include <iostream>
#include <cmath>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

namespace MyGridGenerator
{
    using namespace dealii;
    
    template<int dim>
    void create_coarse_grid(
        Triangulation<dim> &triangulation,
        std::vector<unsigned int> &manifold_ids,
        std::vector<std::string> &manifold_descriptors,
        unsigned int &boundary_count,
        const std::string grid_name,
        const std::vector<double> sizes)
    {
        if (grid_name == "hyper_rectangle")
        {
            GridGenerator::hyper_rectangle(
                triangulation,
                {sizes[0], sizes[1]},
                {sizes[2], sizes[3]},
                true);

            boundary_count = 4;
        }
        else if (grid_name == "hyper_shell")
        {
            GridGenerator::hyper_shell(triangulation, Point<dim>(), sizes[0], sizes[1], 0, 0);
            
            triangulation.set_all_manifold_ids(0);
            
            manifold_ids.push_back(0);
            
            manifold_descriptors.push_back("spherical");
            
            boundary_count = 1;
        }
        else
        {
            throw(ExcNotImplemented());
        }
    }
 
}

#endif
