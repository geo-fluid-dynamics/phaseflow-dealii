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
    void set_all_manifold_ids(Triangulation<dim> &tria, unsigned int id=0)
    {
        auto cell = tria.begin_active();
        auto endc = tria.end();
        for (; cell!=endc; ++cell)
        {
            cell->set_all_manifold_ids(id);
        }
    }
    
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

            boundary_count = pow(2, dim);
        }
        else if (grid_name == "hyper_ball")
        {
            GridGenerator::hyper_ball(triangulation);
            
            MyGridGenerator::set_all_manifold_ids(triangulation, 0);
            
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
