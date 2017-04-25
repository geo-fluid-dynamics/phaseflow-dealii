#include <iostream>

#include <deal.II/grid/manifold_lib.h>

int main(int argc, char* argv[])
{
    const unsigned int dim = 2;
    
    Triangulation<dim> triangulation;
    
    GridGenerator::hyper_ball(triangulation);
    
    const unsigned int manifold_id = 0;
            
    MyGridGenerator::set_all_manifold_ids(triangulation, manifold_id);
   
    SphericalManifold<dim> spherical_manifold;
    
    this->triangulation.set_manifold(manifold_id, spherical_manifold);      
   
    this->triangulation.refine_global(2);
    
    return 0;
}
