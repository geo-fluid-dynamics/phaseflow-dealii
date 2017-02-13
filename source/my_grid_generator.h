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
    void cylinder_with_split_boundaries
        (
        Triangulation<dim> &tria,
        const double L0, const double L1, const double L2
        );
    template<int dim>
    void cylinder_with_split_boundaries
        (
        Triangulation<dim> &tria,
        const double L0, const double L1, const double L2
        )
    {
        // Based on create_coarse_grid in http://dealii.org/8.4.1/doxygen/deal.II/step_14.html
        Assert (dim==2, dealii::ExcNotImplemented());
        // Set the vertices
        static const Point<dim> static_points[] =
        {
            // Top
            Point<dim>(-(L0 + L1), 0.),
            Point<dim>(-L0, 0.),
            Point<dim>(L0, 0.),
            Point<dim>(L0 + L1, 0.),
            // Bottom
            Point<dim>(-(L0 + L1), -L2),
            Point<dim>(-L0, -L2),
            Point<dim>(L0, -L2),
            Point<dim>(L0 + L1, -L2)
        };
        const unsigned int point_count = sizeof(static_points)/sizeof(static_points[0]);
        const std::vector<Point<dim>> points(&static_points[0], &static_points[point_count]);
        const int vertices_per_cell = 4; // ISO C++ forbids variable length array
        const int cell_count = 3;
        // Set the cells
        // Must be careful about the orientation of cells; the deal.II rules are unclear.
        static const int cell_vertices[cell_count][vertices_per_cell] =
        {
            {0, 4, 1, 5}, 
            {1, 5, 2, 6},
            {2, 6, 3, 7}
        };
        std::vector<dealii::CellData<dim> > cells(cell_count, dealii::CellData<dim>());
        for (unsigned int i=0; i< cell_count; ++i)
        {
            for (unsigned int j=0; j < vertices_per_cell; ++j)
            {
                cells[i].vertices[j] = cell_vertices[i][j];
            }
            cells[i].material_id = 0;
        }
        // Assemble the grid
        tria.create_triangulation(points, cells, dealii::SubCellData());
        // Set boundary ID's.
        auto cell = tria.begin_active();
        double epsilon = 1e-8*L0;
        for (; cell != tria.end(); ++cell) {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
                if (cell->face(f)->at_boundary()) {
                    double xf = cell->face(f)->center()[0], yf = cell->face(f)->center()[1];
                    if (yf > -epsilon)
                    { // Top 
                        if (abs(xf) < epsilon)
                        { // Top-Middle
                            cell->face(f)->set_boundary_id(0);
                        }
                        else
                        { // Top-Edge
                            cell->face(f)->set_boundary_id(1);
                        }
                    }
                    else if (abs(yf + L2/2.) < epsilon)
                    { // Side
                        cell->face(f)->set_boundary_id(2);
                    }
                    else
                    { // Bottom
                        cell->face(f)->set_boundary_id(3);
                    }
                }
            }
        }
        //
        set_all_manifold_ids(tria, 0);
    }
  
    /**
    @brief Produce a domain that is the space between two hemisphere-cylinders.
    @detail
    Based on hyper_cube_with_cylindrical_hole
    Boundary ID's follow the vertex numbering, with 
        0 -> right side of outer spherical boundary
        ... counter-clockwise
        4 -> left side of outer spherical boundary
        5 -> right side of inner spherical boundary
        ... counter-clockwise
        9 -> left side of inner spherical boundary
    Manifold ID's:
        0 -> Hemi-spherical manifold radially centered at the origin
        1 -> Cyldrincal manifoldd
    The origin of the coordinate system is at the radial center of the spherical manifold.
    @todo Insert image showing coarse grid with overlaid geometry and labels
    @todo Implement in 3D
    **/
    template<int dim>
    void hemisphere_cylinder_shell (
        Triangulation<dim> &tria,
        const double inner_radius=0.25,
        const double outer_radius=0.5,
        const double inner_length=1.0,
        const double outer_length=1.25);

    template <>
    void hemisphere_cylinder_shell (
        Triangulation<1> &,
        const double,
        const double,
        const double,
        const double)
    {
        Assert (false, ExcNotImplemented()); // There should be no reason to implement this in 1D.
    }
  
    template<>
    void hemisphere_cylinder_shell (
        Triangulation<2> & tria,
        const double inner_radius, // Radius of shell's inside boundary
        const double outer_radius, // Radius of shell's outside boundary
        const double inner_length, // Length of the inside cylinder
        const double outer_length  // Length of the outside cylinder
    )
    {
        /*
        Based on hyper_cube_with_cylindrical_hole
        Boundary ID's follow the vertex numbering, with 
            0 -> right side of outer spherical boundary
            ... counter-clockwise
            4 -> left side of outer spherical boundary
            5 -> right side of inner spherical boundary
            ... counter-clockwise
            9 -> left side of inner spherical boundary
        Manifold ID's:
            0 -> Hemi-spherical manifold radially centered at the origin
            1 -> Cyldrincal manifoldd
        The origin of the coordinate system is at the radial center of the spherical manifold.
        */
        const int dim = 2;
        Assert(outer_radius > inner_radius, ExcInvalidState());
        Assert(outer_length > inner_length, ExcInvalidState());
        // To ensure valid indexing and to simplify extension to 3D:
        //  - Create an hyper_shell in two dimensions.
        //  - Modify it to match our desired shape.
        GridGenerator::hyper_shell(tria, {0, 0}, inner_radius, outer_radius, 5, true);
        // Rotate the grid ninety degrees to obtain symmetry about the y-axis.
        GridTools::rotate(-numbers::PI/2., tria);
        // Modify half of the grid into a cylindrical shell.
        Triangulation<dim>::active_cell_iterator cell = tria.begin_active(), endc = tria.end();
        std::vector<bool> treated_vertices(tria.n_vertices(), false);
        for (; cell != endc; ++cell) {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
                if (cell->face(f)->at_boundary()) {
                    for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
                        unsigned int vv = cell->face(f)->vertex_index(v);
                        if (treated_vertices[vv] == false) {
                            treated_vertices[vv] = true;
                            switch (vv) {
                                case 1:
                                    cell->face(f)->vertex(v) = Point<dim>(outer_radius,0.);
                                    break;
                                case 2:
                                    cell->face(f)->vertex(v) = Point<dim>(outer_radius,outer_length);
                                    break;
                                case 3:
                                    cell->face(f)->vertex(v) = Point<dim>(-outer_radius,outer_length);
                                    break;
                                case 4:
                                    cell->face(f)->vertex(v) = Point<dim>(-outer_radius,0.);
                                    break;
                                case 6:
                                    cell->face(f)->vertex(v) = Point<dim>(inner_radius,0.);
                                    break;
                                case 7:
                                    cell->face(f)->vertex(v) = Point<dim>(inner_radius,inner_length);
                                    break;
                                case 8:
                                    cell->face(f)->vertex(v) = Point<dim>(-inner_radius,inner_length);
                                    break;
                                case 9:
                                    cell->face(f)->vertex(v) = Point<dim>(-inner_radius,0.);
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // Set boundary ID's corresponding to vertex numbers.
        // This uses the same loop structure as above.
        cell = tria.begin_active();
        std::fill(treated_vertices.begin(), treated_vertices.end(), false);
        double eps = 1e-8*inner_radius;
        for (; cell != endc; ++cell) {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f) {
                if (cell->face(f)->at_boundary()) {
                    double xf = cell->face(f)->center()[0], yf = cell->face(f)->center()[1];
                    for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
                        unsigned int vv = cell->face(f)->vertex_index(v);
                        if (treated_vertices[vv] == false) {
                            switch (vv) {
                                case 0:
                                    if (xf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 1:
                                    if (yf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 2:
                                    if (xf < eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 3:
                                    if (yf < inner_length) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 4:
                                    if (yf < -eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 5:
                                    if (xf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 6:
                                    if (yf > eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 7:
                                    if (xf < eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 8:
                                    if (yf < inner_length - eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                                case 9:
                                    if (yf < -eps) {
                                        cell->face(f)->set_boundary_id(vv);
                                        treated_vertices[vv] = true;
                                    }
                                    break;
                            }
                        }
                    }
                }
            }
        }
        // Label the spherical and cylindrical manifolds.
        cell = tria.begin_active();
        for (; cell!=endc; ++cell) {
            if ((cell->center())[1] < eps) {
                cell->set_all_manifold_ids(0); // Spherical
            }
            else {
                cell->set_all_manifold_ids(1); // Cylindrical
            }
        }
    }
  
    template <>
    void hemisphere_cylinder_shell (
        Triangulation<3> &,
        const double,
        const double,
        const double,
        const double)
    {
        Assert (false, ExcNotImplemented()); // A 3D implementation is needed.
    }
 
    using namespace dealii;
    
    template<int dim>
    void create_coarse_grid(
        Triangulation<dim> &triangulation,
        std::vector<unsigned int> &manifold_ids,
        std::vector<std::string> &manifold_descriptors,
        const std::string grid_name,
        const std::vector<double> sizes)
    {
        if (grid_name == "hyper_cube")
        {
            GridGenerator::hyper_cube(
                triangulation,
                sizes[0],
                sizes[1],
                true);
        }
        else if (grid_name == "hyper_rectangle")
        {
            GridGenerator::hyper_rectangle(
                triangulation,
                {sizes[0], sizes[1]},
                {sizes[2], sizes[3]},
                true);
        }
        else if (grid_name == "hyper_shell")
        {
            Point<dim> center; 
            for (unsigned int i = 0; i < dim; i++)
            {
                center[i] = 0.; // The grid will be shifted later for consistency with the other grid methods.
            }
            GridGenerator::hyper_shell
                (
                triangulation,
                center,
                sizes[0], sizes[1],
                0, true 
                );
            MyGridGenerator::set_all_manifold_ids(triangulation, 0);
            manifold_ids.push_back(0);
            manifold_descriptors.push_back("spherical");
        }
        else if (grid_name == "hemisphere_cylinder_shell") 
        {
            MyGridGenerator::hemisphere_cylinder_shell
               (
                triangulation,
                sizes[0], sizes[1], sizes[2], sizes[3]
                );
            manifold_ids.push_back(0);
            manifold_descriptors.push_back("spherical");
            manifold_ids.push_back(1);
            manifold_descriptors.push_back("cylindrical");
        }
        else if (grid_name == "cylinder")
        {
            assert(dim > 1);
            
            if (dim == 2)
            {
                GridGenerator::hyper_rectangle(
                    triangulation,
                    {-sizes[0], 0.},
                    {sizes[0], -sizes[1]},
                    true);
            }
            else if (dim == 3)
            {
                GridGenerator::cylinder(
                    triangulation,
                    sizes[0],
                    sizes[1]/2.);
                Tensor<1,dim> shift_vector;
                shift_vector[1] = -sizes[1]/2.;
                GridTools::shift(shift_vector, triangulation);
            }
            
            manifold_ids.push_back(0);
            manifold_descriptors.push_back("cylindrical");
        }
        else if (grid_name == "cylinder_with_split_boundaries")
        {
            MyGridGenerator::cylinder_with_split_boundaries
                (
                triangulation,
                sizes[0], sizes[1], sizes[2]
                );
            manifold_ids.push_back(0);
            manifold_descriptors.push_back("cylindrical");
        }
        else if (grid_name == "hyper_cube_with_cylindrical_hole") 
        {
            assert(dim==2);
            GridGenerator::hyper_cube_with_cylindrical_hole
               (
                triangulation,
                sizes[0], sizes[1]
                );
            GridTools::copy_boundary_to_manifold_id(triangulation);
            manifold_ids.push_back(1);
            manifold_descriptors.push_back("spherical");  
        }
        else
        {
            throw(ExcNotImplemented());
        }
    }
 
}

#endif
