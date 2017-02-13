#ifndef _extrapolated_field_h_
#define _extrapolated_field_h_

#include <deal.II/base/function.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/fe_field_function.h>

/**
 * @brief Extends the FEFieldFunction with nearest neighbor extrapolation.
 *
 * @detail
 *
 *    This is perhaps useful for interpolating the function onto a similar domain
 *    with non-conforming boundaries, e.g. if the domain has only been slightly shifted or rotated.
 *
 * @author Alexander Zimmerman 2016
*/
namespace MyFunctions
{
    using namespace dealii;

    template<int dim>
    class ExtrapolatedField : public Function<dim>
    {
    public:
        ExtrapolatedField(const DoFHandler<dim> &dof_handler, const Vector<double> &f)
          : Function<dim>(),
            field(dof_handler, f),
            dof_handler_sp(&dof_handler, "ExtrapolatedField")
        {}
        virtual double value(const Point<dim>  &point,
                             const unsigned int component = 0) const;
    private:
        Functions::FEFieldFunction<dim> field;
        SmartPointer<const DoFHandler<dim>,ExtrapolatedField<dim>> dof_handler_sp;
        Point<dim> get_nearest_boundary_vertex(const Point<dim> &point) const;
    };

    template<int dim>
    double ExtrapolatedField<dim>::value(const Point<dim> &point,
                                         const unsigned int component) const
    {
        Assert(component == 0, ExcInternalError());
        double val;
        try 
        {
            val = field.value(point, component);
        }
        catch (GridTools::ExcPointNotFound<dim>)
        {
            val = field.value(get_nearest_boundary_vertex(point));
        }
        return val;
    }

    template <int dim>
    Point<dim> ExtrapolatedField<dim>::get_nearest_boundary_vertex(const Point<dim> &point) const
    {
        double arbitrarily_large_number = 1.e32;
        double nearest_distance = arbitrarily_large_number;
        Point<dim> nearest_vertex;
        for (auto cell : dof_handler_sp->active_cell_iterators())
        {
        if (!cell->at_boundary())
            {
                continue;
            }
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            {
                if (!cell->face(f)->at_boundary())
                {
                    continue;
                }
                for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v)
                {
                    Point<dim> vertex = cell->face(f)->vertex(v);
                    double distance = (point - vertex).norm_square();
                    if (distance < nearest_distance)
                    {
                        nearest_vertex = vertex;
                        nearest_distance = distance;		
                    }
                }
            }
        }
        return nearest_vertex;
    }

} 

#endif
