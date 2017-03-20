#include <deal.II/grid/tria.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>

#include <iostream>
#include <fstream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace FEFieldTools
{
    using namespace dealii;
    
    template<int dim>
    void save_field_parts(
        Triangulation<dim> &field_tria,
        DoFHandler<dim> &field_dof_handler,
        Vector<double> &field_solution)
    {
        {
            std::string file_path = "field_solution";
            std::ofstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_oarchive archive(file_stream);
            archive << field_solution;
        }
        {
            std::string file_path = "field_dof_handler";
            std::ofstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_oarchive archive(file_stream);
            archive << field_dof_handler;
        }
        {
            std::string file_path = "field_triangulation";
            std::ofstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_oarchive archive(file_stream);
            archive << field_tria;
        }
    }

    template<int dim>
    void load_field_parts(
        Triangulation<dim> &field_tria,
        DoFHandler<dim> &field_dof_handler,
        Vector<double> &field_solution,
        FE_Q<dim> &fe)
    {
        {
            std::string file_path = "field_triangulation";
            std::ifstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_iarchive archive(file_stream);
            archive >> field_tria;
        }
        field_dof_handler.distribute_dofs(fe);
        {
            std::string file_path = "field_dof_handler";
            std::ifstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_iarchive archive(file_stream);
            archive >> field_dof_handler;
        }
        {
            std::string file_path = "field_solution";
            std::ifstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_iarchive archive(file_stream);
            archive >> field_solution;
        }
    }

    template<int dim>
    void load_field_parts(
        Triangulation<dim> &field_tria,
        DoFHandler<dim> &field_dof_handler,
        Vector<double> &field_solution,
        FESystem<dim,dim> &fe)
    {
        {
            std::string file_path = "field_triangulation";
            std::ifstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_iarchive archive(file_stream);
            archive >> field_tria;
        }
        field_dof_handler.distribute_dofs(fe);
        {
            std::string file_path = "field_dof_handler";
            std::ifstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_iarchive archive(file_stream);
            archive >> field_dof_handler;
        }
        {
            std::string file_path = "field_solution";
            std::ifstream file_stream(file_path, std::ios::binary);
            if (!file_stream.good())
            {
                throw std::runtime_error("Error while opening the file: " + file_path);
            }
            boost::archive::binary_iarchive archive(file_stream);
            archive >> field_solution;
        }
    }
}
  