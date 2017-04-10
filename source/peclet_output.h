#ifndef peclet_output_h
#define peclet_output_h

template<int dim>
void Peclet<dim>::write_solution()
{
  
    if (this->params.output.write_solution_vtk)
    {
        Output::write_solution_to_vtk(
            "solution-"+Utilities::int_to_string(this->time_step_counter)+".vtk",
            this->dof_handler,
            this->solution);    
    }

}

#endif
