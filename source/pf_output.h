#ifndef pf_output_h
#define pf_output_h

template<int dim>
void Phaseflow<dim>::write_solution()
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
