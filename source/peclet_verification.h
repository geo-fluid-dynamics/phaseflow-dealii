#ifndef _peclet_verification_h_
#define _peclet_verification_h_

template<int dim>
void Peclet<dim>::append_verification_table()
{
    assert(this->params.verification.enabled);

    Vector<float> difference_per_cell(triangulation.n_active_cells());

    VectorTools::integrate_difference(
        this->dof_handler,
        this->solution,
        this->exact_solution_function,
        difference_per_cell,
        QGauss<dim>(dim + 1),
        VectorTools::L2_norm);
        
    double L2_norm_error = difference_per_cell.l2_norm();

    VectorTools::integrate_difference(
        this->dof_handler,
        this->solution,
        this->exact_solution_function,
        difference_per_cell,
        QGauss<dim>(dim + 1),
        VectorTools::L1_norm);
        
    double L1_norm_error = difference_per_cell.l1_norm();

    this->verification_table.add_value("cells", this->triangulation.n_active_cells());
    this->verification_table.add_value("dofs", this->dof_handler.n_dofs());
    this->verification_table.add_value("L1_norm_error", L1_norm_error);
    this->verification_table.add_value("L2_norm_error", L2_norm_error);

}

template<int dim>
void Peclet<dim>::write_verification_table()
{
    const int precision = 14;

    this->verification_table.set_precision("cells", precision);
    this->verification_table.set_scientific("cells", true);

    this->verification_table.set_precision("dofs", precision);
    this->verification_table.set_scientific("dofs", true);

    this->verification_table.set_precision("L2_norm_error", precision);
    this->verification_table.set_scientific("L2_norm_error", true);

    this->verification_table.set_precision("L1_norm_error", precision);
    this->verification_table.set_scientific("L1_norm_error", true);

    std::ofstream out_file(this->verification_table_file_name, std::fstream::app);
    assert(out_file.good());
    this->verification_table.write_text(out_file);
    out_file.close(); 
}
  
#endif
