#ifndef lpcm_post_processor_templates_h
#define lpcm_post_processor_templates_h

#include "../include/post_processor.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

template <int dim, typename VECTOR>
PostProcessor<dim,VECTOR>::PostProcessor (const unsigned int n_mpi, 
					  const unsigned int this_mpi) :
  n_mpi_processes(n_mpi),
  this_mpi_process(this_mpi),
  data_out(this_mpi)
{
  initialized = false;
}

template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::declare_parameters (ParameterHandler &prm, 
						    unsigned int ntables)
{
  prm.enter_subsection ("I/O Options");

  prm.declare_entry ("Problem base name", "out/test", Patterns::Anything());

  prm.leave_subsection();

  prm.enter_subsection("Post Processing");
  prm.declare_entry ("Write error files", "false", Patterns::Bool());
  prm.declare_entry ("Output error tables", "true", Patterns::Bool());
  prm.declare_entry ("Output partitioning", "false", Patterns::Bool());
  prm.declare_entry ("Error file format", "tex", Patterns::Selection("tex|txt"));  
  prm.declare_entry ("Compute error", "true", Patterns::Bool());
  prm.declare_entry ("Table names", "error", Patterns::Anything(),
		     "Comma separated list of table names. ");
  prm.declare_entry ("Solution names", "u", Patterns::Anything(), 
		     "Comma separated list of names for the components. This "
		     "will be used both for error tables in text format and to "
		     "output the solution to a file. Note that in the case "
		     "of a vector function the error name which is used to "
		     "compute the norm (supposing the type of the other "
		     "components is 'AddUp') is the first one.");
  prm.declare_entry ("Solution names for latex", "u", Patterns::Anything(),
		     "Comma separated version of the same thing as above for "
		     "the latex version of the table.");

  // prm.declare_entry ("Ib output format", "msh", Patterns::Selection("raw|msh"));
  // prm.declare_entry ("Ib input file prefix", "ellipse", Patterns::Anything());
  
  prm.enter_subsection("Solution Out Format");
  DataOut<dim>::declare_parameters(prm);
  prm.set("Output format", "vtk");
  prm.leave_subsection();

  for(unsigned int i=0; i<ntables; ++i) {
    char tmp[10];
    sprintf(tmp, "Table %d", i);
    prm.enter_subsection(tmp);

    prm.declare_entry("List of error norms to compute", "Linfty, L2, W1infty, H1",
		      Patterns::Anything(), "Each component is separated by a semicolon, "
		      "and each norm by a comma. Implemented norms are Linfty, L2, W1infty, "
		      "H1 and AddUp, which means that the norm is added to the previous "
		      "component. Useful for vector valued functions.");
    prm.declare_entry("Add convergence rates", "true", Patterns::Bool(),
		      "Evaluate convergence rates and add a column to the table for each "
		      "computed norm. ");
    prm.declare_entry("Latex table caption", "error", Patterns::Anything(),
		      "The caption that will go under the table if we write the file in "
		      "tex format. The default value for this object is the same name "
		      "as the table name.");
    prm.declare_entry("Extra terms", "cells,dofs", 
		      Patterns::Anything(),
		      "The extra columns to add to the table.");

    prm.leave_subsection();
  }
  prm.leave_subsection();
}
template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::parse_parameters (ParameterHandler &prm)
{
  prm.enter_subsection ("I/O Options");

  base_name = prm.get ("Problem base name");

  prm.leave_subsection();

  prm.enter_subsection("Post Processing");

  write_error = prm.get_bool ("Write error files");
  output_error = prm.get_bool ("Output error tables");
  output_partitioning = prm.get_bool ("Output partitioning");

  error_file_format = prm.get ("Error file format");  
  compute_error = prm.get_bool ("Compute error");
  std::string all_names = prm.get ("Table names");
  headers = Utilities::split_string_list(prm.get ("Solution names"));
  latex_headers = Utilities::split_string_list(prm.get ("Solution names for latex"));

  prm.enter_subsection("Solution Out Format");
  data_out.parse_parameters(prm);
  prm.leave_subsection();
    
  if (all_names != "") {
    names = Utilities::split_string_list(all_names);
    types.resize(names.size(), std::vector<NormFlags> (headers.size()));
    add_rates.resize(names.size());
    tables.resize(names.size());
    latex_captions.resize(names.size());
    std::map<std::string, bool> extra;
    extra["dof"] = false;
    extra["cells"] = false;
    extra["dt"] = false;
	
    extras.resize(names.size(), extra);

    for(unsigned int i=0; i<names.size(); ++i) {
      char tmp[10];
      sprintf(tmp, "Table %d", i);
      prm.enter_subsection(tmp);

      all_names = prm.get("List of error norms to compute");
      add_rates[i] = prm.get_bool("Add convergence rates");
      latex_captions[i] = prm.get("Latex table caption");
      std::vector<std::string> all_extras = 
	Utilities::split_string_list(prm.get("Extra terms"));

      for(unsigned int x=0; x< all_extras.size(); ++x)
	extras[i][all_extras[x]] = true;
	    
      prm.leave_subsection();

      std::vector<std::string> all_comps = Utilities::split_string_list(all_names, ';');
      // Check that the input string has all the needed fields
      AssertThrow(all_comps.size() == headers.size(),
		  ExcDimensionMismatch(all_comps.size() , headers.size()));

      for(unsigned int j=0; j<all_comps.size(); ++j) {
	std::vector<std::string> all_types = 
	  Utilities::split_string_list(all_comps[j]);
	for(unsigned int k=0; k<all_types.size(); ++k) {
	  if(all_types[k] == "Linfty") {
	    types[i][j] |= Linfty;
	  } else if(all_types[k] == "L2") {
	    types[i][j] |= L2;
	  } else if(all_types[k] == "W1infty") {
	    types[i][j] |= W1infty;
	  } else if(all_types[k] == "H1") {
	    types[i][j] |= H1;
	  } else if(all_types[k] == "AddUp") {
	    types[i][j] |= AddUp;
	  } else {
	    AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
	  }
	}
      }
    }
  }
  prm.leave_subsection();
  initialized = true;
}


template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::prepare_data_output(const DoFHandler<dim> &dh,
						    const std::string &filename) {
  AssertThrow(initialized, ExcNotInitialized());
  deallog.push("PrepareOutput");
  
  if(data_out.default_suffix() != "") {
    
    std::string fname = base_name + "_" + filename;
    
    // If the output is needed and we have many processes, just output
    // the one we need *in intermediate format*.
    if(n_mpi_processes > 1) {
      fname += ("_" + Utilities::int_to_string(this_mpi_process, 2) + 
		data_out.default_suffix(DataOutBase::deal_II_intermediate)) ;
    } else {
      fname += data_out.default_suffix();
    }
    
    deallog << "Will write on file: " << fname.c_str() << std::endl;
    output_file.open(fname.c_str());
    AssertThrow(output_file, ExcIO());
    data_out.attach_dof_handler (dh);
    
    if(n_mpi_processes > 1) {
      // Output the partitioning
      if(output_partitioning) {
	std::vector<unsigned int> partition_int (dh.get_tria().n_active_cells());
	GridTools::get_subdomain_association (dh.get_tria(), partition_int);
	Vector<double> partitioning(partition_int.begin(),
				    partition_int.end());
	static Vector<double> static_partitioning;
	static_partitioning.swap(partitioning);
	data_out.add_data_vector (static_partitioning, "partitioning");
      }
    }
  }
  deallog.pop();
}


template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::add_data_vector(const VECTOR &data_vector, 
						const std::string &desc)
{
  AssertThrow(initialized, ExcNotInitialized());
  deallog.push("AddingData");
  std::vector<std::string> dd = Utilities::split_string_list(desc);
  if(data_out.default_suffix() != "") {
    if (dd.size() ==1 ) 
      data_out.add_data_vector (data_vector, desc);
    else 
      data_out.add_data_vector (data_vector, dd);
    deallog << "Added data: " << desc << std::endl;
  }
  deallog.pop();
}


template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::write_data_and_clear() {
  AssertThrow(initialized, ExcNotInitialized());
  AssertThrow(output_file, ExcIO());
  deallog.push("WritingData");
  if(data_out.default_suffix() != "") {
    data_out.build_patches();
    
    if(n_mpi_processes > 1) {
      data_out.write_deal_II_intermediate(output_file);
    } else {
      data_out.write(output_file);
    }
    deallog << "Wrote output file." << std::endl;
    data_out.clear();
    output_file.close();
    deallog << "Reset output." << std::endl;
  }
  deallog.pop();
}

template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::output_solution(
						const VectorSpace<dim> &,
						const VECTOR &, 
						const std::string &)
{
  AssertThrow(false, ExcMessage("Deprecated Function"));
}

template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::dump_vector (const VECTOR &, 
					     const std::string & )
{
  Assert(false, ExcNotImplemented());
  // Only specializations exist.
}

template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::output_table (const unsigned int table_no) {
  if (compute_error && this_mpi_process == 0) {
    AssertThrow(initialized, ExcNotInitialized());
    AssertThrow(table_no < names.size(), ExcIndexRange(table_no, 0, names.size()));

    // Add convergence rates
    if(add_rates[table_no]) {
      if(extras[table_no]["dofs"])
	tables[table_no].omit_column_from_convergence_rate_evaluation("dofs");	
      if(extras[table_no]["cells"])
	tables[table_no].omit_column_from_convergence_rate_evaluation("cells");
      if(extras[table_no]["dt"])
	tables[table_no].omit_column_from_convergence_rate_evaluation("dt");
      tables[table_no].evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
      //tables[table_no].evaluate_all_convergence_rates("dofs", ConvergenceTable::reduction_rate_log2);
    }

    if(output_error) tables[table_no].write_text(std::cout);

    if(write_error) {
      std::string filename = base_name + "_" + names[table_no] +
	"." + error_file_format;

      std::ofstream table_file(filename.c_str());

      if(error_file_format != "txt")
	tables[table_no].write_tex(table_file);
      else 
	tables[table_no].write_text(table_file);
      table_file.close();
    }
  }
}

template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::difference(const VectorSpace<dim> & vspace, 
					   const VECTOR &solution1,
					   const VECTOR &solution2, 
					   unsigned int table_no,
					   double dt) {
  AssertThrow(solution1.size() == solution2.size(), ExcDimensionMismatch(
									 solution1.size(), solution2.size()));
  VECTOR solution(solution1);
  solution -= solution2;
  error_from_exact(vspace, solution, 
		   ConstantFunction<dim>(0, headers.size()), table_no, dt);
}


template <int dim, typename VECTOR>
void PostProcessor<dim,VECTOR>::error_from_exact(const VectorSpace<dim> & vspace, 
						 const VECTOR &solution, 
						 const Function<dim> &exact,
						 unsigned int table_no,
						 double dt) 
{
  if (compute_error) {
    AssertThrow(initialized, ExcNotInitialized());
    AssertThrow(table_no < types.size(), ExcIndexRange(table_no, 0, names.size()));
    AssertThrow(exact.n_components == types[table_no].size(), 
		ExcDimensionMismatch(exact.n_components, types[table_no].size()));

    deallog.push("PostProcessing");
    deallog << "Calculating Errors." << std::endl;
    std::vector< std::vector<double> > error( exact.n_components, std::vector<double>(4));
    const unsigned int n_active_cells = vspace.get_tria().n_active_cells();
    const unsigned int n_dofs=vspace.get_dh().n_dofs();

    if(extras[table_no]["cells"]) {
      tables[table_no].add_value("cells", n_active_cells);
      tables[table_no].set_tex_caption("cells", "\\# cells");
      tables[table_no].set_tex_format("cells", "r");
    }
    if(extras[table_no]["dofs"]) {
      tables[table_no].add_value("dofs", n_dofs);
      tables[table_no].set_tex_caption("dofs", "\\# dofs");
      tables[table_no].set_tex_format("dofs", "r");
    }
    if(extras[table_no]["dt"]) {
      tables[table_no].add_value("dt", dt);
      tables[table_no].set_tex_caption("dt", "\\Delta t");
      tables[table_no].set_tex_format("dt", "r");
    }

    bool compute_Linfty = false;
    bool compute_L2 = false;
    bool compute_W1infty = false;
    bool compute_H1 = false;
    bool add_this = false;

    unsigned int last_non_add = 0;

    for(unsigned int component=0; component < exact.n_components; ++component) {
      NormFlags norm = types[table_no][component];

      deallog << "Error flags: " << norm << std::endl;

      // Select one Component
      ComponentSelectFunction<dim> select_component ( component, 1. , exact.n_components);

      Vector<float> difference_per_cell (vspace.get_tria().n_active_cells());

      QGauss<dim> q_gauss((vspace.get_fe().degree+1) * 2);

      // The add bit is set 
      add_this = (norm & AddUp);

      if(!add_this) {
	last_non_add	= component;
	compute_L2	= ( norm & L2 );
	compute_H1	= ( norm & H1 );
	compute_W1infty = ( norm & W1infty ) ;
	compute_Linfty	= ( norm & Linfty );
      }
      // if add is set, we do not modify the previous selection

      if(compute_L2) {
	VectorTools::integrate_difference (//mapping,
					   vspace.get_dh(), //dof_handler,
					   solution, 
					   exact,
					   difference_per_cell,
					   q_gauss,
					   VectorTools::L2_norm,
					   &select_component );
      }

      const double L2_error = difference_per_cell.l2_norm();
      difference_per_cell = 0;
	   
      if(compute_H1) { 
	VectorTools::integrate_difference (//mapping,
					   vspace.get_dh(), //dof_handler,
					   solution,
					   exact,
					   difference_per_cell,
					   q_gauss,
					   VectorTools::H1_norm,
					   &select_component );
      }
      const double H1_error = difference_per_cell.l2_norm();
      difference_per_cell = 0;

      if(compute_W1infty) {
	VectorTools::integrate_difference (//mapping,
					   vspace.get_dh(), //dof_handler,
					   solution,
					   exact,
					   difference_per_cell,
					   q_gauss,
					   VectorTools::W1infty_norm,
					   &select_component );
      }

      const double W1inf_error = difference_per_cell.linfty_norm();

      if(compute_Linfty) {
	VectorTools::integrate_difference (//mapping,
					   vspace.get_dh(), //dof_handler,
					   solution,
					   exact,
					   difference_per_cell,
					   q_gauss,
					   VectorTools::Linfty_norm,
					   &select_component );
      }

      const double Linf_error = difference_per_cell.linfty_norm();

      if(add_this) {
	AssertThrow(component, ExcMessage("Cannot add on first component!"));

	error[last_non_add][0] = std::max(error[last_non_add][0], Linf_error);
	error[last_non_add][1] += L2_error;
	error[last_non_add][2] = std::max(error[last_non_add][2], W1inf_error);
	error[last_non_add][3] += H1_error;

      } else {

	error[component][0] = Linf_error;
	error[component][1] = L2_error;
	error[component][2] = W1inf_error;
	error[component][3] = H1_error;

      }
    }  

    for(unsigned int j=0; j<exact.n_components; ++j) {
      NormFlags norm = types[table_no][j];
      // If this was added, don't do anything
      if( !(norm & AddUp) ) {
	if(norm & Linfty) {
	  std::string name = headers[j] + "_Linfty";
	  std::string latex_name = "$\\| " + 
	    latex_headers[j] + " - " +
	    latex_headers[j] +"_h \\|_\\infty $";
	  double this_error =  error[j][0];

	  tables[table_no].add_value(name, this_error);
	  tables[table_no].set_precision(name, 3);
	  tables[table_no].set_scientific(name, true);
	  tables[table_no].set_tex_caption(name, latex_name);
	}

	if(norm & L2) {
	  std::string name = headers[j] + "_L2";
	  std::string latex_name = "$\\| " + 
	    latex_headers[j] + " - " +
	    latex_headers[j] +"_h \\|_0 $";
	  double this_error =  error[j][1];

	  tables[table_no].add_value(name, this_error);
	  tables[table_no].set_precision(name, 3);
	  tables[table_no].set_scientific(name, true);
	  tables[table_no].set_tex_caption(name, latex_name);
	}
	if(norm & W1infty) {
	  std::string name = headers[j] + "_W1infty";
	  std::string latex_name = "$\\| " + 
	    latex_headers[j] + " - " +
	    latex_headers[j] +"_h \\|_{1,\\infty} $";
	  double this_error =  error[j][2];

	  tables[table_no].add_value(name, this_error);
	  tables[table_no].set_precision(name, 3);
	  tables[table_no].set_scientific(name, true);
	  tables[table_no].set_tex_caption(name, latex_name);
	}
	if(norm & H1){
	  std::string name = headers[j] + "_H1";
	  std::string latex_name = "$\\| " + 
	    latex_headers[j] + " - " +
	    latex_headers[j] +"_h \\|_1 $";
	  double this_error =  error[j][3];

	  tables[table_no].add_value(name, this_error);
	  tables[table_no].set_precision(name, 3);
	  tables[table_no].set_scientific(name, true);
	  tables[table_no].set_tex_caption(name, latex_name);
	}
      }
    }
    deallog.pop();
  }
}

#endif
