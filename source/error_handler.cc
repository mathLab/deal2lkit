#include "../include/error_handler.h"

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

template <int ntables>
ErrorHandler<ntables>::ErrorHandler (const std::string name) :
  ParameterAcceptor(name)
{
  initialized = false;
}

template <int ntables>
void ErrorHandler<ntables>::declare_parameters (ParameterHandler &prm)
{
  prm.declare_entry ("Write error files", "false", Patterns::Bool());
  prm.declare_entry ("Output error tables", "true", Patterns::Bool());
  prm.declare_entry ("Error file format", "tex", Patterns::Selection("tex|txt|gpl|org"));
  prm.declare_entry ("Compute error", "true", Patterns::Bool());
  prm.declare_entry ("Table names", "error", Patterns::List(Patterns::Anything(), 1, ntables),
                     "Comma separated list of table names. ");
  prm.declare_entry ("Solution names", "u", Patterns::Anything(),
                     "Comma separated list of names for the components. This "
                     "will be used both for error tables in text format and to "
                     "output the solution to a file. Note that in the case "
                     "of a vector function the error name which is used to "
                     "compute the norm (supposing the type of the other "
                     "components is 'Add') is the first one.");
  prm.declare_entry ("Solution names for latex", "u", Patterns::Anything(),
                     "Comma separated version of the same thing as above for "
                     "the latex version of the table.");

  // prm.declare_entry ("Ib output format", "msh", Patterns::Selection("raw|msh"));
  // prm.declare_entry ("Ib input file prefix", "ellipse", Patterns::Anything());
  for (unsigned int i=0; i<ntables; ++i)
    {
      char tmp[10];
      sprintf(tmp, "Table %d", i);
      prm.enter_subsection(tmp);

      prm.declare_entry("List of error norms to compute", "Linfty, L2, H1",
                        Patterns::Anything(), "Each component is separated by a semicolon, "
                        "and each norm by a comma. Implemented norms are Linfty, L2, "
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
      prm.declare_entry("Rate key", "",
                        Patterns::Selection("dofs|cells|dt|"),
                        "The key to use to compute the convergence rates.");
      prm.leave_subsection();
    }
}

template <int ntables>
void ErrorHandler<ntables>::parse_parameters (ParameterHandler &prm)
{
  write_error = prm.get_bool ("Write error files");
  output_error = prm.get_bool ("Output error tables");

  error_file_format = prm.get ("Error file format");
  compute_error = prm.get_bool ("Compute error");
  std::string all_names = prm.get ("Table names");
  headers = Utilities::split_string_list(prm.get ("Solution names"));
  latex_headers = Utilities::split_string_list(prm.get ("Solution names for latex"));

  if (all_names != "")
    {
      names = Utilities::split_string_list(all_names);
      Assert(names.size() <= ntables,
             ExcMessage("You tried to construct more tables than you have compiled for."));
      types.resize(names.size(), std::vector<NormFlags> (headers.size()));
      add_rates.resize(names.size());
      tables.resize(names.size());
      latex_captions.resize(names.size());
      std::map<std::string, bool> extra;
      extra["dof"] = false;
      extra["cells"] = false;
      extra["dt"] = false;

      extras.resize(names.size(), extra);
      rate_keys.resize(names.size(), "");

      for (unsigned int i=0; i<names.size(); ++i)
        {
          char tmp[10];
          sprintf(tmp, "Table %d", i);
          prm.enter_subsection(tmp);

          all_names = prm.get("List of error norms to compute");
          add_rates[i] = prm.get_bool("Add convergence rates");
          rate_keys[i] = prm.get("Rate key");
          latex_captions[i] = prm.get("Latex table caption");
          std::vector<std::string> all_extras =
            Utilities::split_string_list(prm.get("Extra terms"));

          for (unsigned int x=0; x< all_extras.size(); ++x)
            extras[i][all_extras[x]] = true;

          prm.leave_subsection();

          std::vector<std::string> all_comps = Utilities::split_string_list(all_names, ';');
          // Check that the input string has all the needed fields
          AssertThrow(all_comps.size() == headers.size(),
                      ExcDimensionMismatch(all_comps.size() , headers.size()));

          for (unsigned int j=0; j<all_comps.size(); ++j)
            {
              std::vector<std::string> all_types =
                Utilities::split_string_list(all_comps[j]);
              for (unsigned int k=0; k<all_types.size(); ++k)
                {
                  if (all_types[k] == "Linfty")
                    {
                      types[i][j] |= Linfty;
                    }
                  else if (all_types[k] == "L2")
                    {
                      types[i][j] |= L2;
                    }
                  else if (all_types[k] == "W1infty")
                    {
                      types[i][j] |= W1infty;
                    }
                  else if (all_types[k] == "H1")
                    {
                      types[i][j] |= H1;
                    }
                  else if (all_types[k] == "AddUp")
                    {
                      types[i][j] |= AddUp;
                    }
                  else
                    {
                      AssertThrow(false, ExcMessage("Didn't recognize a norm type."));
                    }
                }
            }
        }
    }
  initialized = true;
}

template <int ntables>
void ErrorHandler<ntables>::output_table (std::ostream &out, const unsigned int table_no)
{
  if (compute_error)
    {
      AssertThrow(initialized, ExcNotInitialized());
      AssertThrow(table_no < names.size(), ExcIndexRange(table_no, 0, names.size()));

      // Add convergence rates
      if (add_rates[table_no])
        {
          if (extras[table_no]["dofs"])
            tables[table_no].omit_column_from_convergence_rate_evaluation("dofs");
          if (extras[table_no]["cells"])
            tables[table_no].omit_column_from_convergence_rate_evaluation("cells");
          if (extras[table_no]["dt"])
            tables[table_no].omit_column_from_convergence_rate_evaluation("dt");
          if (rate_keys[table_no] == "")
            tables[table_no].evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
          else
            tables[table_no].evaluate_all_convergence_rates(rate_keys[table_no], ConvergenceTable::reduction_rate_log2);
        }

      if (output_error) tables[table_no].write_text(out);

      if (write_error)
        {
          std::string filename = names[table_no] +
                                 "." + error_file_format;

          std::ofstream table_file(filename.c_str());

          if (error_file_format == "tex")
            tables[table_no].write_tex(table_file);
          else if (error_file_format == "txt")
            tables[table_no].write_text(table_file);
          else if (error_file_format == "gpl")
            tables[table_no].write_text(table_file,
                                        TableHandler::table_with_separate_column_description);
          else if (error_file_format == "org")
            tables[table_no].write_text(table_file,
                                        TableHandler::org_mode_table);
          table_file.close();
        }
    }
}

template class ErrorHandler<1>;
template class ErrorHandler<2>;
template class ErrorHandler<3>;
template class ErrorHandler<4>;
