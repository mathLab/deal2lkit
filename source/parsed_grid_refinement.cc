#include <deal2lkit/parsed_grid_refinement.h>


#ifdef DEAL_II_WITH_CXX11

ParsedGridRefinement::ParsedGridRefinement(const std::string &name,
                                           const std::string &strategy,
                                           const double &top_parameter,
                                           const double &bottom_parameter,
                                           const unsigned int &max_cells,
                                           const unsigned int &order) :
  ParameterAcceptor(name),
  strategy(strategy),
  top_parameter(top_parameter),
  bottom_parameter(bottom_parameter),
  max_cells(max_cells),
  order(order)
{}

void ParsedGridRefinement::declare_parameters(ParameterHandler &prm)
{
  add_parameter(prm, &strategy,
                "Refinement strategy", strategy,
                Patterns::Selection("fraction|number"),
                "Refinement strategy to use. fraction|number");

  add_parameter(prm, &top_parameter,
                "Top fraction", std::to_string(top_parameter),
                Patterns::Double(0.0),
                "Refinement fraction.");

  add_parameter(prm, &bottom_parameter,
                "Bottom fraction", std::to_string(bottom_parameter),
                Patterns::Double(0.0),
                "Coarsening fraction.");

  add_parameter(prm, &max_cells,
                "Maximum number of cells (if available)", std::to_string(max_cells),
                Patterns::Integer(0),
                "Maximum number of cells.");


  add_parameter(prm, &order,
                "Order (optimize)", std::to_string(order),
                Patterns::Integer(0),
                "Maximum number of cells.");

}

#endif
