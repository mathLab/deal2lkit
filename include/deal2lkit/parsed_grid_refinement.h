#ifndef __dealii_sak_parsed_grid_refinement_h
#define __dealii_sak_parsed_grid_refinement_h

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>


#include <deal2lkit/config.h>
#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;


// I'm using std::to_string. I didn't bother changing this.
#ifdef DEAL_II_WITH_CXX11

/**
 * A wrapper for refinement strategies.
 */
class ParsedGridRefinement : public ParameterAcceptor
{
public:
  /**
   * Constructor.
   */
  ParsedGridRefinement(const std::string &name="",
                       const std::string &strategy="fraction",
                       const double &top_parameter=.3,
                       const double &bottom_parameter=.1,
                       const unsigned int &max_cells=std::numeric_limits< unsigned int >::max(),
                       const unsigned int &order=2);

  /**
   * Declare local parameters.
   */
  virtual void declare_parameters(ParameterHandler &prm);


  /**
   * Mark cells a the triangulation for refinement or coarsening,
   * according to the given strategy applied to the supplied vector
   * representing local error criteria.
   *
   * Cells are only marked for refinement or coarsening. No refinement
   * is actually performed. You need to call
   * Triangulation::execute_coarsening_and_refinement() yourself.
   */
  template<int dim, class Vector , int spacedim>
  void mark_cells(const Vector &criteria,
                  Triangulation< dim, spacedim > &tria) const;


  /**
   * Mark cells of a distribtued triangulation for refinement or
   * coarsening, according to the given strategy applied to the
   * supplied vector representing local error criteria. If the
   * criterion which is specified in the parameter file is not
   * available, an exception is thrown.
   *
   * Cells are only marked for refinement or coarsening. No refinement
   * is actually performed. You need to call
   * Triangulation::execute_coarsening_and_refinement() yourself.
   */
  template<int dim, class Vector , int spacedim>
  void mark_cells(const Vector &criteria,
                  parallel::distributed::Triangulation< dim, spacedim > &tria) const;

private:
  /**
   * Default expression of this function. "
   */
  std::string strategy;
  double top_parameter;
  double bottom_parameter;
  unsigned int max_cells;
  unsigned int order;
};



template<int dim, class Vector , int spacedim>
void ParsedGridRefinement::mark_cells(const Vector &criteria,
                                      parallel::distributed::Triangulation< dim, spacedim > &tria) const
{
  if (strategy == "number")
    GridRefinement::refine_and_coarsen_fixed_number (tria,
                                                     criteria,
                                                     top_parameter, bottom_parameter,
                                                     max_cells);
  else if (strategy == "fraction")
    GridRefinement::refine_and_coarsen_fixed_fraction (tria,
                                                       criteria,
                                                       top_parameter, bottom_parameter);
  else
    Assert(false, ExcInternalError());

}


template<int dim, class Vector , int spacedim>
void ParsedGridRefinement::mark_cells(const Vector &criteria,
                                      Triangulation< dim, spacedim > &tria) const
{

  if (strategy == "number")
    GridRefinement::refine_and_coarsen_fixed_number (tria,
                                                     criteria,
                                                     top_parameter, bottom_parameter,
                                                     max_cells);
  else if (strategy == "fraction")
    GridRefinement::refine_and_coarsen_fixed_fraction (tria,
                                                       criteria,
                                                       top_parameter, bottom_parameter,
                                                       max_cells);
  // This one does not seem to work properly
  // else if(strategy == "optimize")
  //   GridRefinement::refine_and_coarsen_optimize (tria, order);
  else
    Assert(false, ExcInternalError());
}


#endif

#endif
