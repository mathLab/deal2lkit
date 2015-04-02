#include "parsed_grid_generator.h"
#include <deal.II/grid/grid_generator.h>

template <int dim, int spacedim>
ParsedGridGenerator<dim, spacedim>::ParsedGridGenerator(std::string name) :
    ParameterAcceptor(name)
{}

template <int dim, int spacedim>
std::string ParsedGridGenerator<dim, spacedim>::create_default_value(Point<spacedim> input)
{
  std::ostringstream strs;
  strs << input[0];
  for(unsigned int i=1; i<spacedim; ++i)
    strs<< ","<< input[i];
  std::string def = strs.str();
  return def;

}

template <int dim, int spacedim>
std::string ParsedGridGenerator<dim, spacedim>::create_default_value(std::vector<double> input)
{
  std::ostringstream strs;
  strs << input[0];
  for(unsigned int i=1; i<input.size(); ++i)
    strs<< ","<< input[i];
  std::string def = strs.str();
  return def;

}

template <int dim, int spacedim>
std::string ParsedGridGenerator<dim, spacedim>::create_default_value(std::vector<unsigned int> input)
{
  std::string def = Utilities::int_to_string(input[0]);
  for(unsigned int i=1; i<input.size(); ++i)
    def += "," + Utilities::int_to_string(input[i]);
  return def;
}

template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{

    std::vector<unsigned int> dummy_vec_int(spacedim);
    std::vector<double> dummy_vec_double(spacedim);
    Point<spacedim> dummy_point;
    std::string def_point, def_int, def_double;
    def_point = create_default_value(dummy_point);
    def_int = create_default_value(dummy_vec_int);
    def_double = create_default_value(dummy_vec_double);

    add_parameter(prm, &grid_name,
                  "Grid to generate", "unit_hypercube",
                  Patterns::Selection("file|unit_hypercube|hypercube|subhypercube|hyperrectangle"),//|unit_hyperball|unit_hypershell|subhyperrectangle"),
                  "The grid to generate. You can choose among\n"
                  " file: read grid from a file (use Input grid filename for the filename)\n"
                  " unit_hypercube: create a unit hypercube\n"
                  " unit_hypershell: create a unit hypersphere\n"
                  " hypercube: create a unit hypercube using the additional input\n"
                  );

    add_parameter(prm, &mesh_smoothing,
                  "Mesh smoothing alogrithm", "none",
                  Patterns::MultipleSelection("none|"
                          "limit_level_difference_at_vertices|"
                          "eliminate_unrefined_islands|"
                          "patch_level_1|"
                          "coarsest_level_1|"
                          "allow_anisotropic_smoothing|"
                          "eliminate_refined_inner_islands|"
                          "eliminate_refined_boundary_islands|"
                          "do_not_produce_unrefined_islands|"
                          "smoothing_on_refinement|"
                          "smoothing_on_coarsening|"
                          "maximum_smoothing"));

    add_parameter(prm, &grid_file_name,
                  "Grid file name", "",
                  Patterns::FileName(),
                  "Name of the input grid. All supported deal.II formats. "
                  "The extestion will be used to decide what "
                  "grid format to use.");

    add_parameter(prm, &double_option_one,
                  "First additional double input for the grid", "1.",
                  Patterns::Double(),
                  "First additional double to be used in the generation of the grid. "
                  "The use of it will depend on the specific grid.");

    add_parameter(prm, &double_option_two,
                  "Second additional double input for the grid", "0.5",
                  Patterns::Double(),
                  "Second additional double to be used in the generation of the grid. "
                  "The use of it will depend on the specific grid.");


    add_parameter(prm, &point_option_one,
                  "First additional Point<spacedim> input for the grid", def_point,
                  Patterns::List(Patterns::Double(), spacedim, spacedim),
                  "First additional Point<spacedim> to be used in the generation of the grid. "
                  "The use of it will depend on the specific grid.");

    add_parameter(prm, &point_option_two,
                  "Second additional Point<spacedim> input for the grid", def_point,
                  Patterns::List(Patterns::Double(), spacedim, spacedim),
                  "Second additional Point<spacedim> to be used in the generation of the grid. "
                  "The use of it will depend on the specific grid.");

    add_parameter(prm, &un_int_option_one,
                  "Unsigned int input for the grid", "1",
                  Patterns::Integer(),
                  "Unsigned int to be used in the generation of the grid. "
                  "The use of it will depend on the specific grid.");


}


template <int dim, int spacedim>
parallel::distributed::Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::distributed(MPI_Comm comm) {
    Assert(grid_name != "", ExcNotInitialized());
    auto tria = new parallel::distributed::Triangulation<dim,spacedim>
      (comm);//, get_smoothing());

      create(*tria);
    return tria;
}
// triangulation (mpi_communicator,
//                    typename Triangulation<dim>::MeshSmoothing
//                    (Triangulation<dim>::smoothing_on_refinement |
//                     Triangulation<dim>::smoothing_on_coarsening)),

template <int dim, int spacedim>
Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::serial() {
    Assert(grid_name != "", ExcNotInitialized());
    auto tria = new Triangulation<dim,spacedim>(get_smoothing());
    create(*tria);
    return tria;
}

template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::create(Triangulation<dim,spacedim> &tria)
{
    Assert(grid_name != "", ExcNotInitialized());
    if(grid_name == "unit_hypercube") {
      GridGenerator::hyper_cube(tria);
    }
    else if(grid_name == "hypercube")
    {
      Point<spacedim> center;

      if(double_option_one > double_option_two)
        GridGenerator::hyper_cube (tria,
                                    double_option_two, double_option_one);
      else
        GridGenerator::hyper_cube (tria,
                                    double_option_one, double_option_two);


    }
    else if(grid_name == "subhypercube"){

      if(double_option_one > double_option_two)
      {
        GridGenerator::hyper_cube (tria,
                                    double_option_two, double_option_one);
        tria.refine_global( un_int_option_one);
      }
      else
      {
        GridGenerator::hyper_cube (tria,
                                    double_option_one, double_option_two);
        tria.refine_global( un_int_option_one);

      }
    }
    else if(grid_name == "hyperrectangle")
    {
      GridGenerator::hyper_rectangle (tria,
                                       point_option_two, point_option_one);
    }
    // TO BE DONE WHEN POSSIBLE
    // else if(grid_name == "unit_hyperball")
    // {
    //   Point<spacedim> center;
    //   double radius=1.;
    //
    //   GridGenerator::hyper_ball (tria,
    //                               center, radius);
    //   }
    // else if(grid_name == "subhyperrectangle"){
    //   if(dim == spacedim)
    //   {
    //
    //   std::vector<unsigned int> refinement(spacedim);
    //   for(unsigned int i=0; i<spacedim;++i)
    //     refinement[i] = un_int_option_one;
    //   GridGenerator::subdivided_hyper_rectangle (tria, refinement,
    //                                    point_option_two, point_option_one);
    //   }
    //   else
    //     Assert(true, ExcInternalError("dim != spacedim not supported"))
    //
    // }


}

namespace {
template <int dim, int spacedim>
inline typename Triangulation<dim,spacedim>::MeshSmoothing operator|
(typename Triangulation<dim,spacedim>::MeshSmoothing a,
 typename Triangulation<dim,spacedim>::MeshSmoothing b) {
    return static_cast<typename Triangulation<dim,spacedim>::MeshSmoothing>
           (static_cast<int>(a)|static_cast<int>(b));
}
}

template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::MeshSmoothing
ParsedGridGenerator<dim, spacedim>::get_smoothing() {
    int smoothing = 0;
    std::vector<std::string> ss = Utilities::split_string_list(mesh_smoothing);
    for(auto s : ss) {
        if(s == "none")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::none;
        else if(s == "limit_level_difference_at_vertices")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::limit_level_difference_at_vertices;
        else if(s == "eliminate_unrefined_islands")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::eliminate_unrefined_islands;
        else if(s == "patch_level_1")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::patch_level_1;
        else if(s == "coarsest_level_1")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::coarsest_level_1;
        else if(s == "allow_anisotropic_smoothing")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::allow_anisotropic_smoothing;
        else if(s == "eliminate_refined_inner_islands")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::eliminate_refined_inner_islands;
        else if(s == "eliminate_refined_boundary_islands")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::eliminate_refined_boundary_islands;
        else if(s == "do_not_produce_unrefined_islands")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::do_not_produce_unrefined_islands;
        else if(s == "smoothing_on_refinement")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::smoothing_on_refinement;
        else if(s == "smoothing_on_coarsening")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::smoothing_on_coarsening;
        else if(s == "maximum_smoothing")
            smoothing = smoothing|(int)Triangulation<dim,spacedim>::maximum_smoothing;
        else
            Assert(false, ExcInternalError("Found a non supported smoothing. Should not happen:"+s))
        }
    return static_cast<typename Triangulation<dim,spacedim>::MeshSmoothing>(smoothing);
}



template class ParsedGridGenerator<1,1>;
template class ParsedGridGenerator<1,2>;
template class ParsedGridGenerator<1,3>;
template class ParsedGridGenerator<2,2>;
template class ParsedGridGenerator<2,3>;
template class ParsedGridGenerator<3,3>;
