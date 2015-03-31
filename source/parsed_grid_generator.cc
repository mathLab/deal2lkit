#include "parsed_grid_generator.h"
#include <deal.II/grid/grid_generator.h>

template <int dim, int spacedim>
ParsedGridGenerator<dim, spacedim>::ParsedGridGenerator(std::string name) :
    ParameterAcceptor(name)
{}

template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{
    add_parameter(prm, &grid_name,
                  "Grid to generate", "unit_hypercube",
                  Patterns::Selection("file|unit_hypercube"),
                  "The grid to generate. You can choose among\n"
                  " file: read grid from a file (use Input grid filename for the filename)\n"
                  " unit_hypercube: create a unit hypercube\n");

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
}

template <int dim, int spacedim>
Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::serial() {
    Assert(grid_name != "", ExcNotInitialized());
    auto tria = new Triangulation<dim,spacedim>(get_smoothing());

    if(grid_name == "unit_hypercube") {
      GridGenerator::hyper_cube(*tria);
    }
    
    return tria;
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
