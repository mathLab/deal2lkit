#include <deal.II/base/config.h>
#include "parsed_grid_generator.h"
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <fstream>

std::string extension(const std::string &filename)
{
  std::string::size_type idx;
  idx = filename.rfind('.');
  if (idx != std::string::npos)
    {
      return filename.substr(idx+1);
    }
  else
    {
      return "";
    }
}

template <int dim, int spacedim>
ParsedGridGenerator<dim, spacedim>::ParsedGridGenerator(std::string name) :
  ParameterAcceptor(name)
  //un_int_vec_option_one(dim)
{}

template <int dim, int spacedim>
std::string ParsedGridGenerator<dim, spacedim>::create_default_value(const Point<spacedim> &input)
{
  std::ostringstream strs;
  strs << input[0];
  for (unsigned int i=1; i<spacedim; ++i)
    strs<< ","<< input[i];
  std::string def = strs.str();
  return def;

}

template <int dim, int spacedim>
std::string ParsedGridGenerator<dim, spacedim>::create_default_value(const std::vector<double> &input)
{
  std::ostringstream strs;
  strs << input[0];
  for (unsigned int i=1; i<input.size(); ++i)
    strs<< ","<< input[i];
  std::string def = strs.str();
  return def;

}

template <int dim, int spacedim>
std::string ParsedGridGenerator<dim, spacedim>::create_default_value(const std::vector<unsigned int> &input)
{
  std::string def = Utilities::int_to_string(input[0]);
  for (unsigned int i=1; i<input.size(); ++i)
    def += "," + Utilities::int_to_string(input[i]);
  return def;
}

template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{

  std::vector<unsigned int> dummy_vec_int(dim, 1);
  std::vector<double> dummy_vec_double(spacedim);
  Point<spacedim> dummy_point;
  Point<spacedim> dummy_point_2;
  for(unsigned int d=0; d<dim; ++d)
    dummy_point_2[d]=1.0;
  
  std::string def_point, def_point_2, def_int, def_double;
  def_point = create_default_value(dummy_point);
  def_point_2 = create_default_value(dummy_point_2);
  def_int = create_default_value(dummy_vec_int);
  def_double = create_default_value(dummy_vec_double);

  add_parameter(prm, &grid_name,
                "Grid to generate", "rectangle",
                Patterns::Selection("file|rectangle"), //|unit_hyperball|unit_hypershell|subhyperrectangle"),
                "The grid to generate. You can choose among:\n"
		"- file: read grid from a file using:\n"
		"	- Input grid filename	    : input filename\n\n"
		"- rectangle: create a subdivided hyperrectangle using:\n"
		"	- Optional Point<spacedim> 1: left corner\n"
		"	- Optional Point<spacedim> 2: right corner\n"
		"	- Optional Vector of dim int: subdivisions on each direction\n"
		"	- Optional bool 1	    : colorize grid\n");

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

  add_parameter(prm, &input_grid_file_name,
                "Input grid file name", "",
                Patterns::FileName(),
                "Name of the input grid. All supported deal.II formats. "
                "The extestion will be used to decide what "
                "grid format to use.");

  add_parameter(prm, &double_option_one,
                "Optional double 1", "1.",
                Patterns::Double(),
                "First additional double to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &double_option_two,
                "Optional double 2", "0.5",
                Patterns::Double(),
                "Second additional double to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");


  add_parameter(prm, &point_option_one,
                "Optional Point<spacedim> 1", def_point,
                Patterns::List(Patterns::Double(), spacedim, spacedim),
                "First additional Point<spacedim> to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &point_option_two,
                "Optional Point<spacedim> 2", def_point_2,
                Patterns::List(Patterns::Double(), spacedim, spacedim),
                "Second additional Point<spacedim> to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");


  add_parameter(prm, &un_int_option_one,
                "Optional int 1","1",
                Patterns::Integer(),
                "Unsigned int to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &un_int_vec_option_one,
                "Optional vector of dim int", def_int,
                Patterns::List(Patterns::Integer(1), dim, dim),
                "Vector of positive unsigned int to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &bool_option_one,
                "Optional bool 1","false",
                Patterns::Bool(),
                "Bool be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");


  add_parameter(prm, &output_grid_file_name,
                "Output grid file name", "",
                Patterns::FileName(),
                "Name of the output grid. All supported deal.II formats. "
                "The extestion will be used to decide what "
                "grid format to use. If empty, no grid will be written.");

}

#ifdef DEAL_II_WITH_MPI
template <int dim, int spacedim>
parallel::distributed::Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::distributed(MPI_Comm comm)
{
  Assert(grid_name != "", ExcNotInitialized());
  auto tria = new parallel::distributed::Triangulation<dim,spacedim>
  (comm);//, get_smoothing());

  create(*tria);
  return tria;
}
#endif


template <int dim, int spacedim>
Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::serial()
{
  Assert(grid_name != "", ExcNotInitialized());
  auto tria = new Triangulation<dim,spacedim>(get_smoothing());
  create(*tria);
  return tria;
}

template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::create(Triangulation<dim,spacedim> &tria)
{
  Assert(grid_name != "", ExcNotInitialized());
  if(grid_name == "rectangle")
    {
      GridGenerator::subdivided_hyper_rectangle (tria, un_int_vec_option_one,
						 point_option_two, point_option_one,
						 bool_option_one);
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
  else if (grid_name == "file")
    {
      GridIn<dim, spacedim> gi;
      gi.attach_triangulation(tria);

      std::ifstream in(input_grid_file_name.c_str());
      AssertThrow(in, ExcIO());

      auto ext = extension(input_grid_file_name);
      if (ext == "vtk")
        gi.read_vtk(in);
      else if (ext == "msh")
        gi.read_msh(in);
      else if (ext == "ucd" || ext == "inp")
        gi.read_ucd(in);
      else if (ext == "unv")
        gi.read_unv(in);
      else
        Assert(false, ExcNotImplemented());
    }
  else
    Assert(false, ExcInternalError("Unrecognized grid."));

}

template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::write(const Triangulation<dim,spacedim> &tria) const
{
  if (output_grid_file_name != "")
    {
      GridOut go;
      std::ofstream out(output_grid_file_name.c_str());
      AssertThrow(out, ExcIO());

      go.set_flags(GridOutFlags::Msh(true, true));
      go.set_flags(GridOutFlags::Ucd(false,true, true));

      auto ext = extension(output_grid_file_name);
      if (ext == "vtk")
        go.write_vtk(tria, out);
      else if (ext == "msh")
        go.write_msh(tria, out);
      else if (ext == "ucd" || ext == "inp")
        go.write_ucd(tria, out);
      else if (ext == "vtu")
        go.write_vtu(tria, out);
      else
        Assert(false, ExcNotImplemented());
      out.close();
    }
}


namespace
{
  template <int dim, int spacedim>
  inline typename Triangulation<dim,spacedim>::MeshSmoothing operator|
  (typename Triangulation<dim,spacedim>::MeshSmoothing a,
   typename Triangulation<dim,spacedim>::MeshSmoothing b)
  {
    return static_cast<typename Triangulation<dim,spacedim>::MeshSmoothing>
           (static_cast<int>(a)|static_cast<int>(b));
  }
}

template <int dim, int spacedim>
typename Triangulation<dim,spacedim>::MeshSmoothing
ParsedGridGenerator<dim, spacedim>::get_smoothing()
{
  int smoothing = 0;
  std::vector<std::string> ss = Utilities::split_string_list(mesh_smoothing);
  for (auto s : ss)
    {
      if (s == "none")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::none;
      else if (s == "limit_level_difference_at_vertices")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::limit_level_difference_at_vertices;
      else if (s == "eliminate_unrefined_islands")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::eliminate_unrefined_islands;
      else if (s == "patch_level_1")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::patch_level_1;
      else if (s == "coarsest_level_1")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::coarsest_level_1;
      else if (s == "allow_anisotropic_smoothing")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::allow_anisotropic_smoothing;
      else if (s == "eliminate_refined_inner_islands")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::eliminate_refined_inner_islands;
      else if (s == "eliminate_refined_boundary_islands")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::eliminate_refined_boundary_islands;
      else if (s == "do_not_produce_unrefined_islands")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::do_not_produce_unrefined_islands;
      else if (s == "smoothing_on_refinement")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::smoothing_on_refinement;
      else if (s == "smoothing_on_coarsening")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::smoothing_on_coarsening;
      else if (s == "maximum_smoothing")
        smoothing = smoothing|(int)Triangulation<dim,spacedim>::maximum_smoothing;
      else
        Assert(false, ExcInternalError("Found a non supported smoothing. Should not happen:"+s))
      }
  return static_cast<typename Triangulation<dim,spacedim>::MeshSmoothing>(smoothing);
}



template class ParsedGridGenerator<1,1>;
template class ParsedGridGenerator<1,2>;
template class ParsedGridGenerator<2,2>;
template class ParsedGridGenerator<2,3>;
template class ParsedGridGenerator<3,3>;
