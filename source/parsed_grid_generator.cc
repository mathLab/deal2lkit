//-----------------------------------------------------------
//
//    Copyright (C) 2015 by the deal2lkit authors
//
//    This file is part of the deal2lkit library.
//
//    The deal2lkit library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE at
//    the top level of the deal2lkit distribution.
//
//-----------------------------------------------------------

#include <deal.II/base/config.h>
#include <deal2lkit/parsed_grid_generator.h>
#include <deal2lkit/utilities.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>


namespace
{
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
}

D2K_NAMESPACE_OPEN


template <int dim, int spacedim>
ParsedGridGenerator<dim, spacedim>::ParsedGridGenerator(const std::string _section_name,
                                                        const std::string _grid_type,
                                                        const std::string _input_grid_file,
                                                        const std::string _point_1,
                                                        const std::string _point_2,
                                                        const std::string _colorize,
                                                        const std::string _double_1,
                                                        const std::string _double_2,
                                                        const std::string _double_3,
                                                        const std::string _int_1,
                                                        const std::string _int_2,
                                                        const std::string _vec_of_int,
                                                        const std::string _mesh_smoothing,
                                                        const std::string _output_grid_file)
  :
  ParameterAcceptor(_section_name),
  mesh_smoothing(_mesh_smoothing),
  grid_name(_grid_type),
  input_grid_file_name(_input_grid_file),
  output_grid_file_name(_output_grid_file),
  str_point_1(_point_1),
  str_point_2(_point_2),
  str_colorize(_colorize),
  str_double_1(_double_1),
  str_double_2(_double_2),
  str_double_3(_double_3),
  str_un_int_1(_int_1),
  str_un_int_2(_int_2),
  str_vec_int(_vec_of_int)
{}


template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::declare_parameters(ParameterHandler &prm)
{

  std::vector<unsigned int> dummy_vec_int(dim, 1);
  std::vector<double> dummy_vec_double(spacedim);
  Point<spacedim> dummy_point;
  Point<spacedim> dummy_point_2;
  for (unsigned int d=0; d<dim; ++d)
    dummy_point_2[d]=1.0;

  std::string def_point, def_point_2, def_int, def_double;
  def_point = print(dummy_point);
  def_point_2 = print(dummy_point_2);
  def_int = print(dummy_vec_int);
  def_double = print(dummy_vec_double);

  add_parameter(prm, &grid_name,
                "Grid to generate", grid_name,
                Patterns::Selection("file|rectangle|unit_hyperball|hyper_shell|subdivided_hyper_rectangle|hyper_sphere|hyper_L|half_hyper_ball|cylinder|truncated_cone|hyper_cross|hyper_cube_slit|half_hyper_shell|quarter_hyper_shell|cylinder_shell|torus|hyper_cube_with_cylindrical_hole|moebius|cheese"),
                "The grid to generate. You can choose among:\n"
                "- file: read grid from a file using:\n"
                "	- Input grid filename	    : input filename\n\n"
                "- rectangle: create a subdivided hyperrectangle using:\n"
                "	- Optional Point<spacedim> 1: lower-left corner\n"
                "	- Optional Point<spacedim> 2: upper-right corner\n"
                "	- Optional Vector of dim int: subdivisions on each direction\n"
                "	- Optional bool 1	    : colorize grid\n"
                "- subdivided_hyper_rectangle: create a subdivided hyperrectangle using:\n"
                "	- Optional Point<spacedim> 1: lower-left corner\n"
                "	- Optional Point<spacedim> 2: upper-right corner\n"
                "	- Optional Vector of dim int: subdivisions on each direction\n"
                "	- Optional bool 1	    : colorize grid\n"
                "- hyper_sphere  : generate an hyper sphere with center and radius prescribed:\n"
                "	- Optional Point<spacedim> : center\n"
                "	- Optional double : radius\n"
                "- unit_hyperball  : initialize the given triangulation with a hyperball:\n\n"
                "	- Optional Point<spacedim> : center\n"
                "	- Optional double : radius\n"
                "- subdivided_hyper_rectangle   : create a coordinate-parallel parallelepiped:\n\n"
                "	- Optional std::vector<unsigned int> : number of subdivisions in each coordinate direction\n"
                "	- Optional Point<spacedim> : lower-left corner\n"
                "	- Optional Point<spacedim> : upper-right corner\n"
                "	- Optional bool 1	    : colorize grid\n"
                "- hyper_shell   : create a gird represented by the region between two spheres with fixed center:\n\n"
                "	- Optional Point<spacedim> : center\n"
                "	- Optional double : inner sphere radius\n"
                "	- Optional double : outer sphere radius\n"
                "	- Optional unsigned int : number of cells of the resulting triangulation (In 3d, only 6, 12, and 96 are allowed)\n"
                "	- Optional bool : colorize grid\n"
                "- hyper_L : initialize the given triangulation with a hyper-L. It produces the hypercube with the interval [left,right] without the hypercube made out of the interval [(left+right)/2,right] for each coordinate."
                "	- Optional double : left\n"
                "	- Optional double : right\n"
                "- half_hyper_ball : produce a half hyper-ball around center, which contains four elements in 2d and 6 in 3d. The cut plane is perpendicular to the x-axis: \n"
                "	- Optional Point<spacedim> : center\n"
                "	- Optional double : radius\n"
                "- cylinder: create a cylinder around the x-axis. The cylinder extends from x=-half_length to x=+half_length and its projection into the yz-plane is a circle of radius radius: \n"
                "	- Optional double : radius\n"
                "	- Optional double : half length of the cylinder\n"
                "- truncated_cone : create a cut cone around the x-axis. The cone extends from x=-half_length to x=half_length and its projection into the yz-plane is a circle of radius radius1 at x=-half_length and a circle of radius radius2 at x=+half_length :\n"
                "	- Optional double : radius 1\n"
                "	- Optional double : radius 2\n"
                "	- Optional double : half length\n"
                "- hyper_cross : a center cell with stacks of cell protruding from each surface :\n"
                "	- Optional Vector of dim int: sizes\n"
                "	- Optional bool : colorize grid\n"
                "- hyper_cube_slit : initialize the given Triangulation with a hypercube with a slit. In each coordinate direction, the hypercube extends from left to right :\n"
                "	- Optional double : left\n"
                "	- Optional double : right\n"
                "	- Optional bool : colorize grid\n"
                "- half_hyper_shell : produce a half hyper-shell, i.e. the space between two circles in two space dimensions and the region between two spheres in 3d :\n"
                "	- Optional Point<spacedim> : center\n"
                "	- Optional double : inner radius\n"
                "	- Optional double : outer radius\n"
                "	- Optional unsigned int : number of cells\n"
                "	- Optional bool : colorize grid\n"
                "- quarter_hyper_shell : Produce a domain that is the intersection between a hyper-shell with given inner and outer radius, i.e. the space between two circles in two space dimensions and the region between two spheres in 3d, and the positive quadrant (in 2d) or octant (in 3d). In 2d, this is indeed a quarter of the full annulus, while the function is a misnomer in 3d because there the domain is not a quarter but one eighth of the full shell :\n"
                "	- Optional Point<spacedim> : center\n"
                "	- Optional double : inner radius\n"
                "	- Optional double : outer radius\n"
                "	- Optional unsigned int : number of cells\n"
                "	- Optional bool : colorize grid\n"
                "- cylinder_shell : produce a domain that is the space between two cylinders in 3d, with given length, inner and outer radius and a given number of elements for this initial triangulation. If n_radial_cells is zero (as is the default), then it is computed adaptively such that the resulting elements have the least aspect ratio. The same holds for n_axial_cells :\n"
                "	- Optional double : lenght\n"
                "	- Optional double : inner radius\n"
                "	- Optional double : outer radius\n"
                "	- Optional unsigned int : n_radial_cells\n"
                "	- Optional unsigned int : n_axial_cells\n"
                "- moebius : produce a ring of cells in 3d that is cut open, twisted and glued together again. This results in a kind of moebius-loop :\n"
                "	- Optional unsigned int : number of cells in the loop\n"
                "	- Optional unsigned int : number of rotations (Pi/2 each) to be performed before gluing the loop together\n"
                "	- Optional double : radius of the circle\n"
                "	- Optional double : radius of the cylinder bend together as loop\n"
                "- hyper_cube_with_cylindrical_hole :  produces a square in the xy-plane with a circular hole in the middle :\n"
                "	- Optional double : inner radius\n"
                "	- Optional double : outer radius\n"
                "	- Optional double : lenght\n"
                "	- Optional unsigned int : repetitions (number of subdivisions along the z-direction)\n"
                "	- Optional bool : colorize grid\n"
                "- torus : produce the surface meshing of the torus :\n"
                "	- Optional double : radius of the circle which forms the middle line of the torus containing the loop of cells\n"
                "	- Optional double :  inner radius of the torus\n"
                "- cheese : domain itself is rectangular. The argument holes specifies how many square holes the domain should have in each coordinate direction :\n"
                "	- Optional Vector of dim int: number of holes on each direction\n"
               );

  add_parameter(prm, &mesh_smoothing,
                "Mesh smoothing alogrithm", mesh_smoothing,
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
                "Input grid file name", input_grid_file_name,
                Patterns::FileName(),
                "Name of the input grid. All supported deal.II formats. "
                "The extestion will be used to decide what "
                "grid format to use.");

  add_parameter(prm, &double_option_one,
                "Optional double 1", str_double_1,
                Patterns::Double(),
                "First additional double to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &double_option_two,
                "Optional double 2", str_double_2,
                Patterns::Double(),
                "Second additional double to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &double_option_three,
                "Optional double 3", str_double_3,
                Patterns::Double(),
                "Second additional double to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &point_option_one,
                "Optional Point<spacedim> 1", (str_point_1==""?def_point:str_point_1),
                Patterns::List(Patterns::Double(), spacedim, spacedim),
                "First additional Point<spacedim> to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &point_option_two,
                "Optional Point<spacedim> 2", (str_point_2==""?def_point_2:str_point_2),
                Patterns::List(Patterns::Double(), spacedim, spacedim),
                "Second additional Point<spacedim> to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");


  add_parameter(prm, &un_int_option_one,
                "Optional int 1",str_un_int_1,
                Patterns::Integer(),
                "Unsigned int to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &un_int_option_one,
                "Optional int 2",str_un_int_2,
                Patterns::Integer(),
                "Unsigned int to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &un_int_vec_option_one,
                "Optional vector of dim int", (str_vec_int==""?def_int:str_vec_int),
                Patterns::List(Patterns::Integer(1), dim, dim),
                "Vector of positive unsigned int to be used in the generation of the grid. "
                "The use of it will depend on the specific grid.");

  add_parameter(prm, &colorize,
                "Colorize",str_colorize,
                Patterns::Bool(),
                "Bool be used in the generation of the grid to set colorize. "
                "The use of it will depend on the specific grid.");


  add_parameter(prm, &output_grid_file_name,
                "Output grid file name", output_grid_file_name,
                Patterns::FileName(),
                "Name of the output grid. All supported deal.II formats. "
                "The extestion will be used to decide what "
                "grid format to use. If empty, no grid will be written.");

}

#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_P4EST
template <int dim, int spacedim>
parallel::distributed::Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::distributed(MPI_Comm comm)
{
  Assert(grid_name != "", ExcNotInitialized());
  parallel::distributed::Triangulation<dim,spacedim> *tria = new parallel::distributed::Triangulation<dim,spacedim>
  (comm);//, get_smoothing());

  create(*tria);
  return tria;
}
#endif
#endif


template <int dim, int spacedim>
Triangulation<dim, spacedim> *
ParsedGridGenerator<dim, spacedim>::serial()
{
  Assert(grid_name != "", ExcNotInitialized());
  Triangulation<dim,spacedim> *tria = new Triangulation<dim,spacedim>(get_smoothing());
  create(*tria);
  return tria;
}

namespace
{
  /**
   * This function is used to generate grids when there are no restriction on
   * the template parameters. This is the reason it is called inside the
   * following functions (create_grid).
   */
  template<int dim, int spacedim>
  void
  default_create_grid( Triangulation<dim,spacedim> &tria,
                       std::string grid_name,
                       Point<spacedim> point_option_one,
                       Point<spacedim> point_option_two,
                       double double_option_two,
                       double double_option_one,
                       double double_option_three,
                       bool colorize,
                       unsigned int un_int_option_one,
                       unsigned int un_int_option_two,
                       std::vector<unsigned int> un_int_vec_option_one,
                       std::string input_grid_file_name)
  {
    if (grid_name == "rectangle" || grid_name == "subdivided_hyper_rectangle")
      {
        GridGenerator::subdivided_hyper_rectangle (tria,
                                                   un_int_vec_option_one,
                                                   point_option_two,
                                                   point_option_one,
                                                   colorize);
      }
    else if (grid_name == "file")
      {
        GridIn<dim, spacedim> gi;
        gi.attach_triangulation(tria);

        std::ifstream in(input_grid_file_name.c_str());
        AssertThrow(in, ExcIO());

        std::string ext = extension(input_grid_file_name);
        if (ext == "vtk")
          gi.read_vtk(in);
        else if (ext == "msh")
          gi.read_msh(in);
        else if (ext == "ucd" || ext == "inp")
          gi.read_ucd(in);
        else if (ext == "unv")
          gi.read_unv(in);
        else if (ext == "ar")
          {
            boost::archive::text_iarchive ia(in);
            tria.load(ia, 0);
          }
        else
          Assert(false, ExcNotImplemented());
      }
    // silence warnings of unused variables
    (void)double_option_one;
    (void)double_option_two;
    (void)double_option_three;
    (void)un_int_option_one;
    (void)un_int_option_two;
  }

  /**
   * This function is used to generate grids when spacedim = dim + 1.
   */
  template<int dim>
  void
  create_grid( Triangulation<dim,dim+1> &tria,
               std::string grid_name,
               Point<dim+1> point_option_one,
               Point<dim+1> point_option_two,
               double double_option_two,
               double double_option_one,
               double double_option_three,
               bool colorize,
               unsigned int un_int_option_one,
               unsigned int un_int_option_two,
               std::vector<unsigned int> un_int_vec_option_one,
               std::string input_grid_file_name)
  {
    if (grid_name == "hyper_sphere")
      {
        GridGenerator::hyper_sphere<dim,dim+1> ( tria,
                                                 point_option_one,
                                                 double_option_one);
      }
    else
      {
        default_create_grid<dim, dim+1>( tria,
                                         grid_name,
                                         point_option_one,
                                         point_option_two,
                                         double_option_two,
                                         double_option_one,
                                         double_option_three,
                                         colorize,
                                         un_int_option_one,
                                         un_int_option_two,
                                         un_int_vec_option_one,
                                         input_grid_file_name);
      }
  }

  /**
   * This function is used to generate grids when spacedim = dim.
   */
  template<int dim>
  void
  create_grid( Triangulation<dim,dim> &tria,
               std::string grid_name,
               Point<dim> point_option_one,
               Point<dim> point_option_two,
               double double_option_two,
               double double_option_one,
               double double_option_three,
               bool colorize,
               unsigned int un_int_option_one,
               unsigned int un_int_option_two,
               std::vector<unsigned int> un_int_vec_option_one,
               std::string input_grid_file_name)
  {
    if (grid_name == "unit_hyperball")
      {

        GridGenerator::hyper_ball( tria,
                                   point_option_one,
                                   double_option_one);
      }
    else if (grid_name == "hyper_L")
      {
        GridGenerator::hyper_L  ( tria,
                                  double_option_two,
                                  double_option_one );
      }
    else if (grid_name == "half_hyper_ball")
      {
        GridGenerator::half_hyper_ball (tria,
                                        point_option_one,
                                        double_option_one);
      }
    else if (grid_name == "cylinder")
      {
        GridGenerator::cylinder ( tria,
                                  double_option_one,
                                  double_option_two);
      }
    else if (grid_name == "truncated_cone")
      {
        GridGenerator::truncated_cone  ( tria,
                                         double_option_one,
                                         double_option_two,
                                         double_option_three);
      }
    else if (grid_name == "hyper_cross")
      {
        GridGenerator::hyper_cross  (tria,
                                     un_int_vec_option_one,
                                     colorize);
      }
    else if (grid_name == "hyper_cube_slit")
      {
        GridGenerator::hyper_cube_slit  (tria,
                                         double_option_two,
                                         double_option_one,
                                         colorize);
      }
    else if (grid_name == "half_hyper_shell")
      {
        GridGenerator::half_hyper_shell ( tria,
                                          point_option_one,
                                          double_option_two,
                                          double_option_one,
                                          un_int_option_one,
                                          colorize);
      }
    else if (grid_name == "quarter_hyper_shell")
      {
        GridGenerator::quarter_hyper_shell  ( tria,
                                              point_option_one,
                                              double_option_two,
                                              double_option_one,
                                              un_int_option_one,
                                              colorize);
      }
    else if (grid_name == "cylinder_shell")
      {
        GridGenerator::cylinder_shell ( tria,
                                        double_option_three,
                                        double_option_two,
                                        double_option_one,
                                        un_int_option_one,
                                        un_int_option_two);
      }
    else if (grid_name == "hyper_cube_with_cylindrical_hole")
      {
        GridGenerator::hyper_cube_with_cylindrical_hole ( tria,
                                                          double_option_two,
                                                          double_option_one,
                                                          double_option_three,
                                                          un_int_option_one,
                                                          colorize);
      }
    else if (grid_name == "hyper_shell")
      {
        GridGenerator::hyper_shell( tria,
                                    point_option_one,
                                    double_option_two,
                                    double_option_one,
                                    un_int_option_one,
                                    colorize);
      }
    else if (grid_name == "cheese")
      {
        GridGenerator::cheese(tria,
                              un_int_vec_option_one);
      }
    else
      {
        default_create_grid<dim, dim>( tria,
                                       grid_name,
                                       point_option_one,
                                       point_option_two,
                                       double_option_two,
                                       double_option_one,
                                       double_option_three,
                                       colorize,
                                       un_int_option_one,
                                       un_int_option_two,
                                       un_int_vec_option_one,
                                       input_grid_file_name);
      }
  }

  /**
   * This function is used to generate grids when spacedim = dim = 3.
   */
  void
  create_grid( Triangulation<3,3> &tria,
               std::string grid_name,
               Point<3> point_option_one,
               Point<3> point_option_two,
               double double_option_two,
               double double_option_one,
               double double_option_three,
               bool colorize,
               unsigned int un_int_option_one,
               unsigned int un_int_option_two,
               std::vector<unsigned int> un_int_vec_option_one,
               std::string input_grid_file_name)
  {
    if (grid_name == "moebius")
      {
        GridGenerator::moebius  ( tria,
                                  un_int_option_one,
                                  un_int_option_two,
                                  double_option_two,
                                  double_option_one );
      }
    else
      {
        create_grid<3>( tria,
                        grid_name,
                        point_option_one,
                        point_option_two,
                        double_option_two,
                        double_option_one,
                        double_option_three,
                        colorize,
                        un_int_option_one,
                        un_int_option_two,
                        un_int_vec_option_one,
                        input_grid_file_name);
      }
  }

  /**
   * This function is used to generate grids when dim = 2 and spacedim = 3.
   */
  void
  create_grid( Triangulation<2,3> &tria,
               std::string grid_name,
               Point<3> point_option_one,
               Point<3> point_option_two,
               double double_option_two,
               double double_option_one,
               double double_option_three,
               bool colorize,
               unsigned int un_int_option_one,
               unsigned int un_int_option_two,
               std::vector<unsigned int> un_int_vec_option_one,
               std::string input_grid_file_name)
  {
    if (grid_name == "torus")
      {
        GridGenerator::torus  (tria,
                               double_option_one,
                               double_option_two);
      }
    else
      create_grid<2>( tria,
                      grid_name,
                      point_option_one,
                      point_option_two,
                      double_option_two,
                      double_option_one,
                      double_option_three,
                      colorize,
                      un_int_option_one,
                      un_int_option_two,
                      un_int_vec_option_one,
                      input_grid_file_name);
  }

}






template <int dim, int spacedim>
void ParsedGridGenerator<dim, spacedim>::create(Triangulation<dim,spacedim> &tria)
{
  Assert(grid_name != "", ExcNotInitialized());

  create_grid( tria,
               grid_name,
               point_option_one,
               point_option_two,
               double_option_two,
               double_option_one,
               double_option_three,
               colorize,
               un_int_option_one,
               un_int_option_two,
               un_int_vec_option_one,
               input_grid_file_name);

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

      std::string ext = extension(output_grid_file_name);
      if (ext == "vtk")
        go.write_vtk(tria, out);
      else if (ext == "msh")
        go.write_msh(tria, out);
      else if (ext == "ucd" || ext == "inp")
        go.write_ucd(tria, out);
      else if (ext == "vtu")
        go.write_vtu(tria, out);
      else if (ext == "ar")
        {
          boost::archive::text_oarchive oa(out);
          tria.save(oa, 0);
        }
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
  for (std::string s : ss)
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


D2K_NAMESPACE_CLOSE



template class deal2lkit::ParsedGridGenerator<1,1>;
template class deal2lkit::ParsedGridGenerator<1,2>;
template class deal2lkit::ParsedGridGenerator<2,2>;
template class deal2lkit::ParsedGridGenerator<2,3>;
template class deal2lkit::ParsedGridGenerator<3,3>;
