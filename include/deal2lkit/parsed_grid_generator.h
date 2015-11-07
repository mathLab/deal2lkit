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

#ifndef _d2k_parsed_grid_generator_h
#define _d2k_parsed_grid_generator_h

#include <deal2lkit/config.h>
#include <deal.II/base/config.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

/**
 * Parsed grid generator. Create a grid reading a parameter file. This
 * can be either generated using dealii::GridGenerator functions, or
 * read from a file.
 *
 * All Triangulations in the GridGenerator Namespace are supported, as
 * well as all formats supported by deal.II.
 */
template<int dim, int spacedim=dim>
class ParsedGridGenerator : public ParameterAcceptor
{
public:
  /**
   * Constructor. Takes a name for the section of the Parameter
   * Handler to use.
   *
   * This class is derived from ParameterAcceptor. Once you
   * constructed an object of this class, if you call
   * ParameterAcceptor::parse_all_parameters(prm), also the
   * parameters of this class will be filled with values from the
   * argument ParameterHandler.
   */
  ParsedGridGenerator (const std::string section_name="",
                       const std::string grid_type="rectangle",
                       const std::string input_grid_file="",
                       const std::string opt_point_1="",
                       const std::string opt_point_2="",
                       const std::string opt_colorize="false",
                       const std::string opt_double_1="1.0",
                       const std::string opt_double_2="0.5",
                       const std::string opt_double_3="1.5",
                       const std::string opt_int_1="1",
                       const std::string opt_int_2="2",
                       const std::string opt_vec_of_int="",
                       const std::string mesh_smoothing="none",
                       const std::string output_grid_file="");

  /**
   * Declare possible parameters of this class.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Return a pointer to a newly created Triangulation. It will
   * throw an exception if called before any parsing has
   * occured. It's the user's responsability to destroy the created
   * grid once it is no longer needed.
   */
  Triangulation<dim, spacedim> *serial();

  /**
   * Generate the grid.
   *
   * The following grids are implemented:
   *
   * - **file**-> read grid from a file using:
   *   - *Input grid filename*  : input filename\n
   *
   * - **rectangle**-> create a subdivided hyperrectangle using:
   *   - *Point<spacedim>*: lower-left corner
   *   - *Point<spacedim>*: upper-right corner
   *   - *Vector of dim int*: subdivisions on each direction
   *   - *bool*: colorize grid
   *
   * - **subdivided_hyper_rectangle**-> create a subdivided hyperrectangle using:
   *   - *Point<spacedim>*: lower-left corner
   *   - *Point<spacedim>*: upper-right corner
   *   - *Vector of dim int*: subdivisions on each direction
   *   - *bool*: colorize grid
   *
   * - **hyper_sphere**-> generate an hyper sphere with center and radius prescribed:
   *   - *Point<spacedim>*: center
   *   - *double*: radius
   *
   * - **unit_hyperball**-> initialize the given triangulation with a hyperball:
   *   - *Point<spacedim>*: center
   *   - *double*: radius
   *
   * - **subdivided_hyper_rectangle**-> create a coordinate-parallel parallelepiped:
   *   - std::vector<unsigned int> : number of subdivisions in each coordinate direction
   *   - *Point<spacedim>*: lower-left corner
   *   - *Point<spacedim>*: upper-right corner
   *   - *bool*: colorize grid
   *
   * - **hyper_shell**-> create a gird represented by the region between two spheres with fixed center:
   *   - *Point<spacedim>*: center
   *   - *double*: inner sphere radius
   *   - *double*: outer sphere radius
   *   - *unsigned int*: number of cells of the resulting triangulation (In 3D, only 6, 12, and 96 are allowed)
   *   - *bool*: colorize grid
   *
   * - **hyper_L**-> initialize the given triangulation with a hyper-L. It produces the hypercube with the interval [left,right] without the hypercube made out of the interval [(left+right)/2,right] for each coordinate.:
   *   - *double*: left
   *   - *double*: right
   *
   * - **half_hyper_ball**-> produce a half hyper-ball around center, which contains four elements in 2d and 6 in 3d. The cut plane is perpendicular to the x-axis:
   *   - *Point<spacedim>*: center
   *   - *double*: radius
   *
   * - **cylinder**-> create a cylinder around the x-axis. The cylinder extends from x=-half_length to x=+half_length and its projection into the yz-plane is a circle of radius radius:
   *   - *double*: radius
   *   - *double*: half length of the cylinder
   *
   * - **truncated_cone**-> create a cut cone around the x-axis. The cone extends from x=-half_length to x=half_length and its projection into the yz-plane is a circle of radius radius1 at x=-half_length and a circle of radius radius2 at x=+half_length:
   *   - *double*: radius 1
   *   - *double*: radius 2
   *   - *double*: half length
   *
   * - **hyper_cross**-> a center cell with stacks of cell protruding from each surface:
   *   - *Vector of dim int*: sizes
   *   - *bool*: colorize grid
   *
   * - **hyper_cube_slit**-> initialize the given Triangulation with a hypercube with a slit. In each coordinate direction, the hypercube extends from left to right:
   *   - *double*: left
   *   - *double*: right
   *   - *bool*: colorize grid
   *
   * - **half_hyper_shell**-> produce a half hyper-shell, i.e. the space between two circles in two space dimensions and the region between two spheres in 3D:
   *   - *Point<spacedim>*: center
   *   - *double*: inner radius
   *   - *double*: outer radius
   *   - *unsigned int*: number of cells
   *   - *bool*: colorize grid
   *
   * - **quarter_hyper_shell**-> Produce a domain that is the intersection between a hyper-shell with given inner and outer radius, i.e. the space between two circles in two space dimensions and the region between two spheres in 3D, and the positive quadrant (in 2D) or octant (in 3D). In 2D, this is indeed a quarter of the full annulus, while the function is a misnomer in 3D because there the domain is not a quarter but one eighth of the full shell:
   *   - *Point<spacedim>*: center
   *   - *double*: inner radius
   *   - *double*: outer radius
   *   - *unsigned int*: number of cells
   *   - *bool*: colorize grid
   *
   * - **cylinder_shell**-> produce a domain that is the space between two cylinders in 3D, with given length, inner and outer radius and a given number of elements for this initial triangulation. If n_radial_cells is zero (as is the default), then it is computed adaptively such that the resulting elements have the least aspect ratio. The same holds for n_axial_cells:
   *   - *double*: lenght
   *   - *double*: inner radius
   *   - *double*: outer radius
   *   - *unsigned int*: n_radial_cells
   *   - *unsigned int*: n_axial_cells
   *
   * - **moebius**-> produce a ring of cells in 3d that is cut open, twisted and glued together again. This results in a kind of m\"oebius-loop:
   *   - *unsigned int*: number of cells in the loop
   *   - *unsigned int*: number of rotations (Pi/2 each) to be performed before gluing the loop together
   *   - *double*: radius of the circle
   *   - *double*: radius of the cylinder bend together as loop
   *
   * - **hyper_cube_with_cylindrical_hole**-> produces a square in the xy-plane with a circular hole in the middle:
   *   - *double*: inner radius
   *   - *double*: outer radius
   *   - *double*: length
   *   - *unsigned int*: repetitions (number of subdivisions along the z-direction)
   *   - *bool*: colorize grid
   *
   * - **torus **-> produce the surface meshing of the torus:
   *   - *double*: radius of the circle which forms the middle line of the torus containing the loop of cells
   *   - *double*: inner radius of the torus
   *
   * - **cheese **-> domain itself is rectangular. The argument holes specifies how many square holes the domain should have in each coordinate direction :"
   *   - *Vector of dim int*: number of holes on each direction"
   */
  void create(Triangulation<dim, spacedim> &tria);

  /**
     * Return a pointer to a newly created Triangulation. It will
     * throw an exception if called before any parsing has
     * occured. It's the user's responsability to destroy the created
     * grid once it is no longer needed.
     */
#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_P4EST
  parallel::distributed::Triangulation<dim, spacedim> *distributed(MPI_Comm mpi_communicator);
#endif
#endif

  std::string create_default_value(const std::vector<unsigned int> &input);

  std::string create_default_value(const std::vector<double> &input);

  std::string create_default_value(const Point<spacedim> &input);

  void write(const Triangulation<dim, spacedim> &tria) const;

private:

  typename Triangulation<dim,spacedim>::MeshSmoothing
  get_smoothing();
  /**
   * Mesh smoothing. Parse the type of MeshSmoothing for the
   * generated Triangulation.
   */
  std::string mesh_smoothing;

  /**
   * The grid to generate. Use the name "file" to read from a file.
   */
  std::string grid_name;

  /**
   * The optional argument for the grid generators. We choose two doubles and two points.
   */

  double double_option_one;

  double double_option_two;

  double double_option_three;

  Point<spacedim> point_option_one;

  Point<spacedim> point_option_two;

  unsigned int un_int_option_one;

  unsigned int un_int_option_two;

  bool colorize;

  std::vector<unsigned int> un_int_vec_option_one;

  /**
   * Grid file name.
   */
  std::string input_grid_file_name;

  std::string output_grid_file_name;

  // strings for prm
  std::string str_point_1;
  std::string str_point_2;
  std::string str_colorize;
  std::string str_double_1;
  std::string str_double_2;
  std::string str_double_3;
  std::string str_un_int_1;
  std::string str_un_int_2;
  std::string str_vec_int;
};

D2K_NAMESPACE_CLOSE

#endif
