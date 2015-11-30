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
#include <deal2lkit/utilities.h>

#include <deal.II/base/config.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>

#include <deal2lkit/parameter_acceptor.h>

using namespace dealii;

D2K_NAMESPACE_OPEN

struct PGGHelper;

/**
 * Parsed grid generator. Create a grid by reading a parameter file,
 * either using dealii::GridGenerator functions, or by reading it from
 * a file in supported format.
 *
 * All Triangulations in the GridGenerator Namespace are supported, as
 * well as all file formats supported by GridIn. This class is in all
 * effects a wrapper around the GridGenerator namespace functions,
 * GridIn, and GridOut classes of the \dealii library
 *
 *
 * ParsedGridGenerator can be used both in serial and parallel
 * settings, and a typical usage of this class is
 * \code
 *   // 2D square - serial mesh
 *   // by default it constructs the rectangle whose opposite corner points
 *   // are p1=(10.0,10.0) and p2=(20.0,20.0)
 *   ParsedGridGenerator<2,2> tria_builder_2d("2D mesh",
 *                                            "rectangle",
 *                                            "",
 *                                            "10.0,10.0",
 *                                            "20.0,20.0",
 *                                            "true");
 *
 *   // 3D parallelepiped - parallel distributed mesh
 *   // by default it constructs the parallelepiped whose opposite corner
 *   // points are p1=(7.0,8.0,9.0) and p2=(15.0.16.0,16.7)
 *   ParsedGridGenerator<3,3> tria_builder_3d("3D mesh",
 *                                             "rectangle",
 *                                             "",
 *                                             "7.0,8.0,9.0",
 *                                             "15.0,16.0,16.7");
 *
 *   // call ParameterAcceptor
 *   ParameterAcceptor::initialize("file.prm", "no_descriptions.prm");
 *
 *   // Construct a serial mesh following the indications in "file.prm"
 *   // in the section "2D mesh"
 *   Triangulation<2,2> *tria_serial  = tria_builder_2d.serial();
 *
 *   // Construct a parallel::distributed::Triangulation following the
 *   // indications in "file.prm'', in the section "3D mesh"
 *   parallel::distributed::Triangulation<3,3> *tria_mpi =
 *         tria_builder_3d.distributed(MPI_COMM_WORLD);
 * \endcode
 *
 * Once the function ParameterAcceptor::initialize("file.prm",
 * "no_descriptions.prm") is called, the "no_descriptions.prm" is
 * filled with the following entries:
 *
 * \code{bash}
 * subsection 3D mesh
 *   set Colorize                   = false
 *   set Grid to generate           = rectangle
 *   set Input grid file name       =
 *   set Mesh smoothing alogrithm   = none
 *   set Optional Point<spacedim> 1 = 7.0,8.0,9.0
 *   set Optional Point<spacedim> 2 = 15.0,16.0,16.7
 *   set Optional double 1          = 1.0
 *   set Optional double 2          = 0.5
 *   set Optional int 1             = 1
 *   set Optional vector of dim int = 1,1,1
 *   set Output grid file name      =
 * end
 * subsection 2D mesh
 *   set Colorize                   = true
 *   set Grid to generate           = rectangle
 *   set Input grid file name       =
 *   set Mesh smoothing alogrithm   = none
 *   set Optional Point<spacedim> 1 = 10.0,10.0
 *   set Optional Point<spacedim> 2 = 20.0,20.0
 *   set Optional double 1          = 1.0
 *   set Optional double 2          = 0.5
 *   set Optional int 1             = 1
 *   set Optional vector of dim int = 1,1
 *   set Output grid file name      =
 * end
 * \endcode
 *
 * If the user would then change the parameter file to generate a sphere,
 * or read a file, no change in the code would be necessary, as at run
 * time the new parameters would be used.
 *
 * If the option `Output grid file name` is set to non-empty, when the
 * user calls ParsedGridGenerator::write(), an output grid would be
 * generated using GridOut, choosing the right format according to the
 * extension of the file.
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
   * ParameterAcceptor::parse_all_parameters(), also the
   * parameters of this class will be filled with values from the
   * argument ParameterHandler.
   *
   * This constructor takes optional strings which allow the user to
   * decide what are the **default** values that will be written on
   * the parameter file. The first optional argument specifies the
   * section name within the parameter file. If the section name is
   * empty, by default it is set to ParsedGridGenerator<x,x> where
   * `x,x` is replaced with the actual `dim` and `spacedim` numbers
   * with which the user instantiated the class.
   *
   * Since every GridGenerator function takes a different set of
   * arguments, we provide a list of optional arguments that will be
   * used by the class itself.
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
                       const std::string output_grid_file="",
                       const std::string opt_manifold_descriptors="");


  /**
   * Return a list of implemented grids.
   */
  static std::string get_grid_names();

  /**
   * Declare all parameters of this class.
   */
  virtual void declare_parameters(ParameterHandler &prm);

  /**
   * Return a pointer to a newly created serial Triangulation. It will
   * throw an exception if called before any parsing has occured. It
   * is the user's responsability to destroy the created grid once it
   * is no longer needed.
   */
  Triangulation<dim, spacedim> *serial();

  /**
   * Generate the grid. Fill a user supplied empty Triangulation using
   * the parameter file. If the Triangulation is not empty, an
   * exception is thrown.
   *
   * The following grids are implemented:
   *
   * - **file**-> read grid from a file using:
   *   - *Input grid filename*  : input filename\n
   *
   * - **rectangle**-> calls GridGenerator::SubdividedHyperRectangle() using:
   *   - *Point<spacedim>*: lower-left corner
   *   - *Point<spacedim>*: upper-right corner
   *   - *Vector of dim int*: subdivisions on each direction
   *   - *bool*: colorize grid
   *
   * - **hyper_sphere**-> calls GridGenerator::HyperSphere() using:
   *   - *Point<spacedim>*: center
   *   - *double*: radius
   *
   * - **hyper_ball**-> calls GridGenerator::HyperBall() using:
   *   - *Point<spacedim>*: center
   *   - *double*: radius
   *
   * - **parallelepiped**-> ccalls GridGenerator::Parallelepiped() using:
   *   - std::vector<unsigned int> : number of subdivisions in each coordinate direction
   *   - *Point<spacedim>*: lower-left corner
   *   - *Point<spacedim>*: upper-right corner
   *   - *bool*: colorize grid
   *
   * - **hyper_shell**-> calls GridGenerator::HyperShell() using:
   *   - *Point<spacedim>*: center
   *   - *double*: inner sphere radius
   *   - *double*: outer sphere radius
   *   - *unsigned int*: number of cells of the resulting triangulation (In 3D, only 6, 12, and 96 are allowed)
   *   - *bool*: colorize grid
   *
   * - **hyper_L**-> GridGenerator::HyperL(). It produces the hypercube with the interval [left,right] without the hypercube made out of the interval [(left+right)/2,right] for each coordinate.:
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
   * - **torus**-> produce the surface meshing of the torus:
   *   - *double*: radius of the circle which forms the middle line of the torus containing the loop of cells
   *   - *double*: inner radius of the torus
   *
   * - **cheese**-> domain itself is rectangular. The argument holes specifies how many square holes the domain should have in each coordinate direction :"
   *   - *Vector of dim int*: number of holes on each direction"
   */
  void create(Triangulation<dim, spacedim> &tria);

#ifdef DEAL_II_WITH_MPI
#ifdef DEAL_II_WITH_P4EST
  /**
   * Return a pointer to a newly created parallel Triangulation. It
   * will throw an exception if called before any parsing has
   * occured. It is the user's responsability to destroy the created
   * grid once it is no longer needed.
   */
  parallel::distributed::Triangulation<dim, spacedim> *distributed(MPI_Comm mpi_communicator);
#endif
#endif

  /**
   * Write the given Triangulation to the output file specified in
   * `Output file name`, or in the optional file name.
   *
   * If no `Output file name` is given and filename is the empty string,
   * this function does nothing. If
   * an output file name is provided (either in the input file, or as an argument
   * to this function), then this function will call the
   * appropriate GridOut method according to the extension of the file
   * name.
   */
  void write(const Triangulation<dim, spacedim> &tria,
             const std::string &filename="") const;

private:
  /**
   * Mesh smoothing. Parse the type of MeshSmoothing for the
   * generated Triangulation.
   */
  typename Triangulation<dim,spacedim>::MeshSmoothing
  get_smoothing();

  /**
   * Take @p str_manifold_descriptors and fill @p manifold_descriptors with
   * ids and manifolds.
   */
  void parse_manifold_descriptors(const std::string &str_manifold_descriptors);

  /**
   * Mesh smoothing to apply to the newly created Triangulation. This
   * variable is only used if the method serial() is called. For the
   * method parallel(), mesh smoothing is not yet supported by
   * \dealii.
   */
  std::string mesh_smoothing;

  /**
   * The grid to generate. Use the name "file" to read from a file.
   */
  std::string grid_name;

  /**
   * Manifold descriptors.
   * Pattern to be used: id followed by '=' manifold descriptor
   * each couple of id and  manifold descriptor is separated by '%'
   * Avaible manifold descriptor:
   * - HyperBallBoundary
   * - HyperShellBoundary
   * - CADSurface
   * - CADLine
   */
  std::string optional_manifold_descriptors;


  /**
   * Default Manifold descriptors. This is filled when creating the grid,
   * and is later translated into an actual manifold.
   */
  std::string default_manifold_descriptors;

  /**
   * A map of Manifold associated to the given manifold_ids.
   */
  std::map<types::manifold_id, shared_ptr<Manifold<dim,spacedim> > > manifold_descriptors;

  /**
   * Optional double argument. First option.
   */
  double double_option_one;

  /**
   * Optional double argument. Second option.
   */
  double double_option_two;

  /**
   * Optional double argument. Third option.
   */
  double double_option_three;

  /**
   * Optional Point argument. First Option.
   */
  Point<spacedim> point_option_one;

  /**
   * Optional Point argument. Second Option.
   */
  Point<spacedim> point_option_two;

  /**
   * Optional int argument. First Option.
   */
  unsigned int un_int_option_one;

  /**
   * Optional int argument. First Option.
   */
  unsigned int un_int_option_two;

  /**
   * For all the internally generated meshes that support it, turn on
   * colorizing.
   */
  bool colorize;

  /**
   * Create default manifold descriptors. If set to true, boundary ids will
   * be copied over manifold ids on the newly created
   * triangulation (independently on the value of the variable
   * copy_boundary_to_manifold_ids, and for each triangulation where we know how
   * to create
   * and associate their manifolds, we create them and associate them to
   * the newly created triangulation.
   *
   * This option produces different associations depending on the colorize
   * parameter. The created manifolds will be compatible with the triangulation
   * and the colorize parameter used.
   */
  bool create_default_manifolds;

  /**
   * Copy boundary ids to manifold ids. If set to true, boundary ids will be
   * copied over manifold ids on the newly created
   * triangulation.
   */
  bool copy_boundary_to_manifold_ids;

  /**
   * Copy material ids to manifold ids. If set to true, material ids will be
   * copied over manifold ids on the newly created
   * triangulation.
   */
  bool copy_material_to_manifold_ids;

  /**
   * Optional vector of integers.
   */
  std::vector<unsigned int> un_int_vec_option_one;

  /**
   * Input grid file name.
   */
  std::string input_grid_file_name;

  /**
   * Output grid file name.
   */
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

  /**
   * Helper function to create grids.
   */
  friend struct PGGHelper;
};

D2K_NAMESPACE_CLOSE

#endif
