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
                       const std::string opt_point_3="",
                       const std::string opt_bool="false",
                       const std::string opt_double_1="1.0",
                       const std::string opt_double_2="0.5",
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
   * - **file** -> read grid from a file (*vtk*, *msh*, *ucd*, *inp* and *unv*)
   * using:
   *  - Input grid filename     : input filename
   *
   * - **rectangle** -> create a subdivided hyperrectangle using:
   *  - *Point<spacedim>* : lower-left corner
   *  - *Point<spacedim>* : upper-right corner
   *  - *std::vector<unsigned int>* : subdivisions on each direction
   *  - *bool* : colorize grid
   *
   * - **hyper_sphere** -> generate an hyper sphere with center and radius
   * prescribed:
   *  - *Point<spacedim>* : center
   *  - *double* : radius
   *
   * - **unit_hyperball** -> initialize the given triangulation with a hyperball
   *  - *Point<spacedim>* : center
   *  - *double* : radius
   *
   * - **subdivided_hyper_rectangle** -> create a coordinate-parallel parallelepiped:
   *  - *std::vector<unsigned int>* : number of subdivisions in each coordinate
   *  direction
   *  - *Point<spacedim>* : lower-left corner
   *  - *Point<spacedim>* : upper-right corner
   *  - *bool* : colorize grid
   *
   * - **hyper_shell** -> create a gird represented by the region between
   * two spheres with fixed center:
   *  - *Point<spacedim>* : center
   *  - *double* : inner sphere radius
   *  - *double* : outer sphere radius
   *  - *unsigned int* : number of cells of the resulting triangulation (In 3d,
   *  only 6, 12, and 96 are allowed)
   *  - *bool* : colorize grid
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

  bool bool_option_one;

  std::vector<unsigned int> un_int_vec_option_one;

  /**
   * Grid file name.
   */
  std::string input_grid_file_name;

  std::string output_grid_file_name;

  // strings for prm
  std::string str_point_1;
  std::string str_point_2;
  std::string str_bool;
  std::string str_double_1;
  std::string str_double_2;
  std::string str_double_3;
  std::string str_un_int_1;
  std::string str_un_int_2;
  std::string str_vec_int;
};

D2K_NAMESPACE_CLOSE

#endif
