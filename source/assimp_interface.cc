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

#include <deal2lkit/assimp_interface.h>

#ifdef D2K_WITH_ASSIMP

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#undef AI_CONFIG_PP_RVC_FLAGS
#define AI_CONFIG_PP_RVC_FLAGS           \
  aiComponent_NORMALS |            \
  aiComponent_TANGENTS_AND_BITANGENTS |        \
  aiComponent_COLORS    |        \
  aiComponent_TEXCOORDS   |        \
  aiComponent_BONEWEIGHTS |        \
  aiComponent_ANIMATIONS |           \
  aiComponent_TEXTURES  |          \
  aiComponent_LIGHTS |             \
  aiComponent_CAMERAS |            \
  aiComponent_MATERIALS


using namespace dealii;
using namespace Assimp;

using namespace std;


D2K_NAMESPACE_OPEN

namespace AssimpInterface
{

  template <int dim, int spacedim>
  bool generate_triangulation(const std::string filename,
                              Triangulation<dim,spacedim> &tria,
                              int mesh_index,
                              bool remove_duplicates,
                              double tol)
  {
    // Only good for surface grids.
    AssertThrow(dim<3, ExcImpossibleInDim(dim));

    // Create an istance of the Importer class
    Assimp::Importer importer;

    // And have it read the given file with some  postprocessing
    const aiScene *scene = importer.ReadFile( filename.c_str(),
                                              aiProcess_RemoveComponent   |
                                              aiProcess_JoinIdenticalVertices   |
                                              aiProcess_ImproveCacheLocality  |
                                              aiProcess_SortByPType   |
                                              aiProcess_OptimizeGraph   |
                                              aiProcess_OptimizeMeshes);

    // If the import failed, report it
    if ( !scene)
      {
        cerr << importer.GetErrorString() << endl;
        return false;
      }

    if ( scene->mNumMeshes == 0)
      {
        cout << "Input file contains no meshes." << endl;
        return false;
      }

    AssertThrow((mesh_index == -1) || (mesh_index < (int) scene->mNumMeshes),
                ExcMessage("Too few meshes in the file."));

    cout << "The given file contains "
         << scene->mNumMeshes  << " meshes. ";
    if (mesh_index == -1)
      cout << "Importing all available meshes." << endl;
    else
      cout << "Trying to import mesh number " << mesh_index << endl;

    int start_mesh = (mesh_index == -1 ? 0 : mesh_index);
    int end_mesh = (mesh_index == -1 ? (int) scene->mNumMeshes : mesh_index+1);



    // Deal.II objects are created empty, and then filled with imported file.
    vector<Point<spacedim> > vertices;
    vector<CellData<dim> > cells;
    SubCellData subcelldata;

    // A series of counters to merge cells.
    unsigned int v_offset=0;
    unsigned int c_offset=0;

    // The index of the mesh will be used as a material index.
    for (int m=start_mesh; m<end_mesh; ++m)
      {

        const aiMesh *mesh = scene->mMeshes[m];

        // Check that we know what to do with this mesh, otherwise just
        // ignore it
        if ( (dim == 2) && mesh->mPrimitiveTypes != aiPrimitiveType_POLYGON)
          {
            cout << "Skipping incompatible mesh " << m
                 << "/" << scene->mNumMeshes << "." << endl;
            continue;
          }
        else if ( (dim == 1) && mesh->mPrimitiveTypes != aiPrimitiveType_LINE)
          {
            cout << "Skipping incompatible mesh " << m
                 << "/" << scene->mNumMeshes << "." << endl;
            continue;
          }
        // Vertices
        const unsigned int n_vertices = mesh->mNumVertices;
        const aiVector3D *mVertices = mesh->mVertices;

        // Faces
        const unsigned int n_faces = mesh->mNumFaces;
        const aiFace *mFaces = mesh->mFaces;

        vertices.resize(v_offset+n_vertices);
        cells.resize(c_offset+n_faces);

        cout << "Input mesh has " << n_vertices << " vertices and "
             << n_faces << " faces." << endl;


        for (unsigned int i=0; i<n_vertices; ++i)
          for (unsigned int d=0; d<spacedim; ++d)
            vertices[i+v_offset][d] = mVertices[i][d];

        unsigned int valid_cell = c_offset;
        for (unsigned int i=0; i<n_faces; ++i)
          {
            if (mFaces[i].mNumIndices == GeometryInfo<dim>::vertices_per_cell)
              {
                cout << "Face: " << i << ": ";
                for (unsigned int f=0; f<GeometryInfo<dim>::vertices_per_cell; ++f)
                  {
                    cells[valid_cell].vertices[f] = mFaces[i].mIndices[f]+v_offset;
                    cout << cells[valid_cell].vertices[f] << " ";
                  }
                cout << endl;
                cells[valid_cell].material_id = (types::material_id) m;
                ++valid_cell;
              }
            else
              {
                cout << "Skipping face " << i << " of mesh "
                     << m << ", it has "
                     << mFaces[i].mNumIndices << " vertices. We expect "
                     << GeometryInfo<dim>::vertices_per_cell << endl;
              }
          }
        cells.resize(valid_cell);

        // The vertices are added all at once. Cells are check for
        // validity, so only valid_cells are now present in the deal.II
        // list of cells.
        v_offset += n_vertices;
        c_offset = valid_cell;
      }

    // No cells were read
    if (cells.size() == 0)
      return false;

    if (remove_duplicates)
      {
        vector<unsigned int> considered_vertices;
        GridTools::delete_duplicated_vertices(vertices, cells, subcelldata,
                                              considered_vertices, tol);

        GridTools::delete_duplicated_vertices(vertices, cells, subcelldata,
                                              considered_vertices, tol);

        GridTools::delete_duplicated_vertices(vertices, cells, subcelldata,
                                              considered_vertices, tol);

        GridTools::delete_duplicated_vertices(vertices, cells, subcelldata,
                                              considered_vertices, tol);

        GridTools::delete_duplicated_vertices(vertices, cells, subcelldata,
                                              considered_vertices, tol);
      }
    GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    if (dim == spacedim)
      GridReordering<dim, spacedim>::invert_all_cells_of_negative_grid(vertices,
          cells);

    cout << "Creating tria with " << vertices.size() << " vertices, "
         << cells.size() << "cells. " << std::endl;

    GridReordering<dim, spacedim>::reorder_cells(cells);
    if (dim == 2)
      tria.create_triangulation_compatibility(vertices, cells, subcelldata);
    else
      tria.create_triangulation(vertices, cells, subcelldata);

    cout << "Generated mesh has " << tria.n_vertices() << " vertices and "
         << tria.n_active_cells() << " active cells. " << endl;

    typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active(),
                                                               endc=tria.end();

    unsigned int boundary_faces = 0;
    for (; cell != endc; ++cell)
      if (cell->at_boundary())
        boundary_faces++;

    cout << "Triangulation has " << boundary_faces << " boundary cells." << endl;
    return true;

  }

  // Explicit instantiations
  template bool generate_triangulation(const string, Triangulation<1,1> &, int,
                                       bool, double);
  // Explicit instantiations
  template bool generate_triangulation(const string, Triangulation<1,2> &, int,
                                       bool, double);
  // Explicit instantiations
  template bool generate_triangulation(const string, Triangulation<1,3> &, int,
                                       bool, double);
  // Explicit instantiations
  template bool generate_triangulation(const string, Triangulation<2,2> &, int,
                                       bool, double);
  // Explicit instantiations
  template bool generate_triangulation(const string, Triangulation<2,3> &, int,
                                       bool, double);

}

D2K_NAMESPACE_CLOSE
#endif


