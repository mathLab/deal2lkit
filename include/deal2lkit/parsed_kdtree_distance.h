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

#ifndef _d2k_parsed_kdtree_distance_h
#define _d2k_parsed_kdtree_distance_h

#include <deal2lkit/config.h>
#include <deal2lkit/utilities.h>
#include <deal.II/base/function.h>
#include <deal2lkit/parameter_acceptor.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <deal2lkit/nanoflann.h>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

using namespace dealii;

D2K_NAMESPACE_OPEN

/**
 * A deal2lkit wrapper for the nanoflann library, used to compute the
 * distance from a collection of points. This function uses nanoflann
 * to efficiently partition the space in a tree. The cost of each
 * query is then roughly of order log(n), where n is the number of
 * points stored in this class.
 *
 * Additional methods are provided that give access to other
 * functionalities of the nanoflann library, like searching the n
 * nearest neighbors, or searching the points that fall within a
 * raidius of a target point.
 */
template<int dim>
class ParsedKDTreeDistance : public ParameterAcceptor, public Function<dim>
{
public:
  /**
   * Constructor: takes an optional name for the section. If the
   * optional expression string is given, than it is used to set the
   * expression as soon as the parameters are declared.
   *
   * The max leaf parameter is used to decide how many points per leaf
   * are used in the kdtree algorithm.
   *
   * If the points are not passed to this constructor, then you have
   * to pass them later to this object by calling the set_points()
   * method.
   *
   * Access to any of the methods without first passing a reference to
   * a vector of points will result in an exception. Only a reference
   * to the points is stored, so you should make sure that the life of
   * the the vector you pass is longer than the life of this class, or
   * you'll get undefinite behaviour.
   *
   * If you update the vector of points in someway, remember to call
   * again the set_points() method. The tree and the index are
   * constructed only once, when you pass the points (either at
   * construction time, or when you call set_points()). If you update
   * your points, and do not call again set_points(), your results
   * will likely be wrong.
   */
  ParsedKDTreeDistance(const std::string &name="",
                       const unsigned int &max_leaf_size=10,
                       const std::vector<Point<dim> > &pts=std::vector<Point<dim> >());


  /**
   * Adaptor class used internally by nanoflann. Class actually stores
   * a reference to the vector of points, and generates some helper
   * functions for nanoflann.
   */
  struct PointCloudAdaptor
  {
    /**
     * A typedef used by nanoflann.
     */
    typedef double coord_t;


    /**
     * Reference to the vector of points from which we want to compute
     * the distance.
     */
    const std::vector<Point<dim> > &points; //!< A const ref to the data set origin


    /**
     * The constrcutor needs the data set source.
     */
    PointCloudAdaptor(const std::vector<Point<dim> > &_points) : points(_points) { }


    /**
     * Return number of points in the data set (required by nanoflann).
     */
    inline size_t kdtree_get_point_count() const
    {
      return points.size();
    }


    /**
     * Return the L2 distance between points
     */
    inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2,size_t size) const
    {
      AssertDimension(size, dim);
      coord_t res=0.0;
      for (size_t d=0; d<size; ++d)
        res += (p1[d]-points[idx_p2][d])*(p1[d]-points[idx_p2][d]);
      return std::sqrt(res);
    }


    /**
     * Return the dim'th component of the idx'th point in the class.
     */
    inline coord_t kdtree_get_pt(const size_t idx, int d) const
    {
      AssertIndexRange(d,dim);
      return points[idx][d];
    }


    /**
     * Optional bounding-box computation: return false to default to a
     * standard bbox computation loop.  Return true if the BBOX was
     * already computed by the class and returned in "bb" so it can be
     * avoided to redo it again.  Look at bb.size() to find out the
     * expected dimensionality (e.g. 2 or 3 for point clouds).
     */
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &) const
    {
      return false;
    }
  };


  /**
   * A typedef for the actual KDTree object.
   */
  typedef typename nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloudAdaptor> ,
          PointCloudAdaptor, dim, unsigned int>  KDTree;

  /**
   * Calls the underlying function of ParsedFunction.
   */
  virtual void declare_parameters(ParameterHandler &prm);


  /**
   * Compute the distance between the given point and the cloud of
   * points in the currently stored adaptor.
   *
   * If the adaptor is empty, throw an exception. You should call the
   * function set_points() to initialize a new adaptor, or call the
   * constructor with a valid set of points. If your points change,
   * call again the set_points() method.
   *
   * @param p: the point where we want to evaluate the distance
   *
   * @param component: function component (throws an exception if it
   * is different from zero. This is required by the dealii interface)
   *
   * @return the minimum distance from the collection of points set
   * with the constructor or with the method set_points()
   */
  virtual double value(const Point<dim> &p, const unsigned int component=0) const;


  /**
   * Store a reference to the passed points. After you called this
   * method, you can call the value() method to compute the minimum
   * distance between an evaluation point and the collection of points
   * you passed to this method, or the get_points_within_ball() and
   * the get_closest_points() methods.
   *
   * Notice that the constructor calls this method internally if you
   * pass it a non empty vector of points.
   *
   * Whenever your points change, you should call this method again,
   * since this is the method responsible for building the index and
   * storing the actual tree internally. If you change your points and
   * don't call again this method, any function you call later will
   * happily return wrong values without you noticing.
   *
   * @param pts: a collection of points
   */
  void set_points(const std::vector<Point<dim> > &pts);


  /**
   * A const accessor to the underlying points.
   */
  inline const Point<dim> &operator[](unsigned int i) const
  {
    AssertIndexRange(i, size());
    return adaptor->points[i];
  }


  /**
   * The size of the vector stored by this class.
   */
  inline unsigned int size() const
  {
    if (adaptor)
      return adaptor->points.size();
    else
      return 0;
  };


  /**
   * Fill a vector with the indices and the distance of the points
   * that are at distance less than or equal to the given radius from
   * the target point. Consider preallocating the size of the return
   * vector if you have a wild guess of how many should be there.
   *
   * @param[in] point: the target point
   * @param[in] radius: the radius of the ball
   * @param[out] mathes: indices and distances of the matching points
   * @param[in] sorted: sort the output results in ascending order with respect to distances
   *
   * @return number of points that are within the ball
   */
  unsigned int get_points_within_ball(const Point<dim> &target, const double &radius,
                                      std::vector<std::pair<unsigned int, double> > &matches,
                                      bool sorted=false) const;

  /**
   * Fill two vectors with the indices and distances of the closest
   * points to the given target point. The vectors are filled with
   * indices and distances until there is space in them. You should
   * resize them to the number of closest points you wish to get. An
   * assertion is thrown if the vectors do not have the same size.
   *
   * @param[in] target: the target point
   * @param[out] indices: indices of the matching points
   * @param[out] distances: distances of the matching points
   */
  void get_closest_points(const Point<dim> &target,
                          std::vector<unsigned int> &indices,
                          std::vector<double> &distances) const;

private:
  /**
   * Max number of points per leaf.
   */
  unsigned int max_leaf_size;


  /**
   * A point cloud adaptor, to be filled when set points is called.
   */
  shared_ptr<PointCloudAdaptor> adaptor;


  /**
   * The actual kdtree.
   */
  shared_ptr<KDTree> kdtree;
};

D2K_NAMESPACE_CLOSE

#endif
