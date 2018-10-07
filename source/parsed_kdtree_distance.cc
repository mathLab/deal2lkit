#include <deal2lkit/parsed_kdtree_distance.h>

using namespace dealii;

D2K_NAMESPACE_OPEN


template <int dim>
ParsedKDTreeDistance<dim>::ParsedKDTreeDistance(
  const std::string &            name,
  const unsigned int &           max_leaf_size,
  const std::vector<Point<dim>> &pts)
  : ParameterAcceptor(name)
  , max_leaf_size(max_leaf_size)
{
  if (pts.size() > 0)
    set_points(pts);
}



template <int dim>
void
ParsedKDTreeDistance<dim>::declare_parameters(ParameterHandler &prm)
{
  add_parameter(
    prm,
    &max_leaf_size,
    "Max number of points per leaf",
    std::to_string(max_leaf_size),
    Patterns::Double(),
    "While building the tree, nodes are recursively divided until the number "
    "of points inside is equal or below this value. While doing queries, the "
    "tree algorithm ends by selecting leaf nodes, then performing linear search "
    "(one-by-one) for the closest point to the query within all those in the leaf.");
}



template <int dim>
double
ParsedKDTreeDistance<dim>::value(const Point<dim> & p,
                                 const unsigned int component) const
{
  // Only 1 component in this function
  AssertDimension(component, 0);
  Assert(adaptor, ExcNotInitialized());
  Assert(
    kdtree,
    ExcInternalError(
      "Adaptor is initialized, but kdtree is not. This should not happen."));

  // do a knn search
  const size_t                    num_results = 1;
  size_t                          ret_index;
  double                          dist;
  nanoflann::KNNResultSet<double> resultSet(num_results);

  resultSet.init(&ret_index, &dist);
  kdtree->findNeighbors(resultSet, &p[0], nanoflann::SearchParams());
  // index.knnSearch(query, indices, dists, num_results,
  // mrpt_flann::SearchParams(10));

  // silence a warning about unused variable
  // when compiling in release mode
  (void)component;
  return dist;
}



template <int dim>
unsigned int
ParsedKDTreeDistance<dim>::get_points_within_ball(
  const Point<dim> &                            center,
  const double &                                radius,
  std::vector<std::pair<unsigned int, double>> &matches,
  bool                                          sorted) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());

  Assert(radius > 0, ExcMessage("Radius is expected to be positive."));

  nanoflann::SearchParams params;
  params.sorted = sorted;
  return kdtree->radiusSearch(&center[0], radius, matches, params);
}

template <int dim>
void
ParsedKDTreeDistance<dim>::get_closest_points(
  const Point<dim> &         target,
  std::vector<unsigned int> &indices,
  std::vector<double> &      distances) const
{
  Assert(adaptor, ExcNotInitialized());
  Assert(kdtree, ExcInternalError());
  AssertDimension(indices.size(), distances.size());

  kdtree->knnSearch(&target[0], indices.size(), &indices[0], &distances[0]);
}

template <int dim>
void
ParsedKDTreeDistance<dim>::set_points(const std::vector<Point<dim>> &pts)
{
  Assert(pts.size() > 0, ExcMessage("Expecting a non zero set of points."));
  adaptor = SP(new PointCloudAdaptor(pts));
  kdtree  = SP(new KDTree(
    dim, *adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size)));
  kdtree->buildIndex();
}


template class ParsedKDTreeDistance<1>;
template class ParsedKDTreeDistance<2>;
template class ParsedKDTreeDistance<3>;

D2K_NAMESPACE_CLOSE
