#include <iostream>
#include <cstdlib>

//#define FCL_HAVE_OCTOMAP 0

#include "Eigen/Dense"
#include "fcl/geometry/bvh/detail/BVH_front.h"
#include "fcl/narrowphase/detail/traversal/collision_node.h"
#include "test_fcl_utility.h"
//#include "fcl_resources/config.h"



using namespace std;
using namespace fcl;

template <typename S>
void test_mesh_distance();

int main()
{
    cout << "yetti" << endl;
    Triangle tri;
    test_mesh_distance<double>();
}

bool verbose = false;

template <typename S>
S DELTA() { return 0.001; }


template<typename BV>
void distance_Test(const Transform3<typename BV::S>& tf,
                   const std::vector<Vector3<typename BV::S>>& vertices1, const std::vector<Triangle>& triangles1,
                   const std::vector<Vector3<typename BV::S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method,
                   int qsize,
                   test::DistanceRes<typename BV::S>& distance_result,
                   bool verbose = true);

template <typename S>
bool collide_Test_OBB(const Transform3<S>& tf,
                      const std::vector<Vector3<S>>& vertices1, const std::vector<Triangle>& triangles1,
                      const std::vector<Vector3<S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method, bool verbose);

template<typename BV, typename TraversalNode>
void distance_Test_Oriented(const Transform3<typename BV::S>& tf,
                            const std::vector<Vector3<typename BV::S>>& vertices1, const std::vector<Triangle>& triangles1,
                            const std::vector<Vector3<typename BV::S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method,
                            int qsize,
                            test::DistanceRes<typename BV::S>& distance_result,
                            bool verbose = true);

template <typename S>
bool nearlyEqual(const Vector3<S>& a, const Vector3<S>& b)
{
  if(fabs(a[0] - b[0]) > DELTA<S>()) return false;
  if(fabs(a[1] - b[1]) > DELTA<S>()) return false;
  if(fabs(a[2] - b[2]) > DELTA<S>()) return false;
  return true;
}
template <typename S>
void test_mesh_distance()
{
  std::vector<Vector3<S>> p1, p2;
  std::vector<Triangle> t1, t2;

  test::loadOBJFile("../env.obj", p1, t1);
  test::loadOBJFile("../rob.obj", p2, t2);

  aligned_vector<Transform3<S>> transforms; // t0
  S extents[] = {-3000, -3000, 0, 3000, 3000, 3000};
#ifdef NDEBUG
  std::size_t n = 10;
#else
  std::size_t n = 1;
#endif

  test::generateRandomTransforms(extents, transforms, n);

  double dis_time = 0;
  double col_time = 0;

  test::DistanceRes<S> res, res_now;
  for(std::size_t i = 0; i < transforms.size(); ++i)
  {
    test::Timer timer_col;
    timer_col.start();
    collide_Test_OBB(transforms[i], p1, t1, p2, t2, detail::SPLIT_METHOD_MEAN, verbose);
    timer_col.stop();
    col_time += timer_col.getElapsedTimeInSec();

    test::Timer timer_dist;
    timer_dist.start();
    distance_Test_Oriented<RSS<S>, detail::MeshDistanceTraversalNodeRSS<S>>(transforms[i], p1, t1, p2, t2, detail::SPLIT_METHOD_MEAN, 2, res, verbose);
    timer_dist.stop();
    dis_time += timer_dist.getElapsedTimeInSec();

    distance_Test_Oriented<RSS<S>, detail::MeshDistanceTraversalNodeRSS<S>>(transforms[i], p1, t1, p2, t2, detail::SPLIT_METHOD_BV_CENTER, 2, res_now, verbose);

    
  }

  std::cout << "distance timing: " << dis_time << " sec" << std::endl;
  std::cout << "collision timing: " << col_time << " sec" << std::endl;
}

template <typename S>
bool collide_Test_OBB(const Transform3<S>& tf,
                      const std::vector<Vector3<S>>& vertices1, const std::vector<Triangle>& triangles1,
                      const std::vector<Vector3<S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method, bool verbose)
{
  BVHModel<OBB<S>> m1;
  BVHModel<OBB<S>> m2;
  m1.bv_splitter.reset(new detail::BVSplitter<OBB<S>>(split_method));
  m2.bv_splitter.reset(new detail::BVSplitter<OBB<S>>(split_method));

  m1.beginModel();
  m1.addSubModel(vertices1, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  CollisionResult<S> local_result;
  detail::MeshCollisionTraversalNodeOBB<S> node;
  if(!detail::initialize(node, (const BVHModel<OBB<S>>&)m1, tf, (const BVHModel<OBB<S>>&)m2, Transform3<S>::Identity(),
                 CollisionRequest<S>(), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  collide(&node);

  if(local_result.numContacts() > 0)
    return true;
  else
    return false;
}

template<typename BV, typename TraversalNode>
void distance_Test_Oriented(const Transform3<typename BV::S>& tf,
                            const std::vector<Vector3<typename BV::S>>& vertices1, const std::vector<Triangle>& triangles1,
                            const std::vector<Vector3<typename BV::S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method,
                            int qsize,
                            test::DistanceRes<typename BV::S>& distance_result,
                            bool verbose)
{
  using S = typename BV::S;

  BVHModel<BV> m1;
  BVHModel<BV> m2;
  m1.bv_splitter.reset(new detail::BVSplitter<BV>(split_method));
  m2.bv_splitter.reset(new detail::BVSplitter<BV>(split_method));


  m1.beginModel();
  m1.addSubModel(vertices1, triangles1);
  m1.endModel();

  m2.beginModel();
  m2.addSubModel(vertices2, triangles2);
  m2.endModel();

  DistanceResult<S> local_result;
  TraversalNode node;
  if(!initialize(node, (const BVHModel<BV>&)m1, tf, (const BVHModel<BV>&)m2, Transform3<S>::Identity(), DistanceRequest<S>(true), local_result))
    std::cout << "initialize error" << std::endl;

  node.enable_statistics = verbose;

  distance(&node, nullptr, qsize);

  // points are in local coordinate, to global coordinate
  Vector3<S> p1 = local_result.nearest_points[0];
  Vector3<S> p2 = local_result.nearest_points[1];


  distance_result.distance = local_result.min_distance;
  distance_result.p1 = p1;
  distance_result.p2 = p2;

  if(verbose)
  {
    std::cout << "distance " << local_result.min_distance << std::endl;

    std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
    std::cout << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
    std::cout << node.num_bv_tests << " " << node.num_leaf_tests << std::endl;
  }
}
