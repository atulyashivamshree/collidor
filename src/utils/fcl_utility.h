/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2011-2014, Willow Garage, Inc.
 *  Copyright (c) 2014-2016, Open Source Robotics Foundation
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Open Source Robotics Foundation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** @author Jia Pan

    MODIFIED BY Atulya Shivam Shree
    Date : Dec 12, 2017
    */

#ifndef SRC_UTILS_FCL_UTILITY_H_
#define SRC_UTILS_FCL_UTILITY_H_

#include <array>
#include <fstream>
#include <iostream>

#include "fcl/common/unused.h"

#include "fcl/math/constants.h"
#include "fcl/math/triangle.h"

#include "fcl/geometry/bvh/BVH_model.h"
#include "fcl/geometry/shape/box.h"
#include "fcl/geometry/shape/cylinder.h"
#include "fcl/geometry/shape/sphere.h"

#include "../BVH.h"

namespace fcl {

namespace test {

/// @brief Load an obj mesh file
template <typename S>
HOST_PREFIX void loadOBJFile(const char* filename,
                             std::vector<Vector3<S>>& points,
                             std::vector<Triangle>& triangles);

template <typename S>
HOST_PREFIX void saveOBJFile(const char* filename,
                             std::vector<Vector3<S>>& points,
                             std::vector<Triangle>& triangles);

/// @brief Structure for minimum distance between two meshes and the
/// corresponding nearest point pair
template <typename S>
struct DistanceRes {
  S distance;
  Vector3<S> p1;
  Vector3<S> p2;
};

}  // namespace test
}  // namespace fcl

HOST_PREFIX void convertFCLToBVH(
    const fcl::BVHModel<fcl::RSS<float>>& bvh_fcl,
    const std::vector<fcl::Vector3<float>>& vertices,
    const std::vector<fcl::Triangle>& triangles, BVH* bvh);

HOST_PREFIX void loadOBJToBVH(const std::string filename, BVH* bvh);

#include "fcl_utility-inl.h"

#endif  // SRC_UTILS_FCL_UTILITY_H_
