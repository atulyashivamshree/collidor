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

#ifndef SRC_UTILS_FCL_UTILITY_INL_H_
#define SRC_UTILS_FCL_UTILITY_INL_H_

//============================================================================//
//                                                                            //
//                              Implementations                               //
//                                                                            //
//============================================================================//

namespace fcl {

namespace test {

//==============================================================================
template <typename S>
HOST_PREFIX void loadOBJFile(const char* filename,
                             std::vector<Vector3<S>>& points,
                             std::vector<Triangle>& triangles) {
  FILE* file = fopen(filename, "rb");
  if (!file) {
    std::cerr << "file not exist" << std::endl;
    return;
  }

  bool has_normal = false;
  bool has_texture = false;
  char line_buffer[2000];
  while (fgets(line_buffer, 2000, file)) {
    char* first_token = strtok(line_buffer, "\r\n\t ");
    if (!first_token || first_token[0] == '#' || first_token[0] == 0) continue;

    switch (first_token[0]) {
      case 'v': {
        if (first_token[1] == 'n') {
          strtok(nullptr, "\t ");
          strtok(nullptr, "\t ");
          strtok(nullptr, "\t ");
          has_normal = true;
        } else if (first_token[1] == 't') {
          strtok(nullptr, "\t ");
          strtok(nullptr, "\t ");
          has_texture = true;
        } else {
          S x = (S)atof(strtok(nullptr, "\t "));
          S y = (S)atof(strtok(nullptr, "\t "));
          S z = (S)atof(strtok(nullptr, "\t "));
          points.emplace_back(x, y, z);
        }
      } break;
      case 'f': {
        Triangle tri;
        char* data[30];
        int n = 0;
        while ((data[n] = strtok(nullptr, "\t \r\n")) != nullptr) {
          if (strlen(data[n])) n++;
        }

        for (int t = 0; t < (n - 2); ++t) {
          if ((!has_texture) && (!has_normal)) {
            tri[0] = atoi(data[0]) - 1;
            tri[1] = atoi(data[1]) - 1;
            tri[2] = atoi(data[2]) - 1;
          } else {
            const char* v1;
            for (int i = 0; i < 3; i++) {
              // vertex ID
              if (i == 0)
                v1 = data[0];
              else
                v1 = data[t + i];

              tri[i] = atoi(v1) - 1;
            }
          }
          triangles.push_back(tri);
        }
      }
    }
  }
}

//==============================================================================
template <typename S>
HOST_PREFIX void saveOBJFile(const char* filename,
                             std::vector<Vector3<S>>& points,
                             std::vector<Triangle>& triangles) {
  std::ofstream os(filename);
  if (!os) {
    std::cerr << "file not exist" << std::endl;
    return;
  }

  for (std::size_t i = 0; i < points.size(); ++i) {
    os << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2]
       << std::endl;
  }

  for (std::size_t i = 0; i < triangles.size(); ++i) {
    os << "f " << triangles[i][0] + 1 << " " << triangles[i][1] + 1 << " "
       << triangles[i][2] + 1 << std::endl;
  }

  os.close();
}

}  // namespace test
}  // namespace fcl

HOST_PREFIX Vector3 getVec3Custom(const fcl::Vector3<float>& vec_fcl) {
  return Vector3({vec_fcl(0), vec_fcl(1), vec_fcl(2)});
}

HOST_PREFIX Triangle
getTriangleCustom(const std::vector<fcl::Vector3<float>>& vertices,
                  const fcl::Triangle& tri_fcl) {
  Triangle obj;
  obj.a = getVec3Custom(vertices[tri_fcl[0]]);
  obj.b = getVec3Custom(vertices[tri_fcl[1]]);
  obj.c = getVec3Custom(vertices[tri_fcl[2]]);
  return obj;
}

HOST_PREFIX RSS getRSSCustom(const fcl::RSSf& rss) {
  RSS obj;
  obj.axis.v1.x = rss.axis(0, 0);
  obj.axis.v1.y = rss.axis(1, 0);
  obj.axis.v1.z = rss.axis(2, 0);
  obj.axis.v2.x = rss.axis(0, 1);
  obj.axis.v2.y = rss.axis(1, 1);
  obj.axis.v2.z = rss.axis(2, 1);
  obj.axis.v3.x = rss.axis(0, 2);
  obj.axis.v3.y = rss.axis(1, 2);
  obj.axis.v3.z = rss.axis(2, 2);

  obj.To.x = rss.To(0);
  obj.To.y = rss.To(1);
  obj.To.z = rss.To(2);

  obj.l[0] = rss.l[0];
  obj.l[1] = rss.l[1];

  obj.r = rss.r;

  obj.size = rss.size();

  return obj;
}

HOST_PREFIX void convertFCLToBVH(
    const fcl::BVHModel<fcl::RSS<float>>& bvh_fcl,
    const std::vector<fcl::Vector3<float>>& vertices,
    const std::vector<fcl::Triangle>& triangles, BVH* bvh) {
  bvh->num_bv = bvh_fcl.getNumBVs();
  bvh->num_tri = triangles.size();
  bvh->bv_arr = new BV[bvh->num_bv];
  bvh->tri_arr = new Triangle[bvh->num_tri];

  // copy over all RSS objects and their tree heirarchy
  for (int i = 0; i < bvh_fcl.getNumBVs(); i++) {
    bvh->bv_arr[i].rss = getRSSCustom(bvh_fcl.getBV(i).bv);
    bvh->bv_arr[i].id1 = bvh_fcl.getBV(i).leftChild();
    bvh->bv_arr[i].id2 = bvh_fcl.getBV(i).rightChild();
    bvh->bv_arr[i].idt = bvh_fcl.getBV(i).primitiveId();
  }

  for (int i = 0; i < triangles.size(); i++) {
    bvh->tri_arr[i] = getTriangleCustom(vertices, triangles[i]);
  }
}

HOST_PREFIX void loadOBJToBVH(const std::string filename, BVH* bvh) {
  std::vector<fcl::Vector3f> points;
  std::vector<fcl::Triangle> triangles;

  bool file_exists = std::ifstream(filename.c_str()).good();
  if (!file_exists) {
    std::cerr << "ERROR " << filename << " could not be found " << endl;
    std::exit(EXIT_FAILURE);
  }

  fcl::test::loadOBJFile(filename.c_str(), points, triangles);

  fcl::BVHModel<fcl::RSS<float>> model;
  model.bv_splitter.reset(new fcl::detail::BVSplitter<fcl::RSS<float>>(
      fcl::detail::SPLIT_METHOD_MEAN));

  model.beginModel();
  model.addSubModel(points, triangles);
  model.endModel();

  convertFCLToBVH(model, points, triangles, bvh);
}

#endif  // SRC_UTILS_FCL_UTILITY_INL_H_
