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

#ifndef SRC_UTILS_TEST_FCL_UTILITY_H_
#define SRC_UTILS_TEST_FCL_UTILITY_H_

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

namespace fcl {

namespace test {

class Timer {
 public:
  Timer();
  ~Timer();

  void start();                       ///< start timer
  void stop();                        ///< stop the timer
  double getElapsedTime();            ///< get elapsed time in milli-second
  double getElapsedTimeInSec();       ///< get elapsed time in second (same as
                                      ///< getElapsedTime)
  double getElapsedTimeInMilliSec();  ///< get elapsed time in milli-second
  double getElapsedTimeInMicroSec();  ///< get elapsed time in micro-second

 private:
  double startTimeInMicroSec;  ///< starting time in micro-second
  double endTimeInMicroSec;    ///< ending time in micro-second
  int stopped;                 ///< stop flag
#ifdef _WIN32
  LARGE_INTEGER frequency;  ///< ticks per second
  LARGE_INTEGER startCount;
  LARGE_INTEGER endCount;
#else
  timeval startCount;
  timeval endCount;
#endif
};

struct TStruct {
  std::vector<double> records;
  double overall_time;

  TStruct() { overall_time = 0; }

  void push_back(double t) {
    records.push_back(t);
    overall_time += t;
  }
};

/// @brief Load an obj mesh file
template <typename S>
void loadOBJFile(const char* filename, std::vector<Vector3<S>>& points,
                 std::vector<Triangle>& triangles);

template <typename S>
void saveOBJFile(const char* filename, std::vector<Vector3<S>>& points,
                 std::vector<Triangle>& triangles);

/// @brief Structure for minimum distance between two meshes and the
/// corresponding nearest point pair
template <typename S>
struct DistanceRes {
  S distance;
  Vector3<S> p1;
  Vector3<S> p2;
};

/// @brief Distance data stores the distance request and the result given by
/// distance algorithm.
template <typename S>
struct DistanceData {
  DistanceData() { done = false; }

  /// @brief Distance request
  DistanceRequest<S> request;

  /// @brief Distance result
  DistanceResult<S> result;

  /// @brief Whether the distance iteration can stop
  bool done;
};

//============================================================================//
//                                                                            //
//                              Implementations                               //
//                                                                            //
//============================================================================//

//==============================================================================
template <typename S>
void loadOBJFile(const char* filename, std::vector<Vector3<S>>& points,
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
void saveOBJFile(const char* filename, std::vector<Vector3<S>>& points,
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

#endif  // SRC_UTILS_TEST_FCL_UTILITY_H_
