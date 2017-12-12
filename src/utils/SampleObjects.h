/*
 * SampleObjects.h
 *
 *  Created on: Dec 12, 2017
 *      Author: atulya
 */

#ifndef UTILS_SAMPLEOBJECTS_H_
#define UTILS_SAMPLEOBJECTS_H_

void createSampleObj1(vector<fcl::Vector3<float>>& vertices,
        vector<fcl::Triangle>& triangles)
{
  vertices.push_back(fcl::Vector3<float>(0, 0, 0));
  vertices.push_back(fcl::Vector3<float>(0, 3, 0));
  vertices.push_back(fcl::Vector3<float>(3, 2, 0));
  vertices.push_back(fcl::Vector3<float>(5, 4, 0));
  vertices.push_back(fcl::Vector3<float>(7, 1, 0));

  triangles.push_back(fcl::Triangle(0, 1, 2));
  triangles.push_back(fcl::Triangle(3, 1, 2));
  triangles.push_back(fcl::Triangle(3, 4, 2));
}

void createSampleObj2(vector<fcl::Vector3<float>>& vertices,
        vector<fcl::Triangle>& triangles)
{
  vertices.push_back(fcl::Vector3<float>(8, 2, 0));
  vertices.push_back(fcl::Vector3<float>(10, 5, 0));
  vertices.push_back(fcl::Vector3<float>(12, 2, 0));
  vertices.push_back(fcl::Vector3<float>(11, 0, 0));

  triangles.push_back(fcl::Triangle(0, 1, 2));
  triangles.push_back(fcl::Triangle(0, 3, 2));
}

#endif /* UTILS_SAMPLEOBJECTS_H_ */
