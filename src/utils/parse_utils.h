/*
 *  parse_utils.h
    Description: utility functions for loading and parsing files

    @author Atulya Shivam Shree
    Created on: Dec 12, 2017
    Copyright (c) 2017 Atulya Shivam Shree
 */

#ifndef SRC_UTILS_PARSE_UTILS_H_
#define SRC_UTILS_PARSE_UTILS_H_

using std::map;
using std::string;
using std::vector;

// REQUIRES : filename of the transformations file
// MODIFIES : the vector of transforms
// EFFECTS  : adds in new transformations to the vector
HOST_PREFIX void loadTransformations(vector<Transform3f>& transforms,
                                     const string filename);

// REQUIRES : filename of the config file
// MODIFIES : conf
//            cfg
// EFFECTS  : loads in all parameteres and updates the Config data structure
HOST_PREFIX void loadConfig(map<string, string>& conf, const string filename);

// EFFECTS  : compares if two numbers are equal within a threshold
HOST_PREFIX bool approxEquals(float a, float b, float EPSILON = 1e-4) {
  return (fabs(a - b) < EPSILON);
}

HOST_PREFIX void loadTransformations(vector<Transform3f>& transforms,
                                     const string filename) {
  std::ifstream fin;
  fin.open(filename.c_str());
  if (!fin.is_open()) {
    cout << "TF: " << filename << " could not be opened" << endl;
    exit(EXIT_FAILURE);
  }
  string header_str;
  for (int i = 0; i < 6; i++) fin >> header_str;

  float roll, pitch, yaw, x, y, z;

  // read in RPY, XYZ values from the file and store them as a transform
  while (fin >> roll >> pitch >> yaw >> x >> y >> z) {
    Eigen::Vector3f pos(x, y, z);

    Eigen::AngleAxisf rollAngle(roll, Eigen::Vector3f::UnitZ());
    Eigen::AngleAxisf yawAngle(pitch, Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf pitchAngle(yaw, Eigen::Vector3f::UnitX());
    Eigen::Quaternion<float> q = rollAngle * yawAngle * pitchAngle;
    Eigen::Matrix3f R = q.matrix();
    Transform3f tf;
    tf.linear() = q.matrix();
    tf.translation() = pos;
    transforms.push_back(tf);
  }
}

HOST_PREFIX void loadConfig(map<string, string>& conf, const string filename) {
  std::ifstream fin;
  fin.open(filename.c_str());
  if (!fin.is_open()) {
    cout << "CFG: " << filename << " could not be opened" << endl;
    exit(EXIT_FAILURE);
  }

  string key, content;
  while (fin >> key >> content) {
    if (key[key.size() - 1] == ':') {
      string temp_key = key.substr(0, key.size() - 1);
      conf[temp_key] = content;
    }
  }

  cout << "Parameters:" << endl;
  for (const auto& p : conf) cout << p.first << ": " << p.second << endl;
}

#endif  // SRC_UTILS_PARSE_UTILS_H_
