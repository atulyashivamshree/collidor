# Collidor
Implements collision check on two objects

## File Description
#### CUDA Executables
These are the executables that would do the main computation on CUDA. Look up the Makefile inside **cuda_build** to get an idea of what they are doing

1. **dist_bvh.exe** : Performs the main distance computation and prints out the results. This is independent of the FCL library and can run on any system.
` ./dist_bvh.exe FILE1.bvh FILE2.bvh OUTPUT_PREFIX`

2. **test_triangles.exe** : Performs some basic tests on triangles and runs the code on a large number of triangles to verifiy accuracy. This might need the FCL library to be installed. It would also be useful for running the performance evaluation on a large number of triangles for testing. 
3. **test_rss.exe** : Similar to triangles test but in this case runs the test on RSS. Again this would also need the FCL library

TODO : perform O1 or O2 optimization

#### Non-CUDA executables
These are helper files meant to perform conversions on the different file formats or run tests. These would require the FCL library as a dependency. Use the CMake file inside **src** to generate and run these executables.

1.  **save_to_bvh** : Converts the obj file to a BVH file which is read as an input by the cuda executable.
2.  **compute_dist_bvh** : TODO this is meant to take in two random files as an input and compute the distance in between them over different value of the transformation matrix
3. **verify_dist** : Takes in the two original .obj files and the .outp.csv files as input. It then verifies that the .outp.csv file has correctly measured the distance between the primitive triangles and RSS
4.  **load_files**, **test_fcl_utility** :  sample files to verify the working of the FCL library
5.  **test_rss**, **test_vector** : tests the working of the RSS, triangles and the standard vector library