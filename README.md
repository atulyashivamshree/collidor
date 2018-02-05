# Collidor - Proximty distance query using GPUs
Collidor implements the proximity distance query between two objects using a GPU. 

## Installation
#### Requirements
* libccd
* libeigen3-dev
* FCL [https://github.com/flexible-collision-library/fcl](https://github.com/flexible-collision-library/fcl)

The install script can be used to install the dependencies and FCL directly
```
./install.sh
```

#### Compilation
To compile first set the location of the FCL in the setup.bash file and source it. 
```
export FCL_PATH=/home/atulya/Documents/fcl
```
**Note:** The above needs to be done for every new terminal instance whether for compiling or running the program.
 Then run makefile inside the build folder.

```
cd collidor
cd build
gedit setup.bash
source setup.bash
make
```

## Running
To run the code 
1. Download the CAD model of objects in .obj format
2. Either create a new config file or use one of the existing files in the params folder for example **cessna_cessna.yaml**. Edit the following variables
```
file1: ../CAD/cessna.obj
file2: ../CAD/cessna.obj
transforms: ../CAD/transforms_set2.csv
outp_prefix: cooper_cessna
```
Replace *file1* and *file2* with the.obj files to be used for comparision. The *transforms* variable provides a set of relative transforms in [X,Y,Z,roll,pitch,yaw] format. The file *transforms_set2.csv* contains 50 different 3D transformations. To run it on a smaller set use *transforms_set0.csv* Finally *outp_prefix* is the prefix for output files that are generated. Keep this the same as that of the config filename.
3. Source the environment variables for every new terminal instance
```
source setup.bash
```
This would export the shared libary path of FCL
4. Finally run the query using :
```
cd build
./dist_bvh.exe ../params/cessna_cessna.yaml
```

## Testing 
The program runs a number of test cases to verify the working:

* **Distance test** : This test takes in  a config file and computes the distance between two objects for multiple transformations as specified in the transforms file. It verifies if the computed distance from GPU is the same as that obtained from the FCL library and also shows the timing comparisions. Make sure you are working in a terminal where environment variables have been update (Step 3 of Running).
```
make  test_dist
``` 
If you wish to run tests on another config file you just need to provide in the name of the config file as a parameter. For example if the config file name is cessna_cessna.yaml just run
```
make  test_config=cessna_cessna test_dist
``` 

* **Unit tests** : Unit tests verify the working of the functions that are used inside the gpu algorithm. To run simply type in 
```
make unit_tests
```