# P8: Kidnapped Vehicle Project

The project includes the following files

* CMakeLists.txt
* src folder
* notes.txt (this file)

## Description

The project was addressed the following way:

* Review the [QA session](https://www.youtube.com/watch?v=-3HI3Iw3Z9g).
* Start with samle code provided on the QA.
* Develop and debug the following methods in that order
  * init, change the initial distributions
  * prediction and check on the simulator
  * updateWeights and do the data association on the same method
* Try number of particles to balance speed vs efficiency: 10, 50, 75, 80, 100

## Troubleshooting

* Wrong estimation values
* The best explanation on what to do to update the weights was detailed [here](https://discussions.udacity.com/t/dont-get-the-point-of-transformation/248991/5)

## Test

Follow these steps for validation:

1. Run `./build.sh` to compile the project.
2. Run `./run.sh` to execute the particle_filter program
3. Run the simulator and Start

## Results

Final result detailed on the simulator

Error: x .655 y .401 yaw .006

