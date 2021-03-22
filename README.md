# Isogeometric Analysis with Deal.II 

:octocat:  :electron: 

![extra](https://img.shields.io/badge/status-under%20development-ff69b4)
![extra](https://img.shields.io/badge/maintained-actively-success)

Isogeometric Analysis classes for the [deal.II library](https://www.dealii.org/).  

This repository contains code used to introduce the Isogeometric Analysis (IGA) framework into the standard deal.II Finite Element (FE) library.

<br/>  

---
### 🆕 Recent updates:
Date: **MAR-2021**

* Note: A previous version of this repository has been completely deleted and substituted with this new one.
* Poisson example : :heavy_check_mark: Fixed  (*The code now compiles without errors or warnings*)
* Poisson example : Improved (*Reorganization of the code, improved user interaction*)
* Obstacle example : :heavy_check_mark: Fixed  (*The code now compiles without errors or warnings*)
* Obstacle example : :warning: BUG IN EXECUTION

<br/>  

---
### ℹ️ Usage:

<br/>  

### :triangular_flag_on_post: Poisson example

The code of the **Poisson example** consist of a single file:
* "poisson.cc"

Note: The file "CMakeLists.txt" is used to compile the code.


#### 👉 Instructions to **compile** the code:

```
cd poisson
mkdir build
cd build
cmake .. -DDEAL_II_DIR=/path/to/dealii/installation
make
```
NOTE: The default installation folder of deal.II library should be: `/usr/include/deal.II/`

If everything went well, you should have an executable named `poisson` in the new *build* directory. 


#### 👉 Instructions to **execute** the code:

The executable takes 0 or 5 arguments to run, respectively:
```
./poisson
```
or
```
./poisson  FE_TYPE  QUADRATURE  DEG  CYCLE_START  CYCLE_STOP
```
The accepted arguments are:

* FE_TYPE:
 	* `bernstein`  [default]
 	* `lagrange`
 	* `lobatto`
* QUADRATURE:
 	* `legendre`   [default]
 	* `lobatto`
* DEG:         (es: ` 1 `) degree of the finite element space
* CYCLE_START: (es: ` 0 `) initial refinement of the grid
* CYCLE_STOP:  (es: ` 5 `) final refinement of the grid  

Example: ` ./poisson bernstein legendre 1 0 5 `


<br/>  

### :triangular_flag_on_post: Obstacle example

The code of the **Obstacle example** consist of a series of files incuded in two folders:
* "source" directory
* "include" directory

Note: The file "CMakeLists.txt" is used to compile the code.


#### 👉 Instructions to **compile** the code:

```
cd obstacle
mkdir build
cd build
cmake .. -DDEAL_II_DIR=/path/to/dealii/installation
make
```
NOTE: The default installation folder of deal.II library should be: `/usr/include/deal.II/`

If everything went well, you should have an executable named `exe` in the new *build* directory. 


#### 👉 Instructions to **execute** the code:

The executable takes no arguments. To run use the following command:
```
./exe
```

<br/>  

---
### 🛠️ Installation

Note: The codes have been tested under Xubuntu, version 20.04 LTS.


These examples require **deal.II** version 8.3 or later to work properly.  
Note: The codes have been tested using the version **deal.II 9.2.0**.

In order to setup the environment, follow the following steps:  

1. Install a C++ compiler:  
   $ `gcc -v`                     check the installed version of gcc (9.3.0)  
   $ `clang -v`                   check the installed version of clang (10.0.0)  
   $ `sudo apt install gcc`       install the last version of gcc  
   $ `sudo apt install clang`     install the last version of clang  

2. Update all packages before start the deal.II installation procedure:  
   $ `sudo apt-get update`  

3. Install the deal.II library:  
   $ `sudo apt-get install libdeal.ii-dev`  
  Note: This command will install the main deal.II library, along with all its suggested additional packages:   
  (es: *boost, fftw3, superlu, trilinos, open mpi, muparser, occt, ...* )  
  The required space is about ~ 2 GB for all.  


For any other information, please see the following links:  
[www.dealii.org / Main page](https://www.dealii.org/)  
[www.dealii.org / Installation instructions](https://www.dealii.org/current/readme.html)  

<br/>  

---
### ✍️ Contributors:  

The actual version is the result of an arrangement and improvement of an original version developed by the authors of the article:  
"*Algorithms, data structures and applications for Isogeometric Analysis with the deal.II library*"  
By *Marco Tezzele*, *Nicola Cavallini*, *Luca Heltai*  
SISSA - International School for Advanced Studies  


---
### ©️ License:  
Unless stated otherwise all works are licensed under the:  
*GNU Lesser General Public License v2.1*  

Please see the attached license file.  


