# Isogeometric Analysis with Deal.II 

:octocat:  :electron: 

![extra](https://img.shields.io/badge/status-under%20development-ff69b4)
![extra](https://img.shields.io/badge/maintained-actively-success)

Isogeometric Analysis classes for the [deal.II library](https://www.dealii.org/).  

This repository contains code used to introduce the Isogeometric Analysis (IGA) framework into the standard deal.II Finite Element (FE) library.

<br/>  

---
### 🆕 Major updates:
Last update: **MAR-2021**  

From newest to oldest:  
* Note: Added Doxygen documentation.
* Obstacle example: :warning: BUG IN EXECUTION.
* Obstacle example: :heavy_check_mark: Fixed  (*The code now compiles without errors or warnings*).
* Poisson example : Improved (*Reorganization of the code, improved user interaction*).
* Poisson example : :heavy_check_mark: Fixed  (*The code now compiles without errors or warnings*).
* Note: A previous version of this repository has been completely deleted and substituted with this new one.

<br/>  

---
### 📄 Documentation:

Documentation files can be generated using [Doxygen](https://www.doxygen.nl/).  

To install Doxygen on Linux, please type: `sudo apt-get install doxygen`.  
If you prefere a Graphical User Interface (GUI): `sudo apt-get install doxygen-gui`.  [← suggested]  

If you have doxygen installed on your system, you can generate the documentation files using these commands:  

1. Open the Doxygen GUI interface by using this terminal command: `doxywizard`  
2. A window should open, **click on** "**File**" **>** "**Open**" and **select the file named** "*Doxyfile*" (in the "*doc*" folder).  
3. At this point some fileds in the GUI should have been automatically filled, **go to the** "**Run**" **tab**, and **press on** "**Run doxygen**"  
4. When the program completes, **click on** "**Show HTML output**" to consult the documentation.  
5. At this point you can close the Doxygen GUI.  

Note: To open the documentation at any time, just search for the folder "*html*" into the "*doc*" folder and open the file named "**index.html**".  

Note: For any other question, please refer to the Doxygen official website at this link: https://www.doxygen.nl/index.html.  


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

Note: At this point, if you want to recompile the code after some changes to the source code, 
it is sufficient to go to the *build* directory and use the command $`make`.


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

Note: At this point, if you want to recompile the code after some changes to the source code, 
it is sufficient to go to the *build* directory and use the command $`make`.


#### 👉 Instructions to **execute** the code:

The executable takes no arguments. To run, go to the "build" folder, and use the following command:
```
./exe
```

<br/>  

#### ↘️ Code output:

The main outputs are usually easy-to-consult text files reporting the convergence table and other execution data.  

The **output image files** are saved in **.vtk** format.  
The *Visualization Toolkit* (*VTK*) is an open-source software system for 3D computer graphics, image processing and scientific visualization.  

In order to open *.vtk* files you can use a graphical software like "**VisIt**".  

In order to install VisIt:  
1. Go to the project website: https://wci.llnl.gov/simulation/computer-codes/visit/executables;
2. Download the archive for your system (es: *Linux - x86_64 64 bit Ubuntu 20*) and download also the "*VisIt install script*";
3. Make the install script executable (es: $`chmod 755 visit-install3_1_4`);
4. Use the following command to install the program:  
$`sudo [script_name] [visit version] [System in use] [Installation path]`  
Example:  $`sudo visit-install3_1_4 3.1.4 linux-x86_64-ubuntu20 /usr/local/bin/visit`  
5. During the installation, if asked, select the remote-computing center nearest to your location and press *enter* until the installation has finished.

Steps to correctly see the .vtk results (for beginners):  
1. Open VisIt: `/usr/local/bin/visit/bin/visit`;
2. Click on "Open" > Select the directory > Select the \*.vtk database;
3. Now press on "Add" > "Pseudocolor" > "Active set";
4. Make the plot 3D: "Operators" > "Transform" > "Elevate";
5. Click on "**Draw**" and you should see a 3D movable picture in the Window.

Note: It is possible to add features to the plot as "Contour" or "Mesh".

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


