#/bin/bash

# FILENAME:      COMPILE.sh
# TO RUN DIGIT:  $ ./COMPILE.sh

# NOTE: Probably is necessary for you to change the path to the deal.ii folder
#       to do this, find below "DDEAL_II_DIR" and change "/usr/include/deal.II/".

echo "=============================================================================="
echo "      COMPILE: step-41      "
echo "=============================================================================="
echo " "
echo ">>> Press ENTER to run "
read -n1		# -n1 accept only 1 character
echo " "

# COMPILING COMMANDS:
#---------------------------------------------
mkdir build
cd build
cmake .. -DDEAL_II_DIR=/usr/include/deal.II/
#cmake .. -DDEAL_II_DIR=/usr/local/include/deal.II/
#cmake .. -DDEAL_II_DIR=/home/user100/dealII/dealii-9.2.0/include/deal.II/
printf "\n\n--------------------------------\n  make\n--------------------------------\n\n"
make
#---------------------------------------------
printf "\n\n--------------------------------\n\n"
echo "    COMPILATION ENDED CORRECTLY ! "
echo " "
echo "=============================================================================="

echo " "
echo ">>> Press ENTER to execute the generated code "
echo "    (Close the terminal to exit)"
read -n1		# -n1 accept only 1 character
echo " "

# RUNNING COMMANDS:
#---------------------------------------------
cd build
./exe
#mpiexec -np 2 exe
#---------------------------------------------

echo " "
echo "=============================================================================="
echo " "


