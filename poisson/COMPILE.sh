#/bin/bash

# FILENAME:      COMPILE.sh
# TO RUN DIGIT:  $ ./COMPILE.sh

# NOTE: Probably is necessary for you to change the path to the deal.ii folder
#       to do this, find below "DDEAL_II_DIR" and change "/usr/include/deal.II/".

echo "=============================================================================="
echo "      COMPILE: poisson      "
echo "=============================================================================="
echo " "
echo ">>> Press ENTER to run "
read -n1		# -n1 accept only 1 character
echo " "

# COMPILING COMMANDS:
#---------------------------------------------
mkdir build
cd build
#cmake ..
#cmake .. -DDEAL_II_DIR=/usr/include/deal.II/                    # Default location
cmake .. /u/sw/pkgs/toolchains/gcc-glibc/9/pkgs/dealii/9.2.0/   # PoliMi mkModules

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
./poisson
#---------------------------------------------

echo " "
echo "=============================================================================="
echo " "


