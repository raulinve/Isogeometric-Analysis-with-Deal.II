#/bin/bash

# FILENAME:      COMPILE.sh
# TO RUN DIGIT:  $ ./COMPILE.sh

# NOTE: Probably is necessary for you to change the path to the deal.ii folder
#       to do this, find below "DDEAL_II_DIR" and change "/usr/include/deal.II/".

echo "=============================================================================="
echo "      COMPILE: main-step-41      "
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
./main-step-41
#---------------------------------------------

echo " "
echo "=============================================================================="
echo " "


