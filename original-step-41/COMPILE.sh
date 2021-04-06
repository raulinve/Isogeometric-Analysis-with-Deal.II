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
#--------------------------------------------- [COMMENT/DE-COMMENT SECTION]
#cmake .. -DDEAL_II_DIR=/usr/include/deal.II/                                 # apt-intalled location [COMMENT/DE-COMMENT]
cmake .. -DDEAL_II_DIR=/u/sw/pkgs/toolchains/gcc-glibc/9/pkgs/dealii/9.2.0/   # PoliMi mkModules [COMMENT/DE-COMMENT]
#cmake ..                                                                     # Standard [DO NOT USE]
#cmake .. -DDEAL_II_DIR=/opt/dealii/8.3.0                                     # Other versions [DO NOT USE]
#---------------------------------------------

printf "\n\n--------------------------------\n  make\n--------------------------------\n\n"
make
#---------------------------------------------
printf "\n\n--------------------------------\n\n"
echo "    COMPILATION ENDED [CHECK ABOVE IF ANY ERROR OCCURRED!] "
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


