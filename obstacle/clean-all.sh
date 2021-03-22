#/bin/bash

# This simple script will remove all additional files, 
# maintaining intact only the original project files.
#
# FILENAME:      "clean-all.sh"
# TO RUN: $      ./clean-all.sh

echo "=============================================================================="
echo "      CLEAN ALL      "
echo "=============================================================================="
echo " "
echo "    ARE YOU SURE YOU WANT TO DELETE EVERITHING? "
echo " "
echo "    The only files and folders that will be maintained are: "
echo "      - DIR:  source "
echo "      - DIR:  include "
echo "      - DIR:  doc "
echo "      - FILE: CMakeLists.txt "
echo "      - FILE: clean-all.sh "
echo "      - FILE: COMPILE.sh "
echo " "
echo ">>> Press ENTER to run "
read -n1		# -n1 accept only 1 character
echo " "

find -maxdepth 1 ! -name "source" ! -name "include" ! -name "doc" ! -name "CMakeLists.txt" ! -name "clean-all.sh" ! -name "COMPILE.sh" ! -name . -exec rm -rv {} \;

echo " "
echo "    CODE ENDED CORRECTLY ! "
echo " "
echo "=============================================================================="
echo " "


