##
#  CMake script for the step-41 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "main-step-41")
include_directories(../include)

SET(TARGET_SRC ${TARGET}.cc)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# NOTE: Try to remove "$ENV{mkDealiiPrefix}" from next line if you DO NOT use mkModules and have problems. 
FIND_PACKAGE(deal.II 8.3 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR} $ENV{mkDealiiPrefix}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
