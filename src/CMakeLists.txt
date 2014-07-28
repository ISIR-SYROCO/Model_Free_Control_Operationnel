include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PC_EIGEN_INCLUDE_DIRS}
    ${PC_ORCFRAMEWORK_INCLUDE_DIRS}
	/home/kuka/compil/qpOASES-3.0beta/include
)

link_directories(
	${PC_ORCFRAMEWORK_LIBRARY_DIRS}
	/home/kuka/compil/qpOASES-3.0beta/bin
)

#
# Components, types and plugins.
#
# The CMake 'target' names are identical to the first argument of the
# macros below, except for orocos_typegen_headers, where the target is fully
# controlled by generated code of 'typegen'.
#

# Creates a component library 
# and installs in the directory lib/orocos/
#	
#orocos_library()	
orocos_library(ModelKukaFixed kukafixed.cpp)	
orocos_component(KukaRMControllerRTNET KukaRMC-rtnetcomponent.cpp) # ...you may add multiple source files.
target_link_libraries(ModelKukaFixed ${PC_ORCFRAMEWORK_LIBRARIES})
target_link_libraries(KukaRMControllerRTNET ModelKukaFixed)
target_link_libraries(KukaRMControllerRTNET qpOASES)
#target_link_libraries(KukaSimpleDemo )
#
# You may add multiple orocos_component statements.

#
# Additional headers:
#
# Installs in the include/orocos/kuka_send_joint_positions/ directory
#
orocos_install_headers(KukaRMC-rtnetcomponent.hpp) # ...you may add multiple header files
#
# You may add multiple orocos_install_headers statements.

