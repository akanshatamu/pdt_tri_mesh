#Set up the environment to install OpenMPI
#in command prompt
module load gcc-4.8.2
module load openmpi-1.6.3-gcc-4.7.2

#Clean the object files
make distclean

#Compile the files
make
#Make sure that an object file named Tricall_Pre_Mesh_Gen is created

#Run the code
mpirun -np 8 Tricall_Pre_Mesh_Gen
#-np defines the number of processing threads

#Code executes the main() in Tricall_Pre_Mesh_Gen()
#Tricall_Pre_Mesh_Gen contains code related to parallelization and the code built by Tarek on cut points
#Tricall_Utility contains methods to call the Triangle library
