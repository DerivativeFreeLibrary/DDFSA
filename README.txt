-----------------------------------------------------------
 How to use the derivative-free simulated annealing global
 optimizer DDFSA for bound constrained global optimization 
 problems
-----------------------------------------------------------
 The package provides a FORTRAN90 version of the code.

0- Gunzip and untar the archive in a folder on your computer by
   issuing in a directory of your choice (ex. curdir) the command

   $> tar -xvf DDFSA.tar.gz

2- Edit file curdir/problem.f90 to define your own objective function.
   In particular, modify the subroutines 
   ninit     : which sets problem dimension
   functinit : which sets upper and lower bounds on the variables plus
               a name for the problem and a know global value, if any
   funct     : which defines the objective function

2- At command prompt in curdir execute 

     $> make
 
   which will create the executables 'ddfsa'

4- execute

     $> ./ddfsa

