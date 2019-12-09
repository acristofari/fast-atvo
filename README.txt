-------------------------------------------------------------------------

FAST-ATVO is a software for community detection in an undirected graph
with non-negative weights.

-------------------------------------------------------------------------

Reference paper:
A. Cristofari, F. Rinaldi, F. Tudisco (2019). Total variation based
community detection using a nonlinear optimization approach.
arXiv preprint arXiv:1907.08048.

-------------------------------------------------------------------------

Authors:
Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
Francesco Tudisco (e-mail: f.tudisco@strath.ac.uk)

-------------------------------------------------------------------------


How to use FAST-ATVO
=========================================================================

1 - In this folder, you should see the following files:

    - 'COPYING.txt',
    - 'fast_atvo.cpp',
    - 'fast_atvo.h',
    - 'graph.h',
    - 'main.cpp',
    - 'README.txt',

    plus a subfolder named 'matlab', which includes the following files:

    - 'example_fast_atvo.m',
    - 'example_graph.mat',
    - 'fast_atvo_matlab.cpp',
    - 'fast_atvo_matlab_syntax.txt',
    - 'make.m'.

2 - You can run FAST-ATVO either via command prompt or via Matlab.
    
    To run FAST-ATVO via command prompt do the following steps 2a-2d and
    see the subsequent example, otherwise skip to step 3 to see how to
    run FAST-ATVO via Matlab.

  2a - Prepare a text file with the weight matrix, each line being of the
       following form:

       i_1,j_1 w_1 i_2,j_2 w_2 i_3,j_3 w_3 ... 

       where any tern i_h,j_h w_h represents an edge between nodes i_h
       and j_h with non-negative weight w_h.

       Note that:
         
       - Unspecified weights between two nodes are assumed to be zero (so
       that only positive weights need to be specified).

       - The weight matrix must be upper triangular (it will be
       symmetrized automatically leaving the diagonal unchanged, since the
       graph is assumed to be undirected) and the first nodes of the terns
       must be ordered increasingly. This means that it must hold
       i_1 <= j_1, i_2 <= j_2, i_3 <= j_3, ..., and
       i_1 <= i_2 <= i_3 <= ....

       For instance, consider the following weight matrix:

       [0    0.9   1.5   2     0
       0.9   0     0     0     0
       1.5   0     0     0.8   1.1
       2     0     0.8   0     0.3
       0     0     1.1   0.3   0  ].

       A valid text file will be:

       1,2 0.9 1,3 1.5 1,4 2
       3,4 0.8 3,5 1.1
       4,5 0.3

       Equivalently, lines can even be broken or unified. This means that
       also a text file of the following form will be valid:

       1,2 0.9 1,3 1.5 1,4 2 3,4 0.8
       3,5 1.1 4,5 0.3

  2b - Prepare a text file with the starting point of the algorithm, each
       line containing scalars separated by blank spaces (one value per
       line is also allowed).

       Note that the starting point is a column vector and its length
       must be equal to the number of nodes.

       For instance, consider the following starting point:

       [0
        1
        -0.3
        0
        0.2].

       A valid text file will be:

       0
       1
       -0.3
       0
       0.2

       Equivalently, lines can even be broken or unified. This means that
       also a text file of the following form will be valid:

       0 1 -0.3
       0 0.2

  2c - Compile 'fast_atvo.cpp' and 'main.cpp', then create the executable
       'fast_atvo'. To run FAST-ATVO, you have to type in the command
       prompt:

       fast_atvo GraphFile X0File [options]

       where 'GraphFile' is the name of the file with the weight matrix,
       'X0File' is the name of the file with the starting point of the
       algorithm and '[options]' are non-mandatory arguments that allow
       the user to modify some algorithm parameters and to print the
       final results to files.
       In particular, '[options]' must have the following form:

       -c string
          It is the name of the file where the communities found by
          FAST-ATVO are printed as a 0-1 vector (if the file does not
          exist, then it will be created, whereas existing files with
          the same name will be overwritten).
          If not specified, by default no file is created.
       -m string
          It is the name of the file where the modularity value of the
          communities found by FAST-ATVO is printed (if the file does not
          exist, then it will be created, whereas existing files with the
          same name will be overwritten).
          If not specified, by default no file is created.
       -s string
          It is the name of the file where the solution found by the
          optimization algorithm is printed (if the file does not exist,
          then it will be created, whereas existing files with the same
          name will be overwritten).
          If not specified, by default no file is created.
       -p number greater than 1
          It is the exponent parameter of the objective function.
          If not specified, by default it is equal to 1.4.
       -w number greater than or equal to 1
          It is the maximum size of the working set in the optimization
          algorithm.
          If not specified, by default it is equal to
          max(10,min(1000,0.03*n)), where 'n' is the number of
          non-isolated nodes.
       -i number greater than or equal to 1
          It is the number of outer iterations for the globalization
          strategy.
          If not specified, by default it is equal to 1, i.e., the
          globalization strategy is not activated.
       -l number less than 0
          It is the lower bound on the variables for the optimization
          problem.
          If not specified, by default it is equal to -1.
       -u number greater than 0
          It is the upper bound on the variables for the optimization
          problem.
          If not specified, by default it is equal to 1.
       -r number between 0 and 1
          It is the percentage of negative and positive variables that
          will be set to the lower and upper bound in x0, respectively.
          If not specified, by default it is equal to 1, i.e., all
          non-zero variables in x0 will be set to the bounds.
       -v number between 0 and 2
          It is the verbosity level, to print iteration details of the
          optimization algorithm.
          If not specified, by default it is equal to 0, i.e., there are
          no prints.

  2d - When the algorithm is terminated, final results can be found in
       the files specified in the options (if any). Moreover, if
       verbosity was activated, a file named 'iteration_history.txt' is
       created, where the iteration details of the optimization algorithm
       are repored.

  Example.
       Consider the two example files included in the main folder, i.e.,
       'ExampleGraph.txt' and 'ExampleX0.txt'. They contain a graph and a
       starting point of the algorithm, respectively, according to the
       above described format. To run FAST-ATVO via command prompt, after
       compiling the .cpp files and creating the executable 'fast_atvo',
       you may type

       fast_atvo ExampleGraph ExampleX0 -c ExampleC.txt

       so that the communities found by FAST-ATVO will be printed to the
       file 'ExampleC.txt'.

       Or, if you also wish to print synthetic iteration details of the
       optimization algorithm, you may type

       fast_atvo ExampleGraph ExampleX0 -c ExampleC.txt -v 1

3 - To run FAST-ATVO via Matlab, move to the subfolder 'matlab' and run
    'make.m' to build MEX files.

    See file 'fast_atvo_matlab_syntax.txt' for the required syntax and
    then run 'fast_atvo'.

    See file 'example_fast_atvo.m' for an example.