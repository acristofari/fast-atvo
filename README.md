# Community detection in undirected graphs

*FAST-ATVO* is a solver for community detection problems in undirected graphs
with non-negative weights, using a non-linear optimization approach.

## Reference paper

[A. Cristofari, F. Rinaldi, F. Tudisco (2020). _Total Variation Based
Community Detection Using a Nonlinear Optimization Approach_. SIAM Journal
on Applied Mathematics, 80(3), 1392-1419](https://epubs.siam.org/doi/10.1137/19M1270446).

## Authors

* Andrea Cristofari (e-mail: [andrea.cristofari@unipd.it](mailto:andrea.cristofari@unipd.it))
* Francesco Rinaldi (e-mail: [rinaldi@math.unipd.it](mailto:rinaldi@math.unipd.it))
* Francesco Tudisco (e-mail: [francesco.tudisco@gssi.it](mailto:francesco.tudisco@gssi.it))

## Licensing

This file is part of FAST-ATVO.
FAST-ATVO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
FAST-ATVO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with FAST-ATVO. If not, see <http://www.gnu.org/licenses/>.

Copyright 2019-2022 Andrea Cristofari, Francesco Rinaldi, Francesco
Tudisco.

## How to use FAST-ATVO

1. This directory should contain the following files:

    * `COPYING.txt`,
    * `ExampleGraph.txt`,
    * `ExampleX0.txt`,
    * `fast_atvo.cpp`,
    * `fast_atvo.h`,
    * `graph.h`,
    * `main.cpp`,
    * `README.md`,

    plus a subdirectory named `matlab`, which should contain the following
    files:

    * `example_graph.mat`,
    * `fast_atvo_matlab.cpp`,
    * `main.m`,
    * `make.m`,
    * `syntax.txt`.

2. You can call FAST-ATVO either from the command prompt (see 2a) or from Matlab (see 2b).

   2a. **How to call FAST-ATVO from the command prompt**

      - Prepare a text file with the weight matrix expressed as a square
        upper triangular matrix. Each line of the text file must have the
        following form:

            i_1,j_1 w_1 i_2,j_2 w_2 i_3,j_3 w_3 ...

        where any tern _i\_h,j\_h w\_h_ represents an edge between the nodes _i\_h_
        and _j\_h_ with non-negative weight _w\_h_.

        It must hold _i\_1 <= i\_2 <= i\_3 <= ...._, i.e., the first nodes of the
        terns must be written in a non-decreasing order.

        Note that it must also hold _i\_1 <= j\_1, i\_2 <= j\_2, i\_3 <= j\_3, ...,_
        since the weight matrix must be upper triangular.

        Unspecified weights between two nodes are assumed to be zero (so
        that only positive weights need to be specified).

        For instance, consider the following weight matrix:

            ┏                             ┓
            ┃ 0     0.9   1.5   2     0   ┃
            ┃ 0.9   0     0     0     0   ┃
            ┃ 1.5   0     0     0.8   1.1 ┃
            ┃ 2     0     0.8   0     0.3 ┃
            ┃ 0     0     1.1   0.3   0   ┃
            ┗                             ┛

        A valid text file will be:

            1,2 0.9 1,3 1.5 1,4 2
            3,4 0.8 3,5 1.1
            4,5 0.3

        Equivalently, lines can even be broken or joined together. This means that
        also a text file of the following form will be valid:

            1,2 0.9 1,3 1.5 1,4 2 3,4 0.8
            3,5 1.1 4,5 0.3

      - Prepare a text file with the starting point of the algorithm.
        Each line of the text file must contain scalars separated by blank spaces
        (one value per line is also allowed).

        The starting point must be a vector of length equal to the number of nodes.

        For instance, consider the following starting point:

            ┏      ┓
            ┃  0   ┃
            ┃  1   ┃
            ┃ -0.3 ┃
            ┃  0   ┃
            ┃  0.2 ┃
            ┗      ┛

        A valid text file will be:

            0
            1
            -0.3
            0
            0.2

        Equivalently, lines can even be broken or joined together. This means that
        also a text file of the following form will be valid:

            0 1 -0.3
            0 0.2

      - Compile the files `fast_atvo.cpp` and `main.cpp`, then create the
        executable `fast_atvo`. To run FAST-ATVO, you have to type in the
        command prompt

            fast_atvo GraphFile X0File [options]

        where `GraphFile` is the name of the file with the weight matrix,
        `X0File` is the name of the file with the starting point of the
        algorithm and `[options]` are optional input arguments that allow
        the user to modify some algorithm parameters and to print the
        final results to files.
        In particular, `[options]` must have the following form:

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

      - When the algorithm is terminated, final results can be found in the
        files specified in the options (if any). Moreover, if verbosity was
        activated, a file named `iteration_history.txt` is created, where
        the iteration details of the optimization algorithm are reported.

      - Here is an example.

        Consider the two files `ExampleGraph.txt` and `ExampleX0.txt`
        included in this folder. They contain a weight matrix and a starting point
        of the algorithm, respectively, according to the above described format.
        Create the executable `fast_atvo` by compiling the files `fast_atvo.cpp` and `main.cpp`,
        then type

            fast_atvo ExampleGraph ExampleX0 -c ExampleC.txt

        so that the communities found by FAST-ATVO will be printed to the
        file `ExampleC.txt`.

        Or, if you also wish to print synthetic iteration details of the
        optimization algorithm, you may type

            fast_atvo ExampleGraph ExampleX0 -c ExampleC.txt -v 1 

   2b. **How to call FAST-ATVO from Matlab**

      - Move to the subdirectory `matlab` and run `make.m` to build the MEX file.

      - See the file `usage.txt` to know how to call FAST-ATVO from Matlab, change algorithm parameters and get output values.

      - See the file `main.m` for an example. To run the example, just call `main.m` in Matlab.