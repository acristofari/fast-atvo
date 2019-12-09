% -------------------------------------------------------------------------
%
% This file is part of FAST-ATVO, which is a software for community
% detection in an undirected graph with non-negative weights.
%
% This file provides an example of how to run FAST-ATVO via Matlab.
% For details on the required syntax, type 'help fast_atvo_matlab_syntax'
% in the Matlab prompt.
%
% -------------------------------------------------------------------------
%
% Reference paper:
% A. Cristofari, F. Rinaldi, F. Tudisco (2019). Total variation based
% community detection using a nonlinear optimization approach.
% arXiv preprint arXiv:1907.08048.
%
% -------------------------------------------------------------------------
%
% Authors:
% Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
% Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
% Francesco Tudisco (e-mail: f.tudisco@strath.ac.uk)
%
% Last update:
% December 9th, 2019
%
% Copyright 2019 Andrea Cristofari, Francesco Rinaldi, Francesco Tudisco.
%
% Licensing:
% This file is part of FAST-ATVO.
% FAST-ATVO is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% FAST-ATVO is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with FAST-ATVO. If not, see <http://www.gnu.org/licenses/>.
%
% -------------------------------------------------------------------------

% N.B. In this example, the starting point of FAST-ATVO is the leading
% eigenvector of the modularity matrix (it is the suggested choice in practice)

make(); % build MEX files

t_start = tic;

%-------------------------------
% load graph
%-------------------------------
dataset_name = 'example_graph.mat'; % it contains a non-negative symmetric weight matrix named 'A'
load(dataset_name,'A');
[N,~] = size(A); % get graph dimension
d = A*ones(N,1); % compute vector of (weighted) degrees
vol = sum(d); % compute graph volume

%-------------------------------
% find the leading eigenvector of the modularity matrix to use it as
% starting point of FAST-ATVO
%-------------------------------
t_x0 = toc(t_start);
opts_linear.issym = 1;  
opts_linear.isreal = 1;
Afun = @(x) A*x - d*(d'*x)/vol;
[x0,~] = eigs(Afun,N,1,'la',opts_linear);
t_x0 = toc(t_start) - t_x0;

%-------------------------------
% run FAST-ATVO
%-------------------------------
A_triu = triu(A); % in FAST-ATVO the weight matrix must be passed as an upper triangular matrix
                  % (it will be symmetrized automatically leaving the diagonal unchanged)
t_solver = toc(t_start);
[C,Q,x] = fast_atvo(A_triu,x0);
t_solver = toc(t_start) - t_solver;

%---------------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE FAST-ATVO PARAMETERS ***
%
% In place of the above instruction '[C,Q,x] = fast_atvo(A_triu,x0);',
% do the following:
%
% (1) create a structure having as field names the parameters to be changed
%     and assign them new values, for instance:
%
%       fastatvo_opts.out_it = 5;
%       fastatvo_opts.verb = 1;
%
% (2)  pass this structure to 'fast_atvo' as third input argument, for instance:
%
%       [C,Q,x] = fast_atvo(A_triu,x0,fastatvo_opts);
%---------------------------------------------------------------------------------

fprintf('\n%s\n\n%s\n\n%s%-.4f\n\n%s%i\n%s%.4g\n%s%-.4e%s\n%s%-.4e%s\n%s%-.4e%s\n\n%s\n\n', ...
        '********************** FINAL RESULTS **********************', ...
        'Algorithm: FAST-ATVO', ...
        'community modularity = ', Q, ...
        'number of graph nodes = ', N, ...
        'graph volume = ', vol, ...
        'solving time: ', max(t_solver,0e0), ' seconds', ...
        'time to compute x0: ', max(t_x0,0e0), ' seconds', ...
        'total time: ', max(t_solver+t_x0,0e0), ' seconds', ...
        '***********************************************************');