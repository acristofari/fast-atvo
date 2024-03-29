%  -------------------------------------------------------------------------
% 
%  This file is part of FAST-ATVO, which is a solver for community
%  detection problems in undirected graphs with non-negative weights.
% 
%  -------------------------------------------------------------------------
% 
%  Reference paper:
% 
%  A. Cristofari, F. Rinaldi, F. Tudisco (2020). Total Variation Based
%  Community Detection Using a Nonlinear Optimization Approach. SIAM Journal
%  on Applied Mathematics, 80(3), 1392-1419
% 
%  -------------------------------------------------------------------------
% 
%  Authors:
%  Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
%  Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
%  Francesco Tudisco (e-mail: francesco.tudisco@gssi.it)
% 
%  Last update of this file:
%  April 8th, 2022
% 
%  Licensing:
%  This file is part of FAST-ATVO.
%  FAST-ATVO is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  FAST-ATVO is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.
%  You should have received a copy of the GNU General Public License
%  along with FAST-ATVO. If not, see <http://www.gnu.org/licenses/>.
% 
%  Copyright 2019-2022 Andrea Cristofari, Francesco Rinaldi, Francesco
%  Tudisco.
% 
%  -------------------------------------------------------------------------

% In this file, it is shown how to call FAST-ATVO to solve a user-defined problem.

% Here, the starting point of FAST-ATVO is the leading eigenvector of the
% modularity matrix (it is the suggested choice in practice)

make; % build the MEX file (just the first time)

%-------------------------------
% load graph
%-------------------------------
dataset_name = 'example_graph.mat'; % it contains a non-negative symmetric weight matrix named 'A'
load(dataset_name,'A');
[n,~] = size(A); % get graph dimension
d = A*ones(n,1); % compute vector of (weighted) degrees
vol = sum(d); % compute graph volume

%-------------------------------
% find the leading eigenvector of the modularity matrix to use it as
% starting point of FAST-ATVO
%-------------------------------
opts_linear.issym = 1;  
opts_linear.isreal = 1;
Afun = @(x) A*x - d*(d'*x)/vol;
[x0,~] = eigs(Afun,n,1,'la',opts_linear);

%-------------------------------
% run FAST-ATVO
%-------------------------------
A_triu = triu(A); % in FAST-ATVO the weight matrix must be passed as an upper triangular matrix
                  % (it will be symmetrized automatically leaving the diagonal unchanged)
[C,Q,x] = fast_atvo(A_triu,x0);

%--------------------------------------------------------------------------
% *** EXAMPLE OF HOW TO CHANGE FAST-ATVO PARAMETERS ***
%
% Instead of calling the algorithm by the above instruction
% '[C,Q,x] = fast_atvo(A_triu,x0);', do the following:
%
% (1) create a structure having as field names the names of the parameters
%     to be changed and assign them new values, e.g.,:
%
%       opts.verbosity = 1;
%
% (2) pass the structure to 'fast_atvo' as third input argument, e.g.,:
%
%       [C,Q,x] = fast_atvo(A_triu,x0,opts);
%--------------------------------------------------------------------------

% write statistics to the screen
fprintf('\n%s\n\n%s\n\n%s%-.4f\n\n%s%i\n%s%.4g\n\n%s\n\n', ...
        '*************************************************************', ...
        'Algorithm: FAST-ATVO', ...
        'community modularity = ', Q, ...
        'number of graph nodes = ', n, ...
        'graph volume = ', vol, ...
        '*************************************************************');