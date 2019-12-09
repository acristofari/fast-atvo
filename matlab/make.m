% -------------------------------------------------------------------------
%
% This file is part of FAST-ATVO, which is a software for community
% detection in an undirected graph with non-negative weights.
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

function make()
    mex -largeArrayDims ../fast_atvo.cpp fast_atvo_matlab.cpp
end