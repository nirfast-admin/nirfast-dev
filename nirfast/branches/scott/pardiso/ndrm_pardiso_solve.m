function x=ndrm_pardiso_solve(A,b)
%
% Copyright 2007 Dartmouth College
% Written by Andrea Borsic
% Hanover, 18/07/07
%
% This function calls two me files (real / complex) in order to solve a linear system A*x=b with PARDISO direct parallel solver
%
% Use: x=ndrm_pardiso_solve(A,b);
%
% INPUTS:
% A = sparse system matrix (see below)
% b = RHS (can be a multicolumn vector)
% OUTPUTS:
% x = solution (will be a multicolumn vector is the RHS is)

if (nargin~=2)
    error('the number of input arguments is incorrect');
end % if

warning off;
num_threads=maxNumCompThreads;
warning on;

% we call PARDISO with the same number of computational threads we use in MATLAB (as we would get a crash differently and as this is what we want)

v=ver('matlab');
v=str2num(v.Version); % get the version number as a number

switch v
    case {7.9} % from R2009b the threading library is libomp
            x=mex_pardiso_libomp(A,b,num_threads);
    case {7.5,7.6,7.7,7.8,7.10} % from R2007b to R2009a the threading library is libguide
            x=mex_pardiso_libguide(A,b,num_threads);
    otherwise
        error('Your matlab version should be newer than R2007b for NDRM to work');
end % switch

% Bye !