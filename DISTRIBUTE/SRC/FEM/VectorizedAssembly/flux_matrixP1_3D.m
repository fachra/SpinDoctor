% Copyright (c) 2018, Jan Valdman
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of University of South Bohemia &  Institute of Information Theory and Automation, Czech Republic nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [M_boundary, neumann2area]=flux_matrixP1_3D(neumann,coordinates,coeffs)
%coeffs can be only P0 (elementwise constant) function 
%represented by a collumn vector with size(elements,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally
%Note: P1 coeffs needs a higher integration rule (not implemented yet)

% [I,J,V] = find(X,...) also returns a vector V containing the values
%  that correspond to the row and column indices I and J.

%computes areas on Neumann faces
neumann2area=evaluate_area(neumann,coordinates);
% nargin: Number of function input arguments.
% this will create the Q but just for the boundary nodes
if nargin==2 
    M=mass_matrixP1_2D(neumann,neumann2area);
elseif nargin==3
    M=mass_matrixP1_2D(neumann,neumann2area,coeffs);
end
% this step will create the Q for all nodes
[X,Y,Z]=find(M);
M_boundary=sparse(X,Y,Z,size(coordinates,1),size(coordinates,1));

% M1=M;M2=M_boundary;ind=find(M1>10^-13);ind1=find(M2>10^-13);