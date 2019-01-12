% Copyright (c) 2019, Jan Valdman
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:

% * Redistributions of source code must retain the above copyright notice, this
% list of conditions and the following disclaimer.

% * Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution

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


function [elems2faces, faces2nodes, faces2elems]=get_faces(elems2nodes)
%function: [element2faces, face2nodes]=get_faces(elems2nodes)
%generates faces of a tetrahedral or hexahedral mesh
%requires: deleterepeatedrows
%input:  elems2nodes - matrix, the i-th row contains numbers of nodes corresponding to the i-th elements
%output: elems2faces - matrix, the i-th row contains numbers of faces corresponding to the i-th elements
%output: faces2nodes - matrix, the i-th row contains numbers of nodes corresponding to the i-th face
%example:  [element2faces, face2nodes]=get_faces([1 3 4 5; 7 4 3 5; 5 7 6 4; 6 8 7 4; 2 1 4 5; 2 4 6 5]) 

switch size(elems2nodes,2)
    case {4,10}      %tetrahedral mesh using P1 elements
        is_quad       = 0;               % flag if the surface is rectangle or triangle
        n_surf_nodes  = 3;               % number of main vertex on the surface
        n_surf        = 4;               % number of faces for each element
        %n_main_nodes  = 4;               % number of main vertices on the element
    case {8,20}
        is_quad = 1; 
        n_surf = 6;
        n_surf_nodes  = 4;
        %n_main_nodes = 8;
end

%extracts sets of faces
if ~is_quad
        faces1=elems2nodes(:,[1 2 3]);
        faces2=elems2nodes(:,[1 2 4]);
        faces3=elems2nodes(:,[1 3 4]);
        faces4=elems2nodes(:,[2 3 4]);
else
      faces1=elems2nodes(:,[1 2 3 4]);
      faces2=elems2nodes(:,[1 2 6 5]);
      faces3=elems2nodes(:,[1 4 8 5]);
      faces4=elems2nodes(:,[2 3 7 6]);  
      faces5=elems2nodes(:,[3 4 8 7]);
      faces6=elems2nodes(:,[5 6 7 8]);
end

%as sets of their nodes (vertices)
vertices=zeros(size(elems2nodes,1)*n_surf,n_surf_nodes);
vertices(1:n_surf:end,:)=faces1;
vertices(2:n_surf:end,:)=faces2;
vertices(3:n_surf:end,:)=faces3;
vertices(4:n_surf:end,:)=faces4;
if is_quad
  vertices(5:n_surf:end,:)=faces5;
  vertices(6:n_surf:end,:)=faces6;
end

%repeated sets of nodes (joint faces) are eliminated 
[faces2nodes,elems2faces] = deleteRepeatedRows(vertices);
elems2faces = reshape(elems2faces,n_surf,size(elems2nodes,1))';




