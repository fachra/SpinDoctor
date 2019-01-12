% Copyright (c) 2019, Jan Valdman

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


function boundaryTetgen2element=match_boundary_facets_with_elements(boundaryTetgen,elements)

%this part creates own boundary facets
[elems2faces, faces2nodes]=get_faces(elements);
faces2elems=entryInWhichRows(elems2faces); 
bindex=find(faces2elems(:,2)==0);           %index of boundary faces - belong to one element only
boundary=faces2nodes(bindex,:);             %nodes of boundary faces 
boundary2element=faces2elems(bindex,1);     %element of boundary faces (unique)

%oder all matrices columnwise for comparison
boundary_ordered=sort(boundary,2);
boundaryTetgen_ordered=sort(boundaryTetgen,2);

%find their intersection
[~,ind,indTetgen] = intersect(boundary_ordered,boundaryTetgen_ordered,'rows');

%and assign (unique) elements to boundary facets 
boundaryTetgen2element=zeros(size(boundaryTetgen,1),1);
boundaryTetgen2element(indTetgen)=boundary2element(ind);  %a columnt vector

end

