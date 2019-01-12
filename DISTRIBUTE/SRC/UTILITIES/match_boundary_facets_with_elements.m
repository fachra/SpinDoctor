function boundaryTetgen2element=match_boundary_facets_with_elements(boundaryTetgen,elements)
% Copyright (c) 2019, Jan Valdman

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

