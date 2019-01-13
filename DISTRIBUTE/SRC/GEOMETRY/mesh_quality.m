function qmesh=mesh_quality( node_xyz, tetra_node) 
% Copyright (c) 2019, Van-Dang Nguyen

% % The minimum, over all tetrahedrons, of 3 times the radius of the insphere divided by the radius of the circumsphere.

tetra_v=cell(4,1);
tetra_xyz=cell(4,1);

for idv=1:4
    tetra_v{idv}=tetra_node(idv,:);
    tetra_xyz{idv}=node_xyz(:,tetra_v{idv});
end;

v21 = tetra_xyz{2} - tetra_xyz{1};
v31 = tetra_xyz{3} - tetra_xyz{1};
v41 = tetra_xyz{4} - tetra_xyz{1};
v32 = tetra_xyz{3} - tetra_xyz{2};
v42 = tetra_xyz{4} - tetra_xyz{2};
v43 = tetra_xyz{4} - tetra_xyz{3};

n123 = cross( v21, v31 );
n124 = cross( v41, v21 );
n134 = cross( v31, v41 );
n234 = cross( v42, v32 );    

l123 = vecnorm(n123,2);
l124 = vecnorm(n124,2);
l134 = vecnorm(n134,2);
l234 = vecnorm(n234,2);

pc = (   l234.* tetra_xyz{1} + l134.* tetra_xyz{2}   ...
      + l124.* tetra_xyz{3} + l123.* tetra_xyz{4} )./( l234 + l134 + l124 + l123 );

b = tetra_xyz;
for idv=1:4
    b{idv}(4,:)=1;
end;

det = ...
  b{1}(1,:) .* ( ...
    b{2}(2,:) .* ( b{3}(3,:) .* b{4}(4,:) - b{3}(4,:) .* b{4}(3,:) ) ...
  - b{2}(3,:) .* ( b{3}(2,:) .* b{4}(4,:) - b{3}(4,:) .* b{4}(2,:) ) ...
  + b{2}(4,:) .* ( b{3}(2,:) .* b{4}(3,:) - b{3}(3,:) .* b{4}(2,:) ) ) ...
- b{1}(2,:) .* ( ...
    b{2}(1,:) .* ( b{3}(3,:) .* b{4}(4,:) - b{3}(4,:) .* b{4}(3,:) ) ...
  - b{2}(3,:) .* ( b{3}(1,:) .* b{4}(4,:) - b{3}(4,:) .* b{4}(1,:) ) ...
  + b{2}(4,:) .* ( b{3}(1,:) .* b{4}(3,:) - b{3}(3,:) .* b{4}(1,:) ) ) ...
+ b{1}(3,:) .* ( ...
    b{2}(1,:) .* ( b{3}(2,:) .* b{4}(4,:) - b{3}(4,:) .* b{4}(2,:) ) ...
  - b{2}(2,:) .* ( b{3}(1,:) .* b{4}(4,:) - b{3}(4,:) .* b{4}(1,:) ) ...
  + b{2}(4,:) .* ( b{3}(1,:) .* b{4}(2,:) - b{3}(2,:) .* b{4}(1,:) ) ) ...
- b{1}(4,:) .* ( ...
    b{2}(1,:) .* ( b{3}(2,:) .* b{4}(3,:) - b{3}(3,:) .* b{4}(2,:) ) ...
  - b{2}(2,:) .* ( b{3}(1,:) .* b{4}(3,:) - b{3}(3,:) .* b{4}(1,:) ) ...
  + b{2}(3,:) .* ( b{3}(1,:) .* b{4}(2,:) - b{3}(2,:) .* b{4}(1,:) ) );
gamma = abs ( det );

% InRadius
hin = gamma ./ ( l234 + l134 + l124 + l123 );

% volume of the tetrahedron
vol=1/6*abs(dot(tetra_xyz{1}-tetra_xyz{4},cross(tetra_xyz{2}-tetra_xyz{4},tetra_xyz{3}-tetra_xyz{4})));

l21 = vecnorm(v21,2);
l31 = vecnorm(v31,2);
l41 = vecnorm(v41,2);
l32 = vecnorm(v32,2);
l42 = vecnorm(v42,2);
l43 = vecnorm(v43,2);    

r1 = l21 .*l43  + l31  .*l42 + l41.*l32;
r2 = -l21.*l43 + l31  .*l42 + l41.*l32;
r3 = l21 .*l43 - l31  .*l42 + l41.*l32;
r4 = l21 .*l43  + l31 .*l42 - l41.*l32;

% Circumradius
hout = sqrt(r1.*r2.*r3.*r4)./(24*vol);

tetrahedron_quality = 3.0 * hin ./ hout;
value_max  = max  ( tetrahedron_quality );
value_min  = min  ( tetrahedron_quality );
value_mean = mean ( tetrahedron_quality );
value_var  = var  ( tetrahedron_quality );
    
qmesh.quality=[value_min, value_mean, value_max, value_var];
qmesh.hout = hout;
qmesh.hin = hin;
qmesh.vol = vol;

