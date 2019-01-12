function [FacA,FacC,FacN,EleC,Pts_in] = get_surfacenormal_mesh(Pts,Ele,Fac)


x1 = Pts(1,Fac(1,:)); y1 = Pts(2,Fac(1,:)); z1 = Pts(3,Fac(1,:));
x2 = Pts(1,Fac(2,:)); y2 = Pts(2,Fac(2,:)); z2 = Pts(3,Fac(2,:));
x3 = Pts(1,Fac(3,:)); y3 = Pts(2,Fac(3,:)); z3 = Pts(3,Fac(3,:));

NFac = size(Fac,2);
FacA = zeros(1,NFac);
FacC = zeros(3,NFac);
FacN = zeros(3,NFac);

FacC(1,:) = mean([x1;x2;x3],1);
FacC(2,:) = mean([y1;y2;y3],1);
FacC(3,:) = mean([z1;z2;z3],1);

v1 = [x1-x2;y1-y2;z1-z2];
v2 = [x2-x3;y2-y3;z2-z3];

vnor = -cross(v1,v2,1);
FacA = 1/2*sqrt(sum(vnor.^2,1));
vnor = vnor./sqrt(sum(vnor.^2,1));
FacN = vnor;

NEle = size(Ele,2);    
EleC = zeros(3,NEle);
x1 = Pts(1,Ele(1,:)); y1 = Pts(2,Ele(1,:)); z1 = Pts(3,Ele(1,:));
x2 = Pts(1,Ele(2,:)); y2 = Pts(2,Ele(2,:)); z2 = Pts(3,Ele(2,:));
x3 = Pts(1,Ele(3,:)); y3 = Pts(2,Ele(3,:)); z3 = Pts(3,Ele(3,:));
x4 = Pts(1,Ele(4,:)); y4 = Pts(2,Ele(4,:)); z4 = Pts(3,Ele(4,:));
EleC(1,:) = mean([x1;x2;x3;x4],1);
EleC(2,:) = mean([y1;y2;y3;y4],1);
EleC(3,:) = mean([z1;z2;z3;z4],1);

Pts_in = zeros(3,NFac);

for iFac = 1:NFac
    dist = sqrt(((FacC(1,iFac)-EleC(1,:)).^2+(FacC(2,iFac)-EleC(2,:)).^2+(FacC(3,iFac)-EleC(3,:)).^2));
    [mm,ii] = min(dist);
    node_in = EleC(:,ii);
    Pts_in(1:3,iFac) = node_in(1:3);
    vector_in = node_in - FacC(:,iFac);
    vector_in = vector_in/norm(vector_in);
    veccos = dot(vector_in,FacN(:,iFac))/norm(vector_in,2)/norm(FacN(:,iFac),2);
    if (veccos > 1e-6)
        FacN(:,iFac) = -FacN(:,iFac);
    elseif (abs(veccos) <= 1e-6)
        disp(['normal parallel to inside node, stop']);
        stop
    end 
end

%figure; hold on;
%h=quiver3(FacC(1,:),FacC(2,:),FacC(3,:),FacN(1,:),FacN(2,:),FacN(3,:),0,'k'); 
%h=quiver3(Pts_in(1,:),Pts_in(2,:),Pts_in(3,:),FacN(1,:),FacN(2,:),FacN(3,:),'b');
%h=plot3(EleC(1,:),EleC(2,:),EleC(3,:),'ro');
%view(3); axis equal;

