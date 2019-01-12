% Copyright (c) 2016, Jan Valdman
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


function [A,volumes]=stiffness_matrixP1_3D(elements,coordinates,coeffs)
%coeffs can be either P0 (elementwise constant) or P1 (elementwise nodal) function 
%represented by a collumn vector with size(elements,1) or size(coordinates,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally)

NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!  
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB           
        coord(d,i,:)=coordinates(elements(:,i),d);       
    end   
end
IP=[1/4 1/4 1/4]';    
[dphi,jac] = phider(coord,IP,'P1'); %integration rule, it must be known!  

volumes=abs(squeeze(jac))/factorial(DIM); %%det(J)/3!

dphi = squeeze(dphi); 

if (nargin<3)
    Z=astam(volumes',amtam(dphi,dphi));  
else
    if numel(coeffs)==size(coordinates,1)  %P1->P0 averaging
        coeffs=evaluate_average_point(elements,coeffs);
    end  
    Z=astam((volumes.*coeffs)',amtam(dphi,dphi));
end

Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);
%copy this part for a creation of a new element

X=permute(Y,[2 1 3]); 
A=sparse(X(:),Y(:),Z(:)); 