% Copyright (c) 2017, Jan Valdman
% Copyright (c) 2010, Jan Valdman
% Copyright (c) 2016, Jan Valdman
% Copyright (c) 2015, Charles Robert
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

function entryrows=entryInWhichRows(A)
%function: entryrows=entryInWhichRows(A)
%requires: none
%for every entry of integer matrix A, 
%its rows indices are stored in output matrix,
%zeros entries indicate no more occurence
%example: entryrows=entryInWhichRows([1 2; 1 3; 2 2]) returns
%         entryrows=[1   2   0; 
%                    1   3   3; 
%                    2   0   0]   
%meaning: entry 1 appears in rows 1 and 2
%         entry 2 appears in rows 1 and 3 (twice)
%         entry 3 appears in row  2 only

%size computation; 
r=max(max(A));
repetition=accumarray(A(:),ones(numel(A),1));
c=max(repetition);

%filling rows occurences
%this part should be somehow vectorized!
entryrows=zeros(r,c);
repetition=zeros(r,1);
for i=1:size(A,1)
    for j=1:size(A,2)
       index=A(i,j); 
       repetition(index)=repetition(index)+1; 
       entryrows(index,repetition(index))=i;
    end
end
