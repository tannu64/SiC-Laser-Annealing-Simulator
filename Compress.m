
function [OutputArray]=Compress(InputArray,SkipNumber)
%Compress the array by skip some points
N=size(InputArray,2)-1;
M=size(InputArray,1);
K=floor(N/SkipNumber);
OutputArray=zeros(M,K+1);
for i=1:M
    for j=0:K
        OutputArray(i,j+1)=InputArray(i,j*SkipNumber+1);
    end
end
end
