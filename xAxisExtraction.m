function [xPointSum,xAxis,Dx,DxTrs]=xAxisExtraction(xMesh)
% extract x axis related parameter
% xPointSum: the amount of x point, unit:um
% xAxis:the array of x position,unit:um
% Dx: the size or width at very x position, unit: um
% DxTrs:the length from one x point to next, unit:cm
M=size(xMesh,1);  % number of rows in xMesh
xMeshRenew=zeros(M,4);  % creat a new matrix with original value of 0
xPositionCount=0;  % creat a counting variable to get the x point position
% fill the xMeshRenew matrix in 3rd column with eacg layer's end position
for i=1:M
    xMeshRenew(i,1)=xMesh(i,1);
    xMeshRenew(i,2)=xMesh(i,2);
    xPositionCount=xPositionCount+xMeshRenew(i,1);
    xMeshRenew(i,3)=xPositionCount;
end
% fill the xMeshRenew matrix in 4th column with point number in each layer
xPointCount=0;    % creat a counting variable to get the x point number
xPointSum=0;     %creat a counting variable to get the sum of x point
xPositionCount=0;
for i= 1:M
    xPiontNumber=round((xMeshRenew(i,3)-xPositionCount)/xMeshRenew(i,2));  
    xPointCount=xPointCount+xPiontNumber;
    xPositionCount=xPositionCount+xMeshRenew(i,2)*xPiontNumber;
    xMeshRenew(i,3)=xPositionCount; %renew it
    xMeshRenew(i,4)=xPiontNumber; 
end
%creat xAxis,Dx,DxTrs arrays
xAxis=zeros(1,xPointCount);
Dx=zeros(1,xPointCount);
DxTrs=zeros(1,xPointCount-1);
xPointSum=xPointCount;
xPointCount=0;
xPositionCount=0;  % x position
%fill xAxis,Dx,DxTrs arrays
for i=1:M
    for j=1:xMeshRenew(i,4)
        xPointCount=xPointCount+1;
        xPositionCount=xPositionCount+xMeshRenew(i,2);
        Dx(xPointCount)=xMeshRenew(i,2);
        xAxis(xPointCount)=xPositionCount-0.5*xMeshRenew(i,2);
    end
end 
%extract DxTrs and convert the unit from um to cm
for i=1:xPointSum-1
    DxTrs(i)=10^(-4)*(Dx(i)+Dx(i+1))/2;
end
end
