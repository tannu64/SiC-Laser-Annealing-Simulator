function [tPointSum,tAxis,tDelta]=tAxisExtraction(tMesh)
%time array extraction from tMesh matrix
% tPointSum£ºthe amount of time point
% tAxis: the array of time table, unit is us
% tDelta: the step of time mesh , unit is second
tPointSum=round(tMesh(1)/tMesh(2))+1; %calculate the amount to time point.
tAxis=zeros(1,tPointSum);    %define a array to save all the moment like a time table
tSum=0;
for o = 1:tPointSum
    tAxis(o)=tSum;
    tSum=tSum+tMesh(2);
end
tDelta=(10^(-6))*tMesh(2);  %convert the unit of us to s 
end
function [Capacity,Density,Conductivity,LasPowRatProf,Absorbance,Reflection]=MatLasProtertyExtraction(MatPro,LasPro,RefAbs,xMesh)
%extract material and laser property in arrays or matrixs along depth(x)
nRowInLasPro=size(LasPro,1);%extract the amount of laser beams applied
nRowInRefAbs=size(RefAbs,1);%extract the amount of row in RefAbs 
nColInRefAbs=size(RefAbs,2);%extract the amount of colume in RefAbs
nRowInMatPro=size(MatPro,1);%extract the amount of material applied
%find if there is any unmatched settings
if nRowInRefAbs~=nRowInLasPro   
    error("RefAbs's row number does not match that of LasPro,please correct it");
end
if nColInRefAbs~=(nRowInMatPro+1)
    error("RefAbs's colume number does not match that of MatPro,please correct it");
end
%cite the function block
[xPointSum,xPosition,Dx,~]=xAxisExtraction(xMesh);
W=size(MatPro,1);%the amount of material applied
MatRange=0;%for calculate the sum of all material thickness
%extract tatal material thickness
for k=1:W
    MatRange=MatRange+MatPro(k,1);
end
%extract reflection array
Reflection=zeros(1,nRowInLasPro);
for k=1:nRowInLasPro
    Reflection(k)=RefAbs(k,1);
end
%to tell if th x_mesh is defined exceeding the material boundary
if MatRange<xPosition(xPointSum)
    error("x_mesh exceeds the material boundary!!!" );
end
%creat capacity, density,conduct arrays with original values of zero
Capacity=zeros(1,xPointSum);
Density=zeros(1,xPointSum);
Conduct=zeros(1,xPointSum);
%creat two counting variables
MatRange=0;
m=1;
%fill the capacity, density,conduct arrays with actual value
for k=1:xPointSum
    if (MatRange+MatPro(m,1))<xPosition(k)
        MatRange=MatRange+MatPro(m,1);
        m=m+1;
    end
    Capacity(k)=MatPro(m,2);
    Density(k)=MatPro(m,3);
    Conduct(k)=MatPro(m,4); 
end
%creat conductivity array with original value of zero
Conductivity=zeros(1,xPointSum-1);
%fill the conductivity arrays with actual value
for i=1:xPointSum-1
    Conductivity(i)=(Conduct(i+1)*Dx(i+1)+Conduct(i)*Dx(i))/(Dx(i+1)+Dx(i));
end
%creat the power intensity martix in percentage in material along depth with
%original value of 1
LasPowRatProf=ones(nRowInLasPro,xPointSum);
%creat the Absorbance martix with original value of 0
Absorbance=zeros(nRowInLasPro,xPointSum);%extract absorbance in material along depth
%fill the power intensity martix and absorbance arrays with actual value
for p=1:nRowInLasPro
    MatRange=0;
    m=1;
    for k=1:xPointSum
        if (MatRange+MatPro(m,1))<xPosition(k)
            MatRange=MatRange+MatPro(m,1);
            m=m+1;
        end
        Absorbance(p,k)=RefAbs(p,m+1);
        if k~=1
            LasPowRatProf(p,k)=LasPowRatProf(p,k-1)*exp(-Absorbance(p,k)*(xPosition(k)-xPosition(k-1)));
        end
    end
end
%convert the unit um-1 of Absorbance to cm-1
for i=1:nRowInLasPro
    for j=1:xPointSum
        Absorbance(i,j)=10^4*Absorbance(i,j);
    end
end
end