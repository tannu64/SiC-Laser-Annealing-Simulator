function main
%define material property 
%[thickness-um, Capacity-J/gC, Density-g/cm3, thermal conductivity-W/cmK]
%one support multiple layers
MatPro = [0.1, 0.97, 2.05, 2.05     %C-film
          0.05, 1.3, 3.21, 0.011     %amf-SiC
          100.0, 1.3, 3.21, 2.8];    % crs-SiC
%define laser property
%[wavelength-nm,pulsewidth-us,energydensity-J/cm2,delay-us,mode(1:Gaussian;2:trapezoid)]
%support multiple lasers to work together
LasPro = [527,  0.25,  1.00,  0.00,  1
          527,  0.25,  1.00,  1.00,  1];
%define reflectivity and absorbance
%[reflectivity,absorbance in first layer,absorbance in second layer,.....]
%the number of absorbance must be equal to the number of layers in MatPro
%the number of row must be equal to the number of layers in LasPro
RefAbs =[0.182, 12.3, 1, 0.0001
         0.182, 12.3, 1, 0.0001];
%define x axis mesh [thickness-um, mesh step-um],support multiple mesh scheme     
xMesh = [1.0,   0.01               
         99.0,   0.5]; 
%define t axis mesh 
%[range-us, mesh step-us],
%support only 1 row
tMesh=[10,0.0000001];
%define the temperature profile at various depths we are going to extract
%unit:um; the array must be ordered from small to big.
ProfilAtDepth=[0,MatPro(1,1),MatPro(1,1)+MatPro(2,1)];
%cite the function
[tPointNumber,tAxis,tDelta]=tAxisExtraction(tMesh);
[xPointSum,xAxis,Dx,DxTrs]=xAxisExtraction(xMesh);
[LasNunber,LasPowProf]=LasPowFrofAtTimeExtraction(LasPro,tMesh);
[Capacity,Density,Conductivity,LasPowRatProf,Absorbance,Reflection]=MatLasProtertyExtraction(MatPro,LasPro,RefAbs,xMesh);
[TMax,Tatx]=Calculation(MatPro,LasPro,RefAbs,xMesh,tMesh,ProfilAtDepth);
%compress some arrays if it is too large. 
%NOTICE: etch excel sheet can only store less than 100 000 lines' data.
[tAxis]=Compress(tAxis,10000);
[Tatx]=Compress(Tatx,10000);
[LasPowProf]=Compress(LasPowProf,10000);
%plot the profile
figure(1)
plot(tAxis,LasPowProf);
xlabel('t(us)');
ylabel('PowerDensity(W/cm2)');
figure(2)
plot(xAxis,TMax);
xlabel('x(um)');
ylabel('TMax(C)');
figure(3)
plot(tAxis,Tatx);
xlabel('t(us)');
ylabel('Tatx(C)');
figure(4)
plot(xAxis,LasPowRatProf);
xlabel('x(um)');
ylabel('LasPowRatProf');
%name the location you want to save the data
%NOTICE:change the route to which in your computer
filename = sprintf('C:\\Users\\hp\\Desktop\\Simulation_%s.xls',datestr(now,30))  
%save the data.it is succeed if the putput is 1. 0 means failed.
success=xlswrite(filename,{'MaterialProperty='},'Setting','A1')
success=xlswrite(filename,MatPro,'Setting','B2')
success=xlswrite(filename,{'LaserProperty='},'Setting','A6')
success=xlswrite(filename,LasPro,'Setting','B7')
success=xlswrite(filename,{'R&a='},'Setting','A11')
success=xlswrite(filename,RefAbs,'Setting','B12')
success=xlswrite(filename,{'xMesh='},'Setting','A16')
success=xlswrite(filename,xMesh,'Setting','B17')
success=xlswrite(filename,{'tMesh='},'Setting','A21')
success=xlswrite(filename,tMesh,'Setting','B22')
success=xlswrite(filename,{'Tmax'},'Tmax','B1')
success=xlswrite(filename,xAxis','Tmax','A2')
success=xlswrite(filename,TMax','Tmax','B2')
success=xlswrite(filename,ProfilAtDepth,'Tatx','B1')
success=xlswrite(filename,tAxis','Tatx','A2')
success=xlswrite(filename,Tatx','Tatx','B2')
end 
function [TMax,Tatx]=Calculation(MatPro,LasPro,RefAbs,xMesh,tMesh,ProfilAtDepth)
%extract the Temperature profile
%TMax:the array of highest temperature profile along depth(x)
%Tatx: the matrix of temperature profile along time at various depth(x)

%cite the function block
[Nx,Px,Dx,DxTrs]=xAxisExtraction(xMesh);
[Nt,Pt,Dt]=tAxisExtraction(tMesh);
[Cap,Den,Cond,LPR,a,R]=MatLasProtertyExtraction(MatPro,LasPro,RefAbs,xMesh);
[NL,PL]=LasPowFrofAtTimeExtraction(LasPro,tMesh);
%to tell whether the Dt is set smally enough
KMax=max(Cond);
CMin=min(Cap);
DMin=min(Den);
xMin=min(DxTrs);
r=KMax*Dt/(CMin*DMin*xMin^2);
if r>=0.5
    error("t mesh step is too big!! please compress it");
end
%creat profile matrix to save the temperature data at various depth
ND=size(ProfilAtDepth,2);
Tatx=zeros(ND,Nt);
%to tell whether the ProfileatDepth exceeds the x space
if min(ProfilAtDepth)<0||max(ProfilAtDepth)>max(Px)
    error("ProfilAtDepth exceeds the x space defined by xmesh, please correct it");
end
%creat array to save the hightest temperature along depth
TMax=zeros(1,Nx);  
%creat arrays to save the temperature profile of former moment and current  
T1=27*ones(1,Nx);    
T2=27*ones(1,Nx);
%initialize the rate of progress
CompletePercentage=0; 
for i=1:Nt   %start a loop on the timeline
    %show the rate of progress
    if round(i/Nt*100)~=CompletePercentage   %if the rate of progress changes, output new rate
        CompletePercentage=round(i/Nt*100);       
        fprintf("Complete Percentage: %d %%\n",CompletePercentage);
    end
    % calculate the temperature change of first x point 
    % calculate the powerdensity by thermalconductivity
    DP2=Cond(1)*(T1(2)-T1(1))/(DxTrs(1)*Dx(1)/10^4);   
    % calculate the powerdensity by laser irradiation
    DP3=0;
    for j=1:NL    %start a loop on the multiple laser
        DP3=DP3+PL(j,i)*(1-R(j))*a(j,1)*LPR(j,1);
    end
    % calculate the change of temperature 
    DP1=DP2+DP3;
    T2(1)=T1(1)+DP1*Dt/(Cap(1)*Den(1));
    % calculate the temperature change of inner x points
    for k=2:Nx-1  %x position loop
        DP2=Cond(k)*(T1(k+1)-T1(k))/(DxTrs(k)*Dx(k)/10^4)-Cond(k-1)*(T1(k)-T1(k-1))/(DxTrs(k-1)*Dx(k)/10^4);
        DP3=0;
        for j=1:NL
            DP3=DP3+PL(j,i)*(1-R(j))*a(j,k)*LPR(j,k);
        end
        DP1=DP2+DP3;
        T2(k)=T1(k)+DP1*Dt/(Cap(k)*Den(k));
    end
    % calculate the temperature change of last x point
    DP2=-Cond(Nx-1)*(T1(Nx)-T1(Nx-1))/(DxTrs(Nx-1)*Dx(Nx)/10^4);
    DP3=0;
    for j=1:NL  %laser loop
        DP3=DP3+PL(j,i)*(1-R(j))*a(j,Nx)*LPR(j,Nx);
    end
    DP1=DP2+DP3;
    T2(Nx)=T1(Nx)+DP1*Dt/(Cap(Nx)*Den(Nx));
    T1=T2;  
    %record the data we are going to stroe 
    NDC=1;
    for k=1:Nx
        if TMax(k)<T2(k)    %save the highest temperature profile
            TMax(k)=T2(k);  
        end
        if (-Dx(k)/2)<=(ProfilAtDepth(NDC)-Px(k))&&(ProfilAtDepth(NDC)-Px(k))<(Dx(k)/2)   %save the temperature profile at vaious depth
           Tatx(NDC,i)=T2(k);
           if NDC<ND
                NDC=NDC+1;   
           end
        end     
    end
end
end
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
function [LasNunber,LasPowProf]=LasPowFrofAtTimeExtraction(LasPro,tMesh)
%extract the laser power profile in time line
%LasNunber:the beam number of applied laser
%LasPowProf:laser power profile, unit:W/cm2

%cite function to get x related arrays
[tPointNumber,tAxis,tDelta]=tAxisExtraction(tMesh);
M=size(LasPro,1); %get the beam number of applied laser
LasNunber=M;   %number of laser beam
N=tPointNumber;  %number of x point
LasPowProf=zeros(M,N);  %creat laser power profile matrix with original values of 0
%fill the laser power profile matrix
for i=1:M
    %to tell if the profile is required as Gaussian
    if LasPro(i,5)==1
        tao=LasPro(i,2)/2.355;
        for j=1:N
            LasPowProf(i,j)=10^6*LasPro(i,3)*exp(-(tAxis(j)-3*tao-LasPro(i,4))^2/(2*(tao)^2))/(tao*(2*pi)^0.5);  
        end
    %to tell if the profile is required as trapezoid   
    elseif LasPro(i,5)==2
        %the trapezoid profile is defined with 3us incress and decress
        %process
        if LasPro(i,2)<3
            error("trapezoid laser pulsewidth must be over 3us")
        end
        for j=1:N
            if tAxis(j)<LasPro(i,4)
                LasPowProf(i,j)=0;
            elseif tAxis(j)<(LasPro(i,4)+3)
                LasPowProf(i,j)= (tAxis(j)-LasPro(i,4))/3*(10^6*LasPro(i,3)/LasPro(i,2));
            elseif tAxis(j)<(LasPro(i,4)+LasPro(i,2))
                LasPowProf(i,j)=10^6*LasPro(i,3)/LasPro(i,2);
            elseif tAxis(j)<(LasPro(i,4)+LasPro(i,2)+3)
                LasPowProf(i,j)= (LasPro(i,4)+LasPro(i,2)+3-tAxis(j))/3*(10^6*LasPro(i,3)/LasPro(i,2));
            else
                LasPowProf(i,j)=0;
            end
        end
    else
        error("LasPro's last colume defined wrong, it is only 1 or 2")
    end
end
end
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
