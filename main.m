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