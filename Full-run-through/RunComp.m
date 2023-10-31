clc; clear; close all;
%% Step 1: Time interval definition
% Input data
T=0.4;
N=801; % Use "odd" number to then select half spectrum
start= 0.05;
Dt=T/(N-1);
t=linspace(start,start+T,N);


%% Step 2: Input structural data and mesh

straightThermowell
% Modulus of elasticity, Material density, Poisson's ratio
E= 200e9; nu= 0.3; rho= 7850; 

model= createpde('structural');
nodes= 0.001*msh.POS; elements = msh.TETS;

% Creates geometry and creates model.mesh
geom= geometryFromMesh(model,nodes',elements(:,1:4)');    

% Applying structural properties and assembling object containing FEM matrices
structuralProperties(model,'YoungsModulus',E,'PoissonsRatio',nu,'MassDensity',rho);
FEM= assembleFEMatrices(model);

% Fixing DOFs based on their nodal coordinates
B1= reshape([nodes(:,3)>0 nodes(:,3)>0 nodes(:,3)>0], [length(nodes)*3,1]);
B= diag(double(B1)); B(:, all(~B,1))= []; B= sparse(B);

% Mass Matrix, Stiffness Matrix, Rayleigh damping
d1= 0; d2= 1;
M= B'*FEM.M*B; K= B'*FEM.K*B; C= d1*M+d2*K;


%% Step 3: Reading face centres of fluid mesh

pathF = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\fluid\Re5e5\';
fFaces = ReadingFacesFluid(pathF);
[fAreas, fCentres] = ReadingCentres(pathF,fFaces);


%% Step 4: Reading face centres of solid mesh

pathS = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\solid\';
sFaces = ReadingFacesSolid(pathS);
[sAreas, sCentres] = ReadingCentres(pathS,sFaces);


%% Step 5: Find nearest solid mesh face center to each fluid mesh face

[~,cNearest] = pdist2(fCentres, sCentres, 'euclidean', 'Smallest', 1);


%% Step 6: Read and map forces from fluid mesh to solid mesh

%Initialise F and loop through time files
F =zeros(length(M),length(t));
for ti=1:length(t)
    fForce= ReadingForces(pathF, fFaces, t(ti));
    %Calculate tractions to map
    fTraction = fForce./fAreas;

    %Interpolate
    sTraction = zeros(size(sCentres));
    for i = 1:length(cNearest)
        sTraction(i,:) = fTraction(cNearest(i),:);
    end
    
    %Calculate force on each face of solid mesh
    sForce = sTraction.*sAreas;
    
    %Apply force to nodes and eliminate fixed nodes
    Ft = ApplyF(sForce, nodes, sFaces);
    Ft = B'*Ft;
    F(:,ti) = Ft; 
end


%% Step 7: Frequency solver. Defining nModes and nFreq to solve for

nModes= 20;
nFreq = 30;
U= ReducedOrderT(M, K, C, F, nModes, 2*nFreq, t); 
U= B*U;
U = reshape(U, [length(nodes),3,length(t)]);


%% Step 2: Input structural data and mesh

helixThermowell
% Modulus of elasticity, Material density, Poisson's ratio
E= 200e9; nu= 0.3; rho= 7850; 

model= createpde('structural');
nodes= 0.001*msh.POS; elements = msh.TETS;

% Creates geometry and creates model.mesh
geom= geometryFromMesh(model,nodes',elements(:,1:4)');    

% Applying structural properties and assembling object containing FEM matrices
structuralProperties(model,'YoungsModulus',E,'PoissonsRatio',nu,'MassDensity',rho);
FEM= assembleFEMatrices(model);

% Fixing DOFs based on their nodal coordinates
B1= reshape([nodes(:,3)>0 nodes(:,3)>0 nodes(:,3)>0], [length(nodes)*3,1]);
B= diag(double(B1)); B(:, all(~B,1))= []; B= sparse(B);

% Mass Matrix, Stiffness Matrix, Rayleigh damping
d1= 0; d2= 1;
M= B'*FEM.M*B; K= B'*FEM.K*B; C= d1*M+d2*K;


%% Step 3: Reading face centres of fluid mesh

pathF = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Helical\fluid\Re5e5\';
fFaces = ReadingFacesFluid(pathF);
[fAreas, fCentres] = ReadingCentres(pathF,fFaces);


%% Step 4: Reading face centres of solid mesh

pathS = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Helical\solid\';
sFaces = ReadingFacesSolid(pathS);
[sAreas, sCentres] = ReadingCentres(pathS,sFaces);


%% Step 5: Find nearest solid mesh face center to each fluid mesh face

[~,cNearest] = pdist2(fCentres, sCentres, 'euclidean', 'Smallest', 1);


%% Step 6: Read and map forces from fluid mesh to solid mesh

%Initialise F and loop through time files
F =zeros(length(M),length(t));
for ti=1:length(t)
    fForce= ReadingForces(pathF, fFaces, t(ti));
    %Calculate tractions to map
    fTraction = fForce./fAreas;

    %Interpolate
    sTraction = zeros(size(sCentres));
    for i = 1:length(cNearest)
        sTraction(i,:) = fTraction(cNearest(i),:);
    end
    
    %Calculate force on each face of solid mesh
    sForce = sTraction.*sAreas;
    
    %Apply force to nodes and eliminate fixed nodes
    Ft = ApplyF(sForce, nodes, sFaces);
    Ft = B'*Ft;
    F(:,ti) = Ft; 
end


%% Step 7: Frequency solver. Defining nModes and nFreq to solve for

nModes= 20;
nFreq = 30;
U2= ReducedOrderT(M, K, C, F, nModes, 2*nFreq, t); 
U2= B*U2;

U2 = reshape(U2, [length(nodes),3,length(t)]);

%% Step 9: Plotting the displacement of the Selected DOF 

%plot displacement time history 
figure(1)
hold on

plot(t, squeeze(vecnorm(real(U(6,:,:)))),'k', t, squeeze(vecnorm(real(U2(15,:,:)))),'r')
grid on;
xlabel('Time (s)','Interpreter','latex');
ylabel('Displacement (m)','Interpreter','latex');
saveas(gcf,'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\displacementComp5e5.png')

%% End
