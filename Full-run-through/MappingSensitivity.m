clc; clear; close all;
%% Step 1: Time interval definition
% Input data
T= 1.0;
N= 401; % Use "odd" number to then select half spectrum
start= 3.0;
Dt= T/(N-1);
t= linspace(start,start+T,N);

%% Step 2: Input structural data and mesh
straightAx120A50
% Modulus of elasticity, Material density, Poisson's ratio
E= 200e9; nu= 0.3; rho= 7850; 

model= createpde('structural');
nodes= 0.001*msh.POS; elements = msh.TETS(:,1:4)';

% Creates geometry and creates model.mesh
geom= geometryFromMesh(model,nodes',elements);    

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

pathF = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\fluid\Re1e5\';
fFaces = ReadingFacesFluid(pathF);
[fAreas, fCentres] = ReadingCentres(pathF,fFaces);

%% Step 4: Reading face centres of solid mesh

%pathS = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\solid\';
pathS = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\L2test\Ax120A50\';
sFaces = ReadingFacesSolid(pathS);
[sAreas, sCentres] = ReadingCentres(pathS,sFaces);

%% Step 5: Find nearest solid mesh face center to each fluid mesh face

[~,cNearest] = pdist2(fCentres, sCentres, 'euclidean', 'Smallest', 1);

%% Step 6: Read and map forces from fluid mesh to solid mesh

%Initialise F and loop through time files
F =zeros(length(M),1);

ti=1;
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
%%
Qi = vecnorm(sForce,2,2);
L2 = sqrt(sum(Qi.^2));

error = abs(norm(sum(fForce))-norm(sum(sForce)))/norm(sum(fForce));
[Mode, nat] = eigs(K, M, 2, 'smallestabs');                          % Calculate first nModes mode shapes
%%
answer = sqrt(nat)/(2*pi);

%%
ang=[20 30 40 50];
l2= [2.7706 1.4679 1.4101 1.3732];
plot(ang, l2,'b','LineWidth',1.25)
hold on;
plot(ang, l2,'r.','MarkerSize',18)
xlabel('Angular Divisions','Interpreter','latex','fontname','Times New Roman'); 
ylabel('L2 norm of the Mapped Force','Interpreter','latex','fontname','Times New Roman'); 
grid on;
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 18;
saveas(gcf,append(pathS,'\L2test2'),'epsc')

axial=[20 40 60 80 100 120 140];
l2= [7.5293 3.8413 2.7706 2.1745 1.8332 1.6098 1.452];
figure;
plot(axial, l2,'b','LineWidth',1.25)
hold on;
plot(axial, l2,'r.','MarkerSize',18)
xlabel('Axial Divisions','Interpreter','latex','fontname','Times New Roman'); 
ylabel('L2 norm of the Mapped Force','Interpreter','latex','fontname','Times New Roman'); 
grid on;
ax = gca;
xlim([20,140])
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 18;
saveas(gcf,append(pathS,'\L2test'),'epsc')