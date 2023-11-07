clc; clear; close all;
%% Step 1: Time interval definition

%Time period in seconds
T= 1.0;

%Number of time steps
N= 2001; % Use "odd" number to then select half spectrum

%Starting time
start= 0.2;

%Time step
Dt= T/(N-1);

%Time vector
t= linspace(start,start+T,N);

%% Step 2: Input structural data and mesh
%Import matlab solid mesh
straightThermowell % helixThermowell is used for Helical cases

% Modulus of elasticity, Material density, Poisson's ratio
E= 200e9; nu= 0.3; rho= 7850; 

%Use pde package
model= createpde('structural');

%Import node positions and connectivity
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

%Read fluid mesh faces and face centers
pathF = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\fluid\Re1e6\';
fFaces = ReadingFacesFluid(pathF);
[fAreas, fCentres] = ReadingCentres(pathF,fFaces);


%% Step 4: Reading face centres of solid mesh

%Read solid mesh faces and face centers
pathS = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\solid\';
sFaces = ReadingFacesSolid(pathS);
[sAreas, sCentres] = ReadingCentres(pathS,sFaces);


%% Step 5: Find nearest solid mesh face center to each fluid mesh face

%Find nearest solid mesh face center for each fluid mesh face center
[~,cNearest] = pdist2(fCentres, sCentres, 'euclidean', 'Smallest', 1);


%% Step 6: Read and map forces from fluid mesh to solid mesh

%Initialise F and loop through time files
F =zeros(length(M),length(t));

for ti=1:length(t)
    % Read Force vector field from Openfoam file
    fForce= ReadingForces(pathF, fFaces, t(ti));

    % Calculate tractions to map
    fTraction = fForce./fAreas;

    % Interpolate
    sTraction = zeros(size(sCentres));
    for i = 1:length(cNearest)
        sTraction(i,:) = fTraction(cNearest(i),:);
    end
    
    % Calculate force on each face of solid mesh
    sForce = sTraction.*sAreas;    
    
    % Apply force to nodes and eliminate fixed nodes
    Ft = ApplyF(sForce, nodes, sFaces);
    Ft = B'*Ft;
    F(:,ti) = Ft; 

end


%% Step 7: Frequency solver. Defining nModes and nFreq to solve for

%Defining number of modes and frequencies solve for
nModes= 50;
nFreq = 30;

%Solve for displacement 
U= ReducedOrderT(M, K, C, F, nModes, 2*nFreq, t);

%Include fixed nodes
U= B*U;

%% Step 8: Solve structural system using Newmark

%Solve for displacement using benchmarking method
gamma=1/2;beta=1/4;
u= NewMark(M, K, C, F, gamma, beta, Dt, N);

%Include fixed nodes
u= B*u;

%% Step 9: Stress Calc and visualisation

%Calculate stress for each element
[eqStress, p] = stress(nodes, elements, U, E, nu, t);


%% End
