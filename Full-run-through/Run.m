clc; clear; close all;
%% Step 1: Time interval definition
% Input data
T= 1.0;
N= 2001; % Use "odd" number to then select half spectrum
start= 0.2;
Dt= T/(N-1);
t= linspace(start,start+T,N);
%%
fontsize=18;

%% Step 2: Input structural data and mesh
straightThermowell
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

pathF = 'C:\Users\callu\OpenFoam-cases\Thermowell\finalRuns\MPhil\3Dfinite\Straight\fluid\Re1e6\';
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

%Lift and drag check
fDrag = 0*t; fLift = 0*t; sDrag = 0*t; sLift = 0*t;
%Total Force check
fTotF = 0*t;
sTotF = 0*t;
error = 0*t;

constF= 0.5*998*9e-3*33.5^2;
%constF= 0.5*998*9e-3*16.7^2;
%constF= 0.5*998*9e-3*6.7^2;
%constF= 0.5*998*9e-3*3.35^2;

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
    
    fDrag(ti) = sum(fForce(:,1));
    fLift(ti) = sum(fForce(:,2));    
    sDrag(ti) = sum(sForce(:,1));
    sLift(ti) = sum(sForce(:,2));
    fTotF(ti) = norm(sum(fForce));
    sTotF(ti) = norm(sum(sForce));
    error(ti) = abs(norm(sum(fForce))-norm(sum(sForce)))/norm(sum(fForce));
%     if ti==1
%         writeTractionFile(pathS, sTraction, t(ti))
%         writeTractionFileF(pathF, fTraction, t(ti))
%     end
end

%%
data = dlmread(append(pathF,'postProcessing\forceCoeffs\0\coefficient.dat'), "\t", 13, 0);

time = data(:,1); 
Cd = data(:,2); 
Cl = data(:,4); 

t1=time;
time(t1<start) = []; time(t1>t(end)) =[];
Cl(t1<start) = []; Cl(t1>t(end)) =[];
Cd(t1<start) = []; Cd(t1>t(end)) =[];

%Cd(time <=0.4) = Cd(time<=0.4)*(3.35/16.7)^2;
%Cd(time >0.4) = Cd(time>0.4)*(6.7/16.7)^2;
%Cl(time <=0.4) = Cl(time<=0.4)*(3.35/16.7)^2;
%Cl(time >0.4) = Cl(time>0.4)*(6.7/16.7)^2;
%%
figure(3)
plot(time, Cd, 'b','LineWidth',1.25);
hold on
plot(time, Cl, 'r','LineWidth',1.25);
grid on;
% xticks([0.4 0.6 0.8 1.0 1.2 1.4])
ax = gca;
ax.FontSize = fontsize;
xlabel('Time (s)','Color','k') 
ylabel('Force Coefficient','Color','k') 
ylim([-0.5 1.5])
% xlim([0.4,1.4])
% saveas(gcf,append(pathF,'\ClCd'),'epsc')
%%
Fs = 1/(t(end)-t(1)); 
L = length(t)-rem(length(t),2);
L2 = length(time)-rem(length(time),2);

f = Fs*(0:L/2); 
f2 = Fs*(0:L2/2); 

fftL = fft(sLift);
P2 = abs(fftL/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

fftL2 = fft(Cl*constF);
P4 = abs(fftL2/L2);
P3 = P4(1:L2/2+1);
P3(2:end-1) = 2*P3(2:end-1);

%%
figure(5)
%plot(f2(2:end), P3(2:end), 'b', f(2:end), P1(2:end), 'r')
plot(f2(2:end), P3(2:end), 'b', 'LineWidth',1.25)
hold on;
plot(f(2:end), P1(2:end), 'r', 'LineWidth',1.25)
grid on;
ax = gca;
ax.FontSize = fontsize;
xlim([0,400])
%ylim([0,2.5])

xlabel('Frequency (Hz)','Color','k'); 
ylabel('Lift Amplitude (N)','Color','k');
%saveas(gcf,append(pathF,'\fftLift'),'epsc')

figure(6)
%plot(t, fTotF, 'b', 'LineWidth',1.25)
hold on;
plot(t, sTotF, 'b','LineWidth',1.25)
grid on; 
% xticks([0.4 0.6 0.8 1.0 1.2 1.4])
% xlim([0.4,1.4])

ax = gca;
ax.FontSize = fontsize;
xlabel('Time (s)','Color','k');
ylabel('Force (N)','Color','k');
saveas(gcf,append(pathF,'\totF'),'epsc')

figure(7)
plot(t, error,'b', 'LineWidth', 1.25)
% xticks([0.4 0.6 0.8 1.0 1.2 1.4])
% xlim([0.4,1.4])

grid on; 
ax = gca;
ax.FontSize = fontsize;
xlabel('Time (s)','Color','k');
ylabel('Relative Error','Color','k');
%saveas(gcf,append(pathF,'\errorTotF'),'epsc')

%% Step 7: Frequency solver. Defining nModes and nFreq to solve for

nModes= 50;
nFreq = 30;
U= ReducedOrderT(M, K, C, F, nModes, 2*nFreq, t); 
U= B*U;


%% Step 8: With the singular value decomposition

% rank= 10;
% U2= ReducedOrderSVD(M, K, C, F, nModes, 2*nFreq, rank, t); 
% U2= B*U2;


%% Step 9: Solve structural system using Newmark

% gamma=1/2;beta=1/4;
% u= NewMark(M, K, C, F, gamma, beta, Dt, N);
% u= B*u;


%% Step 10: Plotting the displacement of the Selected DOF 
 
%plot displacement time history 
figure(1)

plot(t, real(U(15,:)) - mean(real(U(15,:))),'b' , 'LineWidth', 1.25);
grid on;
xlabel('Time (s)','Color','k')
ylabel('Displacement (m)','Color','k')
% xlim([0.4,1.4])

ax = gca;
ax.FontSize = fontsize;
% xticks([0.4 0.6 0.8 1.0 1.2 1.4])
%xticks([0.15 0.25 0.35 0.45 0.55]);
%xlim([0.05, .45]);
%saveas(gcf,append(pathF,'\displacement'),'epsc')

%% Stress Calc and visualisation

[eqStress, p] = stress(nodes, elements, U, E, nu, t);
centre = elCentre(nodes,elements);

%%
% figure('Position', [10 10 600 900])
% scatter3(centre(:,1),centre(:,2),centre(:,3), 6, p(:,200), 'filled')
% patch('Faces',msh.TRIANGLES(:,1:3),'Vertices',nodes,'FaceColor','r','EdgeColor','k','FaceAlpha',0.)
% colormap jet
% c = colorbar;
% w = c.LineWidth;
% c.LineWidth = 0.5;
% c.FontSize=20;
% c.Location = 'southoutside';
% grid off;
% axis equal
% axis off
% saveas(gcf,append(pathF,'\elementp'),'epsc')

figure('Position', [10 10 1000 666])
scatter3(centre(:,1),centre(:,2),centre(:,3), 8, eqStress(:,200), 'filled')
patch('Faces',msh.TRIANGLES(:,1:3),'Vertices',nodes,'FaceColor','r','EdgeColor','k','FaceAlpha',0.)
colormap jet
c = colorbar;
w = c.LineWidth;
c.LineWidth = 0.5;
c.FontSize=18;
c.Location = 'southoutside';
grid off;
axis equal
axis off
%rotate(gcf,45,[0 0 0],[0 1 0]);
view([-90,200,-200])
camroll(-50)

%saveas(gcf,append(pathF,'\eqStress'),'epsc')

%%
figure;
plot(t,max(abs(p))-mean(max(abs(p))), 'b', 'LineWidth', 1.25);
grid on;
xlabel('Time (s)','Color','k');
ylabel('Pressure (Pa)','Color','k');
%xticks([0.4 0.6 0.8 1.0 1.2 1.4])
% xlim([0.4,1.4])

ax = gca;
ax.FontSize = fontsize;
%saveas(gcf,append(pathF,'\tvsp'),'epsc')

figure;
plot(t,max(abs(eqStress))-mean(max(abs(eqStress))),'b', 'LineWidth', 1.25);
grid on;
xlabel('Time (s)','Color','k');
ylabel('von Mises Stress (Pa)','Color','k');
%xlim([0.4,1.4])
%xticks([0.4 0.6 0.8 1.0 1.2 1.4])
ax = gca;
ax.FontSize = fontsize;
%saveas(gcf,append(pathF,'\tvseqStress'),'epsc')

fprintf('Mean displacement = %3f \n', mean(real(U(6,:))))
fprintf('Mean pressure = %3f \n', mean(max(abs(p))))
fprintf('Mean stress = %3f \n', mean(max(abs(eqStress))))

%% End
