clear;clc;
mrstModule add ad-core ad-blackoil ad-props mrst-gui ad-fi
% Define model ------------------------------------------------------------
nx = 25; ny = 20; nz = 5;
G = cartGrid([nx, ny, nz]);
G = computeGeometry(G);
%rock.perm  = repmat(300*milli*darcy, [G.cells.num, 1]);
% Set up permeability based on K-indices
%rock.poro  = repmat(0.3, [G.cells.num, 1]);
[I, J, K] = gridLogicalIndices(G);

px       = 200*milli*darcy*ones(G.cells.num,1);
px(K==2) = 400*milli*darcy;
px(K==3) = 600*milli*darcy;
px(K==4) = 800*milli*darcy;
px(K==5) = 300*milli*darcy;

% Introduce anisotropy by setting K_x = 5*K_z.
perm = [px, px, 0.1*px];
rock = makeRock(G, perm, 0.3);

                 

% Wells and initial rates -------------------------------------------------

% one injection wells
c1 = 4*(nx*ny) + (2:15)'; c2= 1002:nx:1275; c3 = 3*(nx*ny) + (2:15)';
%nInj = 1;
%cellsWell1 =  1 : nx*ny : nx*ny*nz;
%radius     = .1;
W = struct([]);
W = addWell(W, G, rock, (500:nx*ny:2500),'Name', 'I', 'radius', 5*inch,...
'Type', 'rate', 'Val', 1.0/day(), 'comp_i', [1, 0], 'Sign', 1);

%cellsWell2 =  nx : nx*ny : nx*ny*nz;
%W = addWell(W, G, rock, cellsWell2, 'Type', 'rate', ...
 %           'Val', 1.0/day(), 'Radius', radius, 'comp_i', [1,0], 'sign', 1, 'name', 'I2');
disp('Well #1: '); display(W(1));

% Two production wells
%nProd = 2;
W = addWell(W, G, rock, c1, 'Name', 'P1', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'comp_i', [0, 1], 'Sign', -1);
W = addWell(W, G, rock, c2, 'Name', 'P2', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'comp_i', [0, 1],'Sign', -1);
W = addWell(W, G, rock, c3, 'Name', 'P3', 'radius', 5*inch, ...
'Type', 'bhp', 'Val', 1000*psia, 'comp_i', [0, 1], 'Sign', -1);

 %       disp('Well #4: '); display(W(4));
plotGrid(G, 'FaceColor', 'none'); view(3);
[ht, htxt, hs] = plotWell(G, W, 'radius', 0.1, 'height', 2, 'Color', 'r');
set(htxt, 'FontSize', 16);

figure,            
plotCellData(G, log(rock.perm(:,1)));
plotWell(G, W), view([1 1 1])
%%
pv = poreVolume(G, rock);
time = 720;
irate = sum(pv)/(time); %inject one pore volume over simulation period
% Prod schedule
n = 100; % number of time steps
dt = time/n;
timesteps = repmat(dt, n, 1);
% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W);
%%
fluid = initSimpleADIFluid('phases', 'WO', ...
                           'rho', [1000, 700], ...
                           'n', [2, 2], ...
                           'mu', [1, 5]*centi*poise ...
                           );
                    
 
% Set up twophase, immiscible model with fully implicit discretization
model = TwoPhaseOilWaterModel(G, rock, fluid);

% Set up initial reservoir at 5000 psi pressure and 90% oil saturation.
State0 = initResSol(G, 5000*psia, [0.1, 0.9]);
Time = cumsum(schedule.step.val);

[wellSols, state, report]= simulateScheduleAD(State0, model, schedule);
plotWellSols(wellSols, Time);


qOs     =  -getWellOutput(wellSols, 'qOs', 2);
qWs     =  -getWellOutput(wellSols, 'qWs', 2);
% Compute NPV
d   = 0.05;    % yearly discount factor
ro  = 26000;      % oil revenue/price ($/stb)
rwp =  2154;      % water production handling costs ($/stb)
rwi =  2154;      % water injection cost ($/stb) 

npvopts = {'OilPrice',             ro , ...
           'WaterProductionCost', rwp , ...
           'WaterInjectionCost',  rwi , ...
           'DiscountFactor',        d};
vals     = cell2mat(NPVOW(G, wellSols, schedule, npvopts{:}));

% Plot evolution of NPV and indicate peak value:
npv = cumsum(vals);
figure,  plot(time, cumsum(vals), '--b', 'LineWidth', 2);
title('Evolution of NPV [naira]')
legend('base')
xlabel('Time')
ylabel('NPV')