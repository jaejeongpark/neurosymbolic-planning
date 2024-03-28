%% Main file to use for a single call of the tool
% This tool computes an interval over-approximation of the finite-time
% reachable set for a continuous-time system dx/dt = f(t,x,p) or
% discrete-time system x^+ = f(t,x,p) with input p, based on the intervals 
% of initial states and of input values.

% Authors:  
%   Pierre-Jean Meyer, <pjmeyer -AT- berkeley.edu>, EECS, UC Berkeley
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
%   Murat Arcak, <arcak -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 13th of October 2018

%% Initialization
close all
clear
% Folder containing the various over-approximation methods
addpath('OA_methods')   
% Folder containing useful tools and functions
addpath('Utilities')    
% Folder containing user-provided files required by some methods 
% (eg: signs, bounds or functions of the Jacobian or Sensitivity matrices)
addpath('Input_files')  

%% Choice of the example system
% 1: Unicycle (discrete-time)
% 2: 4-state vehicle (discrete-time)
global system_choice
system_choice = 2;  

global dt
global v
global alpha
        
%% Definition of the reachability problem (for each example system)
%global u        % Some systems need an external control input
switch system_choice
    case 1
        %% Unicycle 
        dt = 0.1;
        v = 2;
        
        x_low = [0.; 0.; -pi/4];
        x_up =  [1.; 1.; pi/4];
        n_x = length(x_low);
        
        p_low = [-pi/2];
        p_up =  [pi/2];
        n_p = length(p_low);
        
        t_init = 0;        
        bool_discrete_time = 1;
    
    case 2
        %% Unicycle (K,b as control input)
        % state [x, y, theta]
        % input [K1, K2, K3, b]
        
        dt = 0.1;
        v = 2;
        
        x_low = [1.; 1.; 5*pi/4];
        x_up =  [2.; 2.; 7*pi/4];
        n_x = length(x_low);
        
        p_low = [-1.; -1.; -1.; -1.];
        p_up  = [0.; 0.; 0.; 0.];
        n_p = length(p_low);
        
        t_init = 0;        
        bool_discrete_time = 1;
        
    case 3
        %% 4-state vehicle 
        dt = 0.1;
        alpha = 0.2;
        
        x_low = [0.; 0.; -pi/4; -1];
        x_up =  [1.; 1.;  pi/4; 1];
        n_x = length(x_low);
        
        p_low = [-pi/2; -1.];
        p_up =  [ pi/2;  1.];
        n_p = length(p_low);
        
        t_init = 0;        
        bool_discrete_time = 1;
end       

%% Interval of initial states (defined by two column vectors)
%x_low = [0;0];      % Lower bound
%x_up = [0;0];       % Upper bound

%% Interval of allowed input values (defined by two column vectors)
% If the system has no input, define: p_low = 0; p_up = 0;
%p_low = [0;0];      % Lower bound
%p_up = [0;0];       % Upper bound

%% Time interval
%t_init = 0;         % Initial time
%t_final = 1;        % Final time (for continuous-time systems)

%% Call of the main over-approximation function
% succ_low is the lower bound of the over-approximation interval
% succ_up is the upper bound of the over-approximation interval

% Call for continuous-time systems over the time range [t_init,t_final]
% [succ_low,succ_up] = TIRA([t_init,t_final],x_low,x_up,p_low,p_up);

% Call for discrete-time systems for one step starting from time t_init
% [succ_low,succ_up] = TIRA(t_init,x_low,x_up,p_low,p_up);

if bool_discrete_time
    % Call for discrete-time systems for one step starting from time t_init
    [succ_low,succ_up] = TIRA(t_init,x_low,x_up,p_low,p_up);
else
    % Call for continuous-time systems over the time range [t_init,t_final]
    [succ_low,succ_up] = TIRA([t_init,t_final],x_low,x_up,p_low,p_up);
end

disp(succ_low)
disp(succ_up)
assert(false, 'done')

%% Optional choice of the over-approximation method to be used

% The optional choice of the over-approximation method can either be done
%   in the file Solver_parameters.m, 
%       then calling the function TIRA as above
%
%   or here, defining the integer 'OA_method' to one of the methods below:
%           1 - Continuous-time monotonicity
%           2 - Continuous-time componentwise contraction and growth bound
%           3 - Continuous-time mixed-monotonicity
%           4 - Continuous-time sampled-data mixed-monotonicity
%           5 - Discrete-time mixed-monotonicity
%       then adding 'OA_method' as the last argument of function TIRA:
%       [succ_low,succ_up] = TIRA([t_init,t_final],x_low,x_up,p_low,p_up,OA_method);
%       [succ_low,succ_up] = TIRA(t_init,x_low,x_up,p_low,p_up,OA_method);

% If it is not provided, function TIRA.m will pick the most
% suitable method that can be used with the additional system descriptions
% provided by the user in the './Input_files' folder.

%% Compute successors from random initial states and disturbances
sample_succ_number = 10000;
fprintf('\nCompute successors from %d random initial states and parameters ...\n', sample_succ_number)
tic
rand_succ = NaN(n_x,sample_succ_number);
for i = 1:sample_succ_number
    x0 = x_low + rand(n_x,1).*(x_up-x_low);
    p = p_low + rand(n_p,1).*(p_up-p_low);
    if bool_discrete_time
        rand_succ(:,i) = System_description(t_init,x0,p);
    else
        [~,x_traj] = ode45(@(t,x) System_description(t,x,p),[t_init t_final],x0);
        rand_succ(:,i) = x_traj(end,:)';
    end
end
toc

%% Plot the reachable set and over-approximations
% Only plot if the system has between 2 and 20 dimensions (ie max 10 plots)
if n_x > 1 && n_x <= 20         
    % Set indices of the subplots
    if ~mod(n_x,2)  
        % Even number of states
        plot_dimensions = [(1:2:n_x)' (2:2:n_x)'];
    else                
        % Odd number of states
        plot_dimensions = [(1:2:n_x-1)' (2:2:n_x)';n_x-1 n_x];
    end
    
    % Subplots
    for i = 1:size(plot_dimensions,1)
        figure
        hold on
        grid on

        % Plot the successors from random initial states
        for j = 1:sample_succ_number
            plot(rand_succ(plot_dimensions(i,1),j),rand_succ(plot_dimensions(i,2),j),'k.');
        end

        % Over-approximation interval
        handle_OA_IA = rectangle('Position',[succ_low(plot_dimensions(i,:))' (succ_up(plot_dimensions(i,:))-succ_low(plot_dimensions(i,:)))']);
        set(handle_OA_IA,'edgecolor',[0 0.5 0],'linewidth',2)

        % Legend
        xlabel(['$x_',num2str(plot_dimensions(i,1)),'$'],'Interpreter','Latex','FontSize',20)
        ylabel(['$x_',num2str(plot_dimensions(i,2)),'$'],'Interpreter','Latex','FontSize',20)
    end
end