%% User-provided function returning bounds on the Jacobian matrices
% The user can either provider global bounds, or a function of the input
% arguments (t_init,t_final,x_low,x_up,p_low,p_up) returning local bounds
% Jacobian definitions:
%   to states:  J_x(t) = d(System_description(t,x,p))/dx
%   to inputs:  J_p(t) = d(System_description(t,x,p))/dp

% List of inputs
%   t_init: initial time
%   t_final: time at which the reachable set is approximated (for continuous-time system only)
%       for discrete-time system, a dummy value can be provided
%   [x_low,x_up]: interval of initial states (at time t_init)
%   [p_low,p_up]: interval of allowed input values

% List of outputs
%   [J_x_low,J_x_up]: bounds of the Jacobian with respect to the state
%   [J_p_low,J_p_up]: bounds of the Jacobian with respect to the input

% Authors:  
%   Pierre-Jean Meyer, <pjmeyer -AT- berkeley.edu>, EECS, UC Berkeley
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 13th of October 2018

function [J_x_low,J_x_up,J_p_low,J_p_up] = UP_Jacobian_Bounds(t_init,t_final,x_low,x_up,p_low,p_up)
n_x = length(x_low);
n_p  = length(p_low);

%% Default values as NaN (not a number)
J_x_low = NaN(n_x);
J_x_up = NaN(n_x);
J_p_low = NaN(n_x,n_p);
J_p_up = NaN(n_x,n_p);

%% User-provided Jacobian bounds
% Can be either global bounds
% or local bounds depending on inputs: t_init,t_final,x_low,x_up,p_low,p_up

% If System_description.m has no input variable 'p', uncomment below:
% J_p_low = zeros(n_x,n_p);
% J_p_up = zeros(n_x,n_p);

global system_choice
global dt
global v
global alpha

switch system_choice
    case {1, 2, 3}
        theta_low = x_low(3);
        theta_up = x_up(3);
        % Min/Max for cosine
        theta_min_cos = pi;
        if theta_low>-pi && theta_up<pi
            theta_min_cos = max(abs(theta_low),abs(theta_up));
        end
        theta_max_cos = 0;
        if (theta_low>0 || theta_up<0) && theta_up<2*pi && theta_low>-2*pi
            theta_max_cos = min(abs(mod(theta_low+pi,2*pi)-pi),abs(mod(theta_up+pi,2*pi)-pi));
        end
        % Min/Max for sine
        theta_min_sin = 3*pi/2;
        if theta_low>-pi/2 && theta_up<3*pi/2
            theta_min_sin = max(abs(theta_low-pi/2),abs(theta_up-pi/2)) + pi/2;
        end
        theta_max_sin = pi/2;
        if (theta_low>pi/2 || theta_up<pi/2) && theta_up<5*pi/2 && theta_low>-3*pi/2
            theta_max_sin = min(abs(mod(theta_low-pi/2+pi,2*pi)-pi),abs(mod(theta_up-pi/2+pi,2*pi)-pi)) + pi/2;
        end
end
        

%% Jacobian bounds
switch system_choice
    case 1
        %% Unicycle 
        J_x_low = eye(n_x);
        J_x_up = eye(n_x);
        J_p_low = [0; 0; dt];
        J_p_up = [0; 0; dt];

        % Bounds of the partial derivatives
        if v > 0
            J_x_low(1,3) = -dt*v*sin(theta_max_sin);
            J_x_up(1,3) = -dt*v*sin(theta_min_sin);
            J_x_low(2,3) = dt*v*cos(theta_min_cos);
            J_x_up(2,3) = dt*v*cos(theta_max_cos);
        else
            J_x_low(1,3) = -dt*v*sin(theta_min_sin);
            J_x_up(1,3) = -dt*v*sin(theta_max_sin);
            J_x_low(2,3) = dt*v*cos(theta_max_cos);
            J_x_up(2,3) = dt*v*cos(theta_min_cos);
        end
        
    case 2
        %% Unicycle (K,b as control input)
        J_x_low = eye(n_x);
        J_x_up = eye(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up = zeros(n_x,n_p);
        
        % Bounds on J_p
        J_p_low(3,1) = dt*x_low(1);
        J_p_up(3,1) = dt*x_up(1);
        J_p_low(3,2) = dt*x_low(2);
        J_p_up(3,2) = dt*x_up(2);
        J_p_low(3,3) = dt*x_low(3);
        J_p_up(3,3) = dt*x_up(3);
        J_p_low(3,4) = dt;
        J_p_up(3,4) = dt;
        
         % Bounds of the partial derivatives
        if v > 0
            J_x_low(1,3) = -dt*v*sin(theta_max_sin);
            J_x_up(1,3) = -dt*v*sin(theta_min_sin);
            J_x_low(2,3) = dt*v*cos(theta_min_cos);
            J_x_up(2,3) = dt*v*cos(theta_max_cos);
        else
            J_x_low(1,3) = -dt*v*sin(theta_min_sin);
            J_x_up(1,3) = -dt*v*sin(theta_max_sin);
            J_x_low(2,3) = dt*v*cos(theta_max_cos);
            J_x_up(2,3) = dt*v*cos(theta_min_cos);
        end
        J_x_low(3,1) = dt*p_low(1);
        J_x_up(3,1) = dt*p_up(1);
        J_x_low(3,2) = dt*p_low(2);
        J_x_up(3,2) = dt*p_up(2);
        J_x_low(3,3) = 1+dt*p_low(3);
        J_x_up(3,3) = 1+dt*p_up(3);
        
    case 3
        %% 4-state vehicle 
        J_x_low = eye(n_x);
        J_x_up = eye(n_x);
        J_p_low = zeros(n_x, n_p);
        J_p_up = zeros(n_x, n_p);
        
        % Bounds on J_p
        J_p_low(3,1) = dt*x_low(4);
        J_p_up(3,1) = dt*x_up(4);
        J_p_low(4,2) = dt;
        J_p_up(4,2) = dt;
       
        % Bounds of the partial derivatives
        J_x_low(1,3) = min([-dt*x_low(4)*sin(theta_min_sin), -dt*x_low(4)*sin(theta_max_sin), -dt*x_up(4)*sin(theta_min_sin), -dt*x_up(4)*sin(theta_max_sin)]);
        J_x_up(1,3)  = max([-dt*x_low(4)*sin(theta_min_sin), -dt*x_low(4)*sin(theta_max_sin), -dt*x_up(4)*sin(theta_min_sin), -dt*x_up(4)*sin(theta_max_sin)]);
        
        J_x_low(2,3) = min([dt*x_low(4)*cos(theta_min_cos), dt*x_low(4)*cos(theta_max_cos), dt*x_up(4)*cos(theta_min_cos), dt*x_up(4)*cos(theta_max_cos)]);
        J_x_up(2,3)  = max([dt*x_low(4)*cos(theta_min_cos), dt*x_low(4)*cos(theta_max_cos), dt*x_up(4)*cos(theta_min_cos), dt*x_up(4)*cos(theta_max_cos)]);
        
        J_x_low(1,4) = dt*cos(theta_min_cos);
        J_x_up(1,4)  = dt*cos(theta_max_cos);
        
        J_x_low(2,4) = dt*sin(theta_min_sin);
        J_x_up(2,4)  = dt*sin(theta_max_sin);
        
        J_x_low(3,4) = dt*p_low(1);
        J_x_up(3,4)  = dt*p_up(1);
        
        J_x_low(4,4) = 1-dt*alpha;
        J_x_up(4,4)  = 1-dt*alpha;
        
    case 102
        %% Linear system 2 states
        J_x_low = zeros(n_x);
        J_x_up  = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up  = zeros(n_x,n_p);
        
        J_x_low(1,1) = 1+dt*p_low(1);
        J_x_up(1,1)  = 1+dt*p_up(1);
        J_x_low(1,2) = dt*p_low(2);
        J_x_up(1,2)  = dt*p_up(2);
        J_x_low(2,1) = dt*p_low(3);
        J_x_up(2,1)  = dt*p_up(3);
        J_x_low(2,2) = 1+dt*p_low(4);
        J_x_up(2,2)  = 1+dt*p_up(4);
        
        J_p_low(1,1) = dt*x_low(1);
        J_p_up(1,1)  = dt*x_up(1);
        J_p_low(1,2) = dt*x_low(2);
        J_p_up(1,2)  = dt*x_up(2);
        J_p_low(1,5) = dt;
        J_p_up(1,5) = dt;
        J_p_low(2,3) = dt*x_low(1);
        J_p_up(2,3)  = dt*x_up(1);
        J_p_low(2,4) = dt*x_low(2);
        J_p_up(2,4)  = dt*x_up(2);
        J_p_low(2,6) = dt;
        J_p_up(2,6) = dt;
    
    case 104
        %% Linear system 4 states
        J_x_low = zeros(n_x);
        J_x_up  = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up  = zeros(n_x,n_p);
        
        J_x_low(1,1) = 1;
        J_x_up(1,1)  = 1;
        J_x_low(1,3) = dt;
        J_x_up(1,3)  = dt;
        J_x_low(2,2) = 1;
        J_x_up(2,2)  = 1;
        J_x_low(2,4) = dt;
        J_x_up(2,4)  = dt;
        
        J_x_low(3,1) = dt*p_low(1);
        J_x_up(3,1)  = dt*p_up(1);
        J_x_low(3,2) = dt*p_low(2);
        J_x_up(3,2)  = dt*p_up(2);
        J_x_low(3,3) = 1+dt*p_low(3);
        J_x_up(3,3)  = 1+dt*p_up(3);
        J_x_low(3,4) = dt*p_low(4);
        J_x_up(3,4)  = dt*p_up(4);
        
        J_x_low(4,1) = dt*p_low(5);
        J_x_up(4,1)  = dt*p_up(5);
        J_x_low(4,2) = dt*p_low(6);
        J_x_up(4,2)  = dt*p_up(6);
        J_x_low(4,3) = dt*p_low(7);
        J_x_up(4,3)  = dt*p_up(7);
        J_x_low(4,4) = 1+dt*p_low(8);
        J_x_up(4,4)  = 1+dt*p_up(8);
        
        J_p_low(3,1) = dt*x_low(1);
        J_p_up(3,1)  = dt*x_up(1);
        J_p_low(3,2) = dt*x_low(2);
        J_p_up(3,2)  = dt*x_up(2);
        J_p_low(3,3) = dt*x_low(3);
        J_p_up(3,3)  = dt*x_up(3);
        J_p_low(3,4) = dt*x_low(4);
        J_p_up(3,4)  = dt*x_up(4);
        J_p_low(3,9) = dt;
        J_p_up(3,9)  = dt;
        
        J_p_low(4,5) = dt*x_low(1);
        J_p_up(4,5)  = dt*x_up(1);
        J_p_low(4,6) = dt*x_low(2);
        J_p_up(4,6)  = dt*x_up(2);
        J_p_low(4,7) = dt*x_low(3);
        J_p_up(4,7)  = dt*x_up(3);
        J_p_low(4,8) = dt*x_low(4);
        J_p_up(4,8)  = dt*x_up(4);
        J_p_low(4,10) = dt;
        J_p_up(4,10)  = dt;
        
    case 106
        %% Linear system 6 states
        J_x_low = zeros(n_x);
        J_x_up  = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up  = zeros(n_x,n_p);
        
        J_x_low(1,1) = 1;
        J_x_up(1,1)  = 1;
        J_x_low(1,3) = dt;
        J_x_up(1,3)  = dt;
        J_x_low(2,2) = 1;
        J_x_up(2,2)  = 1;
        J_x_low(2,5) = dt;
        J_x_up(2,5)  = dt;
        J_x_low(3,3) = 1;
        J_x_up(3,3)  = 1;
        J_x_low(3,4) = dt;
        J_x_up(3,4)  = dt;
        J_x_low(5,5) = 1;
        J_x_up(5,5)  = 1;
        J_x_low(5,6) = dt;
        J_x_up(5,6)  = dt;
        
        J_x_low(4,1) = dt*p_low(1);
        J_x_up(4,1)  = dt*p_up(1);
        J_x_low(4,2) = dt*p_low(2);
        J_x_up(4,2)  = dt*p_up(2);
        J_x_low(4,3) = dt*p_low(3);
        J_x_up(4,3)  = dt*p_up(3);
        J_x_low(4,4) = 1+dt*p_low(4);
        J_x_up(4,4)  = 1+dt*p_up(4);
        J_x_low(4,5) = dt*p_low(5);
        J_x_up(4,5)  = dt*p_up(5);
        J_x_low(4,6) = dt*p_low(6);
        J_x_up(4,6)  = dt*p_up(6);
        
        J_x_low(6,1) = dt*p_low(7);
        J_x_up(6,1)  = dt*p_up(7);
        J_x_low(6,2) = dt*p_low(8);
        J_x_up(6,2)  = dt*p_up(8);
        J_x_low(6,3) = dt*p_low(9);
        J_x_up(6,3)  = dt*p_up(9);
        J_x_low(6,4) = dt*p_low(10);
        J_x_up(6,4)  = dt*p_up(10);
        J_x_low(6,5) = dt*p_low(11);
        J_x_up(6,5)  = dt*p_up(11);
        J_x_low(6,6) = 1+dt*p_low(12);
        J_x_up(6,6)  = 1+dt*p_up(12);
        
        J_p_low(4,1) = dt*x_low(1);
        J_p_up(4,1)  = dt*x_up(1);
        J_p_low(4,2) = dt*x_low(2);
        J_p_up(4,2)  = dt*x_up(2);
        J_p_low(4,3) = dt*x_low(3);
        J_p_up(4,3)  = dt*x_up(3);
        J_p_low(4,4) = dt*x_low(4);
        J_p_up(4,4)  = dt*x_up(4);
        J_p_low(4,5) = dt*x_low(5);
        J_p_up(4,5)  = dt*x_up(5);
        J_p_low(4,6) = dt*x_low(6);
        J_p_up(4,6)  = dt*x_up(6);
        J_p_low(5,13) = dt;
        J_p_up(5,13)  = dt;
        
        J_p_low(6,7) = dt*x_low(1);
        J_p_up(6,7)  = dt*x_up(1);
        J_p_low(6,8) = dt*x_low(2);
        J_p_up(6,8)  = dt*x_up(2);
        J_p_low(6,9) = dt*x_low(3);
        J_p_up(6,9)  = dt*x_up(3);
        J_p_low(6,10) = dt*x_low(4);
        J_p_up(6,10)  = dt*x_up(4);
        J_p_low(6,11) = dt*x_low(5);
        J_p_up(6,11)  = dt*x_up(5);
        J_p_low(6,12) = dt*x_low(6);
        J_p_up(6,12)  = dt*x_up(6);
        J_p_low(6,14) = dt;
        J_p_up(6,14)  = dt;
        
    case 108
        %% Linear system 8 states
        J_x_low = zeros(n_x);
        J_x_up  = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up  = zeros(n_x,n_p);
        
        J_x_low(1,1) = 1;
        J_x_up(1,1)  = 1;
        J_x_low(1,3) = dt;
        J_x_up(1,3)  = dt;
        J_x_low(2,2) = 1;
        J_x_up(2,2)  = 1;
        J_x_low(2,6) = dt;
        J_x_up(2,6)  = dt;
        J_x_low(3,3) = 1;
        J_x_up(3,3)  = 1;
        J_x_low(3,4) = dt;
        J_x_up(3,4)  = dt;
        J_x_low(4,4) = 1;
        J_x_up(4,4)  = 1;
        J_x_low(4,5) = dt;
        J_x_up(4,5)  = dt;
        J_x_low(6,6) = 1;
        J_x_up(6,6)  = 1;
        J_x_low(6,7) = dt;
        J_x_up(6,7)  = dt;
        J_x_low(7,7) = 1;
        J_x_up(7,7)  = 1;
        J_x_low(7,8) = dt;
        J_x_up(7,8)  = dt;
        
        J_x_low(5,1) = dt*p_low(1);
        J_x_up(5,1)  = dt*p_up(1);
        J_x_low(5,2) = dt*p_low(2);
        J_x_up(5,2)  = dt*p_up(2);
        J_x_low(5,3) = dt*p_low(3);
        J_x_up(5,3)  = dt*p_up(3);
        J_x_low(5,4) = dt*p_low(4);
        J_x_up(5,4)  = dt*p_up(4);
        J_x_low(5,5) = 1+dt*p_low(5);
        J_x_up(5,5)  = 1+dt*p_up(5);
        J_x_low(5,6) = dt*p_low(6);
        J_x_up(5,6)  = dt*p_up(6);
        J_x_low(5,7) = dt*p_low(7);
        J_x_up(5,7)  = dt*p_up(7);
        J_x_low(5,8) = dt*p_low(8);
        J_x_up(5,8)  = dt*p_up(8);
        
        J_x_low(8,1) = dt*p_low(9);
        J_x_up(8,1)  = dt*p_up(9);
        J_x_low(8,2) = dt*p_low(10);
        J_x_up(8,2)  = dt*p_up(10);
        J_x_low(8,3) = dt*p_low(11);
        J_x_up(8,3)  = dt*p_up(11);
        J_x_low(8,4) = dt*p_low(12);
        J_x_up(8,4)  = dt*p_up(12);
        J_x_low(8,5) = dt*p_low(13);
        J_x_up(8,5)  = dt*p_up(13);
        J_x_low(8,6) = dt*p_low(14);
        J_x_up(8,6)  = dt*p_up(14);
        J_x_low(8,7) = dt*p_low(15);
        J_x_up(8,7)  = dt*p_up(15);
        J_x_low(8,8) = 1+dt*p_low(16);
        J_x_up(8,8)  = 1+dt*p_up(16);
        
        J_p_low(5,1) = dt*x_low(1);
        J_p_up(5,1)  = dt*x_up(1);
        J_p_low(5,2) = dt*x_low(2);
        J_p_up(5,2)  = dt*x_up(2);
        J_p_low(5,3) = dt*x_low(3);
        J_p_up(5,3)  = dt*x_up(3);
        J_p_low(5,4) = dt*x_low(4);
        J_p_up(5,4)  = dt*x_up(4);
        J_p_low(5,5) = dt*x_low(5);
        J_p_up(5,5)  = dt*x_up(5);
        J_p_low(5,6) = dt*x_low(6);
        J_p_up(5,6)  = dt*x_up(6);
        J_p_low(5,7) = dt*x_low(7);
        J_p_up(5,7)  = dt*x_up(7);
        J_p_low(5,8) = dt*x_low(8);
        J_p_up(5,8)  = dt*x_up(8);
        J_p_low(5,17) = dt;
        J_p_up(5,17)  = dt;
        
        J_p_low(8,9) = dt*x_low(1);
        J_p_up(8,9)  = dt*x_up(1);
        J_p_low(8,10) = dt*x_low(2);
        J_p_up(8,10)  = dt*x_up(2);
        J_p_low(8,11) = dt*x_low(3);
        J_p_up(8,11)  = dt*x_up(3);
        J_p_low(8,12) = dt*x_low(4);
        J_p_up(8,12)  = dt*x_up(4);
        J_p_low(8,13) = dt*x_low(5);
        J_p_up(8,13)  = dt*x_up(5);
        J_p_low(8,14) = dt*x_low(6);
        J_p_up(8,14)  = dt*x_up(6);
        J_p_low(8,15) = dt*x_low(7);
        J_p_up(8,15)  = dt*x_up(7);
        J_p_low(8,16) = dt*x_low(8);
        J_p_up(8,16)  = dt*x_up(8);
        J_p_low(8,18) = dt;
        J_p_up(8,18)  = dt;
        
    case 110
        %% Linear system 10 states
        J_x_low = zeros(n_x);
        J_x_up  = zeros(n_x);
        J_p_low = zeros(n_x,n_p);
        J_p_up  = zeros(n_x,n_p);
        
        J_x_low(1,1) = 1;
        J_x_up(1,1)  = 1;
        J_x_low(1,3) = dt;
        J_x_up(1,3)  = dt;
        J_x_low(2,2) = 1;
        J_x_up(2,2)  = 1;
        J_x_low(2,7) = dt;
        J_x_up(2,7)  = dt;
        J_x_low(3,3) = 1;
        J_x_up(3,3)  = 1;
        J_x_low(3,4) = dt;
        J_x_up(3,4)  = dt;
        J_x_low(4,4) = 1;
        J_x_up(4,4)  = 1;
        J_x_low(4,5) = dt;
        J_x_up(4,5)  = dt;
        J_x_low(5,5) = 1;
        J_x_up(5,5)  = 1;
        J_x_low(5,6) = dt;
        J_x_up(5,6)  = dt;
        J_x_low(7,7) = 1;
        J_x_up(7,7)  = 1;
        J_x_low(7,8) = dt;
        J_x_up(7,8)  = dt;
        J_x_low(8,8) = 1;
        J_x_up(8,8)  = 1;
        J_x_low(8,9) = dt;
        J_x_up(8,9)  = dt;
        J_x_low(9,9) = 1;
        J_x_up(9,9)  = 1;
        J_x_low(9,10) = dt;
        J_x_up(9,10)  = dt;
        
        J_x_low(6,1) = dt*p_low(1);
        J_x_up(6,1)  = dt*p_up(1);
        J_x_low(6,2) = dt*p_low(2);
        J_x_up(6,2)  = dt*p_up(2);
        J_x_low(6,3) = dt*p_low(3);
        J_x_up(6,3)  = dt*p_up(3);
        J_x_low(6,4) = dt*p_low(4);
        J_x_up(6,4)  = dt*p_up(4);
        J_x_low(6,5) = dt*p_low(5);
        J_x_up(6,5)  = dt*p_up(5);
        J_x_low(6,6) = 1+dt*p_low(6);
        J_x_up(6,6)  = 1+dt*p_up(6);
        J_x_low(6,7) = dt*p_low(7);
        J_x_up(6,7)  = dt*p_up(7);
        J_x_low(6,8) = dt*p_low(8);
        J_x_up(6,8)  = dt*p_up(8);
        J_x_low(6,9) = dt*p_low(9);
        J_x_up(6,9)  = dt*p_up(9);
        J_x_low(6,10) = dt*p_low(10);
        J_x_up(6,10)  = dt*p_up(10);
        
        J_x_low(10,1) = dt*p_low(11);
        J_x_up(10,1)  = dt*p_up(11);
        J_x_low(10,2) = dt*p_low(12);
        J_x_up(10,2)  = dt*p_up(12);
        J_x_low(10,3) = dt*p_low(13);
        J_x_up(10,3)  = dt*p_up(13);
        J_x_low(10,4) = dt*p_low(14);
        J_x_up(10,4)  = dt*p_up(14);
        J_x_low(10,5) = dt*p_low(15);
        J_x_up(10,5)  = dt*p_up(15);
        J_x_low(10,6) = dt*p_low(16);
        J_x_up(10,6)  = dt*p_up(16);
        J_x_low(10,7) = dt*p_low(17);
        J_x_up(10,7)  = dt*p_up(17);
        J_x_low(10,8) = dt*p_low(18);
        J_x_up(10,8)  = dt*p_up(18);
        J_x_low(10,9) = dt*p_low(19);
        J_x_up(10,9)  = dt*p_up(19);
        J_x_low(10,10) = 1+dt*p_low(20);
        J_x_up(10,10)  = 1+dt*p_up(20);
        
        J_p_low(6,1) = dt*x_low(1);
        J_p_up(6,1)  = dt*x_up(1);
        J_p_low(6,2) = dt*x_low(2);
        J_p_up(6,2)  = dt*x_up(2);
        J_p_low(6,3) = dt*x_low(3);
        J_p_up(6,3)  = dt*x_up(3);
        J_p_low(6,4) = dt*x_low(4);
        J_p_up(6,4)  = dt*x_up(4);
        J_p_low(6,5) = dt*x_low(5);
        J_p_up(6,5)  = dt*x_up(5);
        J_p_low(6,6) = dt*x_low(6);
        J_p_up(6,6)  = dt*x_up(6);
        J_p_low(6,7) = dt*x_low(7);
        J_p_up(6,7)  = dt*x_up(7);
        J_p_low(6,8) = dt*x_low(8);
        J_p_up(6,8)  = dt*x_up(8);
        J_p_low(6,9) = dt*x_low(9);
        J_p_up(6,9)  = dt*x_up(9);
        J_p_low(6,10) = dt*x_low(10);
        J_p_up(6,10)  = dt*x_up(10);
        J_p_low(6,21) = dt;
        J_p_up(6,21)  = dt;
        
        J_p_low(10,11) = dt*x_low(1);
        J_p_up(10,11)  = dt*x_up(1);
        J_p_low(10,12) = dt*x_low(2);
        J_p_up(10,12)  = dt*x_up(2);
        J_p_low(10,13) = dt*x_low(3);
        J_p_up(10,13)  = dt*x_up(3);
        J_p_low(10,14) = dt*x_low(4);
        J_p_up(10,14)  = dt*x_up(4);
        J_p_low(10,15) = dt*x_low(5);
        J_p_up(10,15)  = dt*x_up(5);
        J_p_low(10,16) = dt*x_low(6);
        J_p_up(10,16)  = dt*x_up(6);
        J_p_low(10,17) = dt*x_low(7);
        J_p_up(10,17)  = dt*x_up(7);
        J_p_low(10,18) = dt*x_low(8);
        J_p_up(10,18)  = dt*x_up(8);
        J_p_low(10,19) = dt*x_low(9);
        J_p_up(10,19)  = dt*x_up(9);
        J_p_low(10,20) = dt*x_low(10);
        J_p_up(10,20)  = dt*x_up(10);
        J_p_low(10,22) = dt;
        J_p_up(10,22)  = dt;
end




















