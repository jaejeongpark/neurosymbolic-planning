%% Descripion of the time-varying system with input
% either continuous-time dynamics: dx/dt = f(t,x,p)
% or discrete-time successor: x^+ = f(t,x,p)

% List of inputs
%   t: time
%   x: state
%   p: input value

% List of outputs
%   dx: 
%    continuous-time: vector field evaluated at time t, state x and input p
%    discrete-time: one-step successor from state x at time t with input p  

% Authors:  
%   Pierre-Jean Meyer, <pjmeyer -AT- berkeley.edu>, EECS, UC Berkeley
%   Alex Devonport, <alex_devonport -AT- berkeley.edu>, EECS, UC Berkeley
% Date: 13th of October 2018

function dx = System_description(t,x,p)

%% Default value as NaN vector (not a number)
n_x = length(x);
dx = NaN(n_x,1);

%% User-provided system description
% For continuous-time system: dx is the time derivative dx/dt
% For discrete-time system: dx is the one-step successor x^+

global system_choice
global dt
global v
global alpha

switch system_choice
    case 1
        %% Unicycle 
        % states:
        %     x(1), x(2): planar (x,y) coordinates of the unicycle
        %     x(3): unicycle facing angle
        % inputs:
        %     p: angular velocity
        
        dx = x + [dt*v*cos(x(3)); ...
                  dt*v*sin(x(3)); ...
                  dt*p(1)];
              
    case 2
        %% Unicycle (K,b as control input)
        % State: [x, y, theta]
        % Input: [K1, K2, K3, b]
        
        dx = x + [dt*v*cos(x(3)); ...
                  dt*v*sin(x(3)); ...
                  dt*(p(1)*x(1)+p(2)*x(2)+p(3)*x(3)+p(4))];
              
    case 3
        %% 4-state vehicle    
        dx = x + [dt*x(4)*cos(x(3)); ...
                  dt*x(4)*sin(x(3)); ...
                  dt*x(4)*p(1); ...
                  -alpha*dt*x(4)+dt*p(2)];
    
    case 102
        %% Linear system 2 states
        dx = x + [dt*p(1)*x(1)+dt*p(2)*x(2)+dt*p(5); ...
                  dt*p(3)*x(1)+dt*p(4)*x(2)+dt*p(6)];
              
    case 104
        %% Linear system 4 states
        dx = x + [dt*x(3); ...
                  dt*x(4); ...
                  dt*p(1)*x(1)+dt*p(2)*x(2)+dt*p(3)*x(3)+dt*p(4)*x(4)+dt*p(9); ...
                  dt*p(5)*x(1)+dt*p(6)*x(2)+dt*p(7)*x(3)+dt*p(8)*x(4)+dt*p(10)];
              
    case 106
        %% Linear system 6 states
        dx = x + [dt*x(3); ...
                  dt*x(5); ...
                  dt*x(4); ...
                  dt*p(1)*x(1)+dt*p(2)*x(2)+dt*p(3)*x(3)+dt*p(4)*x(4)+dt*p(5)*x(5)+dt*p(6)*x(6)+dt*p(13); ...
                  dt*x(6); ...
                  dt*p(7)*x(1)+dt*p(8)*x(2)+dt*p(9)*x(3)+dt*p(10)*x(4)+dt*p(11)*x(5)+dt*p(12)*x(6)+dt*p(14)];
              
    case 108
        %% Linear system 8 states
        dx = x + [dt*x(3); ...
                  dt*x(6); ...
                  dt*x(4); ...
                  dt*x(5); ...
                  dt*p(1)*x(1)+dt*p(2)*x(2)+dt*p(3)*x(3)+dt*p(4)*x(4)+dt*p(5)*x(5)+dt*p(6)*x(6)+dt*p(7)*x(7)+dt*p(8)*x(8)+dt*p(17); ...
                  dt*x(7); ...
                  dt*x(8); ...
                  dt*p(9)*x(1)+dt*p(10)*x(2)+dt*p(11)*x(3)+dt*p(12)*x(4)+dt*p(13)*x(5)+dt*p(14)*x(6)+dt*p(15)*x(7)+dt*p(16)*x(8)+dt*p(18)];
    
    case 110
        %% Linear system 10 states
        dx = x + [dt*x(3); ...
                  dt*x(7); ...
                  dt*x(4); ...
                  dt*x(5); ...
                  dt*x(6); ...
                  dt*p(1)*x(1)+dt*p(2)*x(2)+dt*p(3)*x(3)+dt*p(4)*x(4)+dt*p(5)*x(5)+dt*p(6)*x(6)+dt*p(7)*x(7)+dt*p(8)*x(8)+dt*p(9)*x(9)+dt*p(10)*x(10)+dt*p(21); ...
                  dt*x(8); ...
                  dt*x(9); ...
                  dt*x(10); ...
                  dt*p(11)*x(1)+dt*p(12)*x(2)+dt*p(13)*x(3)+dt*p(14)*x(4)+dt*p(15)*x(5)+dt*p(16)*x(6)+dt*p(17)*x(7)+dt*p(18)*x(8)+dt*p(19)*x(9)+dt*p(20)*x(10)+dt*p(22)];
end

