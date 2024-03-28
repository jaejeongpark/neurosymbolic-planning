%% Entry file to use TIRA for computing posteriors
% Load abstract states and controller partitions, return posteriors

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


%% Load abstract states and controller partitions
load_file = '../main/data/cells.mat';
save_file = '../main/data/posteriors.mat';

%load_file = '../main/cav_backup/scalability_linear/data_n10/cells.mat'

load(load_file)

save_data = true;

%% Set system specific parameters
global system_choice
global dt
global v 
global alpha

system_choice = system_dict.sys_choice;
dt = system_dict.dt;
t_init = 0;        

switch system_choice
    case 1 
        %% Unicycle (u as control input)
        v = system_dict.v;
        
    case 2
        %% Unicycle (K,b as control input)
        % state [x, y, theta]
        % input [K1, K2, K3, b]
        v  = system_dict.v;
        
     case 3
        %% 4-state vehicle 
        alpha = system_dict.alpha;
end

%disp(dt);
%disp(v);


%% Call of the main over-approximation function
num_states = size(states_low, 1);
n_x = size(states_low(1,:), 2);

num_partitions = size(partitions_low, 1);
n_p = size(partitions_low(1,:), 2);

%posteriors_low = [];
%posteriors_up  = [];
%tic
%for i = 1:num_states
%    if mod(i, 50) == 0
%        fprintf('\nComputing state %d\n', i)
%    end
%    x_low = states_low(i,:);
%    x_low = reshape(x_low, [n_x, 1]);
%    x_up  = states_up(i,:);
%    x_up  = reshape(x_up, [n_x, 1]);
%    for j = 1:num_partitions
%        %fprintf('\nComputing partition %d', j)
%        p_low = partitions_low(j,:);
%        p_low = reshape(p_low, [n_p, 1]);
%        p_up  = partitions_up(j,:);
%        p_up = reshape(p_up, [n_p, 1]);
%        
%        [succ_low, succ_up] = TIRA(t_init, x_low, x_up, p_low, p_up);
%        
%        succ_low = reshape(succ_low, [1, n_x]);
%        succ_up  = reshape(succ_up, [1, n_x]);
%        posteriors_low = [posteriors_low; succ_low];
%        posteriors_up  = [posteriors_up; succ_up];
%    end
%end
%toc

% Save data in batches to avoid the data containers grow too big.
% TODO: Not save in batch, but combine batches to one array before save. 
% TODO: Save data in multiple batches is not tested, only timed. 
num_batch = 1;
num_succ = num_states * num_partitions;
batch_size = num_states/num_batch;
fprintf('\nnum_batch: %d\n', num_batch);
fprintf('\nbatch_size: %d\n', batch_size);
fprintf('\nnum_succ: %d\n', num_succ);

tic       
for batch=0:num_batch-1  
    fprintf('\nStart batch %d\n', batch);
    posteriors_low = zeros(batch_size: n_x);
    posteriors_up  = zeros(batch_size: n_x);
   
    count = 0;
    for i = batch_size*batch+1:batch_size*(batch+1)
        if mod(i, 50) == 0
            fprintf('\nComputing state %d\n', i)
        end
        x_low = states_low(i,:);
        x_low = reshape(x_low, [n_x, 1]);
        x_up  = states_up(i,:);
        x_up  = reshape(x_up, [n_x, 1]);
        for j = 1:num_partitions
            %fprintf('\nComputing partition %d', j)
            p_low = partitions_low(j,:);
            p_low = reshape(p_low, [n_p, 1]);
            p_up  = partitions_up(j,:);
            p_up = reshape(p_up, [n_p, 1]);
        
            [succ_low, succ_up] = TIRA(t_init, x_low, x_up, p_low, p_up);
            succ_low = reshape(succ_low, [1, n_x]);
            succ_up  = reshape(succ_up, [1, n_x]);
            count = count+1;
        
            posteriors_low(count,:) = succ_low;
            posteriors_up(count,:)  = succ_up;
        end
    end
    %disp(posteriors_low)
    %disp(posteriors_up)
    disp(size(posteriors_low))
    disp(size(posteriors_up))
    if save_data
        save(save_file, 'posteriors_low', 'posteriors_up')
        fprintf('\nData saved\n')
    end  
end
toc






