% Run experiments to plot results for Section 5.2

clc
clear
rng(0)

%-- Add paths --
addpath('results')
addpath('codes')

%-- Set parameters --
% Set the oversampling parameter
p = 0;
p = 5;

tic
for matrix = 1:4
    
    matrix
    
    %-- Select matrix --
    if matrix == 1
       
        filename = 'results/exponential_integrator';
        
        
    elseif matrix == 2
        
        filename = 'results/estrada';
        
    elseif matrix == 3
        
        filename = 'results/quantum_spin';
        
    elseif matrix == 4
        
        filename = 'results/synthetic_log';
        
    end
    
    filename = append(filename,'_p=',num2str(p));
    %-- Plot results --
    error_plotter(filename)
    
end
toc