% Run experiments to plot results for Section 5.3

clc
clear
rng(0)

%-- Add paths --
addpath('results')
addpath('codes')

%-- Set parameters --
% Determine what to keep constant in the two methods
constant = 'basis'; % The single vector method produces a basis that has the same dimension as the block method. Extra effort spent on approximating quadratic form. 
%constant = 'quadratic_form'; % The single vector method produces a basis that has a larger dimension than the block method. No extra effort spent on approximating quadratic form.

tic
for matrix = 1
    
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

    filename = append(filename,'_svk',constant);
    
    %-- Plot results --
    svk_error_plotter(filename)
    
end
toc