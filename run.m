function varargout = run(numSteps, pauseLen, makeVideo)
% RUN PS3 EKF Feature-Based SLAM
%   RUN(ARG)
%   RUN(ARG, PAUSELEN)
%   RUN(ARG, PAUSELEN, MAKEVIDEO)
%      ARG - is either the number of time steps, (e.g. 100 is a complete
%            circuit) or a data structure from a previous run.
%      PAUSELEN - set to `inf`, to manually pause, o/w # of seconds to wait
%                 (e.g., 0.1 is the default)
%      MAKEVIDEO - set to 1 to enable make a video, default is 0 (disabled)
%
%   Note: more parameters can be controlled in the run.m file itself via
%   fields of the Param structure.
close all;
addpath('./util');
addpath('./vicpark');
% Check whether we have input of set to default value
if ~exist('numSteps', 'var') || isempty(numSteps)
    numSteps = 500;
end
if ~exist('pauseLen', 'var') || isempty(pauseLen)
    pauseLen = 0.1;
end
if ~exist('makeVideo','var') || isempty(makeVideo)
    makeVideo = false;
end

clear global Param State Data;
global Param;
global State;
% select which data association method to use in ekfupdate.m, choices are:
%   nn    - incremental maximum likelhood nearest neighbor

% select which update method to use in ekfupdate.m, choices are:
%   batch  - batch updates
%   seq    - sequential updates
Param.updateMethod = 'batch';

% size of bounding box for VP data set plotting
Param.bbox = 0; % bbox = 20 [m] speeds up graphics

% Structure of global State variable
%===================================================
State.Ekf.t     = 0;          % time
State.Ekf.mu    = zeros(3,1); % robot initial pose
State.Ekf.Sigma = zeros(3,3); % robot initial covariance
State.Ekf.iR    = 1:3;        % 3 vector containing robot indices
State.Ekf.iM    = [];         % 2*nL vector containing map indices
State.Ekf.iL    = {};         % nL cell array containing indices of landmark i
State.Ekf.sL    = [];         % nL vector containing signatures of landmarks
State.Ekf.nL    = 0;          % scalar number of landmarks
State.Ekf.mu_results = [];    % 3*timestep vector containing the ekf mu result at each time step
%===================================================
% !!!Run vp slam model 
runvp(numSteps,pauseLen, makeVideo);


% Plot results
figure(2)
% Calculate Determinant of Landmark covariance matrix
det_array = zeros(1,State.Ekf.nL);
for i = 1:State.Ekf.nL    
    det_array(i) = det(State.Ekf.Sigma(2*i+2:2*i+3,2*i+2:2*i+3));    
end

plot(1:State.Ekf.nL, det_array, '-x');
xlabel('feature id');
ylabel('Determinant');
% Set measurement of x axis
set(gca, 'XTick', [1:1:State.Ekf.nL]);

figure(3)
% plot correlate matrix
subplot(1,2,1);
correlate_matrix = corrcov(State.Ekf.Sigma(4:2*State.Ekf.nL+3,4:2*State.Ekf.nL+3));
image = (ones(size(correlate_matrix))-correlate_matrix)*255;
imshow(image, [0 255]);
% Plot correlation value of each landmark
subplot(1,2,2);
hold on;
for i = 1:State.Ekf.nL
    cor_array = [];
    for j = 1:State.Ekf.nL
        cor_array = [cor_array, correlate_matrix(2*i-1,2*j-1)];
    end
    plot(1:State.Ekf.nL, cor_array);
    set(gca,'XTick',[1:1:State.Ekf.nL]);
end


