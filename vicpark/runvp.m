function runvp(nSteps,pauseLen,makeVideo)

global Param;
global State;
global Data;

if makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end
Data = load_vp_si();

% Initalize Params
%===================================================
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = zeros(3);
State.Ekf.mu_results = [State.Ekf.mu_results, State.Ekf.mu];
global AAr;
AAr = [0:360]*pi/360;


figure(1); clf;
axis equal;

t_correction = zeros(1, min(nSteps, length(Data.Laser.time)));
num_landmark = zeros(1, min(nSteps, length(Data.Laser.time)));
t_predict = zeros(1,min(nSteps, length(Data.Laser.time)));
ci = 1; % control index
t = min(Data.Laser.time(1), Data.Control.time(1));
for k=1:min(nSteps, length(Data.Laser.time))
    tic;
    while (Data.Control.time(ci) < Data.Laser.time(k))
       % Control availabe
       % Get control value
       dt = Data.Control.time(ci) - t;
       t = Data.Control.time(ci);
       u = [Data.Control.ve(ci), Data.Control.alpha(ci)]';
       
       % !!!
       % prediction step for predicting motion in next time step
       ekfpredict_vp(u, dt);
       ci = ci+1;
    end
    t_predict(k) = toc;
    % observation available
    % Get observation value
    t = Data.Laser.time(k);
    z = detectTreesI16(Data.Laser.ranges(k,:));
    z_ = z-repmat([0;pi/2;0],1,size(z,2));
    Param.maxObs = size(z,2);
    tic;
    
    % !!!
    % Correction step for correcting prediction step with observation
    num_landmark(k) = ekfupdate(z_);
    State.Ekf.mu_results = [State.Ekf.mu_results, State.Ekf.mu(State.Ekf.iR)];
    State.Ekf.t = State.Ekf.t+1;
    t_correction(k) = toc;
    doGraphics(z);
    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end
    if makeVideo
        F = getframe(gcf);
        switch votype
          case 'avifile'
            vo = addframe(vo, F);
          case 'VideoWriter'
            writeVideo(vo, F);
          otherwise
            error('unrecognized votype');
        end
    end
end
if makeVideo
    fprintf('Writing video...');
    switch votype
      case 'avifile'
        vo = close(vo);
      case 'VideoWriter'
        close(vo);
      otherwise
        error('unrecognized votype');
    end
    fprintf('done\n');
end
figure(2)
% Plot prediction cpu time
subplot(3,1,1);
plot(1:min(nSteps, length(Data.Laser.time)), t_predict);
title('Prediction CPU time');
xlabel('iteration step');
ylabel('CPU time');
% Plot correction cpu time
subplot(3,1,2);
plot(1:min(nSteps, length(Data.Laser.time)), t_correction);
title('Correction CPU time');
xlabel('iteration step');
ylabel('CPU time');
subplot(3,1,3);
plot(1:min(nSteps, length(Data.Laser.time)), num_landmark);
title('number of landmarks at each iteration');
xlabel('iteration steps');
ylabel('number of landmarks');

%==========================================================================
function doGraphics(z)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!
global State;
% plot the robot and 3-sigma covariance ellipsoid
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on;

plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:3, 1:3), 'm', 0, 0, 0);

% restrict view to a bounding box around the current pose
% BB=50;
% axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);
% plot GPS path and Ekf reuslts path
h1 = plot( State.Ekf.mu_results(1,:), State.Ekf.mu_results(2,:), 'r' );
%plot( Data.Gps.x(2:2+State.Ekf.t), Data.Gps.y(2:2+State.Ekf.t), 'g');


% project raw sensor detections in global frame using estimate pose
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr-pi/2);
    yl = yr + r*sin(b+tr-pi/2);
    h3 = plot([xr; xl], [yr; yl],'g',xl,yl,'b*');
end
for i=1:State.Ekf.nL
    
    h2 = plotcov2d( State.Ekf.mu(2*i+2), State.Ekf.mu(2*i+3), State.Ekf.Sigma(2*i+2:2*i+3,2*i+2:2*i+3), 'm', 0);

end
legend([h1, h2, h3(2)], 'Ekf Correction', 'Landmark', 'Data Association');
title('EKF SLAM VISUALIZATION');
xlabel('x');
ylabel('y');
hold off;

