function ekfpredict_vp(u, dt)
% EKF-SLAM prediction for Victoria Park process model

    global Param;
    global State;
    map_num = State.Ekf.nL*2;
    % Get robot pose vector from the robot state vector
    robot_pose = State.Ekf.mu(State.Ekf.iR);
    a = Param.a;
    b = Param.b;
    L = Param.L;
    % Calculate how much we have to add to the currect mean value of robot
    % pose state vector
    delta_mu = [u(1)*cos(robot_pose(3))-(u(1)/L)*tan(u(2))*(a*sin(robot_pose(3))+b*cos(robot_pose(3)));
                u(1)*sin(robot_pose(3))+(u(1)/L)*tan(u(2))*(a*cos(robot_pose(3))-b*sin(robot_pose(3)));
                (u(1)/L)*tan(u(2))].*dt;
    % Predict next mean value of robot pose state vector
    State.Ekf.mu(State.Ekf.iR) = robot_pose + delta_mu;
    State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
    
    % Predict Covariance matrix of the entire robot state vector
    % Calculate Jacobian matrix for robot pose state
    G = eye(3+map_num);
    G(1,3) = -dt*(u(1)*sin(State.Ekf.mu(3))-(u(1)/L)*tan(u(2))*(a*cos(State.Ekf.mu(3))-b*sin(State.Ekf.mu(3))));
    G(2,3) = dt*(u(1)*cos(State.Ekf.mu(3))+(u(1)/L)*tan(u(2))*(-a*sin(State.Ekf.mu(3))-b*cos(State.Ekf.mu(3))));
    
    R = zeros(3+map_num);
    % Calculate Jacobian matrix for control value
    V = [cos(State.Ekf.mu(3))-(1/L)*tan(u(2))*(a*sin(State.Ekf.mu(3))+b*cos(State.Ekf.mu(3))), -(u(1)/L)*sec(u(2)).^2*(a*sin(State.Ekf.mu(3))+b*cos(State.Ekf.mu(3)));
         sin(State.Ekf.mu(3))+(1/L)*tan(u(2))*(a*cos(State.Ekf.mu(3))-b*sin(State.Ekf.mu(3))), (u(1)/L)*sec(u(2)).^2*(a*cos(State.Ekf.mu(3))-b*sin(State.Ekf.mu(3)));
         (1/L)*tan(u(2))                                       , (u(1)/L)*sec(u(2)).^2].*dt;
    R(1:3, 1:3) = V*Param.Qu*V';
    QF = zeros(3+map_num);
    QF(1:3, 1:3) = Param.Qf;
    % Prediect Sigma
    State.Ekf.Sigma = G*State.Ekf.Sigma*G' + G*QF*G' + R;
end
