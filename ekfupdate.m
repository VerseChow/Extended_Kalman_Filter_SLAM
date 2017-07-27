function count_landmark = ekfupdate(z)
    global Param;
    global State;
    % !!!!!!!!!
    % Data association step
    % Returns state vector indices pairing observations with landmarks Li
    % For value of Li:
    % 0  -- For this observation, we have to initialize a new landmark
    % -1 -- For this observation, it is meanless, we have to eliminate that
    % other number -- landmark index
    Li = da_nn(z(1:2,:)); 
    num = size(z,2);
    count = 0;
    if strcmp(Param.updateMethod,'batch')
        for i=1:num
            if Li(i)==0
                % Initialize new landmark if Li index equals 0
                index = initialize_new_landmark(z(:,i), Param.R);
                Li(i) = index;% update the new landmark index
            end
            if Li(i)==-1
                count = count+1;% update the total number of meanless observation
            end
        end
        count = num-count;% Calculate valid landmark counts
    end
    % For sequential correction method
    if strcmp(Param.updateMethod,'seq') 
        for i = 1:num
            if Li(i) == -1
                continue;
            end
            index = Li(i);
            if Li(i) == 0
                % Initualize new landmark if Li index equals 0
                index = initialize_new_landmark(z(:,i), Param.R);
            end
%%%%%%%%%%%%%Correction Step by Step via Seq method
            [Ht, zt] = getHt(index);
            St = Ht*State.Ekf.Sigma*Ht'+Param.R;
            % Calculate Kalman Gain
            Kt = State.Ekf.Sigma*Ht'/St;
            % Correct robot pose state vector
            State.Ekf.mu = State.Ekf.mu+Kt*( minimizedAngle(z(1:2,i)-zt));
            State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
            temp = eye(size(State.Ekf.Sigma))-Kt*Ht;
            % Correct covariance matrix
            State.Ekf.Sigma = temp*State.Ekf.Sigma*temp'+Kt*Param.R*Kt';
        end
    % For batch correction method    
    elseif strcmp(Param.updateMethod,'batch')
        zt = zeros(2*count,1);
        dz = zt;
        Ht = zeros(2*count,3+2*State.Ekf.nL);
        j = 1;
        R = [];
        for i=1:num
            if Li(i)~=-1
%%%%%%%%%%%%%Stack the observation to batch type
                deltaX = State.Ekf.mu(2+2*Li(i))-State.Ekf.mu(1);
                deltaY = State.Ekf.mu(3+2*Li(i))-State.Ekf.mu(2); 
                q = deltaX^2+deltaY^2;
                zt(2*j-1) = sqrt(q);
                zt(2*j) = atan2(deltaY, deltaX)-State.Ekf.mu(3);
                dz(2*j-1) = z(1,i)-zt(2*j-1);
                dz(2*j) = minimizedAngle(z(2,i)-zt(2*j));
                Ht(2*j-1:2*j,1:3) = [-deltaX/sqrt(q), -deltaY/sqrt(q), 0; deltaY/q, -deltaX/q, -1];
                Ht(2*j-1:2*j,2*Li(i)+2:3+2*Li(i)) = [deltaX/sqrt(q), deltaY/sqrt(q);-deltaY/q, deltaX/q];
                R = blkdiag(R, Param.R);
                j = j+1;
            end
        end
%%%%%%%%%%%%%Correction Step of Batch
        S = Ht*State.Ekf.Sigma*Ht'+R;
        % Calculate Kalman Gain
        Kt = State.Ekf.Sigma*Ht'/S;
        % Correct robot pose state vector
        State.Ekf.mu = State.Ekf.mu+Kt*dz;
        State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
        temp = eye(size(State.Ekf.Sigma))-Kt*Ht;
        % Correct covariance matrix
        State.Ekf.Sigma = temp*State.Ekf.Sigma*temp'+Kt*R*Kt';
    end
    count_landmark = count;
end
%%%%%%%%%%%%%%Construct Ht matrix
function [Ht, zt] = getHt(index)
    global State;
    Ht = zeros(2,3+2*size(State.Ekf.sL,2));
    zt = zeros(2,1);
    deltaX = State.Ekf.mu(2+2*index)-State.Ekf.mu(1);
    deltaY = State.Ekf.mu(3+2*index)-State.Ekf.mu(2);
    q = deltaX^2+deltaY^2;
    Ht(1:2,1:3) = [-deltaX/sqrt(q), -deltaY/sqrt(q), 0; deltaY/q, -deltaX/q, -1];
    Ht(1:2,2*index+2:3+2*index) = [deltaX/sqrt(q), deltaY/sqrt(q);-deltaY/q, deltaX/q];
    zt(1) = sqrt(q);
    zt(2) = minimizedAngle(atan2(deltaY, deltaX)-State.Ekf.mu(3));
end
