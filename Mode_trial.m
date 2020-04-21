%% Part B: Evaluation of Handling Qualities
clear 
clc
%% Initialisation
% Define 2 flight missions
for fm = 1:2
    
    switch fm
        case 1
            % Flight data for martime mission
            run('albatross_parameters_maritime');
            
            V = v_cruise;
            Cd = cdo;
            CL = cl_climb; % beginning of cruise
            m = TOW;
            
            
            
        case 2
            % Flight data for airfield mission
            run('albatross_parameters_airfield');
            
            V = v_cruise;
            Cd = cdo;
            CL = cl_climb; % beginning of cruise
            m = TOW;
            
    end

    % Longitudinal Matrices
    [A_Lon,B_Lon] = LonMatrix_Trim_fold2(q_bar,V,Cd,CL,m);
    % Obtain Eigenvalue and eigenvector
    [Evec ,Eval] = eig(A_Lon) ;
    Eval = diag(Eval) ;
    [wn1,dp1,p1] = damp(A_Lon);
    
    % time constant
    tau1 = 1./(dp1.*wn1);
    % Plot poles to analyse response
%     
%     figure (1)
%     plot(real(Eval),imag(Eval),'*','MarkerSize',12)
%     xlabel('Re')
%     ylabel('Im')
% %     xlim([-3.5 1])
%     legend('maritime', 'airfield')
%     legend boxoff
%     hold on
%     grid on
%     grid minor
%     
  %     % Lateral-directional dynamics
    [A_Lat,B_Lat] = LatMatrix_Trim_fold2(q_bar,V,m);
    
    % Obtain Eigenvalue and eigenvector
    [Evec2 ,Eval2] = eig(A_Lat) ;
    Eval2 = diag(Eval2) ;
    [wn2,dp2,p2] = damp(A_Lat);
    % time constant
    tau2 = 1./(dp2.*wn2);
    % Plot poles to analyse response
    
%     figure (2)
%     plot(real(Eval2),imag(Eval2),'*','MarkerSize',12)
%     xlabel('Re')
%     ylabel('Im')
% %     xlim([-3.5 1])
%  	legend('maritime', 'airfield')
%     legend boxoff
%     
%     hold on
%     grid on
%     grid minor
% %     
%% Begin controllability analysis
    % Begin Euler integration loop
    
    % Define Time Period
    T_max = 100;
    % Step
    dt = 0.01 ;
    % Number of Steps
    N = T_max/dt ;
    
    % Initialise time vector
    T = zeros(1,N);
    % Initialise state vector
    n = length(A_Lon);
    X_Lon = zeros(n,N);
    X_Lat = zeros(n,N);
    % Initialise control vector
    m = length(B_Lon);
    U_Lon = zeros(2,N);
    U_Lat = zeros(2,N);
            %% Analyse dynamics for each control impulse
            for mnvre = 1:3
    
                % Euler integration
                for i = 2:N
                    T(i) = i*dt;
                    [U_Lon(:,i-1),U_Lat(:,i-1)] = Controls(T(i),mnvre);
                    X_Lon_dot = A_Lon*X_Lon(:,i-1) + B_Lon*U_Lon(:,i-1);
                    X_Lon(:,i) = X_Lon(:,i-1) + X_Lon_dot*dt;
    
                    X_Lat_dot = A_Lat*X_Lat(:,i-1) + B_Lat*U_Lat(:,i-1);
                    X_Lat(:,i) = X_Lat(:,i-1) + X_Lat_dot.*dt;
                end
    
                for mode = 1:2
                    if mode == 1
                        xl = [0 , 20] ;
    
                        if mnvre == 1
                            % Elevator
                            xl = [0 , 10];
%                             figure (3)
%                             SPM_cases(X_Lon,X_Lat,T,xl,V)
                            % Aileron
                        elseif mnvre ==  2
%                             figure (4)
%                             Roll_cases(X_Lon,X_Lat,T,xl,V)
                            % Rudder
                        elseif mnvre == 3
%                             figure (5)
%                             Dutch_cases(X_Lon,X_Lat,T,xl,V)
                        end
                    elseif mode == 2
                        xl = [0 , 100] ;
                        if mnvre == 1
                            % Elevator
%                             figure (6)
%                             LPM_cases(X_Lon,X_Lat,T,xl,V)
                            % Aileron
    
                            % Rudder
                        elseif mnvre == 3
                            xl = [0 , 40];
%                             figure (7)
%                             Spiral_cases(X_Lon,X_Lat,T,xl,V)
                        end
                    end
                end
            end
end
