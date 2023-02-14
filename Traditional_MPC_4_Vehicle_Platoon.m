

% point stabilization + Single shooting
% clear all
% close all
% clc


Ts = 0.1;
A = [0 ,1, 0 ; 0 , 0 , 1; 0 0 ,-1/0.51];
B = [0;0;1/0.51];

[Ad1,Bd1] = c2d(A,B,Ts);


A = [0 ,1, 0 ; 0 , 0 , 1; 0 0 ,-1/0.51];
B = [0;0;1/0.51];

[Ad1,Bd1] = c2d(A,B,Ts);


A2 = [0 ,1, 0 ; 0 , 0 , 1; 0 0 ,-1/0.6];
B2 = [0;0;1/0.6];

[Ad2,Bd2] = c2d(A2,B2,Ts);

A3 = [0 ,1, 0 ; 0 , 0 , 1; 0 0 ,-1/0.55];
B3 = [0;0;1/0.55];

[Ad3,Bd3] = c2d(A3,B3,Ts);

A4 = [0 ,1, 0 ; 0 , 0 , 1; 0 0 ,-1/0.49];
B4 = [0;0;1/0.49];

[Ad4,Bd4] = c2d(A4,B4,Ts);


%%%% Define Platoon Agents
%%%% Define Platoon Agents
% Leader
Leader.Ad = Ad1;
Leader.Bd = Bd1;
Leader.x0 = [200;22.2;0];
Leader.x1 = [0;0;0];
Leader.States=[];

% Follower1 
Follower1.Ad = Ad2;
Follower1.Bd = Bd2;
Follower1.x0 = [195;20;0];
Follower1.x1=[0;0;0];
Follower1.States=[];

% Follower2 
Follower2.Ad = Ad3;
Follower2.Bd = Bd3;
Follower2.x0 = [175;25;0];
Follower2.x1=[0;0;0];
Follower2.States=[];

% Follower3 
Follower3.Ad = Ad4;
Follower3.Bd = Bd4;
Follower3.x0 = [162;21;0];
Follower3.x1=[0;0;0];
Follower3.States=[];


%%%%%%% Number of Agents
Number_of_agent = 4;

N = 20; % prediction horizo


% CasADi v3.4.5
addpath('C:/Users/Alireza/Desktop/Casadi/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
%%% Define States of Agents
x1 = SX.sym('x1'); x_dot1 = SX.sym('x_dot1'); x_ddot1 = SX.sym('x_ddot1');
x2 = SX.sym('x2'); x_dot2 = SX.sym('x_dot2'); x_ddot2 = SX.sym('x_ddot2');
x3 = SX.sym('x3'); x_dot3 = SX.sym('x_dot3'); x_ddot3 = SX.sym('x_ddot3');
x4 = SX.sym('x4'); x_dot4 = SX.sym('x_dot4'); x_ddot4 = SX.sym('x_ddot4');
states = [x1;x_dot1;x_ddot1;x2;x_dot2;x_ddot2;x3;x_dot3;x_ddot3;x4;x_dot4;x_ddot4]; n_states = length(states);

tau1 = SX.sym('tau1');
tau2 = SX.sym('tau2');
tau3 = SX.sym('tau3');
tau4 = SX.sym('tau4');
controls = [tau1;tau2;tau3;tau4]; n_controls = length(controls);

%%% Define Reference 
Ref_states = [Leader.x0];
x__ = Leader.x0;
for i=1:N
    x_ = Leader.Ad*x__ + Leader.Bd*0;
    x__ = x_;
    Ref_states = [Ref_states,x__];
end
Temp = Ref_states;
for i=1:Number_of_agent-1  
    Temp2 = Temp - [10*i;0;0];
    Ref_states= [Ref_states;Temp2];
end


%%% Define A AND B for whole of System 
A_main = blkdiag(Leader.Ad,Follower1.Ad,Follower2.Ad,Follower3.Ad);
B_main = blkdiag(Leader.Bd,Follower1.Bd,Follower2.Bd,Follower3.Bd);


%%% Define Parameters for Solver 
U = SX.sym('U',n_controls,N); 
P = SX.sym('P',n_states + N*(n_states)+ Number_of_agent*2 - 2);
X = SX.sym('X',n_states,(N+1));

% compute solution symbolically
X(:,1) = P(1:n_states); % initial state
for k = 1:N
    st = X(:,k);  con = U(:,k);
    st_next  = A_main*st+ B_main*con;
    X(:,k+1) = st_next;
end

ff=Function('ff',{U,P},{X});


obj = 0; 
g = [];  

Q_ = [10,0,0;0,4,0;0,0,1];
R_ =  1;
Q = blkdiag(Q_,Q_,Q_,Q_);
R = blkdiag(R_,R_,R_,R_);

D = [10;10;10];
L = [10;20;30];

st  = X(:,1); 
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P((Number_of_agent*3*k)+1:((Number_of_agent*3)*(k+1))))'*Q*(st-P(((Number_of_agent*3)*k)+1:((Number_of_agent*3)*(k+1)))) + (con)'* R *(con) ; % calculate obj
    for i = 0:(Number_of_agent-2)
        obj = obj + 80*(st((3*i)+1)-st((3*(i+1)+1))-P(n_states + N*(n_states)+1+i))^2;
    end
    for i = 0:(Number_of_agent-2)
        obj = obj + (st((3*i)+2)-st((3*(i+1)+2)))^2;
    end
    for i = 0:(Number_of_agent-2)
        obj = obj + (st((3*i)+3)-st((3*(i+1)+3)))^2;
    end
    for i= 1:(Number_of_agent-1)
        obj = obj + (st(1)-st((3*i)+1)-P(n_states + N*(n_states)+1 + Number_of_agent-2+i))^2;
    end
    for i= 1:(Number_of_agent-1)
        obj = obj + (st(2)-st((3*i)+2))^2;
    end
    for i= 1:(Number_of_agent-1)
        obj = obj + (st(3)-st((3*i)+3))^2;
    end 
end



% compute constraints
for k = 1:N+1   % box constraints due to the map margins
    g = [g ; X(2,k)];   
    g = [g ; X(3,k)]; 
    g = [g ; X(5,k)];   
    g = [g ; X(6,k)]; 
    g = [g ; X(8,k)];   
    g = [g ; X(9,k)]; 
    g = [g ; X(11,k)];   
    g = [g ; X(12,k)]; 
end



% make the decision variables one column vector
OPT_variables = reshape(U,N*Number_of_agent,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;


Temp3 = [0;-10];
args.lbg = repmat(Temp3,(Number_of_agent)*(N+1),1);
Temp3 = [27.7;10];
args.ubg = repmat(Temp3,(Number_of_agent)*(N+1),1);

Temp3 = [-15;-15;-15;-15];
Temper_A = repmat(Temp3,N,1);
args.lbx(1:Number_of_agent*N) = Temper_A;
Temp3 = [10;10;10;10];
Temper_A = repmat(Temp3,N,1);
args.ubx(1:Number_of_agent*N) = Temper_A;

t0 = 0;
Leadr.x0 = [ 200 ; 22.2 ; 0];
x0 = [Leader.x0;Follower1.x0;Follower2.x0;Follower3.x0];    % initial condition.


xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,Number_of_agent);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 10; % Maximum simulation time

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];
main_loop = tic;
U_final_ =[];
delta_u = [];
u_back = 0 ;

while(mpciter < sim_tim / Ts) % new - condition for ending the loop
    % args.p   = [x0;xs]; % set the values of the parameters vector
    %----------------------------------------------------------------------
    args.p(1:(Number_of_agent*3)) = x0; % initial condition of the robot posture
    for k = 1:N %new - set the reference to track
        args.p(((Number_of_agent*3)*k)+1:((Number_of_agent*3)*(k+1))) = Ref_states(:,k);   
    end
   % D = [8;8;8];
   % L = [ 8 ; 16 ; 24];
    args.p((Number_of_agent*3)*(N+1)+1:((Number_of_agent*3)*(N+1)+1)+((Number_of_agent-1)*2 -1))=[D',L'];
    %----------------------------------------------------------------------    
    % initial value of the optimization variables
    args.x0  = [reshape(u0',(Number_of_agent)*N,1)];
    
    
             v1_lim = x0(2);
             v2_lim = x0(5);
             v3_lim = x0(8);
             v4_lim = x0(11);
             
             a1_lim = x0(3);
             a2_lim = x0(6);
             a3_lim = x0(9);
             a4_lim = x0(12);
             
             lim_u_1 = (0.5*1*2.2*0.35*v1_lim^2 + 150 + 0.5 *1*2.2*0.35*v1_lim*a1_lim)/1500;
             lim_u_2 = (0.5*1*2.2*0.35*v2_lim^2 + 150 + 0.5 *1*2.2*0.35*v2_lim*a2_lim)/1500;
             lim_u_3 = (0.5*1*2.2*0.35*v3_lim^2 + 150 + 0.5 *1*2.2*0.35*v3_lim*a3_lim)/1500;
             lim_u_4 = (0.5*1*2.2*0.35*v4_lim^2 + 150 + 0.5 *1*2.2*0.35*v4_lim*a4_lim)/1500;
             Temp33 = [ -10-lim_u_1 ; -10-lim_u_2; -10-lim_u_3; -10-lim_u_4];
             Temp34 = [5-lim_u_1;5-lim_u_2;5-lim_u_3;5-lim_u_4];
    
    
    
             Temper_A = repmat(Temp33,N,1);
             args.lbx(1:Number_of_agent*N) = Temper_A;
             Temper_A = repmat(Temp34,N,1);
             args.ubx(1:Number_of_agent*N) = Temper_A;

    
    
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x)',Number_of_agent,N); % get controls only from the solution
    
    t(mpciter+1) = t0;
    
    % Apply the control and shift the solution
    
    t0 = t0 + Ts;
    U_final = u(:,1);
    
    
    u_cl= [u_cl, U_final];
    X___ = A_main*x0 + B_main*U_final;
    u0 = [u(2:size(u,1),:);u(size(u,1),:)];
    x0 = full(X___);
    %%%%%%%
    
    
    %%%%%%%
    
    
    xx(:,mpciter+2) = x0;
    mpciter
    mpciter = mpciter + 1;
    
    Ref_states = [Leader.x0];
    x__ = Leader.x0;
    for i=1:N
        x_ = Leader.Ad*x__ + Leader.Bd*0;
        if i ==1
            Leader.x0 = x_;
        end
        x__ = x_;
        Ref_states = [Ref_states,x__];
    end
    Temp = Ref_states;
    for i=1:Number_of_agent-1 
        Temp2 = Temp - [10*i;0;0];
        Ref_states= [Ref_states;Temp2];
    end
    
end
main_loop_time = toc(main_loop);
average_mpc_time = main_loop_time/(mpciter-1)


Time = 0:0.1:sim_tim;
figure 
p1 = plot(Time,xx(1,:));
hold on 
p2 = plot(Time,xx(4,:));
hold on 
p3 = plot(Time,xx(7,:));
hold on 
p4 = plot(Time,xx(10,:));
legend('Leader Position','Follower1 Position','Follower2 Position','Follower3 Position');
title('position')
ylabel('Position (m)');
xlabel('Time (S)');
title('Platoon of 4 Vehicle');
p1.LineWidth = 2;
p2.LineWidth = 2;
p3.LineWidth = 2;
p4.LineWidth = 2;
set(gca,'FontSize',12,'FontName','helvetica');
grid on

figure 
p1 = plot(Time,xx(1,:)-xx(4,:));
hold on 
p2 = plot(Time,xx(4,:)-xx(7,:));
hold on 
p3 = plot(Time,xx(7,:)-xx(10,:));
hold on 
legend('Distance Leader and Foolower1 ','Distance Follower1 and Follower2','Distance Follower2 and Follower3');

ylabel('Distance (m)');
xlabel('Time (S)');
title('Platoon of 4 Vehicle');
p1.LineWidth = 2;
p2.LineWidth = 2;
p3.LineWidth = 2;
set(gca,'FontSize',12,'FontName','helvetica');
grid on



figure 
p1 = plot(Time,xx(2,:));
hold on 
p2 = plot(Time,xx(5,:));
hold on 
p3 = plot(Time,xx(8,:));
hold on 
p4 = plot(Time,xx(11,:));
hold on 
legend('Leader Velocity','Follower1 Velocity','Follower2 Velocity','Follower3 Velocity');
ylabel('Velocity (m/S)');
xlabel('Time (S)');
title('Platoon of 4 Vehicle');
p1.LineWidth = 2;
p2.LineWidth = 2;
p3.LineWidth = 2;
p4.LineWidth = 2;
set(gca,'FontSize',12,'FontName','helvetica');
grid on

figure 
p1 = plot(Time,xx(3,:));
hold on 
p2 = plot(Time,xx(6,:));
hold on 
p3= plot(Time,xx(9,:));
hold on 
p4 = plot(Time,xx(12,:));
hold on 
legend('Leader Acceleration','Follower1 Acceleration','Follower2 Acceleration','Follower3 Acceleration');
ylabel('Acceleration (m/S^2)');
xlabel('Time (S)');
title('Platoon of 4 Vehicle');
p1.LineWidth = 2;
p2.LineWidth = 2;
p3.LineWidth = 2;
p4.LineWidth = 2;
set(gca,'FontSize',12,'FontName','helvetica');
grid on


Time = 0:0.1:sim_tim-0.1;


figure 
p1 = plot(Time,u_cl(1,:));
hold on 
p2 = plot(Time,u_cl(2,:));
hold on 
p3 = plot(Time,u_cl(3,:));
hold on 
p4 = plot(Time,u_cl(4,:));
hold on 
legend('Leader Input Signal','Follower1 Input Signal','Follower2 Input Signal','Follower3 Input Signal');
ylabel('u (m/S^2)');
xlabel('Time (S)');
title('Platoon of 4 Vehicle');
p1.LineWidth = 2;
p2.LineWidth = 2;
p3.LineWidth = 2;
p4.LineWidth = 2;
set(gca,'FontSize',12,'FontName','helvetica');
grid on


