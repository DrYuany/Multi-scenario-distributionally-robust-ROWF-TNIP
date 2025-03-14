clear,clc
mpc = (loadcase('case30_20250223'));
% 41 existing lines
% 41*2 = 82 candidate lines
Ae = 1:41; 
Ac = 42:82;
disp(['existing lines: ', num2str(length(Ae))]);
disp(['alternative lines: ', num2str(length(Ac))]);
% --- set hyper parameters of optimizer ---
bigM = 10^6;
MAXLOOP = 10;
option = sdpsettings('solver','gurobi','verbose', 0);

%% -----  initialize -----
%% read data
define_constants;
mpc.bus(mpc.bus(:, PD) < 0, PD) = 0;
mpc = ext2int(mpc);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);

% w_gen(:,1):wind farm bus;
% w_gen(:,2):wind farm capacity
% build 1200MW*3=3600MW wind farm out of 7 candidate wind farm
w_gen(:,1) = [3,12,13,14,15,18,23]';
w_gen(:,2) = 1 *12*ones(length(w_gen(:,1)),1); 

[ref, pv, pq] = bustypes(bus, gen);
%% system info
nb      = size(bus, 1);
ng      = size(gen, 1);
nw      = size(w_gen, 1);
nt      = 24;
ns      = 10; % number of scenarios            
disp(['existing lines: ', num2str(length(Ae))]);
disp(['alternative lines: ', num2str(length(Ac))]);
ne      = length(Ae);
nc      = length(Ac);

%---unit dynamic matrix preparing for unit commitment constraints
UD = zeros(nt - 1, nt);     
for i = 1:nt-1
    UD(i, i) = -1;
    UD(i, i+1) = 1;
end

%% set network parameters
%---build addmitance matrix Bf and connection matrix Cft
%---B = Cft'*Bf, Bf * Va is the vector of real branch powers injected
branch(Ae, BR_STATUS) = 1;
branch(Ac, BR_STATUS) = 0;
[Bf1, CCft] = makeBf(baseMVA, bus, branch);
Bf1 = Bf1(Ae, :);
Cft = CCft';
Cft1 = Cft(:, Ae);
Cft2 = Cft(:, Ac);
%---existing branches
branch(Ae, BR_STATUS) = 0;
branch(Ac, BR_STATUS) = 1;
Bf2 = makeBf(baseMVA, bus, branch);
Bf2 = Bf2(Ac, :);
Fmax1 = 3004*ones(length(Ae),1) / baseMVA;
Fmax1 = Fmax1*1;
Fmin1 = -Fmax1;
Fmax2 = 3004*ones(length(Ac),1) / baseMVA;
Fmax2 = Fmax2*1;
Fmin2 = -Fmax2;
%---bus angle limit
Vamax = 2*ones(nb, 1);
Vamin = -2*ones(nb, 1);
Vamax(ref) = 0;
Vamin(ref) = 0;
%---generation correlation matrix Cg
Cg = sparse(gen(:, GEN_BUS), 1:ng, 1, nb, ng);
Cw = sparse(w_gen(:, 1), 1:nw, 1, nb, nw);
%---generator capacity
% 2	6	13	16	22	23	27	28, the buses of 8 generators 
gen_mat = [2000	6000	2000	2500	5500	1500	3000	1500
0.31	0.29	0.31	0.31	0.3	0.32	0.31	0.32
0.02	0.3	0.02	0.02	0.3	0.01	0.06	0.01
2	5	2	3	4	2	3	2];

Pgmax = gen_mat(1,:)';
Pgmax = Pgmax / baseMVA;
% Pgmin([6:8]) = 0; % 6-8 gas
Pgmin = 0.3*Pgmax; % 6-8 gas
tag_on = [2,5,7];

UT = gen_mat(4,:)';
DT = UT;

% warm start cost,M￥
CU_temp = gen_mat(3,:)';
% generation cost ￥
CP_temp= gen_mat(2,:)';

% number of switchable lines
NV = 0;

%---unit ramp
ramp   = 1;
U_ramp = ramp * Pgmax;
D_ramp = -ramp * Pgmax;

%---load curtailment limit
Rmin = zeros(nb, 1);
Rmax = bus(:, PD) / baseMVA;

%% load,wind, solar curve
load_per = [0.779330007	0.744810331	0.716477974	0.692875903	0.679070542	0.688258921	0.729626972	0.840101379	0.923576146	0.968888194	0.976540808	0.966828235	1	0.99727195	0.99227202	0.989742939	0.958209636	0.920798023	0.913111244	0.941306943	0.956806858	0.947510274	0.917345484	0.82437904];
wind_per = [0.1926	0.2281	0.1652	0.1383	0.1203	0.1258	0.0848	0.0667	0.0427	0.0291	0.0215	0.0191	0.0386	0.0569	0.0569	0.0895	0.148	0.1766	0.1977	0.2466	0.2319	0.2513	0.2391	0.1926];
solar_per = [0,0,0,0,0,0,0.0100000000000000,0.0600000000000000,0.170000000000000,0.290000000000000,0.350000000000000,0.400000000000000,0.420000000000000,0.410000000000000,0.360000000000000,0.280000000000000,0.210000000000000,0.100000000000000,0.0300000000000000,0,0,0,0,0];
load_ratio = 1;
load_basic = load_ratio*bus(:,3) / baseMVA;
for i = 1:ns
    max_pd(:,:,i) = repmat(load_basic,1,nt).*repmat(load_per(1,1:nt),nb,1);
end

% 24000MW generation capacity, 18900MW max pd, 12.5% reserve (without wind power)
pd0 = (210/1.892)*max_pd; 

%% set cost parameters
% 20 typhoon days
Time = (24/nt)*20
% length of candidate lines, km
line_len = [33	25	25	33	30	67	62	32	25	44	52	70	8	42	53	16	20	25	40	20	25	30	30	20	50	35	30	30	16	58	30	120	32	12	46	72	20	25	33	26	37]';
% 700w/km,30years, 换算亿元
CI = [line_len]*70*1e-4;
CV = 0.0001*ones(ne,1);

CP   = CP_temp * 1000 * baseMVA * 365 * 1e-8;
CP = repmat(CP,1,ns);
% CP(:,1:9) = 0;

% M￥换算亿元
CU   = CU_temp * 1e-2;
CU = repmat(CU,1,ns);
CU(:,1:ns-1) = CU(:,1:ns-1) * (Time/(ns-1)); 
CU(:,ns) = CU(:,ns)*(365-Time);
CD   = CU;

% load-curtailment cost, CR1 for important load, CR2 for normal load
CR   = 1 * 1000 * baseMVA * 365 * 1e-8 * ones(nb, 1); %
CR = repmat(CR,1,ns);
CR1   = 10*CR; %10￥/kWh
CR2   = 1*CR; % 1￥/kWh, 5 times of the commercial load price

% CDR for demand response cost, which is ignored in this paper
CDR = 2.5 * 1000 * baseMVA * 365 * 1e-8 * ones(nb, 1);
CDR = repmat(CDR,1,ns);

% wind curtailment cost
CWR   = 0.4 * 1000 * baseMVA * 365 * 1e-8 * ones(nw, 1); 
CWR = repmat(CWR,1,ns);

% wind farm investment cost,7500￥/kW，20 years depreciation，0.08*7500=601￥/kW，
CW = 600 * 1000 * baseMVA * 1e-8 * ones(nw, 1);
% LCOE of offshore wind power 0.33￥/kWh
CW1 = 0.33 * 1000 * baseMVA  * 365 * 1e-8;

%% define variables
% primal variables
m = MAXLOOP;
x   = binvar(nc, 1);               % decision of lines
v   = binvar(ne, 1);               % decision of lines
uw  = binvar(nw, 1);
w   = binvar(ne, ns);               % status of existing lines
y   = binvar(nc, 1);               % status of candidate lines
u   = binvar(ng, nt, ns);               % status of generators
u_st   = binvar(ng, nt, ns);               % status of generators
u_sd   = binvar(ng, nt, ns);               % status of generators

va   = sdpvar(nb, nt, ns, m);         
pg  = sdpvar(ng, nt, ns, m);  
pwr  = sdpvar(nw, nt, ns, m);  
pr   = sdpvar(nb, nt, ns, m);
pr1   = sdpvar(nb, nt, ns, m); % load-shedding for important load
pr2   = sdpvar(nb, nt, ns, m); 
dr1   = sdpvar(nb, nt, ns, m); % day-ahead demand response 
dr2   = sdpvar(nb, nt, ns, m); % in-day demand response, dr1,dr2 are not considered 
pf1  = sdpvar(ne, nt, ns, m);  
pf2  = sdpvar(nc, nt, ns, m);  

% dual variables
u1  = sdpvar(nb, nt, ns);
u2  = sdpvar(ne, nt, ns);
u3  = sdpvar(ne, nt, ns);
u4  = sdpvar(ne, nt, ns);
u5  = sdpvar(ne, nt, ns);
u6  = sdpvar(ne, nt, ns);
u7  = sdpvar(ne, nt, ns);
u8  = sdpvar(nc, nt, ns);
u9  = sdpvar(nc, nt, ns);
u10 = sdpvar(nc, nt, ns);
u11 = sdpvar(nc, nt, ns);
u12 = sdpvar(ng, nt, ns);
u13 = sdpvar(ng, nt, ns);
u14  = sdpvar(nb, nt, ns);
u15  = sdpvar(nb, nt, ns);
u16  = sdpvar(nb, nt, ns);
u17  = sdpvar(nb, nt, ns);
u18  = sdpvar(nb, nt, ns);
u19  = sdpvar(nb, nt, ns);
u20  = sdpvar(nb, nt, ns);
u21  = sdpvar(ng, nt-1, ns);
u22  = sdpvar(ng, nt-1, ns);
u23  = sdpvar(nw, nt, ns);
u24  = sdpvar(nw, nt, ns);
u25  = sdpvar(nb, nt, ns);
u26  = sdpvar(nb, nt, ns);

% auxiliary variables
zw_u = binvar(nw,nt,ns);
zw_d = binvar(nw,nt,ns);
bt_zw_u = sdpvar(nb,nt,ns);
bt_zw_d = sdpvar(nb,nt,ns);
bt_zw_u_24 = sdpvar(nw,nt,ns);
bt_zw_d_24 = sdpvar(nw,nt,ns);

alfa = sdpvar(1, 1);

%% fault distributiona probability setting
p_s  = sdpvar(ns,1); % the actual probability
p_k  = sdpvar(ns,1); % auxiliary variable

% fault_tag1: the 9 fault scenarios
% fault_tag=[22,22,23,22,22,21,23,8,8,19,21,21,22,24,26;23,24,24,23,23,23,24,9,9,21,23,24,19,19,16;24,30,30,30,24,24,27,19,10,35,25,27,10,16,10;32,32,32,32,30,25,29,21,40,36,26,29,8,35,36];
% fault_tag1 = fault_tag(:,[2,3,7,8,9,10,11,12,13]);
fault_tag1 = [22	22	7	7	6	15	8	21	22
30	23	18	9	10	21	9	23	19
32	24	23	15	11	22	10	25	10
41	26	24	23	41	27	40	26	8];
% fault_tag1 = [...
% 32	23	7	9	6	21	8	21	22
% 30	23	7	9	6	21	8	21	19
% 32	24	18	15	11	27	10	25	19
% 32	24	18	15	11	27	10	25	22];
% % fault_tag1([1:2],:)=[];

% A: the fault lines in all fault scenarios
A= unique(fault_tag1(:));

% p_0 the theratical probability distribution
% p_0 =  [0.18;0.17;0.16;0.1;0.07;0.06;0.11;0.08;0.07];
p_0 = [0.16 	0.17 	0.15 	0.07 	0.10 	0.04 	0.08 	0.11 	0.12]';
p_0 = [p_0*Time./365;1-Time/365];
% w0:the fault status of each line and initialization of w0 
w0 = ones(ne,ns);
    for jjj = 1:ns-1
        w0(fault_tag1(:,jjj),jjj) = 0;
    end
w0(:,ns) = 1;


yita = 0.9; % the total deviation of 
theta1 = 0.8; % theta1 theta 2 the deviation of fault probability distribution 
theta2 = 1.2;
k1 = 0.2; % k1,k2 the percentage of improtant and normal load
k2 = 0.8;
kr_1 = 0.0; % kr_1,_2 are demand response parameters
kr_2 = 0.0; % 

%% wind power CVaR-uncertainty set
load('pw_30.mat'); % the wind speed and output ration of wind farms at each bus 
load('v_w_30.mat');
% % the CVAR computation model
% v_in = 3; % cut in cut out rated wind speed of turbine
% v_rated = 12;
% v_out = 25;
% pro_w = 0;
% com_pro = 0;
% lam_norm = 0.03;
% gama = 0.19;% 
% for vv = 1:size(v_w,1)
%     for vvv = 1:nt
%         if v_w(vv,vvv) <= v_rated
%             lam_w(vv,vvv) = lam_norm;
%         elseif v_w(vv,vvv) >= v_rated && v_w(vv,vvv) <= v_out
%             lam_w(vv,vvv) = lam_norm*exp(-gama*(v_w(vv,vvv)-v_rated))^-1;
%         else
%            lam_w(vv,vvv) = lam_norm*exp(-gama*(v_out-v_rated))^-1; 
%         end
%     end
% end
% w_rate = lam_w(:,1);
% 
% for www = 1:nw
% for i = 1:120
%     pro_w(www,i) = nchoosek(120,i)*power(w_rate(www),i).*power(1-w_rate(www),720/6-i);
% end
% end
% com_pro = zeros(nw,120);
% for www = 1:720/6
%   com_pro(:,www) = sum(pro_w(:,1:www),2);
% end

N_normal = zeros(nw,1);
% N_99,90,80: turbine fault number of 7 wind farms with confidential level 99,90,80
N_99 = [20,50,48,55,55,55,24]'; 
N_90 = [16,45,42,49,49,49,19]';
N_80 = [14,42,39,46,46,46,17]';
w_fault_per = (120-[N_normal,N_99,N_90,N_80])/120; % percenatge of fault tubines
wind_per_normal = repmat(wind_per,nw,1);
pw_normal = repmat(w_gen(:,2),1,nt).*wind_per_normal; % output under normal scenario
% tag_N select confidential level, tag_beishu controls the forecasting error
% benchmark error is 20% under typhoon and 15% under normal weather
tag_N = 2;
tag_beishu = 0.2;
for i = 1:ns - 1
    w_u(:,:,i) = ((1+tag_beishu)/1.2)*repmat(w_fault_per(1:nw,tag_N),1,nt).*repmat(w_gen(:,2),1,nt).*pw_up_revise(1:nw,1:nt);
    w_d(:,:,i) = ((1-tag_beishu)/0.8)*repmat(w_fault_per(1:nw,tag_N),1,nt).*repmat(w_gen(:,2),1,nt).*pw_down_revise(1:nw,1:nt);
 end
    w_u(:,:,ns) = (1+tag_beishu-0.05)*pw_normal;
    w_d(:,:,ns) = (1-tag_beishu+0.05)*pw_normal;

% the theoretical wind farm power is the average of upper and lower bound under typhoon
z_aver = 0.5*(w_u+w_d);
w_u = w_u - z_aver;
w_d = z_aver - w_d;

[w_a,w_b] = find((pw_up_revise(:,24)-pw_down_revise(:,24))==0);
w_equal = 0;
w_equal = [w_a,w_b];

TAO = 50; % the budget of wind power uncertainty

%% initialization of the first iteration
idx_x = [];
z0=ones(ng,1);
y0=ones(nc,1);
x0 = zeros(nc, MAXLOOP+1);
uw0 = zeros(nw, MAXLOOP+1);
v0 = ones(ne,MAXLOOP+1);
u0 = ones(ng,nt,ns,MAXLOOP+1);
u0(1:end,:,:,1) = 0;
u_st0 = ones(ng,nt,ns,MAXLOOP+1);
u_sd0 = zeros(ng,nt,ns,MAXLOOP+1);
dr1_0(:,:,:,1) = zeros(nb,nt,ns);
cost_CU_0(1) = 0;
LB=0;
UB=inf;
%% computation start
tic
master_cons = [];
master_new_con = [];
master_new_con1 = [];
for i = 1:MAXLOOP
    master_new_con = [];
    dual_cons = [];
    dual_new_con = [];
    dual_new_con1 = [];
    dual_new_con2 = [];
    dual_new_con3 = [];
    dual_new_con4 = [];
    dual_new_con5 = [];
    dual_new_con6 = [];
    dual_new_con7 = [];
    dual_obj = 0;
    dual_obj1 = 0;
    dual_obj2 = 0;
   %% ----------------------------------- Hrdening module
    % the hardening is implemented before the start of planning to improve the theoretical fault probability distribution    
    % NH select how much lines to be hardened 
    % Nstra select which harden strategy, total 3 strategies
    NH = 15;
    Nstra = 1;
    cost_celue_NH = 0; % total cost of hardening
    if i==1       
        [cost_celue_NH1, p_s_0_2] =  harden_NEW_30(fault_tag1, NH, Nstra, p_0, i,Time);
         p_0 = p_s_0_2;       
    end
    cost_celue_NH = cost_celue_NH1;
    %----------------------------------- solve slave problem    
    for k = 1:ns
        for j = 1:nt      
            dual_obj = - (pd0(:,j)-dr1_0(:,j,k,i))'*u1(:,j,k)...
               + (Cw*(uw0(:,i).*w_u(:,j,k)))'*bt_zw_u(:,j,k) - (Cw*(uw0(:,i).*w_d(:,j,k)))'*bt_zw_d(:,j,k) + (Cw*(uw0(:,i).*z_aver(:,j,k)))'*u1(:,j,k)...
               + (-bigM*(1 - v0(:,i)))'*u2(:,j,k) -(bigM*(1 - v0(:,i)))'*u3(:,j,k)...
               + (-bigM*v0(:,i))'*u4(:,j,k) - (bigM*v0(:,i))'*u5(:,j,k)...
               + Fmin1'*u6(:,j,k) - Fmax1'*u7(:,j,k)...
               + (-bigM*(1-x0(:,i)))'*u8(:,j,k) - (bigM*(1-x0(:,i)))'*u9(:,j,k)...
               + (diag(x0(:,i))*Fmin2)'*u10(:,j,k) - (diag(x0(:,i))*Fmax2)'*u11(:,j,k) ... 
               + (u0(:,j,k,i).*Pgmin)'*u12(:,j,k)- (u0(:,j,k,i).*Pgmax)'*u13(:,j,k)...
               - (k1*(pd0(:,j)-dr1_0(:,j,k,i)))'*u15(:,j,k) - (k2*(pd0(:,j)-dr1_0(:,j,k,i)))'*u17(:,j,k)...
               + Vamin'*u19(:,j,k) - Vamax'*u20(:,j,k)...
               - (uw0(:,i).*w_u(:,j,k))'*bt_zw_u_24(:,j,k) + (uw0(:,i).*w_d(:,j,k))'*bt_zw_d_24(:,j,k) - (uw0(:,i).*z_aver(:,j,k))'*u24(:,j,k)... 
               - (kr_2*pd0(:,j))'*u26(:,j,k);
            
            dual_new_con = [   
              -Cft1'*u1(:,j,k) + u3(:,j,k) - u2(:,j,k) + u5(:,j,k) - u4(:,j,k) + u7(:,j,k) - u6(:,j,k) == 0;                          % pf1
              -Cft2'*u1(:,j,k) + u9(:,j,k) - u8(:,j,k) + u11(:,j,k) - u10(:,j,k) == 0;                          % pf2
              (diag(w0(:,k))*Bf1)'*u2(:,j,k) - (diag(w0(:,k))*Bf1)'*u3(:,j,k) + Bf2'*u8(:,j,k) - Bf2'*u9(:,j,k) + u20(:,j,k) - u19(:,j,k) == 0;    % va
       
               - u1(:,j,k) - u18(:,j,k) == 0;                % pr

               Cw'*u1(:,j,k) + u23(:,j,k) - u24(:,j,k) == p_s(k)*CWR(:,k);                % pwr
               u14(:,j,k) - u15(:,j,k) + u18(:,j,k) == p_s(k)*CR1(:,k);                % pr1
               u16(:,j,k) - u17(:,j,k) + u18(:,j,k) == p_s(k)*CR2(:,k);                % pr2  
               -Cg'*u1(:,j,k) + u12(:,j,k) - u13(:,j,k) == p_s(k)*CP(:,k);  % pg  
               -u1(:,j,k) + k1*(u15(:,j,k) - u14(:,j,k)) + k2*(u17(:,j,k) - u16(:,j,k)) + u26(:,j,k) - u25(:,j,k) <= p_s(k)*CDR(:,k);  % DR2
     
                - bigM*(1-Cw*zw_u(:,j,k)) <= bt_zw_u(:,j,k) - u1(:,j,k) <= bigM*(1-Cw*zw_u(:,j,k));
                - bigM*Cw*zw_u(:,j,k) <= bt_zw_u(:,j,k) <= bigM*Cw*zw_u(:,j,k);
                - bigM*(1-Cw*zw_d(:,j,k)) <= bt_zw_d(:,j,k) - u1(:,j,k) <= bigM*(1-Cw*zw_d(:,j,k));
                - bigM*Cw*zw_d(:,j,k) <= bt_zw_d(:,j,k) <= bigM*Cw*zw_d(:,j,k);
        
                - bigM*(1-zw_u(:,j,k)) <= bt_zw_u_24(:,j,k) - u24(:,j,k) <= bigM*(1-zw_u(:,j,k));
                - bigM*zw_u(:,j,k) <= bt_zw_u_24(:,j,k) <= bigM*zw_u(:,j,k);
                - bigM*(1-zw_d(:,j,k)) <= bt_zw_d_24(:,j,k) - u24(:,j,k) <= bigM*(1-zw_d(:,j,k));
                - bigM*zw_d(:,j,k) <= bt_zw_d_24(:,j,k) <= bigM*zw_d(:,j,k);
        
                u2(:,j,k) >= 0; u3(:,j,k) >= 0;u4(:,j,k) >= 0; u5(:,j,k) >= 0; u6(:,j,k) >= 0; u7(:,j,k) >= 0; u8(:,j,k) >= 0;
                u9(:,j,k) >= 0; u10(:,j,k) >= 0; u11(:,j,k) >= 0; u12(:,j,k) >= 0; u13(:,j,k) >= 0; 
                u14(:,j,k) >= 0; u15(:,j,k) >= 0; u16(:,j,k) >= 0; u17(:,j,k) >= 0;
                u19(:,j,k) >= 0; u20(:,j,k) >= 0; u23(:,j,k) >= 0; u24(:,j,k) >= 0;
                u25(:,j,k) >= 0; u26(:,j,k) >= 0;
                ];
        dual_obj1 = dual_obj1 + dual_obj;
        dual_cons = [dual_cons,dual_new_con];
        end
    end
    
    % the fault probaility distribution uncertainty set
    dual_new_con4 = [p_k(1:ns-1) >= p_s(1:ns-1) - p_0(1:ns-1);
    p_k(1:ns-1) >= p_0(1:ns-1) - p_s(1:ns-1);
    sum(p_s(1:ns-1)) == sum(p_0(1:ns-1));
    sum(p_k(1:ns-1)) <= yita*sum(p_0(1:ns-1));
    theta1*p_0(1:ns-1) <= p_s(1:ns-1) <= theta2*p_0(1:ns-1);
    p_s(1:ns-1) >= 0; p_k(1:ns-1) >= 0;
    p_s(ns) == p_0(ns);
    p_k(ns) == 0];    
    dual_cons = [dual_cons,dual_new_con4];  

    % the wind power uncertainty set
    for ss = 1:ns-2 % the wind farm output under 9 possible fault scenarios are the same
        dual_new_con5 = [dual_new_con5,zw_u(:,:,ss) == zw_u(:,:,ss+1);];
        dual_new_con5 = [dual_new_con5,zw_d(:,:,ss) == zw_d(:,:,ss+1);];
    end
    dual_cons = [dual_cons,dual_new_con5]; 

    dual_cons = [dual_cons,zw_u + zw_d == 1;]; 
      
    if i >= 2 % the tag_wind records which wind farm is built
        dual_new_con7 = [...
                % the typhoon effect lies in 12 level typhoon wind speed
               zw_u(1,:,:) == zeros(nt,ns); 
               sum(sum(sum(zw_u(2,1:9,:)) + sum(zw_u(3,1:8,:)) + sum(zw_u(4,1:7,:)) + sum(zw_u(5,1:4,:)) + sum(zw_u(6,1:7,:)) + sum(zw_u(7,1:5,:)))) == TAO;
               sum(sum(sum(sum(zw_u(:,:,:))))) == TAO;
         ];
        dual_cons = [dual_cons,dual_new_con7]; 
    end

    if i == 1 % there is no wind farm built in the first iteration
        dual_new_con7 = [...
             sum(sum(sum(zw_u(:,:,1)+zw_u(:,:,10)))) == TAO;
             sum(sum(sum(sum(zw_u(:,:,:))))) == TAO;
            ];
        dual_cons = [dual_cons,dual_new_con7]; 
     end

    dual_cons = [dual_cons,dual_new_con6];        
    result_dual = optimize(dual_cons, -dual_obj1, option);
    op_sub=value(dual_obj1)
    p_s_0_1(:,i) = value(p_s(:));

   
%% ---------------------------------------------------------------------
    for iii = 1:ns
        cost_dr_temp(iii) = sum(CDR(:,iii)'*dr1_0(:,:,iii,i));
        cost_dr(iii) = cost_dr_temp(iii)*p_s_0_1(iii,i);
    end
    UB = CI'*x0(:,i)+cost_celue_NH+CW'*(uw0(:,i).*w_gen(:,2))+CV'*(1-v0(:,i))+cost_CU_0(i)+sum(cost_dr)+op_sub
    if (UB-LB)<=10^-1
        [idx] = find(x0(:,i)==1);
        AAA = branch(Ae,1:2);AAA(idx,:)
        [iduw] = find(uw0(:,i)==1);
        [iduw]

        p_s_0_1
        APR=value(pr);
        APW = value(pwr);
        AAPR2 = zeros(nb,nt,ns);
        AAPR2 = zeros(nb,nt,ns);
        AAPR = zeros(nt,ns);
        AAPW = zeros(nt,ns);   

        for jj=1:ns
            for kk=1:nt
                AAPR(kk,jj) = sum(APR(:,kk,jj,i-1));
                AAPR_RATE(kk,jj) = AAPR(kk,jj)/sum(pd0(:,kk,jj));
                AAPW(kk,jj) = sum(APW(:,kk,jj,i-1));
                AAPW_RATE(kk,jj) = AAPW(kk,jj)/sum(uw0(:,i).*pw_0(:,kk,jj,i-1));

                end
            end
        total_pr = 0;
        for jj = 1:9
           total_pr = total_pr + AAPR(:,jj) *  p_s_0_1(jj,i-1);
        end
        left_pr = (sum(pd0(:,:,1))' - total_pr)./sum(pd0(:,:,1))';

        curtail_pd_rate = 0;
        curtail_pw_rate = 0;
        curtail_pd_cost = 0;
        curtail_pw_cost = 0;
 
        for jj = 1:ns
            curtail_pd_rate =  curtail_pd_rate + p_s_0_1(jj,i-1)*sum(AAPR(:,jj))/sum(sum(pd0(:,:,jj)));
            curtail_pw_rate =  curtail_pw_rate + p_s_0_1(jj,i-1)*sum(AAPW(:,jj))/sum(sum(uw0(:,i)'*pw_0(:,:,jj,i-1)));
        end
        for jj = 1:ns
%             curtail_pd_cost = curtail_pd_cost + p_s_0_1(jj,i-1)*(CR(1,jj)*sum(AAPR(:,jj)));
            curtail_pw_cost = curtail_pw_cost + p_s_0_1(jj,i-1)*(CWR(1,jj)*sum(AAPW(:,jj)));
        end
        cost_CU=0;
          for k = 1:ns
            cost_CU =  cost_CU + sum(CU(:,k)'*u_st0(:,:,k,i) + CD(:,k)'*u_sd0(:,:,k,i));
        end
   
        cost_pg_TDS = 0;
        pr1_TDS = 0;
        pr2_TDS = 0;
        wind_TDS = 0;
        cost_pg_NOS = 0;
        pr1_NOS = 0;
        pr2_NOS = 0;
        wind_NOS = 0;
        cost_pg_TDS = 0;
        cost_CU_TDS = 0;
        cost_CU_NOS = 0;
        CP   = CP_temp * 1000 * baseMVA * 365 * 1e-8;
        CP = repmat(CP,1,ns);
  
      for ii = 1:ns-1
            cost_pg_TDS = cost_pg_TDS+sum(CP(:,ii)'*pg0(:,:,ii,i))*p_s_0_1(ii,i-1);
            pr1_TDS = pr1_TDS+sum(sum(pr10(:,:,ii,i)))*p_s_0_1(ii,i-1);
            pr2_TDS = pr2_TDS+sum(sum(pr20(:,:,ii,i)))*p_s_0_1(ii,i-1);
            wind_TDS = wind_TDS+sum(sum(pw_0(iduw,:,ii,i-1)-APW(iduw,:,ii,i-1)))*p_s_0_1(ii,i-1);
            cost_CU_TDS =  cost_CU_TDS + sum(CU(:,ii)'*u_st0(:,:,ii,i) + CD(:,ii)'*u_sd0(:,:,ii,i));
      end

      ii = ns;
        cost_pg_NOS = sum(CP(:,ii)'*pg0(:,:,ii,i))*p_s_0_1(ii,i-1);
        pr1_NOS = sum(sum(pr10(:,:,ii,i)))*p_s_0_1(ii,i-1);
        pr2_NOS = sum(sum(pr20(:,:,ii,i)))*p_s_0_1(ii,i-1);
        wind_NOS = sum(sum(pw_0(iduw,:,ii,i-1)-APW(iduw,:,ii,i-1)))*p_s_0_1(ii,i-1);
        cost_CU_NOS =  cost_CU_NOS + sum(CU(:,ii)'*u_st0(:,:,ii,i) + CD(:,ii)'*u_sd0(:,:,ii,i));      

      pr1_NOS_cost =  CR1(1) * pr1_NOS/1;
      pr2_NOS_cost =  CR2(1) * pr2_NOS/1;
      pr1_TDS_cost =  CR1(1) * pr1_TDS/1;
      pr2_TDS_cost =  CR2(1) * pr2_TDS/1;
      pr1_cost = pr1_NOS_cost + pr1_TDS_cost;
      pr2_cost = pr2_NOS_cost + pr2_TDS_cost;
      wind_TDS_cost = CW1*wind_TDS;
      wind_NOS_cost = CW1*wind_NOS;
      wind_cost = wind_NOS_cost + wind_TDS_cost;
      pr1_sum = pr1_NOS + pr1_TDS;
      pr2_sum = pr2_NOS + pr2_TDS;

      curtail_pd_cost  =  pr1_cost + pr2_cost;
      curtail_cost = curtail_pd_cost + curtail_pw_cost;   
      operation_cost = value(alfa) - curtail_cost; 
 

     hard_r = [NH,Nstra]
     wind_r = [tag_N,tag_beishu]
     Pg_r = Pgmax(1)/Pgmin(1)
     CR_r = [k1,k2]
     CR_r2 = [CR1(1)/CR(1),CR2(1)/CR(1)]
     Time
     TAO
     MAXLOOP
     [idx]'  
    disp('规划成本、加固成本，总运行费用, 总切负荷费用、总弃风费用、切负荷率、弃风率')
    disp( [CI'*x0(:,i), cost_celue_NH, operation_cost+cost_CU,curtail_pd_cost,curtail_pw_cost,curtail_pd_rate,curtail_pw_rate])
    disp('火电费用、TDS发电成本、NOS发电成本，重要切负荷成本, 一般切符合成本') 
    disp([operation_cost+cost_CU,cost_CU_TDS + cost_pg_TDS, cost_CU_NOS + cost_pg_NOS, pr1_cost,pr2_cost])
            break
    end
        operate_cost = 0;
        PG_cost =  0;
        PR1_cost = 0;
        PR2_cost = 0;
        PWR_cost = 0; 
        operate_cost1 = 0;
        PG_cost1 =  0;
        PR1_cost1 = 0;
        PR2_cost1 = 0;
        PWR_cost1 = 0; 

        zw_u_0(:,:,:,i) = value(zw_u(:,:,:));
        zw_d_0(:,:,:,i) = value(zw_d(:,:,:));
        
        p_s_0(:,i) = p_s_0_1(:,i);
        pw_0(:,:,:,i) = z_aver(:,:,:) + zw_u_0(:,:,:,i).*w_u - zw_d_0(:,:,:,i).*w_d;
        pw_0 = round(pw_0*100)/100; % this might sovle the large matrix coefificient range problem 
    %% solve master problem
%     master_cons = [];
    for k = 1:ns
        for j = 1:nt
        master_new_con1 = [
        Cg*pg(:,j,k,i) + Cw*(uw.*pw_0(:,j,k,i) - pwr(:,j,k,i)) + pr(:,j,k,i) - Cft1*pf1(:,j,k,i) - Cft2*pf2(:,j,k,i) == pd0(:,j) - dr1(:,j,k,i) - dr2(:,j,k,i); 
        -bigM*(1-v(:,1)) <= pf1(:,j,k,i) - diag(w0(:,k))*Bf1*va(:,j,k,i);
        pf1(:,j,k,i) - diag(w0(:,k))*Bf1*va(:,j,k,i) <= bigM*(1-v(:,1));
        -bigM*v(:,1) <= pf1(:,j,k,i) <=  bigM*v(:,1);  
        Fmin1 <= pf1(:,j,k,i) <= Fmax1;      
        -bigM*(1-x) <= pf2(:,j,k,i) - diag(y0)*Bf2*va(:,j,k,i);
        pf2(:,j,k,i) - diag(y0)*Bf2*va(:,j,k,i) <= bigM*(1-x);
        diag(x)*Fmin2 <= pf2(:,j,k,i) <= diag(x)*Fmax2;
        u(:,j,k).*Pgmin <= pg(:,j,k,i) <= u(:,j,k).*Pgmax; % u12 u13
        0 <= pr1(:,j,k,i) <= k1*(pd0(:,j) - dr1(:,j,k,i) - dr2(:,j,k,i)); % u14~u15
        0 <= pr2(:,j,k,i) <= k2*(pd0(:,j) - dr1(:,j,k,i) - dr2(:,j,k,i)); % u16~u17
        pr(:,j,k,i) == pr1(:,j,k,i) + pr2(:,j,k,i); % u18
        Vamin <= va(:,j,k,i) <= Vamax; % u19~u20
        0 <= pwr(:,j,k,i) <= uw.*pw_0(:,j,k,i); % u23~u24
        0 <= dr1(:,j,k,i) <= kr_1*pd0(:,j);
        0 <= dr2(:,j,k,i) <= kr_2*pd0(:,j);% u25~u26 
        ];
        master_new_con = [master_new_con,master_new_con1];
        operate_cost = operate_cost + CP(:,k)'*pg(:,j,k,i) + CR1(:,k)'*pr1(:,j,k,i) + CR2(:,k)'*pr2(:,j,k,i) +  CWR(:,k)'*pwr(:,j,k,i) + CDR(:,k)'*(dr1(:,j,k,i)+dr2(:,j,k,i));
        PG_cost =  PG_cost + CP(:,k)'*pg(:,j,k,i);
        PR1_cost = PR1_cost + CR1(:,k)'*pr1(:,j,k,i);
        PR2_cost = PR2_cost + CR2(:,k)'*pr2(:,j,k,i);
        PWR_cost = PWR_cost +  CWR(:,k)'*pwr(:,j,k,i); 
        end
           
        master_cons = [master_cons, master_new_con];
        master_new_con = [];
        operate_cost1 = operate_cost1 + p_s_0(k,i)*operate_cost;
        PG_cost1 =  PG_cost1 + p_s_0(k,i)*PG_cost;
        PR1_cost1 = PR1_cost1 + p_s_0(k,i)*PR1_cost;
        PR2_cost1 = PR2_cost1 + p_s_0(k,i)*PR2_cost;
        PWR_cost1 = PWR_cost1 + p_s_0(k,i)*PWR_cost; 
        operate_cost = 0;
         PG_cost =  0;
        PR1_cost = 0;
        PR2_cost = 0;
        PWR_cost = 0; 
    end

%% unit commitment and other constriants
 % without wind and fault uncertainty,[16    28    31    35    38] should be bulid
 % only nos,[16    24  28    31    35    38];
        base_x = [16    28    31    35    38];
%           base_x = [16    24  28    31    35    38];
%         base_x = [ 9    10    16    21    23    24    25    28    31    35    38];
        for k = 1:ns   
        for j = 1:nt-1
              master_new_con2 = [
                dr1(:,j+1,k,i) - dr1(:,j,k,i) <= 0.5*kr_1*pd0(:,j);
                dr1 == 0; dr2 == 0;
                sum(x(setdiff(Ae,[A;base_x']))) <= 0.01;
%                 sum(x) == 6;
                 x(base_x) == 1; 
              ];
                master_new_con = [master_new_con,master_new_con2];
        end
        end
    

        master_new_con = [master_new_con,x([idx_x]) == 1;];
        master_new_con = [master_new_con,sum(v) >= ne - NV;];
        master_new_con = [master_new_con, sum(uw) == 3];

    
        for j = 2:nt
            for kk = j:min(UT+j-1,nt)
                master_new_con = [master_new_con,-u(:,j-1,:)+u(:,j,:)-u(:,kk,:)<=0];
            end
            for kk = j:min(DT+j-1,nt)
                master_new_con = [master_new_con,u(:,j-1,:)-u(:,j,:)+u(:,kk,:)<=1;];
            end
        end
        for j = 2:nt
            master_new_con3 = [
                -u(:,j-1,:)+u(:,j,:)-u_st(:,j,:)<=0;
                u(:,j-1,:)-u(:,j,:)-u_sd(:,j,:)<=0;
                ];
            master_new_con = [master_new_con,master_new_con3];
        end
%                 load('u0_no_typhoon_30.mat');
%                 load('u0_only_TDS.mat');
        for k = 1:ns-2
            master_new_con4 = [
                          
                u(tag_on,:,:) == ones(length(tag_on),nt,ns);
%                 u(:,:,k) == u0_no_typhoon_30;
%                 u(:,:,ns) == u0_no_typhoon_30;
                u(:,:,k) == u(:,:,k+1);

                ];
            master_new_con = [master_new_con,master_new_con4];
        end
        cost_CU = 0;     
        for k = 1:ns
            cost_CU = cost_CU + sum(CU(:,k)'*u_st(:,:,k) + CD(:,k)'*u_sd(:,:,k));
        end
        master_cons = [master_cons, master_new_con];
        master_new_con4 = [alfa >= operate_cost1];
        master_cons = [master_cons, master_new_con4];
        master_obj = CI'*x + cost_celue_NH + CW'*(uw.*w_gen(:,2)) + cost_CU + alfa;
        result_master = optimize(master_cons, master_obj, option);
        op_mas=value(master_obj);
        LB = op_mas
        x0(:, i+1) = value(x);
        [idx_x] = find(x0(:,i+1)==1);
    %     if i<=1
        uw0(:,i+1) = value(uw);
    %     else
    %         uw0(:,i+1) = uw0(:,i);
    %     end
        v0(:, i+1) = value(v);
        tag_wind_1 = find(uw0(:,i+1) == 1);
        tag_wind_0 = find(uw0(:,i+1) == 0);
        u0(:, :, :, i+1) = value(u(:,:,:));
        u_st0(:,:,:,i+1) = value(u_st(:,:,:));
        u_sd0(:,:,:,i+1) = value(u_sd(:,:,:)); 
        dr1_0(:,:,:,i+1) = value(dr1(:,:,:,i));
        dr2_0(:,:,:,i+1) = value(dr2(:,:,:,i));  
        pr0(:,:,:,i+1) = value(pr(:,:,:,i));
        pr10(:,:,:,i+1) = value(pr1(:,:,:,i));
        pr20(:,:,:,i+1) = value(pr2(:,:,:,i));
        pg0(:,:,:,i+1) = value(pg(:,:,:,i));
        cost_CU_0(i+1) = value(cost_CU);
        cost(:, i) = [value(CI'*x),value(CW'*(uw.*w_gen(:,2))),value(alfa),value(master_obj)];
        disp(cost(:,i))

     if (UB-LB)<=10^-1 || i == MAXLOOP
 
        [idx] = find(x0(:,i+1)>=0.9);
        AAA = branch(Ae,1:2);AAA(idx,:)
        [iduw] = find(uw0(:,i+1)==1);
        [iduw]'
        p_s_0_1
        APR=value(pr);
        APW = value(pwr);
        AAPR2 = zeros(nb,nt,ns);
        AAPR2 = zeros(nb,nt,ns);
        AAPR = zeros(nt,ns);
        AAPW = zeros(nt,ns);   

        for jj=1:ns
            for kk=1:nt
                AAPR(kk,jj) = sum(APR(:,kk,jj,i));
                AAPR_RATE(kk,jj) = AAPR(kk,jj)/sum(pd0(:,kk,jj));
                AAPW(kk,jj) = sum(APW(:,kk,jj,i));
                AAPW_RATE(kk,jj) = AAPW(kk,jj)/sum(uw0(:,i+1).*pw_0(:,kk,jj,i));

                end
            end
        total_pr = 0;
        for jj = 1:ns-1
           total_pr = total_pr + AAPR(:,jj) *  p_s_0_1(jj,i);
        end
        left_pr = (sum(pd0(:,:,1))' - total_pr)./sum(pd0(:,:,1))';

        curtail_pd_rate = 0;
        curtail_pw_rate = 0;
        curtail_pd_cost = 0;
        curtail_pw_cost = 0;
                cost_CU = 0; 
        for jj = 1:ns
            curtail_pd_rate =  curtail_pd_rate + p_s_0_1(jj,i)*sum(AAPR(:,jj))/sum(sum(pd0(:,:,jj)));
            curtail_pw_rate =  curtail_pw_rate + p_s_0_1(jj,i)*sum(AAPW(:,jj))/sum(sum(uw0(:,i+1)'*pw_0(:,:,jj,i)));
        end
        for jj = 1:ns
%             curtail_pd_cost = curtail_pd_cost + p_s_0_1(jj,i-1)*(CR(1,jj)*sum(AAPR(:,jj)));
            curtail_pw_cost = curtail_pw_cost + p_s_0_1(jj,i)*(CWR(1,jj)*sum(AAPW(:,jj)));
        end
           for k = 1:ns
            cost_CU = cost_CU + sum(CU(:,k)'*u_st0(:,:,k,i+1) + CD(:,k)'*u_sd0(:,:,k,i+1));
        end
   
        cost_pg_TDS = 0;
        pr1_TDS = 0;
        pr2_TDS = 0;
        wind_TDS = 0;
        cost_pg_NOS = 0;
        pr1_NOS = 0;
        pr2_NOS = 0;
        wind_NOS = 0;
        cost_CU_TDS = 0;
        cost_CU_NOS = 0;
        CP   = CP_temp * 1000 * baseMVA * 365 * 1e-8;
        CP = repmat(CP,1,ns);

      for ii = 1:ns-1
        cost_pg_TDS = cost_pg_TDS+sum(CP(:,ii)'*pg0(:,:,ii,i+1))*p_s_0_1(ii,i);
        pr1_TDS = pr1_TDS+sum(sum(pr10(:,:,ii,i+1)))*p_s_0_1(ii,i);
        pr2_TDS = pr2_TDS+sum(sum(pr20(:,:,ii,i+1)))*p_s_0_1(ii,i);
        wind_TDS = wind_TDS+sum(sum(pw_0(iduw,:,ii,i)-APW(iduw,:,ii,i)))*p_s_0_1(ii,i);
        cost_CU_TDS = cost_CU_TDS + sum(CU(:,ii)'*u_st0(:,:,ii,i+1) + CD(:,ii)'*u_sd0(:,:,ii,i+1)); 
      end
   
        ii = ns;
        cost_pg_NOS = sum(CP(:,ii)'*pg0(:,:,ii,i+1))*p_s_0_1(ii,i);
        pr1_NOS = sum(sum(pr10(:,:,ii,i+1)))*p_s_0_1(ii,i);
        pr2_NOS = sum(sum(pr20(:,:,ii,i+1)))*p_s_0_1(ii,i);
        wind_NOS = sum(sum(pw_0(iduw,:,ii,i)-APW(iduw,:,ii,i)))*p_s_0_1(ii,i);
        cost_CU_NOS = cost_CU_NOS + sum(CU(:,ii)'*u_st0(:,:,ii,i+1) + CD(:,ii)'*u_sd0(:,:,ii,i+1)); 
  
      
      pr1_NOS_cost =  CR1(1) * pr1_NOS/1;
      pr2_NOS_cost =  CR2(1) * pr2_NOS/1;
      pr1_TDS_cost =  CR1(1) * pr1_TDS/1;
      pr2_TDS_cost =  CR2(1) * pr2_TDS/1;
      pr1_cost = pr1_NOS_cost + pr1_TDS_cost;
      pr2_cost = pr2_NOS_cost + pr2_TDS_cost;
      wind_TDS_cost = CW1*wind_TDS;
      wind_NOS_cost = CW1*wind_NOS;
      wind_cost = wind_NOS_cost + wind_TDS_cost;
      pr1_sum = pr1_NOS + pr1_TDS;
      pr2_sum = pr2_NOS + pr2_TDS;

      curtail_pd_cost  =  pr1_cost + pr2_cost;
      curtail_cost = curtail_pd_cost + curtail_pw_cost;   
      operation_cost = value(alfa) - curtail_cost; 
     hard_r = [NH,Nstra]
     wind_r = [tag_N,tag_beishu]
     Pg_r = Pgmax(1)/Pgmin(1)
     CR_r = [k1,k2]
     CR_r2 = [CR1(1)/CR(1),CR2(1)/CR(1)]
     TAO
     MAXLOOP
     Time
    disp('规划成本、加固成本，总运行费用, 总切负荷费用、总弃风费用、切负荷率、弃风率')
    disp( [CI'*x0(:,i+1), cost_celue_NH, operation_cost+cost_CU,curtail_pd_cost,curtail_pw_cost,curtail_pd_rate,curtail_pw_rate])
    disp('火电费用、TDS发电成本、NOS发电成本，重要切负荷成本, 一般切符合成本') 
    disp([operation_cost+cost_CU,cost_CU_TDS + cost_pg_TDS, cost_CU_NOS + cost_pg_NOS, pr1_cost,pr2_cost])
    
    [idx]'  
     break
    end
end
toc
% save workspace
% curDate = datestr(now,'mm-dd');
% diy = 'N-2';
% filename = ['workspace_',diy,'_',curDate,'.mat']
% save(filename);

  


