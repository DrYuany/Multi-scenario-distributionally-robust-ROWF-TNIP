# Multi-scenario-distributionally-robust-ROWF-TNIP

**Illustration of all files:**
(1) Multi_scenario_distributionally_robust_OWF_TNIP.m: the main program in MATLAB version: 
(2) case30_20250223.m: the data of the modified IEEE-30 bus system: 
(3) makeBf.m and harden_NEW_30.m: two self-defined functions used in main program
(4) all .mat files: the data needed in the main program

**Requirements:**
(1) Matpower
(2) Yalmip
(3) GUROBI (If CPLEX, revise the code option = sdpsettings('solver','gurobi','verbose', 0);

**Illustration of important parameters in the code:**
(1) Parameters of uncertiaty set of OWF output: 
    TAO is the uncertianty budget; 
    tag_N select the CVaR number of failed turbines with different confidential level;
    tag_beishu is the forecasting error;
(2) Parameters of uncertiaty set of fault: 
    yita: the total deviation of probability distribution
    theta1, theta2: the deviation bound of fault probability distribution 
(3) Parameters of short-term measures:  
    k1,k2: the percentage of improtant and normal load
    NH: number of llines to be hardened
    Nstra: select which hardening strategy
    




