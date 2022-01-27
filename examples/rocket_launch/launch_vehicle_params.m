[g0,w,Re,rho0,h0,Af,cd,mu,m_srbt,m_srb0,T_srb,t_srbbt,I_srb,m_1tot,m_1p,m_10,T_1,t_1bt,I_1,m_2tot,m_2p,m_20,T_2,t_2bt,I_2] 
g0 = 9.80665;
w = [0; 0; 7.29211585e-5]; % Angular velocity of earth relative to inertial space
Re = 6378145;
rho0 = 1.225;
h0 = 7200;
Af = 4*pi;
cd = 0.5;
mu = 3.986012e14;


srb.m  = 19290;  % SRB Total mass
srb.mp = 17010;  % SRB propellant mass
srb.m0 = m_srbt - m_srbp; % SRB dry mass
srb.T  = 628500; % SRB Thrust 
srb.tb = 75.2;   % SRB Burn time
srb.I  = T_srb/(g0 * m_srbp/t_srbbt); % SRB Specific impulse

e1.m  = 104380;  % Main engine total mass
e1.mp = 95550;   % Main engine propellant mass
e1.m0 = m_1tot - m_1tot; % Main engine dry mass
e1.T  = 1083100; % Main engine thrust
e1.tb = 261;     % Main engine burn time
e1.I  = T_1/(g0 * m_1p/t_1bt); % Main engine specific impulse

e2.m  = 19300;  % Secondary engine total mass
e2.mp = 16820;  % Secondary engine propellant mass
e2.m0 = m_2tot - m_2p; % Secondary engine dry weight
e2.T  = 110094; % Secondary engine thrust
e2.tb = 700;    % Secondary engine burn time
e2.I  = T_2/(g0 * m_1p/t_2bt); % Secondary engine specific impulse