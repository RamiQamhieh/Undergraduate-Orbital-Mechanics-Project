function [transfer] = planetJump(p1,p2,t0,theta_arrive,pf)

%% Brief
%       Analayzes the potential transfer orbit between two planets 'p1' and
%     'p2' at a given departure time based on the location of the arrival
%     planet when the transfer orbit and it's own orbit intersect. See
%     outputs section for explicit analyisis criteria.

%       Has two modes, based on a variable in the struct 'pf' called 'yn'.
%       If set to y, then the dv's calculated from PlanetJump.m will
%       consider the relative velocity of the spacecraft orbit around the
%       heliocentric velocity of the planet, as well as the dv needed to
%       escape this orbit. If pf.yn is set to 'n', will evaluta the
%       transfer as if the spacecraft starts at the planets heliocentric
%       velocity with no sphere of influence to escape. Spacecraft's
%       initial orbit assumed to be circular.

%     * Note * To avoid issues with imaginary results from the application
%     of LambertSolverND.m and unrealistic orbits, the transfer orbit is
%     only evaluated if the difference in Mean anomoly for the arrival
%     planet 'p2' from the start of the transfer to the end of the transfer
%     is greater than 60 degrees. This should ensure that all inputs for
%     time of transfer are reasonable.

%% Inputs

% p1: id of planet to depart
% p2: id of planet to arrive at
% t0: Time of departure (days)
% Theta_arrive: The true anomoly of the second planet at which you desire
% to arrive at
% pf:   data struct regarding the orbit of the spacecraft around the planet
%       'p1'
%     - yn : 'y' will consider the data in pf, 'n' will run simulation
%     ignoring planetary orbit
%     - rmag : radius of spacecraft orbit
%     - t0 : start of orbit (days)
%     - th0 : true anomoly of spacecraft @ t0

%% Outputs
%    transfer: Struct with the following analysis 


%    transfer.time : Time of transfer (days)
%    transfer.wait_time : N/A to this function, but needed for similarity
%    transfer.net_time : Same as above
% 
%    transfer.d_Th :  true anomoly of planet 1 at the beggining and end of
%                       tranfer, radians
%    transfer.dv : 1x3 vector of dv's. if yn is 'n', one dv and two
%                   zeros. If yn is 'y', both dv's followed by the sum
%    transfer.choice : 'retrograde' or 'prograde' Lambert's solution
%                       choosen
%    transfer.v_leave : Heliocentric velocity leaving p1
%    transfer.v_arrive : Heliocentric velocity arriving at p2
%    transfer.region : for yn = 'y', detrmines if spacecraft is in good
%                       position to exit sphere of influence
%    transfer.cost : N/A for this function, here for similarity
%    transfer.leave_date : day of departure
%    transfer.arrive_date : day of arrival
%    transfer.evaluate : whether or not transfer was evaluated, based on
%                       the criteria of the difference in mean anomoly of
%                       planet p2 throughout the transfer

global uSun
global planet_mu
% t0 in dayyyyyyssss?? sure

% p1 and p2 are respective planet i.d.'s
% t0 in days
% v_esc relative to velocity of planet 1


% time
date = datevec(t0); % So t0 in days then?

[r, v, jd, coe, coe2] = PlanetData2(p1, date(1), date(2), date(3), date(4),date(5),date(6), uSun);
planet1.r=r ; planet1.v=v; 
planet1.coe=[coe,coe2]; planet1.mu = planet_mu(p1);

%% Calc delta-t transfer

% arrival true anomoly --> Mean anomoly
[~, ~, ~, coe, coe2] = PlanetData2(p2, date(1), date(2), date(3), date(4),date(5),date(6), uSun);
planet2.coe=[coe,coe2];
planet2.start_m = planet2.coe(10);
planet2.rate = sqrt(uSun/planet2.coe(1)^3);
planet2.rate_days = planet2.rate*24*60*60; % days
e = planet2.coe(5);
E = 2*atan( sqrt((1-e)/(1+e)) * tan(theta_arrive/2) );
if(E<0)
    E = E + 2*pi;
end
Me_arrive = E - e*sin(E);
% New true anomoly --> Mean anomoly
% dt  = dMean_anomoly/mean_orbital_rate
dMe = Me_arrive - planet2.coe(10);
if(dMe < 0)
    dMe = dMe + 2*pi;
end
transfer.time = dMe/planet2.rate_days;
time_sec = transfer.time*60*60*24;
arrive_time = t0 + transfer.time;
date = datevec(arrive_time);
[r, v, ~, coe, coe2] = PlanetData2(p2, date(1), date(2), date(3), date(4),date(5),date(6), uSun);

planet2.r=r ;
planet2.v=v; 
transfer.d_me = [dMe,planet2.start_m,coe2(4)];
[r, v, ~, coe, coe2] = PlanetData2(p2, date(1), date(2), date(3), date(4),date(5),date(6), uSun);
transfer.d_Th = [planet1.coe(6),coe(6)];

transfer.Consider_SOI = pf.yn;

transfer.Theta_arrive = theta_arrive;
%% Worth evaluating criteria
if transfer.d_me(1) > (60*pi/180)
    transfer.evaluate = 'y';

%% Calc necessary dv 

% Lambert's solver
[v1_pro,v2_pro] = LambertSolverND(planet1.r,planet2.r,time_sec,uSun,'prograde');
[v1_retro,v2_retro] = LambertSolverND(planet1.r,planet2.r,time_sec,uSun,'retrograde');
% imag(v1_pro)
% if (abs(imag(v1_pro))>0 | abs(imag(v1_retro))>0)
%     v1_pro = [1E30,0,0]
%     v1_retro = [1E30,0,0]
% end
% dv = v_new - v_old    
    % v_old:
    v_old = planet1.v;
    if strcmp(pf.yn,'y')
            pf.n = sqrt(planet1.mu/pf.r^3)*60*60*24;
            pf.th = vpa(mod(1000*2*pi+(t0-pf.t0)*pf.n+pf.th0,2*pi),4);
            Vc = sqrt(planet1.mu/pf.r);
            pf.v = Vc.*[-sin(pf.th);cos(pf.th);0];
            v_old = v_old + pf.v;
    end       
    % v_new:  (hard part)
    v_new = [v1_pro, v1_retro];
    transfer.v_want_he = v_new(:,1);
    dv_pro = norm(v_new(:,1)-v_old);
    dv_retro = norm(v_new(:,2)-v_old);
    if strcmp(pf.yn,'y')
        % Change desired v2's to pf frame 
        % Call funtion ToBeNamed to find what the v1 would be needed to
        % match the direction of desired v2 from interplanetary lambert
        % solution(s) (within pf frame from step above).
        % v_new is v1 output. v2 output is v_inf.
        % find magnitude of desired v2 in the pf frame. 
        % redefine to be the v required to exit circular orbit
        v2_pf(:,1) = v_new(:,1) - planet1.v;
        v2_pf_dot(1) = mod(2*pi+atan2(v2_pf(2,1),v2_pf(1,1)),2*pi);
        v_desired_mag_pf(1) = norm(v2_pf(:,1));
        v2_pf(:,2) = v_new(:,2) - planet1.v;      
        v2_pf_dot(2) = mod(2*pi+atan2(v2_pf(2,2),v2_pf(1,2)),2*pi);
        v_desired_mag_pf(2) = norm(v2_pf(:,2));
        [yesno(1),v_esc(1),v_inf(1)]  = LeaveOrbit(v2_pf_dot(1),pf.r,pf.th,planet1.mu);
        [yesno(2),v_esc(2),v_inf(2)] = LeaveOrbit(v2_pf_dot(2),pf.r,pf.th,planet1.mu);
        dv_pro(1) = abs(v_esc(1)-Vc) ;
        dv_pro(2) = abs(v_desired_mag_pf(1)-v_inf(1)) ;
        dv_retro(1) = abs(v_esc(2)-Vc) ;
        dv_retro(2) =  abs(v_desired_mag_pf(2)-v_inf(2)) ;
        dv_pro(3) = sum(dv_pro(1:2));
        dv_retro(3) = sum(dv_retro(1:2));
    end   
if dv_retro<dv_pro

    transfer.choice = 'retro';
    transfer.v_leave = v1_retro;
    transfer.v_arrive = v2_retro;
    if strcmp(pf.yn,'y')
        transfer.dv = vpa(dv_retro,4);
        transfer.region = yesno(2);
    else
        transfer.region = 'N/A';
        transfer.dv = [vpa(dv_retro,4),0,0];
    end
else
    transfer.choice = 'pro';
    transfer.v_leave = v1_pro;
    transfer.v_arrive = v2_pro;
    if strcmp(pf.yn,'y')
        transfer.dv = vpa(dv_pro,4);
        transfer.region = yesno(1);
    else
        transfer.region = 'N/A';
        transfer.dv = [vpa(dv_pro,4),0,0];
    end
end
transfer.cost = 0; % To be used in other function, have to include here so structs have same elements
transfer.wait_time = 0;
transfer.net_time = 0;
transfer.leave_date = round(datevec(t0));
transfer.arrive_date = round(datevec(t0+transfer.net_time));
else    
   % Not to be tried such as to avoid problems with imaginary results in
   % Lambert's solutions
   transfer.time = 0;
   transfer.wait_time = 0;
   transfer.net_time = 0;

   transfer.d_Th = [0 0];
   transfer.dv = [0 0 0];
   transfer.choice = 'N/A';
   transfer.v_leave = [0 0 0];
   transfer.v_arrive = [0 0 0];
   transfer.region = 'N/A';
   transfer.cost = 0;
   transfer.leave_date = [0 0 0 0 0 0];
   transfer.arrive_date = [0 0 0 0 0 0];
   transfer.evaluate = 'n';
end
% Note that this method chooses retrograde vs. prograde solely on the
% initial delta-v value and does not consider the direction or magnitude of
% the velocity upon reaching the sphere of influence of planet two.


% change above s.t. retrograde vs prograde is determined by the cost function 
