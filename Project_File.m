%% Declarations
global uSun planet_mu planet_rad
uSun = 1.327E11;
earthmass = 5.9722E24;
planet_mass = earthmass*[.055,.815,1,.107,318];
G = 6.6740831313131E-20;
planet_mu = G*planet_mass;  
planet_rad = [2439 6051 6378 3389 69911];

% Establish time zero
date0 = [2019 4 1 0 0 0]; % April the 4th, 2019
day0 = datenum(date0);

% Other criteria:
k = 50; % Experimentally found to be a useful value for demonstration
d_theta = pi/8;

%% Earth orbit and transfer to Mars

p1 = 3; p2 = 4; % Planet I.D.'s for reference
% Describe orbit:
pf_Ea.r = planet_rad(3)+750;
pf_Ea.t0 = day0;
pf_Ea.th0 = 0; % (Implies that the radius between Earth and the spacecraft 
% lies paralel to the x-axis of the Helio-centric frame and the beggnining
% of the orbit)

% Evauating transfer options using function guidemanuever.m (see foreword):
[easy,fast,cheap,results] = guidemanuever(p1,p2,day0,k,d_theta) ;

% Analysis:
fprintf('Earth To Mars transfer analysis: \n')
fprintf('Minimizing dv: net_dv = %.4f \t time in space = %.1f days  \n',easy.dv(1),easy.net_time)
fprintf('Minimizing time in space: net_dv = %.4f \t time in space = %.1f days  \n',fast.dv(1),fast.net_time)
fprintf('Minimizing cost: net_dv = %.4f \t time in space = %.1f days  \n',cheap.dv(1),cheap.net_time)

% Using chepest route from guidemanuever.m to select the desired specific
% time and theta of arrival (at mars), considering the escape of the
% Earth's sphere of influence

[easy,fast,cheap,~] = chooseManuever(p1,p2,cheap.Theta_arrive,datenum(cheap.leave_date),k,pf_Ea,d_theta,day0);

% The cheapest route has now been chosen 
transfer(1).Name = 'Earth to Mars'
transfer(1).dv = cheap.dv(1:2);
transfer(1).date = cheap.leave_date;
transfer(1).time = cheap.net_time;
v_arrival_mars = cheap.v_arrive;
arrive_day_mars = day0 + cheap.net_time;


%% Time on Mars
% global uSun planet_mu planet_rad
% uSun = 1.327E11;
% earthmass = 5.9722E24;
% planet_mass = earthmass*[.055,.815,1,.107,318];
% G = 6.6740831313131E-20;
% planet_mu = G*planet_mass;  
% planet_rad = [2439 6051 6378 3389 69911];
% uSun = 1.327E11;
% earthmass = 5.9722E24;
% planet_mass = earthmass*[.055,.815,1,.107,318];
% G = 6.6740831313131E-20;
% planet_mu = vpa(G*planet_mass,8) ;
% planet_rad = [2439 6051 6378 3389 69911];
% k = 50; % Experimentally found to be a useful value for demonstration
% d_theta = pi/8;
% arrive_day_mars = 737516;
% Need aiming radius
	% The mars parking orbit radius was chosen to be large such as to minimize
% the amount of energy needed to eventually escape the orbit. Additionally,
% the outputs of guidemanuever.m are selected based on lambert's results
% from the velocity of thefirst planet. When we call this function to leave
% mars, the smaller reltive velocity of a larger parking radius will
% provide more accurate guidings for the next interplanetary jump
    % Choosing parking altitude of 1500 km:
p1 = 4;
pf_Ma.r = planet_rad(p1)+1500;

Mars_Vc = sqrt(planet_mu(p1)/pf_Ma.r^3);

% Start orbit:

t = datevec(arrive_day_mars);
[~, mars_v, ~, ~, ~] = PlanetData2(p1, t(1), t(2), t(3), t(4), t(5), t(6), uSun);
v_arrive_pf = v_arrival_mars - mars_v;
v_pf_dot = mod(2*pi + atan2(v_arrive_pf(2),v_arrive_pf(1)),2*pi);

% Can now use Leave_Orbit.m "in reverse" of its intended application to
% calculate how fast the spacecraft should be entering the sphere of
% influence and then what speed will be left at radius of perigee of the
% hyperpoblic orbit.
theta_mars = mod(v_pf_dot+pi,2*pi)+.001;
 [~,v1,v2] = LeaveOrbit(v_pf_dot,pf_Ma.r,theta_mars,planet_mu(p1));
dv_enter_mars(1) = abs(norm(v_arrive_pf)-v2);
dv_enter_mars(2) = abs(v1 - Mars_Vc);
transfer(2).Name = 'Enter mars'
transfer(2).dv = dv_enter_mars
transfer(2).date = arrive_day_mars;
transfer(2).time = 0;
fprintf('Arrival date Mars  : \t \t ')
datevec(arrive_day_mars)
pf_Ma.t0 = arrive_day_mars;
pf_Ma.th0 = theta_mars;

% Transfer to Jupiter
% This process can simulated identically to the exit from Earth to mars.
% First, need to let two weeks pass by to collect datat around Mars.
p2 =5;
leave_day_mars = arrive_day_mars + 14;
fprintf('Leaving Mars on \n')
datevec(leave_day_mars)
clear cheap;
[easy,fast,cheap,~] = guidemanuever(p1,p2,leave_day_mars,k,d_theta);

[easy,fast,cheap,~] = chooseManuever(p1,p2,cheap.Theta_arrive,datenum(cheap.leave_date),k,pf_Ma,d_theta,leave_day_mars);
transfer(3).Name = 'Mars to Jupiter'
transfer(3).dv = cheap.dv(1:2);
transfer(3).date = cheap.leave_date;
transfer(3).time = cheap.net_time

cheap.time
v_arrival_Jup = cheap.v_arrive
arrive_day_Jup = leave_day_mars + cheap.net_time;
fprintf('Arriving at Jup:\n')
datevec(arrive_day_Jup)

%% Flyby @ Jupiter

% First, will employ first rwo of results of guidemanuever.m, which employ
% a zero wait-time. This will tell us where we want our velocity to be
% headed after the flyby such as to be in good position to hit venus.
p1 = 5;
p2 = 2;
% For flyby, can run guidemanuever to find a good heliocentric velocity to
% head towards Venus, ignoring any results with a wait time.
[~,~,~,results] = guidemanuever(p1,p2,arrive_day_Jup,k,d_theta);

results = results(:,1) % This will analyze results of zero wait time, applicable to a flyby
mindv = 1E30;
for i = 1:length(results)
    DV(i) = results(i).dv(1);
    cost(i) = results(i).dv(1);
    time(i) = results(i).time;
    if DV(i) < mindv
        mindv = DV(i);
        minimum = i;
    end
end


v_leave_Jupiter = results(i).v_leave;
v_arrive_Venus = results(i).v_arrive;
Jup_to_Venus = results(i).net_time;
t = datevec(arrive_day_Jup);
[~, Jup_v, ~, ~, ~] = PlanetData2(p1, t(1), t(2), t(3), t(4), t(5), t(6), uSun);
v_arrive_pf = v_arrival_Jup - Jup_v;
v_leave_pf = v_leave_Jupiter - Jup_v;
mags = [norm(v_arrive_pf),norm(v_leave_pf)];
dots = [atan2(v_arrive_pf(2),v_arrive_pf(1)),atan2(v_leave_pf(2),v_leave_pf(1))];


% Finding appropriate radius of perigee for desired flyby:
    % Approaching with unchanged arrival velocity magnitude and direction.
    % Will design turn angle 'Beta' such that the direction of the leaving
    % velocity is the direction of the chose lambert's solution to Venus.
    % This Beta value can then be used to find the desired eccentricity of
    % the planetary flyby. Using the magnitude of the arrival velocity, we
    % can find establish how the radius of perigee will determine the
    % velocity at perigee from [v_inf^2 = sqrt(Vrp^2 - 2*u/r)].
    % Now can use the orbit equation, subbing (rp*Vrp) for the angular
    % mometnum h, and solve for the radius of perigee of our desired flyby.
v1_opp = dots(1)+pi; % Have already seen the direction of 'v_inf 1' and 
% know ahead of time that I don't need to account for the singularity of
% crossing theta = 2*pi
B = (v1_opp-dots(2))/2; % Beta is half the angular difference of the two assymptotes
e = sec(B);
% Storing info
syms rp_sym;
Vrp_sym = sqrt(mags(1)^2 - 2*planet_mu(5)/rp_sym);
h_sym = Vrp_sym * rp_sym;
rp = solve(rp_sym == h_sym^2/planet_mu(5)/(1+e),rp_sym);
rp = rp(rp>0);

aim = @(rp,u,v_inf) sqrt(rp^2+2*u*rp/(norm(v_inf))^2);
Flyby_Jup.date = arrive_day_Jup;
Flyby_Jup.v1_pf = v_arrive_pf;
Flyby_Jup.v2_pf = v_leave_pf;
Flyby_Jup.rp = rp;
Flyby_Jup.B = B;
Flyby_Jup.e = e;
Flyby_Jup.aim = aim(rp,planet_mu(5),v_arrive_pf);
Flyby_Jup.dv = abs(norm(v_leave_pf)-norm(v_arrive_pf));
Flyby_Jup.pv = Jup_v

transfer(4).Name = "Leaving Jupiter Flyby towards Venus";
transfer(4).dv = Flyby_Jup.dv;
transfer(4).date = arrive_day_Jup;
transfer(4).time = 0;

%% Flyby @ Venus

p1 = 2;
p2 = 3;
arrive_day_Venus = arrive_day_Jup+Jup_to_Venus;
fprintf('arriving on Venus \n')
datevec(arrive_day_Venus)
[~,~,~,results] = guidemanuever(p1,p2,arrive_day_Venus,k,d_theta);

results = results(:,1); % This will analyze results of zero wait time, applicable to a flyby
mindv = 1E30;
for i = 1:length(results)
    DV(i) = results(i).dv(1);
    cost(i) = results(i).dv(1);
    time(i) = results(i).time;
    if DV(i) < mindv
        mindv = DV(i);
        minimum = i;
    end
end

v_leave_Venus = results(i).v_leave;
v_arrive_Earth = results(i).v_arrive;
Venus_to_Earth = results(i).net_time;
t = datevec(arrive_day_Venus);
[~, Venus_v, ~, ~, ~] = PlanetData2(p1, t(1), t(2), t(3), t(4), t(5), t(6), uSun);
v_arrive_pf = v_arrive_Venus - Venus_v;
v_leave_pf = v_leave_Venus - Venus_v;
mags = [norm(v_arrive_pf),norm(v_leave_pf)];
dots = [atan2(v_arrive_pf(2),v_arrive_pf(1)),atan2(v_leave_pf(2),v_leave_pf(1))];

B = (v1_opp-dots(2))/2; % Beta is half the angular difference of the two assymptotes
e = -sec(B);
% Storing info
syms rp_sym;
Vrp_sym = sqrt(mags(1)^2 - 2*planet_mu(p1)/rp_sym);
h_sym = Vrp_sym * rp_sym;
rp = solve(rp_sym == h_sym^2/planet_mu(p1)/(1+e),rp_sym);
rp = rp(rp>0)

Flyby_Venus.date = arrive_day_Venus;
Flyby_Venus.v1_pf = v_arrive_pf;
Flyby_Venus.v2_pf = v_leave_pf;
Flyby_Venus.rp = rp;
Flyby_Venus.B = B;
Flyby_Venus.e = e;
Flyby_Venus.aim = aim(rp,planet_mu(p1),v_arrive_pf);
Flyby_Venus.dv = abs(norm(v_leave_pf)-norm(v_arrive_pf)) % Will leave 
% sphere of influence with same magnitude that it left. Need dv to match solution from laberts to reach next planet

transfer(5).Name = "Leaving Venus Flyby towards Earth";
transfer(5).dv = Flyby_Venus.dv;
transfer(5).date = arrive_day_Venus;
transfer(5).time = 0;

%% Back to Earth
% Similar to process of landing at Mars. Back to 750 km Parking
p1 = 3;
park = pf_Ea.r;

Orbit_Vc = sqrt(planet_mu(p1)/park^3);

% Start orbit:
arrive_day_Earth = arrive_day_Venus+Venus_to_Earth;
t = datevec(arrive_day_Earth);
[~, Earth_v, ~, ~, ~] = PlanetData2(p1, t(1), t(2), t(3), t(4), t(5), t(6), uSun);
v_arrive_pf = v_arrive_Earth - Earth_v;
v_pf_dot = mod(2*pi + atan2(v_arrive_pf(2),v_arrive_pf(1)),2*pi);

% Can now use Leave_Orbit.m "in reverse" of its intended application to
% calculate how fast the spacecraft should be entering the sphere of
% influence and then what speed will be left at radius of perigee of the
% hyperpoblic orbit.
theta_dv = mod(v_pf_dot+pi,2*pi)+.001;
[~,v1,v2] = LeaveOrbit(v_pf_dot,park,theta_dv,planet_mu(p1));
dv_enter_Earth(1) = abs(norm(v_arrive_pf)-v2);
dv_enter_Earth(2) = abs(v1 - Mars_Vc);
transfer(6).dv = dv_enter_Earth;
transfer(6).leavedate = arrive_day_Earth;
transfer(6).Name = 'Return to Earth orbit';
transfer(6).time = 0;
fprintf('Arrival date Earth  : \t \t ')
datevec(arrive_day_Earth)

%% Print Answer

disp('/t/t/t Answers \n\n\n')
for i =1:6
    disp(transfer(i).Name)
    fprintf('Date of dv')
    disp(datevec(transfer(i).date))
    fprintf('dv or series of dv')
    disp(vpa(transfer(i).dv,4))
end
