function [yn,v1,v2] = LeaveOrbit(v2dot,rmag,theta,u)
%% Brief:

%   Designed for leaving a circular parking orbit via a hyperbolic escape
%   path, but can be used in the reverse for arriving at a parking orbit.
%   radius of perigee of escaping orbit matches of parking radius. 
%   ** Note ** : Only evaluates escape if spacecraft is in a
%   suitable position to esacpe. This is done based on creating four
%   quadrants based on the desired escape velocity direction and evaluting
%   if the spacecraft is in the right one.

%% Inputs:

% v2dot: unit vector of the velocity escaping the orbit ( or entering the
%       orbit)
% rmag: The radius of parking orbit
% theta: the true anomoly of the spacecraft in the planetary frame when the
%       escape begins or ends
% u:    gravitational constant of planet

%% Outputs:

% yn:  Whether or not spacecraft is in good region to leave its orbit
% v1:  Velocity of escape orbit at perigee
% v2:  Velocity of escape orbit upon leaving sphere of influence, or
% 'v_inf'
yn = 'n';
% Find desired B from pf.th and v2dot

% Find vmag for such B based on rmag and u This is v1

% Find v_inf for such hyperbolic trajectory. This is v2
theta = vpa(mod(2*pi+theta,2*pi),4);
v2_opp = vpa(mod(v2dot + pi,2*pi),4);
theta_opp = vpa(mod(theta + pi,2*pi),4);
v2n_ccw = vpa(mod(v2_opp + 3/2*pi,2*pi),4);
v2n_cw = vpa(mod(v2_opp+pi/2,2*pi),4);
quad_v2_opp =0;
if(cos(v2_opp)>0 & sin(v2_opp)<0)
    quad_v2_opp = 4;
end

if(theta > v2_opp)
    if(theta < v2n_cw)
        yn = 'y';
    elseif quad_v2_opp == 4
        yn = 'y';
    end
elseif quad_v2_opp == 4
    if theta < v2n_cw
        yn = 'y';
    end
end
if yn == 'y'    
    % Now find B, e, v1, v2
    %'Region true'
    B = (theta_opp - v2dot);
    if B < 0
        B = B + 2*pi;
    end
    e = vpa(sec(B),4);
    % B = acos(1/e) or cos(-1/e) or somethin'
    % 
    % % Save B for analyiss or somethin' if not effective use in running program
    %     % Alwasy better to have more data to analyze :\
    %         % I'm never gonna finish this project 
    syms hsym;
    h = solve(rmag == hsym^2 / u /(1+e), hsym);
    h = real(h);
    h = vpa(h(h>0),3);
    v1 = vpa(u/h*(1+e),3);
    v2 = vpa(sqrt(v1^2 - 2*u/rmag),3);
else
    v1 = 1E30;
    v2 = 1E30; % Hopefully will avoid problems


% Establish regions  from x axis ccw

% Assign each region a case string such as ('no-go','ccw',and'cw')

% if strcmpr(string described above, 'ccw') then assign proper Beta turn
% angle value




% if strcmpr(strng^,'cw') same thing but different result
% if strcmpr(string^,'nogo') make v'v such that dv is ridicculously large
% and won't be selected ? Unless can find better way to not consider this
% ccase, arbitrarily large numbers seems like bad idea when working with
% space units which can get really large numbers. Also cleaner coding
% technique to just analyze the data I want to directly 

    % Good idea: Call this function (ToBeNamed.m) again within itself with
    % a different theta if theta is in nogo region

end
