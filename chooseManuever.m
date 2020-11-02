function [save_energy, save_time, save_money, data] = chooseManuever(p1,p2,Theta_arrive,leave_day,k,pf,d_theta,start_day)

%% Brief
%    Uses output from guidemanuever.m 'Theta_arrive' and 'wait_time' and
%    examines a more specific range of these independent variables centered
%    around the results from guideManuever.m. The data struct 'pf' (see
%    below section) is given the string 'y/n' and assigned a value of 'y',
%    meaning the relative velocity of the spacecraft to planet 'p1' as well
%    as the exit of the spacecraft from said planet's sphere of influence
%    are considered.
%    Similar to guidemanuever.m, all combonations of specific time of exit
%    and true anomoly of planet of id 'p2' at the spacecraft's arrival are
%    tested for total "dv", time of transfer, and cost. 
%% Inputs
% pf: contains info for orbit of spacecraft around planet of id 'p1' wihthin the Sphere of Influence 
%     - r
%     - t0 start of orbit
%     - th0 @ start of orbit
%% Starting code
global planet_mu
d_cost = @(dv,time) k*dv+time;

% time range centered at input^ and range based on period of orbit in pf
T = 2*pi*pf.r^(3/2)/sqrt(planet_mu(p1));
T_days = T/60/60/24
[time_range] = double(leave_day + [-1,1]*T_days/2);
% [time_range_ampl] = vpa(time_range,10)
time_arr = linspace(time_range(1),time_range(2),41);
time_arr = time_arr(1:end-1);

% Theta arrive should be range within smallest range from guidemanuever

[Theta_range] = Theta_arrive +[-1,1].*d_theta;
Theta_arr = linspace(Theta_range(1),Theta_range(2),6);
Theta_arr = mod(Theta_arr(1:end-1) + 2*pi, 2*pi);

mindv = 1E30;
mintime = 1E30;
mincost = 1E30;

count = 1; 
count2 = 1;
count3 = 1;
pf.yn = 'y';
pf.n = double(sqrt(planet_mu(p1)/pf.r^3)*60*60*24);
for i = time_arr 
       
        pf.th = double(mod(2*pi*1000+(i-pf.t0)*pf.n+pf.th0,2*pi));
        count2 = 1;
    for j = Theta_arr
        trial = planetJump(p1,p2,i,j,pf);
        str = trial.evaluate;
        if strcmp(str,'y')
            unfiltered_data = planetJump(p1,p2,i,j,pf);
                %Regions(count,count2) = 'n';
            if strcmp(unfiltered_data.region,'y')
                Regions(count,count2) = 'y';
                data(count3) = unfiltered_data;
                data(count3).wait_time = i-start_day;
                data(count3).net_time = data(count3).time + i-start_day;
                data(count3).cost = d_cost(data(count3).dv,data(count3).time);
                data(count3).leave_date = i;
                if data(count3).dv(3) < mindv
                    mindv = data(count3).dv;
                    Least_dv = count3;
                end
                if data(count3).time < mintime
                    mintime = data(count3).time;
                    Least_time = count3;
                end
                if data(count3).cost < mincost
                    mincost = data(count3).cost;
                    Least_cost = count3;
                end
                count3 = count3+1;
            else
                %'&&&&&&&&&'
            end
        end
        count2 = count2+1;
    end
    count = count+1;
end

save_energy = data(Least_dv);
save_time = data(Least_time);
save_money = data(Least_cost);
