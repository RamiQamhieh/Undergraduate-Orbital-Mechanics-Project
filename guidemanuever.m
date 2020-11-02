function [save_energy,save_time,save_money,data] = guidemanuever(p1,p2,time0,k,d_theta)

%% Brief

% Considers a wide range of departure times and anomolies of arrival for 
% potential interplanetary transfer orbits 

global uSun

%% Inputs 

% p1
% p2
% time0
% k SHOULDNT NEED K CHANGE THIS


%% Outputs

% save_energy
% save_time
% save_money
% data 

%% Declarations

d_cost = @(dv,time) k*dv+time;
%time in days
date = datevec(time0);
[~, ~, ~, coe, coe2] = PlanetData2(p1, date(1), date(2), date(3), date(4),date(5),date(6), uSun);
planet1.coe=[coe,coe2];
[~, ~, ~, coe, coe2] = PlanetData2(p2, date(1), date(2), date(3), date(4),date(5),date(6), uSun);
planet2.coe=[coe,coe2];
Tp1 = 2*pi*planet1.coe(1)^(3/2)/sqrt(uSun);
Tp2 = 2*pi*planet2.coe(1)^(3/2)/sqrt(uSun);
count = 0;
mindv = 1E30;
mintime = 1E30;
mincost = 1E30;

if Tp2 > Tp1
    wait_seg = Tp2/5;
else 
    wait_seg = Tp1/5;
end
wait_seg = wait_seg /60/60/24; % days
%% Creating data struct matrix
pf.yn = 'n';
for i=0:wait_seg:9*wait_seg % Goes from 0 to 9/5 of the larger of the two orbit periods
% 'Running new wait time...' % Make sure code is goin while waiting for function to complete
    start_time_days = time0 + i;
    count = count+1;
    count2 = 0;
    date = datevec(start_time_days);
    for j = 0: d_theta : (2*pi - d_theta)       
            unfiltered_data = planetJump(p1,p2,start_time_days,j,pf);
                   if unfiltered_data.evaluate == 'y'
                        count2 = count2+1;
                        data(count,count2) = unfiltered_data;
                        data(count,count2).wait_time= i; % Account for wait time + transfer time
                        data(count,count2).net_time = data(count,count2).time+i;
                        data(count,count2).cost = d_cost(data(count,count2).dv(1),data(count,count2).time); 
                        if data(count,count2).dv(1)<mindv
                            mindv = data(count,count2).dv(1);
                            Least_dv = [count,count2];
                        end
                         if data(count,count2).net_time<mintime
                            mintime = data(count,count2).time;
                            Least_time = [count,count2];
                         end
                        if data(count,count2).cost<mincost
                            mincost = data(count,count2).cost;
                            Least_expensive = [count,count2];
                        end
                   else
                    end
    end
end
save_energy = data(Least_dv(1),Least_dv(2));
save_time = data(Least_time(1),Least_time(2));
save_money = data(Least_expensive(1),Least_expensive(2));
