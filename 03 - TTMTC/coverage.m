%%

clear, clc
close all

addpath(genpath('Functions'));

n_days = 1;
deltaTimeAq = 0;                        % used for "locking" the S/C

Canberra_gmat = readmatrix("Canberra.csv");
Canberra_gmat = Canberra_gmat(1:n_days+1,:);
Goldstone_gmat = readmatrix("Goldstone.csv");
Goldstone_gmat = Goldstone_gmat(1:n_days+1,:);
Madrid_gmat = readmatrix("Madrid.csv");
Madrid_gmat = Madrid_gmat(1:n_days+1,:);

% Test delta times
% (date2mjd2000(Canberra_gmat(1,7:12)) - date2mjd2000(Canberra_gmat(1,1:6))) * 24 * 3600
% Canberra_gmat(1, 13)

Canberra_times = zeros(size(Canberra_gmat,1),2);
Goldstone_times = zeros(size(Goldstone_gmat,1),2);
Madrid_times = zeros(size(Madrid_gmat,1),2);

for i=1:size(Canberra_times,1)
    Canberra_times(i,1) = date2mjd2000(Canberra_gmat(i,1:6));
    Canberra_times(i,2) = date2mjd2000(Canberra_gmat(i,7:12));
    Goldstone_times(i,1) = date2mjd2000(Goldstone_gmat(i,1:6));
    Goldstone_times(i,2) = date2mjd2000(Goldstone_gmat(i,7:12));
    Madrid_times(i,1) = date2mjd2000(Madrid_gmat(i,1:6));
    Madrid_times(i,2) = date2mjd2000(Madrid_gmat(i,7:12));
end

timeSpan = linspace(min(min(Canberra_times(1,1),Goldstone_times(1,1)),Madrid_times(1,1)), max(max(Canberra_times(end,2),Goldstone_times(end,2)),Madrid_times(end,2)), 1000);

Canberra_times_tot = onesZeros(Canberra_times, timeSpan, deltaTimeAq);
Goldstone_times_tot = onesZeros(Goldstone_times, timeSpan, deltaTimeAq);
Madrid_times_tot = onesZeros(Madrid_times, timeSpan, deltaTimeAq);

[Canberra_times_tot, Goldstone_times_tot, Madrid_times_tot] = stripVect(Canberra_times_tot, Goldstone_times_tot, Madrid_times_tot);

figure
plot(Canberra_times_tot, LineWidth=2, HandleVisibility="off")
hold on
plot(Goldstone_times_tot, LineWidth=2, HandleVisibility="off")
plot(Madrid_times_tot, LineWidth=2, HandleVisibility="off")

scatter(0,0,nan,"filled",MarkerFaceColor="#0072BD",MarkerEdgeColor="#0072BD")
scatter(0,0,nan,"filled",MarkerFaceColor="#D95319",MarkerEdgeColor="#D95319")
scatter(0,0,nan,"filled",MarkerFaceColor="#EDB120",MarkerEdgeColor="#EDB120")

% legend("Canberra", "Goldstone", "Madrid", Orientation="horizontal",Location="northoutside")
legend("Canberra", "Goldstone", "Madrid")
set(gcf,'position',[600,300,420,140])

% new_vect = Canberra_times_tot + Goldstone_times_tot + Madrid_times_tot;
% figure
% plot(new_vect)

totLen = length(Canberra_times_tot);
set(gca,'XTick',0:(totLen/(3*n_days)):totLen);
set(gca,'XTickLabel',["7:30" "15:30" "23:30"]);
xlabel("Hour UTC", Interpreter='latex', FontSize=10)
set(gca,'YTick', [0 1]);
set(gca,'YTickLabel', [0 1]);
ylabel("Visibility", Interpreter='latex', FontSize=10)

%%

plot2pdf ('TTMTC')

%% check (run first the other section)

close all

Canberra_times_tot1 = onesZeros(Canberra_times, timeSpan, 0);
Canberra_times_tot2 = onesZeros(Canberra_times, timeSpan, 5);

plot(Canberra_times_tot1)
hold on
plot(Canberra_times_tot2)


%% functions

function new_array = onesZeros(array, timeSpan, deltaTimeAq)

new_array = zeros(size(timeSpan));

k = 1;
i = 1;
acquisition = 0;
for f=1:length(timeSpan)
    if timeSpan(i)>array(k,1)
        if acquisition == 0
            i = i + deltaTimeAq;
        end
        acquisition = 1;
    end
    if acquisition
        if timeSpan(i)>array(k,2)
            acquisition = 0;
            k = k + 1;
        end
    end

    if acquisition
        new_array(i) = 1;
    end
    if k == size(array,1)
        break
    end

    i = i + 1;  
end

if deltaTimeAq ~= 0
    i = length(new_array);
    stop = 0;
    for f=length(new_array):-1:1
        if i == 1
            stop = 1;
        end
        if ~stop
            if new_array(i) == 0
                if new_array(i-1) == 1
                    new_array = [new_array(1:i-deltaTimeAq-1) zeros(1,deltaTimeAq) new_array(i:end)];
                    i = i - deltaTimeAq;
                end
            end
            
        end
        i = i - 1;
    end
end

end


function [vect1, vect2, vect3] = stripVect(vect1, vect2, vect3)

if length(vect1) ~= length(vect2) || length(vect1) ~= length(vect3)
    error("Not the same length")
end

initial1 = vect1(end);
initial2 = vect2(end);
initial3 = vect3(end);

for i=length(vect1):-1:1
    if vect1(i) ~= initial1 || vect2(i) ~= initial2 || vect3(i) ~= initial3
        break
    end
end

vect1 = vect1(1:i+1);
vect2 = vect2(1:i+1);
vect3 = vect3(1:i+1);

end
