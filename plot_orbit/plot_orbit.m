%% plot orbit

clear, clc
close all

format longG
addpath('time')

[t, r_soho, v_soho] = extract_horizons('./horizon/soho_2005.txt');
[~, r_sun, ~] = extract_horizons('./horizon/sun_2005.txt');
[~, r_moon, ~] = extract_horizons('./horizon/moon_2005.txt');

% return

% tspan = t0 + dt / (24 * 3600);
% plot(tspan, vecnorm(v_soho));
% figure
% plot(tspan, v_soho(1, :));
% figure
% plot(tspan, v_soho(2, :));
% figure
% plot(tspan, v_soho(3, :));
% return

fig = figure;
xlabel('\textbf{Earth-Sun dir. (x) [km]}', Interpreter='latex', FontSize=10)
ylabel('\textbf{y [km]}', Interpreter='latex', FontSize=10)
zlabel('\textbf{h vect. (z) [km]}', Interpreter='latex', FontSize=10)
grid on, hold on, axis equal
view(320, 10)

xlim([0, 1800000]);
ylim([-1000000, 1000000]);
zlim([-200000, 120000]);

plot3(1.5e6, 0, 0, '.', MarkerSize=20)
plot3(0, 0, 0, '.', MarkerSize=20)

soho_rot = zeros(3, size(r_soho, 2)-2);
moon_rot = soho_rot;
% an = animatedline(0, 0, 0);
for i=2:size(r_soho, 2)-1
	x = r_sun(:, i);
	x = x / norm(x);
	diff_sol = r_sun(:, i+1) - r_sun(:, i-1);
	diff_sol = diff_sol / norm(diff_sol);
	z = cross(x, diff_sol);
	z = z / norm(z);
	y = cross(z, x);

	% plot3([0, z(1)], [0, z(2)], [0, z(3)]);
	% grid on
	% drawnow

	M = [x, y, z]';
	soho_rot(:, i-1) = M * r_soho(:, i);
	moon_rot(:, i-1) = M * r_moon(:, i);

	% t = t0 + dt(i) / (24 * 3600);
	% date = jd2date(t);

	% datestr = strcat(num2str(date(1)),'-',num2str(date(2)),'-',num2str(date(3)));

	% addpoints(an, soho_rot(1, i), soho_rot(2, i), soho_rot(3, i))
	% title(datestr)
	% drawnow
end
plot3(soho_rot(1, :), soho_rot(2, :), soho_rot(3, :))
% plot3(moon_rot(1, :), moon_rot(2, :), moon_rot(3, :))

% x, y, z plot
% x = M * x;
% y = M * y;
% z = M * z;
% mult = 1000000;
% plot3([x(1)*mult; 0], [x(2)*mult; 0], [x(3)*mult; 0], 'r')
% plot3([y(1)*mult; 0], [y(2)*mult; 0], [y(3)*mult; 0], 'g')
% plot3([z(1)*mult; 0], [z(2)*mult; 0], [z(3)*mult; 0], 'b')


%% animazione

clc
close all

figure
hold on, axis equal, grid on

max_frame = max(soho_rot(1,:));
min_frame = min(soho_rot(1,:));
set(gca,'XLim',[min(min_frame, 0) max_frame])
max_frame = max(soho_rot(2,:));
min_frame = min(soho_rot(2,:));
set(gca,'YLim',[min(min_frame, 0) max_frame])
max_frame = max(soho_rot(3,:));
min_frame = min(soho_rot(3,:));
set(gca,'ZLim',[min(min_frame, 0) max_frame])

plot3(0, 0, 0, '.g', markerSize = 20)

for i = 1:length(soho_rot)
    view(i/10, 10)
    a = plot3(soho_rot(1, 1:i), soho_rot(2, 1:i), soho_rot(3, 1:i), 'b');
    b = plot3(soho_rot(1, i), soho_rot(2, i), soho_rot(3, i), 'or');
    drawnow

    delete(a)
    delete(b)
end
