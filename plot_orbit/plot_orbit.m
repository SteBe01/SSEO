clear
close all

format longG

addpath('time')

[t, r_soho, v_soho] = extract_horizons('./soho_2005.txt');
[~, r_sun, ~] = extract_horizons('./sun_2005.txt');
[~, r_moon, ~] = extract_horizons('./moon_2005.txt');

%t = date2jd([1996, 2, 14, 0, 0, 0])

%tspan = t0 + dt / (24 * 3600);
%plot(tspan, vecnorm(v_soho));
%figure
%plot(tspan, v_soho(1, :));
%figure
%plot(tspan, v_soho(2, :));
%figure
%plot(tspan, v_soho(3, :));
%return

xlabel('earth-sun(x)')
ylabel('sun dir(y)')
zlabel('z')
plot3(1.5e6, 0, 0, '.', MarkerSize=20)
grid on;
hold on;
plot3(0, 0, 0, '.', MarkerSize=10);
view(60, 30);
xlim([0, 1800000]);
ylim([-1000000, 1000000]);
zlim([-200000, 120000]);
%
%an = animatedline(0, 0, 0);
for i=2:size(r_soho, 2)-1
	x = r_sun(:, i);
	x = x / norm(x);
	diff_sol = r_sun(:, i+1) - r_sun(:, i-1);
	diff_sol = diff_sol / norm(diff_sol);
	z = cross(x, diff_sol);
	z = z / norm(z);
	y = cross(z, x);

	%plot3([0, z(1)], [0, z(2)], [0, z(3)]);
	%grid on
	%drawnow

	M = [x, y, z]';
	soho_rot(:, i) = M * r_soho(:, i);
	moon_rot(:, i) = M * r_moon(:, i);

	%t = t0 + dt(i) / (24 * 3600);
	%date = jd2date(t);

	%datestr = strcat(num2str(date(1)),'-',num2str(date(2)),'-',num2str(date(3)));

	%addpoints(an, soho_rot(1, i), soho_rot(2, i), soho_rot(3, i))
	%title(datestr)
	%drawnow
end
plot3(soho_rot(1, :), soho_rot(2, :), soho_rot(3, :))
plot3(moon_rot(1, :), moon_rot(2, :), moon_rot(3, :))
