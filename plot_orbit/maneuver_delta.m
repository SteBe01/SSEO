clear
close all

format longG

addpath('time')

%for 14 feb 1996 insertion manuever
[t, r_soho, v_soho] = extract_horizons('./horizon/soho_ins_fine.txt');
[~, r_sun, ~] = extract_horizons('./horizon/sun_ins_fine.txt');
[~, r_moon, ~] = extract_horizons('./horizon/moon_ins_fine.txt');

%for 20 march 1996 trim maneuver
%[t, r_soho, v_soho] = extract_horizons('./horizon/trim_soho.txt');
%[~, r_sun, ~] = extract_horizons('./horizon/trim_sun.txt');
%[~, r_moon, ~] = extract_horizons('./horizon/trim_moon.txt');

sun_mu = 1.327124400e11;
earth_mu = 3.986004418e5;
moon_mu = 4.9048695e3;

dt = diff(t);
dt = dt(1) * 24 * 3600; %convert to seconds

pred_error = [];
for i=1:length(t)-1
	earth_vec = -r_soho(:, i);
	sun_vec = r_sun(:, i) - r_soho(:, i);
	moon_vec = r_moon(:, i) - r_soho(:, i);

	a_earth = earth_mu / norm(earth_vec)^3 * earth_vec;
	a_sun = sun_mu / norm(sun_vec)^3 * sun_vec;
	a_moon = moon_mu / norm(moon_vec)^3 * moon_vec;

	gravity_a = a_earth + a_sun + a_moon;

	predicted_next_dv = v_soho(:, i) + gravity_a * dt;

	pred_error = [pred_error, v_soho(:, i+1) - predicted_next_dv];
end

%print maneuver date
[~, index] = max(vecnorm(pred_error));
jd2date(t(index))

figure
plot(t(1:end-1), pred_error)
grid on
legend('x', 'y', 'z')
