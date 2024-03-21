clear
close all

format longG

addpath('./utils/');

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);

%xcc burn
[t, r_soho, v_soho] = extract_horizons('../horizons_data/xcc_suncentre/soho.txt');
[~, r_earth, v_earth] = extract_horizons('../horizons_data/xcc_suncentre/earth.txt');
[~, r_moon, v_moon] = extract_horizons('../horizons_data/xcc_suncentre/moon.txt');

%insertion in halo orbit
%[t, r_soho, v_soho] = extract_horizons('../horizons_data/haloins_suncentre/soho.txt');
%[~, r_earth, v_earth] = extract_horizons('../horizons_data/haloins_suncentre/earth.txt');
%[~, r_moon, v_moon] = extract_horizons('../horizons_data/haloins_suncentre/moon.txt');

%trim maneuver
%[t, r_soho, v_soho] = extract_horizons('../horizons_data/trim_suncentre/soho.txt');
%[~, r_earth, v_earth] = extract_horizons('../horizons_data/trim_suncentre/earth.txt');
%[~, r_moon, v_moon] = extract_horizons('../horizons_data/trim_suncentre/moon.txt');

%around 1 year of orbit mainteinance
%[t, r_soho, v_soho] = extract_horizons('../horizons_data/maintain_suncentre/soho.txt');
%[~, r_earth, v_earth] = extract_horizons('../horizons_data/maintain_suncentre/earth.txt');
%[~, r_moon, v_moon] = extract_horizons('../horizons_data/maintain_suncentre/moon.txt');

error = [];
for i = 1:length(t)-1
	y0 = [r_soho(:, i); v_soho(:, i)];

	%convert time to seconds for ode45, as velocity is in km/s
	dt = (t(i+1) - t(i)) * 24 * 3600;

	[t_slice, y_slice] = ode45(@dy, [0, dt], y0, options, r_earth(:, i), v_earth(:, i), r_moon(:, i), v_moon(:, i));
	error = [error, y_slice(end, 4:6)' - v_soho(:, i+1)];
end

dv = sum(vecnorm(error))

figure
plot(t(1:end-1), error)
grid on

function dy = dy(t, y, earth_pos0, earth_v, moon_pos0, moon_v)
	sun_mu = 1.327124400e11;
	earth_mu = 3.986004418e5;
	moon_mu = 4.9048695e3;

	r_soho = y(1:3);

	%for small enough deltaT this is a good approximation
	earth_pos = earth_pos0 + earth_v * t;
	moon_pos = moon_pos0 + moon_v * t;

	earth_vec = earth_pos - r_soho;
	moon_vec = moon_pos - r_soho;
	sun_vec = -r_soho;

	a_earth = earth_mu / norm(earth_vec)^3 * earth_vec;
	a_moon = moon_mu / norm(moon_vec)^3 * moon_vec;
	a_sun = sun_mu / norm(sun_vec)^3 * sun_vec;

	gravity_a = a_earth + a_sun + a_moon;

	dy = zeros(6, 1);
	dy(1:3) = y(4:6);
	dy(4:6) = gravity_a;
end
