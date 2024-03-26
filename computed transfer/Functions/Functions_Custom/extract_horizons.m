function[t, r, v] = extract_horizons(filename)
	%% EXTRACT NASA HORIZON VECTOR TABLE {r, v}
	%% ENABLE CSV IN HORIZONS OPTION
	% t0 start time in Julian Day Number, Barycentric Dynamical Time
	% dt time elapsed from t0 [s]
	% r  position vector
	% v  velocity vector

	t = [];
	r = [];
	v = [];
	fid = fopen(filename);
	tline = fgetl(fid);
	start_extract = false;
	while ischar(tline)
		if contains(tline, 'SOE')
			start_extract = true;
			tline = fgetl(fid);
		end
		if contains(tline, 'EOE')
			break
		end

		if start_extract
			data = textscan(tline, '%f %s %f %f %f %f %f %f', 'Delimiter', ',');
			t = [t, data{1}];
			r = [r, [data{3}; data{4}; data{5}]];
			v = [v, [data{6}; data{7}; data{8}]];
		end

		tline = fgetl(fid);
	end
	fclose(fid);
end
