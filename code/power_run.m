

%%%%%%%%%% things that needs to be changed %%%%%%%%%%
%% 1. a list of all voxel names, in the format of "voxel_x_y_z"; data should be a single column
% voxel = importdata('/data/brutus_data32/Hongjie/ACP_R199_voxel_name.txt', "\n");
%% 2. sample frequency for rsfmri data; 0.78 for ACP and ?? for UKBB
fs = 1/0.78; % Sample frequency (samples per unit time or space)
%% 3. input path "file_name"; should be voxel-wise; 
%%%%% a matrix with subject on the row and time point on the column;
%% 4. output path "out_path"

%%%%%%%%%% --- %%%%%%%%%%

region = importdata('/data/brutus_data34/Hongjie/ACP_all_region_name.txt', "\n");

% idx_reg = evalin('base', argv(){1}); % Get the argument
idx_reg = str2double(getenv('ARG1'));

fprintf('This is region %s.\n', char(region(idx_reg)));

temp = strcat("/data/brutus_data34/Hongjie/ACP_rsTS_by_voxel/", region(idx_reg), "/voxel_name.txt");
voxel = importdata(temp, "\n");

out_dir = char(strcat("/data/brutus_data34/Hongjie/ACP_power/", region(idx_reg), "/"));
% mkdir(out_dir);

%%%%%%%%%% --- %%%%%%%%%%


for v = 1:length(voxel)
	out_path = strcat("/data/brutus_data34/Hongjie/ACP_power/", region(idx_reg), "/", voxel(v), ".txt");

	if exist(out_path, 'file') == 0
		tic
		fprintf('This is voxel No. %d.\n', v);

		file_name = strcat("/data/brutus_data34/Hongjie/ACP_rsTS_by_voxel/", region(idx_reg), "/rsTS_", voxel(v), '.txt');

		rts = readtable(file_name);
		rts = table2array(rts);
		d1 = transpose(rts);

		m = length(rts(1,:)); % original sample length
		ss = length(rts(:,1));
		n = pow2(nextpow2(m));  % transform length
		f = (0:n-1)*(fs/n);

		pow_all = zeros(ss,n/2);
		for i = 1:ss
			x = rts(i,:);
			if isnumeric(x)==0
				x = str2double(x);
			end
			y = fft(x,n);        % DFT of signal
			power = abs(y).^2/n; 
			pow_all(i,:) = power(1:floor(n/2));
		end
	
		t2 = toc;
		% fprintf('Elapsed time: %.4f seconds\n', t2);
		dlmwrite(out_path, pow_all, 'delimiter', '\t');
	end
end

