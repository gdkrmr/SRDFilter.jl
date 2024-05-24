% smooth all csv files in 'data/' and save the smoothed version in 'data/smooth'

% smooth every column in csv files
% csv files must be header less and can only contain numerical columns

% save one file per smoother [smoothMS, smoothMS1] and parameter combination degs x ms



source("sincsmoother.m")

% parameters, will iterate over all combinations
degs = [2, 4, 6, 8, 10];
ms = 1:10;


% get csv files in "data/" directory
files = readdir("data");
idxs = strfind(files, ".csv");
filescsvidx = cellfun(@(x) length(x) > 0, idxs);
filescsv = files(filescsvidx);

% iterate over files
for fi = 1:numel(filescsv)
  filename = filescsv{fi};
  csv_data = csvread(strcat("data/", filename));
  %% iterate over degrees
  for d = degs
    % iterate over halfkernel widths
    for m = ms
      % data to be saved
      csv_data_smooth = [];
      csv_data_smooth1 = [];
      % iterate over all columns of csv data
      for ci = 1:size(csv_data, 2)
        data = csv_data(:, ci)';
        % smooth both kernels
        out = smoothMS(data, d, m)';
        out1 = smoothMS1(data, d, m)';
        % add smoothed data
        csv_data_smooth = [csv_data_smooth, out];
        csv_data_smooth1 = [csv_data_smooth1, out1];
      end
      % save data
      csv_smooth_filename = strcat("data/smooth/MS_", filename(1:end - 4), "_", "deg_", num2str(d), "_", "m_", num2str(m), ".csv");
      csv_smooth1_filename = strcat("data/smooth/MS1_", filename(1:end - 4), "_", "deg_", num2str(d), "_", "m_", num2str(m), ".csv");
      csvwrite(csv_smooth_filename, csv_data_smooth);
      csvwrite(csv_smooth1_filename, csv_data_smooth1);
    end
  end
end
