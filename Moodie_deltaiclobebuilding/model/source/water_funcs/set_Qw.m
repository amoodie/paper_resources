function [Qw] = set_Qw(def, num)
    % this function returns one year's worth (365 d) of discharge data for modeling
    % the definition of num is dependent on def
    %
    switch def
        case 'constant'
        % formulation for constant Q at each t
            yr = repmat(num, 365, 1);
            Qw = yr;
        case 'mean_historical'
        % formulation for repeating daily mean Q at each t, from YRIHR field data
            discharge_data = csvread('LIJ_daily_mean_data.csv', 1, 1); % read from row 1, column 1 (exclude header and ids)
            yr = discharge_data(:,15); % load the moving average into Qw
            Qw = yr;
        case 'engineered'
        % formulation for engineered discharge
            yr = repmat(400, 365, 1);
            yr(180:255) = num; % flood
            yr(165:180) = linspace(yr(165), yr(180), length(165:180));
            yr(255:270) = linspace(yr(255), yr(270), length(255:270));
            Qw = yr;
        case 'historical'
        % formulation for choosing a start point in the hist data
            discharge_data = csvread('LIJ_historical_dailyQw.csv', 0, 0); % read from row 0, column 0 (no header and ids)
            smooth = conv(discharge_data(1:end), ones(20, 1) / 20, 'same');
            index = 366 + (365 * (num - 2));
            Qw = smooth(index:index+365-1);
    end
    Qw(Qw < 1) = 1;
end