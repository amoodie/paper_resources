function [s, traj] = virtualdelta_wrapper_CSDMS(preavulThresh, preavulTrigg, mouSwitch, Qwcall, Qwnum)
%% mother function
    
    con.g = 9.81; % gravitational constant
    con.rho_f = 1000;
    con.rho_s = 2650;
    con.Rr = (con.rho_s / con.rho_f) - 1;
    con.phi = 0.4; % did lab experiment to justify
    con.nu = 1.004 * 1e-6; % kinematic viscosity, m^2/s
    
    L = 500e3; % total model length, meters
    nAvulSeek = 11; % number of avulsion cycles to run for
    nx = 400; % number of spatial steps
    dx = L/nx; % spatial step
    x = 0:dx:L; % define x-coordinates
    
    au = 0.5; % winding coefficient
    dt = 0.25; % timestep, explicit to fraction of day
    dtsec = dt * 86400; % convert timestep to seconds
    
    D50 = 90e-6; % D50 in meters
    
    [Qw] = set_Qw(Qwcall, Qwnum); % discharge, m3/s
    [p1, p2] = set_GEH('HM'); % coefficients to Generalized Engelund Hansen formulation
    preavulSwitch = 'any-setup+over+topset'; % avulsion setup/trigger method

    if isnan(preavulTrigg) % convert nan input to ignore overbank flow
        preavulSwitch = 'any-setup+topset';
    end
    
    sigma = 0; % subsidence rate, mm/yr
    
    Bc0 = 400; % channel width
    Bf0 = Bc0*6; % floodplain width
    Bo0 = 9e3; % overbank width
    thet = 2; % plume spreading angle
    
    gamm = 90; % fan delta angle, degrees
    S0 = 6.4e-5; % initial fan slope
    zed0 = 15; % initial topset apex elevation
    
    H0 = 0; % downstream boundary for backwater calc
    Cf = 0.001; % used for backwater, need to come up with another solution I think?
    Qwbf = 3000; % bankfull discharge, m3/s
    HnfunDW = (@(Q, B, Cf, S0, g) ((Q./B).^(2/3).*(Cf./(g.*S0)).^(1/3))); % generalized Darcy Weisbach, from Ganti et al., 2014 supp info
    Hnbf = HnfunDW(Qwbf, Bc0, Cf, S0, 9.81); % bankfull normal flow depth
    Qwform = 1300; % formative discharge, m3/s
    Hform = HnfunDW(Qwform, Bc0, Cf, S0, 9.81); % formative normal flow depth
    d.Hform = Hform;
    d.Lblow = Hform/S0;
    d.Hnbf = Hnbf;
    
    [zed] = set_zedi(zed0, S0, nx, dx); % set initial topset surface
    [eta] = set_etai(zed, Hnbf, S0, dx); % set initial channel bed
    
    radi = find(zed <= 0, 1)*dx; % initial delta radius (from model boundary), meters, intersection of delta with sea level
    moui = radi; % inital mouth location, meters 
    gammapex = radi * 0.85; % apex of delta spreading angle
    mou = moui; % initial mouth
    rad = radi; % initial radius of delta
    
    [Bc, Be] = set_B(get_idx(x, mou), Bc0, Bf0, Bo0, thet, get_idx(x, rad), nx, dx); % Bc is channel width, Be is Exner width
    
    cnt.l = 1; % time loop counter
    cnt.N = 2; % print counter
    cnt.y = 0; % year counter
    cnt.d = 1; % day counter
    cnt.dfrac = 0; % day fraction counter
    cnt.delap = 0; % count days elapsed, no resets
    
    % set some initial conditions
    [qs] = set_ICs(eta, Qw(1), Bc, H0, Cf, p1, p2, D50, nx, dx, con); % necessary initial condition
    d.lij_loc = moui - 56e3; % location of lijin (a datum)
    [d.lij_idx] = get_idx(x, moui - 56e3);
    Hbf = get_backwater_dBdx(eta, get_slope(eta, nx, dx), Bc, H0, Cf, Qwbf, nx, dx); % bankful depth along domain
    flood_surf = (Hbf+eta); % flood surface elevation
    Hbasin = H0 - eta; % depth in basin
    zed = [zed(zed>flood_surf) flood_surf(~(zed>flood_surf))]; % update delta dtopset to be 0 in basin
    
    etai = eta; % save initial eta profile
    etalastavul = etai; % eta of just after last avulsion
    s.mou(1) = mou;
    
    nAvul = 0;
    while nAvul < nAvulSeek
        mou0 = mou; % old mou location
        
        % bed CALCS
        [S] = get_slope(eta, nx, dx); % slope
        [H] = get_backwater_dBdx(eta, S, Bc, H0, Cf, Qw(cnt.d), nx, dx); % water surface
        [U] = get_U(Qw(cnt.d), H, Bc); % velocity
        [qs] = qtGEH(H, U, S, D50, p1, p2, con); % sediment transport
        qsBC = qs(1)*1.0005; % adjustment to upstream boundary condition
        [eta] = new_eta(eta, H, au, qs, qsBC, Bc, Be, nx, dx, dtsec, con); % exner equation (sed conservation)
        [eta] = apply_subsidence(eta, sigma, dtsec); % subsidence
        [dep] = get_thick(eta, etalastavul); % calc deposit thickness (for avulsion check)
        [Bc, Be, mou] = update_long(H, Hbasin, dep, Hform, thet, rad, mou0, Bc0, Bf0, Bo0, x, nx, dx, mouSwitch); % update params of long profile from lobe prog.
        
        cnt.l = cnt.l + 1;
        
        % handle counting for the daily varying Qw
        if cnt.dfrac+eps < 1 % if less than one day has passed
            cnt.dfrac = cnt.dfrac + dt; % add a timestep
        else % if equal to one day has passed
            cnt.d = cnt.d + 1; % add a daystep
            cnt.delap = cnt.delap + 1; % add an elapsed day
            cnt.dfrac = dt; % reset day fraction counter
            [activate, avulloc] = check_avul(eta, zed, ...
                H, Hnbf, rad, preavulSwitch, preavulThresh, preavulTrigg, x, nx); % check for avulsion at this daystep
            switch activate
                case 1 % if true = avulsion occurs
                    [Vol_delta, Vol_lobe] = calc_volume(eta, Bc, Be, rad, etalastavul, L, x, dx, con.phi); % volume of sed in delta portions
                    [rad, zed] = avulsion_plan(rad, zed, gamm, gammapex, Bf0, Hbasin, Vol_delta, Vol_lobe, x, dx, con); % redistribute sediment in planform
                    [eta, Bc, Be, mou] = avulsion_long(etai, eta, zed, Hnbf, Bc0, Bf0, Bo0, avulloc, thet, rad, x, nx, dx); % redistribute sediment in long-prof
                    etalastavul = eta; % update 
                    Hbasin = H0 - eta; % update
                    nAvul = nAvul + 1; % update
                    disp(['nAvul = ' num2str(nAvul) ', iter = ' num2str(cnt.delap)]) % display for user
            end
            [s] = to_storage_all(s, cnt, mou, avulloc, rad, eta, H, zed); % store information for plotting later
%             [s] = to_storage_lite(s, cnt, mou, avulloc, rad)
            if cnt.d == 366 % if a year has passed
                cnt.y = cnt.y + 1; % add a yearstep
                cnt.d = 1; % reset day counter
            end
        end
    end
    
    virtFIG = figure();
    plot_virt(x, s, L, gamm, gammapex, Bf0, dx, virtFIG)
    
    
    [s, d] = do_maths(s, d, dx); % do some useful calculations for exploring data and plotting
    s.d = d; % for output data
end

%% setup/initialization functions
function [zed] = set_zedi(zed0, S0, nx, dx)
    zed = linspace(zed0, zed0 - S0*(nx*dx), nx+1);
end

function [eta] = set_etai(zed, Hbf, S0, dx)
    eta = zed - Hbf;
    basin = eta < -16;
    Nbasin = length(find(basin));
    brk = eta(find(~basin, 1, 'last'));
    eta(basin) = linspace(brk, brk - (S0*1e-1)*(dx*Nbasin), Nbasin);
end

function [Bc, Be] = set_B(mou, Bc0, Bf0, Bo0, thet, rad, nx, dx)
    Bc = zeros(1, nx+1);
    Bc(1:mou) = Bc0;
    Bc(mou+1:end) = Bc0 + 2 * (linspace(1, (nx+1-mou)*dx, (nx+1-mou)) * tand(thet));
    Be = zeros(1, nx+1);
    Be(1:rad) = Bc0 + Bf0;
    Be(rad+1:mou) = Bc0 + Bo0; % constant width of lobe width
    Be(mou+1:end) = Bc0 + Bf0 + 2 * (linspace(1, (nx+1-mou)*dx, (nx+1-mou)) * tand(thet));
end

function [qs] = set_ICs(eta, Qw, B, H0, Cf, p1, p2, D50, nx, dx, con)
    [S] = get_slope(eta, nx, dx);
    [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx);
    [U] = get_U(Qw, H, B);
    [qs] = qtGEH(H, U, S, D50, p1, p2, con);
end

function [Qw] = set_Qw(def, num)
    switch def
        case 'con'
        % formulation for constant Q at each t
            Qw = repmat(1300, 365, 1);
        case 'mean'
        % formulation for repeating daily mean Q at each t, from YRIHR field data
            discharge_data = csvread('./LIJ_daily_mean_data.csv', 1, 1); % read from row 1, column 1 (exclude header and ids)
            yr = discharge_data(:,15); % load the moving average into Qw
%             Qw = repmat(yr, 1, ceil(Tyears));
            Qw = yr;
        case 'hist'
        % formulation for historical daily mean Q at each t, from YRIHR field data
            discharge_data = csvread('./LIJ_historical_dailyQw1950-2000.csv', 0, 0); % read from row 0, column 0 (no header and ids)
            smooth = conv(discharge_data(1:end), ones(20, 1) / 20, 'same');
            Qw = smooth;
        case 'eng'
        % formulation for engineered discharge
            yr = repmat(1000, 365, 1);
            yr(180:255) = num; % flood
            yr(165:180) = linspace(yr(165), yr(180), length(165:180));
            yr(255:270) = linspace(yr(255), yr(270), length(255:270));
%             Qw = repmat(yr, 1, ceil(Tyears));
            Qw = yr;
        case 'num'
        % formulation for choosing a start point in the hist data
            discharge_data = csvread('./LIJ_historical_dailyQw1950-2000.csv', 0, 0); % read from row 0, column 0 (no header and ids)
            smooth = conv(discharge_data(1:end), ones(20, 1) / 20, 'same');
            index = 366 + (365 * (num - 2));
            Qw = horzcat(smooth(index:end), smooth(1:index-1));
    end
    Qw(Qw < 1) = 1;
end

function [p1, p2] = set_GEH(def)
    switch def
        case 'EH'
        % use original Engelund and Hansen, 1967 relation for bed matrial load
            p1 = 0.05; % EH original alpha
            p2 = 2.5; % EH original n
        case 'HM'
        % use modfied EH from Ma et al., 201? relation for YR
            p1 = 0.895; % HMqt alpha, Lijin
            p2 = 1.678; % HMqt n, Lijin
        case 'adj'
        % set to whatever, just exploring parameter space here
            p1 = 0.5; % adjust for reasonable
            p2 = 2; % adjust for reasonable
    end
end

%% general functions
function [S] = get_slope(eta, nx, dx)
    % return slope of input bed (eta)
    S = zeros(1,nx+1);
    S(1) = (eta(1) - eta(2)) / dx;
    S(2:nx) = (eta(1:nx-1) - eta(3:nx+1)) ./ (2*dx);
    S(nx + 1) = (eta(nx) - eta(nx + 1)) / dx;
end

function [idx] =  get_idx(x, loc)
    % get index of some location 'loc' *downstream* from upper boundary
    diffs = abs(x - loc);
    [~, idx] = min(diffs);
end

function [s] = to_storage_all(s, cnt, mou, avul, rad, eta, H, zed)
    s.eta(:, cnt.delap) = eta;
    s.H(:, cnt.delap) = H;
    s.zed(:, cnt.delap) = zed;
    s.mou(cnt.delap) = mou;
    s.avul(cnt.delap) = avul;
    s.rad(cnt.delap) = rad;
%     s.Qw(cnt.delap) = Qw;
%     s.Lb(cnt.delap) = Lb;
end

function [s] = to_storage_lite(s, cnt, mou, avul, rad)
    s.mou(cnt.delap) = mou;
    s.avul(cnt.delap) = avul;
    s.rad(cnt.delap) = rad;
end

function [dep] = get_thick(eta, baseline)
    dep = (eta - baseline);
end

%% bed evolution functions
function [eta] = new_eta(eta, H, au, qs, qs_u, Bc, Be, nx, dx, dt, con)
    eta0 = eta;
    Q = qs .* Bc; % changed this to only account for the width of the channel for Q calc, but use Be for exner calc
    Qu = qs_u * Bc(1);
    dQdx(nx+1) = (Q(nx+1) - Q(nx)) / dx; % error possibly on typo edit??
    dQdx(1) = au*(Q(1)-Qu)/dx + (1-au)*(Q(2)-Q(1))/dx;
    for i=2:nx % transport gradient, positive means erosion 
        dQdx(i) = au*(Q(i)-Q(i-1))/dx + (1-au)*(Q(i+1)-Q(i))/dx;
    end
    eta = eta0 - ((1/(1-con.phi)) .* dQdx .* dt .* (1 ./ Be));
end

function [eta] = apply_subsidence(eta, sigma, dtsec)
    eta0 = eta;
    dz = (-(sigma / 1000) / 31536e3) * dtsec;
    eta = eta0 + dz;
end

function [Bc, Be, mou] = update_long(H, Hbasin, depthick, Hform, thet, rad, mou0, Bc, Bf, Bo, x, nx, dx, mouSwitch)
    mou0_idx = get_idx(x, mou0);
    rad_idx = get_idx(x, rad);
    switch mouSwitch
        case 'na'
            mou = mou0;
        case 'rem'
            if Hbasin(mou0_idx) - depthick(mou0_idx) < Hform % if critical remaining
                mou = mou0 + dx;
            else
                mou = mou0;
            end
    end
    [Bc, Be] = set_B(get_idx(x, mou), Bc, Bf, Bo, thet, get_idx(x, rad), nx, dx); % mou, Bc0, Bf0, Bo0, nx, dx
end

function [activate, loc] = check_avul(eta, zed, H,  Hnbf, rad, def, thresh, over, x, nx)
    % check for an avulsion based on the inpute conditions
    switch def
        case {'na', 'ac1'}
            activate = false;
            loc = NaN;
        case 'any-setup+topset'
            % deposit thickness exceeds threshold of the flow depth *and* is within topset
            setup = zed - eta < (Hnbf * (1-thresh));
            topset = (1:nx+1) < get_idx(x, rad)-1;
            check = and(setup, topset);
            activate = any(check);
            loc = find(check, 1);
        case 'any-setup+over+topset'
            % deposit thickness exceeds threshold of the flow depth *and* levee is exceeded by threshold *and* is within topset
            setup = zed - eta < (Hnbf * (1-thresh));
            curr = eta + H;
            trigger = (curr - zed) > over .* H;
            topset = (1:nx+1) < get_idx(x, rad)-1;
            check = and(and(setup, trigger), topset); % if ALL setup and trigger and topset == 1 at a node, check = 1 at that node
            activate = any(check);
            loc = find(check, 1);
    end
    if isempty(loc)
        loc = NaN;
    end
end

function [eta, Bc, Be, mou] = avulsion_long(etai, eta, zed, Hnbf, Bc0, Bf0, Bo0, loc, thet, rad, x, nx, dx)
    mou = rad;
    mou_idx = get_idx(x, mou);
    rad_idx = get_idx(x, rad);
    win = 11; % width of taper, i.e. number of nodes from loc to taper
    hwin = round((win-1)/2);
        eta(1:loc-hwin) = eta(1:loc-hwin);
        eta(loc-hwin:loc+hwin) = linspace(eta(loc-hwin), etai(loc+hwin), ...
            length(eta(loc-hwin:loc+hwin)));
        eta(loc+hwin+1:end) = etai(loc+hwin+1:end);
    [Bc, Be] = set_B(mou_idx, Bc0, Bf0, Bo0, thet, rad_idx, nx, dx);
end

function [rad, zed] = avulsion_plan(rad, zed, gamm, gammapex, Bf0, Hbasin, Vol_delta, Vol_lobe, x, dx, con)
    zed0 = zed; % old topset surface
    rad0 = rad; % old radius location
    
    rad0_idx = get_idx(x, rad0); % old radius index
    shorelen = (rad0 - gammapex) * deg2rad(gamm); % length of the shoreline at old radius
    % numerical integration below
    idxlim = 50000;
    int = cumtrapz((0:1:idxlim).*(6.4e-5)) .* ...
        (((rad0-gammapex):1:(rad0-gammapex+idxlim)) .* deg2rad(gamm)) + ...
        (Hbasin(rad0_idx).*(0:1:idxlim)) .* (((rad0-gammapex):1:(rad0-gammapex+idxlim)) .* deg2rad(gamm));
    drad = find((sum(Vol_lobe)/con.phi) < int, 1, 'first');
    if isempty(drad)
        error('increase idx limit')
    end
    
    rad = rad0 + drad; % new radius
    rad_idx = get_idx(x, rad);
    gammapex_idx = get_idx(x, gammapex); % index of delta apex
    Adelfun = (@(rad, gamm) (pi .* rad.^2) .* (deg2rad(gamm)/(2*pi)));
    A_deltacum = Adelfun(x, gamm); % pretend as if this is calculated FROM gammapex_idx
    A_delta = [0 (A_deltacum(2:end) - A_deltacum(1:end-1))]; % area is diff between two cells
    
    Bneck = Bf0*1; % this multiplier is like the channel belt width
    A_delta(A_delta < Bneck*dx) = Bneck*dx;
    
    neck_idx = false(size(x));
    neck_idx(1:gammapex_idx) = true;
    top_idx = false(size(x));
    top_idx(gammapex_idx+1:rad0_idx) = true;

    zed(neck_idx) = zed0(neck_idx) + (sum(Vol_delta(neck_idx)) / (Bneck.*x(gammapex_idx)) / (1-con.phi));
    zed(top_idx) = zed0(top_idx) + (Vol_delta(top_idx) ./ A_delta(2:length(find(top_idx))+1) ./ (1-con.phi));
    
    Ch = (12*sum(Vol_delta(top_idx))) / (pi * ((rad_idx-gammapex_idx)*dx)^2);
    zed(top_idx) = zed0(top_idx) + linspace(Ch, 0, length(find(top_idx)));
end

function [Vol_delta, Vol_lobe] = calc_volume(eta, Bc, Be, rad, ETAlastavul, L, x, dx, phi)
    deta_avul = eta - ETAlastavul;
    rad_idx = get_idx(x, rad);
    Vol_delta = (deta_avul(1:rad_idx) .* [dx/2 repmat(dx, 1, rad_idx-1)] .* (Be(1:rad_idx)-Bc(1:rad_idx))) .* 1-phi; 
    Vol_lobe = (deta_avul(rad_idx+1:end) .* [repmat(dx, 1, length(deta_avul(rad_idx+1:end))-1) dx/2] .* Be(rad_idx+1:end)) .* 1-phi;
end

%% water functions 
function [H] = get_backwater_dBdx(eta, S, B, H0, Cf, Qw, nx, dx)
    % backwater formulated for changing width
    H = NaN(1,nx+1); % preallocate depth 
    H(nx+1) = abs(H0 - eta(nx+1)); % water depth at downstream boundary
    for i = nx:-1:1
        % predictor step: computation of a first estimation of the water depth Hp
        [Frsqp] = get_froude(Qw, H(i+1), B(i)); % calculate froude from conditions at i+1
        [dBdx] = get_dBdx(B, nx, dx);
        [dHdxp] = get_dHdx_dBdx(S(i+1), Cf, Frsqp, H(i+1), B(i), dBdx(i)); % get dHdx width changing form
        Hp = H(i+1) - dHdxp * dx; % solve for H prediction
        % corrector step: computation of H
        [Frsqc] = get_froude(Qw, Hp, B(i)); % calculate froude at i with prediction depth
        [dHdxc] = get_dHdx_dBdx(S(i), Cf, Frsqc, Hp, B(i), dBdx(i)); % doaa
        % convolution of prediction and correction, trapezoidal rule
        H(i) = H(i+1) - ( (0.5) * (dHdxp + dHdxc) * dx );
    end
    if Qw < 75
        [~, idx] = find(eta <= 0, 1);
        H(eta >= -0.2) = H(idx + 1);
    end
end

function [dBdx] = get_dBdx(B, nx, dx)
    dBdx = (B(2:nx+1) - B(1:nx)) ./ dx;
end

function [dHdx] = get_dHdx_dBdx(S_loc, Cf, Frsq, H_loc, B_loc, dBdx)
    % formulation to get dHdX for a changing width backwater formulation
    dHdx = ((S_loc-(Cf*Frsq))/(1-Frsq)) + (Frsq/(1-Frsq)) * (H_loc/B_loc) * (dBdx);
end

function [Frsq] = get_froude(Qw, H, B)
    g = 9.81; % gravitational acceleration constant
    Frsq = ( Qw^2 / (g * B^2* H^3) );
end

function [U] = get_U(Qw, H, B)
    U = Qw ./ (H .* B);
end


%% sediment transport functions
function [qt] = qtGEH(H, U, S, D50, p1, p2, con)
    % function for qt following the nondimensionalized formulation by Ma et al., 2017, some relations from Wright and Parker, 2004
    taua = H .* S ./ (con.Rr*D50);
    tau = taua .* (con.Rr) * con.g * D50 * con.rho_f; % use d50met as proxy for Dgeom, tau == total stress
    Cf = (tau) ./ (con.rho_f .* U.^2);
    qs_n = (p1 .* (taua .^p2)) ./ Cf; % qs_n == qs*, i.e. non-dimensionalized
    qs = qs_n .* sqrt(con.Rr * con.g * (D50 ^3)); % dimensionalize
    qt = qs;
end

%% output functions
function [s, d] = do_maths(s, d, dx)
    d.avulidx = find(~isnan(s.avul));
    d.avullen = (s.mou(d.avulidx-1) - s.avul(d.avulidx)*dx);
    d.lobelen = (s.mou(d.avulidx-1) - s.rad(d.avulidx-1));
    d.avullenmean = mean(d.avullen)/1000;
    d.avultime = [NaN (d.avulidx(2:end) - d.avulidx(1:end-1))];
    d.avultimemean = mean(d.avulidx(2:end) - d.avulidx(1:end-1))/365;
end

function [data] = calc_agg(s, idx)
    data.detadt = (s.eta(idx, 2:end) - s.eta(idx, 1:end-1)) * 1000; % mm/d
    data.dzeddt = (s.zed(idx, 2:end) - s.zed(idx, 1:end-1)) * 1000;
    data.thick = s.eta(idx, 2:end) - s.eta(idx, 1);
end

function [plan] = make_planform(x, gamm, gammapex, rad)
    plan = zeros(100, 2);
    pts = linspace( deg2rad(90), deg2rad(90-gamm), 100);
    plan(:, 1) = cosd(gamm/2)*gammapex + ((rad-gammapex) * cos(pts)); % x coords
    plan(:, 2) = sind(gamm/2)*gammapex + ((rad-gammapex) * sin(pts)); % y coords
end

%% plotting functions
function [] = plot_virt(x, s, L, gamm, gammapex, Bf0, dx, virtFIG)
    mdl = fitlm(find(~isnan(s.avul)), s.avul(~isnan(s.avul)).*dx);
    relcoords = [0.5, 0.2];
    CI = coefCI(mdl, 0.05);
    err = mean(abs(CI(2,:) - mdl.Coefficients.Estimate(2)));

    nPrint = 20;
    toPrint = round( linspace(1, size(s.eta, 2), nPrint) );
    colmap = parula(nPrint);
    
    figure(virtFIG);
    subplot(2, 2, [1 2]);
    cla
    hold on
        [t0f] = plot(x/1000, s.eta(:, 1), 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '-'); % initial
        for pp = 1:nPrint
            p = toPrint(pp);
            plot(x/1000, s.eta(:, p), 'LineWidth', 1, 'Color', colmap(pp, :), 'LineStyle', '-'); % time lines
            plot(x/1000, s.eta(:, p) + s.H(:, p), 'LineWidth', 1, 'Color', colmap(pp, :), 'LineStyle', ':'); % water
%             plot(x, s.zed(:, p), 'LineWidth', 1, 'Color', colmap(pp, :), 'LineStyle', ':'); % zed
            line = plot([s.mou(p)/1000 s.mou(p)/1000], ylim, 'LineWidth', 1, 'Color', colmap(pp, :), 'LineStyle', '--');
        end
        legend(line, {'mouth'});
        xlim([0*L/1000, L/1000]);
        ylim([-20 20])
        title( 'long model', 'FontSize', 14);
        xlabel('x distance (km)','FontSize', 14);
        ylabel('z distance (m)','FontSize', 14);
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
        box on
    hold off
    subplot(2, 2, 3); % planform
    cla
    hold on
        [plan] = make_planform(x, gamm, gammapex, s.rad(1));
        
        [bgun] = plot([0 (cosd(gamm/2)*gammapex)/1000], [Bf0/2/1000 (Bf0/2 + sind(gamm/2)*gammapex)/1000], 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '-'); % background upper
        [bgu] = plot([(cosd(gamm/2)*gammapex)/1000 plan(1, 1)/1000], [(sind(gamm/2)*gammapex)/1000 plan(1, 2)/1000], 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '-'); % background upper
        
        [bgln] = plot([Bf0/2/1000 (Bf0/2 + cosd(gamm/2)*gammapex)/1000], [0 (sind(gamm/2)*gammapex)/1000], 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '-'); % background lower
        [bgl] = plot([(cosd(gamm/2)*gammapex)/1000 plan(end, 1)/1000], [(sind(gamm/2)*gammapex)/1000 plan(end, 2)/1000], 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '-'); % background lower
        
        [t0f] = plot(plan(:, 1)/1000, plan(:, 2)/1000, 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '-'); % initial
        for pp = 1:nPrint
            p = toPrint(pp);
            [plan] = make_planform(x, gamm, gammapex, s.rad(p));
            plot(plan(:, 1)/1000, plan(:, 2)/1000, 'LineWidth', 1, 'Color', colmap(pp, :), 'LineStyle', '-'); % time lines
        end
        axis equal
        xlim([0 L/1000]);
        ylim([0 L/1000])
        title( 'planform model','FontSize', 14);
        xlabel('x distance (km)','FontSize', 14);
        ylabel('y distance (km)','FontSize', 14);
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
        box on
    hold off
    s4 = subplot(2, 2, 4);
    cla
    hold on
        plot(1:size(s.eta, 2), s.avul*(x(2)-x(1))/1000, '*k')
        plot(1:size(s.eta, 2), s.rad/1000, 'LineStyle', '-', 'Color', [0 0 0], 'LineWidth', 2)
        plot(1:size(s.eta, 2), s.mou/1000, 'LineStyle', ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
        plot([find(s.mou>s.rad(1), 1, 'first'):size(s.eta, 2)], ...
            (mdl.Coefficients.Estimate(2).*[find(s.mou>s.rad(1), 1, 'first'):size(s.eta, 2)] + mdl.Coefficients.Estimate(1))/1000, ...
            'LineStyle', '--', 'Color', [0 0 0], 'LineWidth', 1.2)
            xlim([1 size(s.eta, 2)+1])
            set(gca, 'XTick', (1:(size(s.eta, 2)/5):size(s.eta, 2)+1), 'XTickLabel', floor((1:(size(s.eta, 2)/5):size(s.eta, 2)+1)/(365)))
            ylim([min([s.mou s.rad])/1000-5 max(s.mou)/1000+5])
            title('Avulsion timing and location', 'FontSize', 14)
            xlabel('elapsed time (years)')
            ylabel('location from datum (km)')
        params = {['rate = ', num2str(round(mdl.Coefficients.Estimate(2)*365.25/1000, 2)), '\pm', num2str(round(err*365.25/1000, 2)), ' km/yr'], ...
            ['r^2 = ', num2str(round(mdl.Rsquared.Ordinary, 2))]};
        format = sprintf('%s\n', params{:});
        annot = text(relcoords(1), relcoords(2), format(1:end-1), ...
            'Color', [0 0 0], 'Parent', s4, 'units', 'normalized', 'BackgroundColor', [1 1 1]);
        legend('avulsion', 'rad', 'mou', 'Location', 'NorthWest')
        set(gca, 'FontSize', 10, 'LineWidth', 1.5)
        box on
    set(virtFIG,'Visible','on');
    set(virtFIG, 'Pos', [100 100 1000 700]);
    set(virtFIG, 'PaperPositionMode', 'auto')
end
