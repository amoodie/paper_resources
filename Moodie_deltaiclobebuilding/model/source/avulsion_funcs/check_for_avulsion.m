function [activate, loc] = check_for_avulsion(eta, zed, H, Hnbf, rad_idx, gammapex_idx, preavulSwitch, preavulThresh, preavulTrigg, closeThresh, nx)
    switch preavulSwitch
        case 'any-setup+over+topset'
            % deposit thickness exceeds threshold of the flow depth *and* levee is exceeded by threshold *and* is within topset
            setup = zed - eta < (Hnbf * (1-preavulThresh));
            curr = eta + H;
            trigger = (curr - zed) > preavulTrigg .* H;
            topset = and((1:nx+1) < rad_idx-1, (1:nx+1) > gammapex_idx-1);
            check = and(and(setup, trigger), topset); % if ALL setup and trigger and topset == 1 at a node, check = 1 at that node
            activate = any(check);
            if activate
                close_setup = (zed - eta) / (Hnbf * (1-preavulThresh)); % remake setup for closeness
                new_setup = close_setup < (1 + closeThresh);
                new_check = and(new_setup, topset);
                loc = find(new_check, 1, 'first');
            else
                loc = [];
            end
        case 'any-setup+topset'
            % deposit thickness exceeds threshold of the flow depth *and* is within topset
            setup = zed - eta < (Hnbf * (1-preavulThresh));
            topset = and((1:nx+1) < rad_idx-1, (1:nx+1) > gammapex_idx-1);
            check = and(setup, topset);
            activate = any(check);
            if activate
                close_setup = (zed - eta) / (Hnbf * (1-preavulThresh)); % remake setup for closeness
                new_setup = close_setup < (1 + closeThresh);
                new_check = and(new_setup, topset);
                loc = find(new_check, 1, 'first');
            else
                loc = [];
            end
        case {'qingshuigou', 'na', 'ac1'}
            activate = false;
            loc = NaN;
        case 'any-setup'
            % deposit thickness exceeds threshold of the flow depth
            setup = zed - eta < (Hnbf * (1-preavulThresh));
            check = setup;
            activate = any(check);
            loc = find(check, 1);
        case 'any-over'
            % water surface exceeds the topset height (Hbf+etai) (i.e. fixed levees)
            curr = eta + H; % current water surface
            trigger = (curr - zed) > preavulTrigg .* H;
            check = trigger;
            activate = any(check);
            loc = find(check, 1);
        case 'any-setup+over'
            % deposit thickness exceeds threshold of the flow depth *and* levee is exceeded by threshold
            setup = zed - eta < (Hnbf * (1-preavulThresh));
            curr = eta + H;
            trigger = (curr - zed) > preavulTrigg .* H;
            check = and(setup, trigger); % if both setup and trigger == 1 at a node, check = 1 at that node
            activate = any(check);
            loc = find(check, 1);
        case 'testing'
            % define when to trigger an avulsion (for 'sustainability gained')
        case 'imposed'
            % use for condition at Xihekou (sp?) at 10000 m3/s
    end
    if isempty(loc)
        loc = NaN;
    end
end