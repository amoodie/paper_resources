function [analysis] = convert_to_analysis(s, d)
    % convert the full data frames down to small storage of needed info

    analysis.Lblow = d.Lblow;
    
    analysis.Hnform = d.Hnform;
    analysis.Hnbf = d.Hnbf;
    
    analysis.lij_loc = d.lij_loc;
    analysis.lij_idx = d.lij_idx;
    
    analysis.avul_idx = d.avul_idx;
    
    analysis.avul_len = d.avul_len;
    analysis.avul_len_mean = d.avul_len_mean;
    
    analysis.avul_time = d.avul_time;
    analysis.avul_time_mean = d.avul_time_mean;
    
    analysis.lobe_len = d.lobe_len;
    analysis.lobe_len_mean = d.lobe_len_mean;

end
