function spm_remove_nans(f_name)
     display(['Removing NaNs for: ' f_name ])
     V = spm_vol(f_name);
     Y = spm_read_vols(V);
     Y(isnan(Y))  = 0 ;
     spm_write_vol(V,Y);
     fprintf(' ...done\n');
