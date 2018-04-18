function fix_SNR_bbrain(nii_fname,xyz_cut,nii_fname_out)
%function fix_SNR_bbrain(nii_fname,xyz_cut,nii_fname_out)
%Goal: To zero every voxel within the xyz_limits denoted in xyz_cut and
%save the output in nii_fname_out. 
%E.g. fix_SNR_bbrain('T1.nii',[XMINLIM XMAXLIM YMINLIM YMAXLIM ZLMINIM ZMAXLIM,'fixed_T1.nii' ]
%E.g. fix_SNR_bbrain('T1.nii',[18201 107 255 33 178 ,'fixed_T1.nii' ] -->
%                                
% ** FSLVIEW viewer yield xyz-1 locations (as it starts in 0 instead of 1)
% Created by Rodrigo Perea

if nargin <2
    error('Two argument are mandatory. Type help');
end

if nargin < 3
    nii_fname_out = ['fixed_SNR_' nii_fname ];
end

%Opening the file:
[vvol ,  hh ] = openIMG(nii_fname);

%Verify that it doesn't exist:
if exist(nii_fname_out,'file') ~= 0 
    error(['Filename: ' nii_fname_out ' exists. Exiting now to avoid replacing the file'])
end

AA=1;

for ix=1:hh.dim(1)
    if (ix < xyz_cut(1) || ix > xyz_cut(2))
        vvol(ix,:,:) = 0 ;
    end
end

for iy=1:hh.dim(2)
    if (iy < xyz_cut(3) || iy > xyz_cut(4))
        vvol(:,iy,:) = 0 ;
    end
end

for iz=1:hh.dim(3)
    if (iz < xyz_cut(5) || iz > xyz_cut(6))
        vvol(:,:,iz) = 0 ;
    end
end

hh.fname = nii_fname_out;
fprintf(['Writing out to: ' hh.fname  ]);
spm_write_vol(hh,vvol);
remove_nans_spm(nii_fname_out);
fprintf('...done\n')

end