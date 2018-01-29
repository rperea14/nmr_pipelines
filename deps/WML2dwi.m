function WML2dwi(SessionID,projectID,OUTPUT_dir, WML_dir)
%output_DIR = function WML2dwi(SessionID,projectID)
%The goal is to coregister white matter lesions (WMLs) developed by the LST
%algorithm (http://www.applied-statistics.de/lst.html) into diffusion
%imaging space. 
%
% INPUT:    SessionID: ID of the participants (e.g. '150401_8CSAD00009')
%           projectID: ID of the project (e.g. 'ADRC')
%           output_DIR: location where the data will be outputted. 
%           WML_dir: (optional): directory where the WML information exist 
%
% EXAMPLE TO RUN: WML2dwi('150401_8CSAD00009','ADRC','/autofs/eris/bang/ADRC/PROJECTS/IMGS_WMLs_Tracts')
if nargin <3
    error(' Not enough arguments. Type ''help WML2dwi'' ');
end


if nargin <4
    if strcmp(projectID,'ADRC')
        WML_dir='/eris/bang/ADRC/PROJECTS/WMLs_LST_LGA/FLAIRS/'
    elseif strcmp(projectID,'HAB')
        WML_dir='/cluster/hab/HAB/Project1/wmh_LST_LPA/'
        display(['WML_dir is now: ' WML_dir 'but need to respecify...']);
        error('HAB project have crosssectional and longitudinal WML processes separately. Please specify')
    else
        error(['projectID (2nd argument): ' projectID ' has not been implemented'])
    end
else
    display('Make sure you append the las ''/ in the directory ');
    pause(1);
end



%INITIALIZE VARIABLES
if strcmp(OUTPUT_dir(end),'/') %check if the slash is present...
    OUTPUT_dir=[OUTPUT_dir SessionID filesep ] ;
else
    OUTPUT_dir=[OUTPUT_dir filesep SessionID filesep ] ;
end
system(['mkdir -p ' OUTPUT_dir ]);


in_T1_fname=[WML_dir 'm' SessionID '_FLAIR.nii']; % ',' stands for bias correction image
in_WML_probmap=[WML_dir 'ples_lpa_m' SessionID '_FLAIR.nii'];

out_T1_fname=[OUTPUT_dir 'm' SessionID '_FLAIR.nii']; % ',' stands for bias correction image
r_T1_fname=[OUTPUT_dir 'r_m' SessionID '_FLAIR.nii'];
out_WML_probmap=[OUTPUT_dir 'ples_lpa_m' SessionID '_FLAIR.nii'];
r_WML_probmap=[OUTPUT_dir 'r_ples_lpa_m' SessionID '_FLAIR.nii'];


if strcmp(projectID,'ADRC')
    in_b0=['/eris/bang/ADRC/Sessions/' SessionID '/DWIs/06_CoRegDWIs/combined_preproc_b0.nii.gz'];
    tmp_out_b0=[OUTPUT_dir 'b0_' SessionID '.nii.gz'];
    out_b0=[OUTPUT_dir 'b0_' SessionID '.nii'];
elseif strcmp(projectID,'HAB')
    error('in_bo has not been implemented to HAB dataset yet...')    
end


%Checking if 'in' files exist:
fprintf('\nChecking ''in'' file...');
if exist(in_T1_fname,'file')
    display('in_T1_fname exists!')
end

if exist(in_WML_probmap,'file')
    display('in_WML_probmap exists!')
end

if exist(in_b0,'file')
    display('in_b0 exists!')
end
fprintf('done checking\n');

%Checking if 'out' files exist:
fprintf('\nChecking ''out'' file...');
if exist(out_T1_fname,'file')
    display('out_T1_fname exists! Skipping copy!');
else
    display('copying T1_fname...');
    system(['cp ' in_T1_fname ' ' out_T1_fname]);
    fprintf('done.\n');
end

if exist(out_WML_probmap,'file')
    display('out_WML_probmap exists! Skipping copy');
else
    display('copying WML_probmap...');
    system(['cp ' in_WML_probmap ' ' out_WML_probmap]);
    fprintf('done.\n');
end

if exist(out_b0,'file')
    display('out_b0 exists! Skipping copy');
else
    display('copying b0...');
    system(['cp ' in_b0 ' ' tmp_out_b0]);
    system(['gunzip ' tmp_out_b0  ]);
    fprintf('done.\n');
end
fprintf('done checking \n');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IMPLEMENTATION STARTING NOW...
%QuickCoReg_rdp(out_b0,out_T1_fname,{out_WML_probmap})
if ~exist(r_WML_probmap)
    QuickCoReg_rdp(out_T1_fname,out_b0,{out_WML_probmap})
    %Removing nans:
    rm_spm_nans(r_T1_fname);
    rm_spm_nans(r_WML_probmap);
else
   display([ r_WML_probmap ' exists. Skipping...']) 
end




