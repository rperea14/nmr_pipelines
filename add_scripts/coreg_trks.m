
fileID=fopen('./thelist.list');
SUBJID=textscan(fileID,'%s');
SUBJID=SUBJID{1};

for ii=1:1 %numel(SUBJID)
    display (['In: ' SUBJID{ii}]);
    if exist(['./' SUBJID{ii}   '/coreg_DTI.trk'], 'file') == 0
        rotrk_copyheader(['./' SUBJID{ii}   '/DTI.trk'],['./' SUBJID{ii}  '/dwi_aparc.a2009+aseg.nii'],['./' SUBJID{ii}   '/coreg_DTI.trk']);
    end
    if exist(['./' SUBJID{ii}   '/coreg_GQI.trk'], 'file') == 0
        rotrk_copyheader(['./' SUBJID{ii}   '/QBI.trk'],['./' SUBJID{ii}  '/dwi_aparc.a2009+aseg.nii'],['./' SUBJID{ii}   '/coreg_GQI.trk']);
    end
    if exist(['./' SUBJID{ii}   '/coreg_DTI.trk'], 'file') == 0
        rotrk_copyheader(['./' SUBJID{ii}   '/QBI.trk'],['./' SUBJID{ii}  '/dwi_aparc.a2009+aseg.nii'],['./' SUBJID{ii}   '/coreg_QBI.trk']);
    end
    
end