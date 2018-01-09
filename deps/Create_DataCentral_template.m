%this is a template to create Tables in data central


FX_fields = fields(obj.Trkland.fx.data);
%add fx to all FX_fields
for ii=1:numel(FX_fields)
    FX_fields{ii} =  ['fx_' FX_fields{ii} ] ; 
end
CING_fields = fields(obj.Trkland.cingulum.data);
%add 'cing' ...
for ii=1:numel(CING_fields)
    CING_fields{ii} =  ['cing_' CING_fields{ii} ] ; 
end
HIPPOCING_fields=fields(obj.Trkland.hippocing.data);
%add 'hippocing' ...
for ii=1:numel(HIPPOCING_fields)
    HIPPOCING_fields{ii} =  ['hippocing_' HIPPOCING_fields{ii} ] ; 
end


all_fields = vertcat(FX_fields,CING_fields, HIPPOCING_fields);
%%


%Creating the DataCentral Query
clear DCTL_CMD
DC_INIT = ['''CREATE TABLE rdp20.TRKLAND ( MRI_Session_ID INT PRIMARY KEY, '];
DC_COLS = ''; 
for ii=1:numel(all_fields)
    %Volumes should be INT
    if ~isempty(strfind(all_fields{ii},'vol'))
        DC_COLS   = [ DC_COLS ' ' all_fields{ii} ' INT, '] ;
        %FA and non_FA values are FLOATS
    else
        DC_COLS = [ DC_COLS ' ' all_fields{ii} ' FLOAT, '] ;
%     else 
%         error(['I am not sure what field should be: ' all_fields{ii} ] );
    end
end
%Getting rid of the last column
DC_COLS = DC_COLS(1:end-2);

%DC_COLS = [ DC_COLS ' PRIMARY KEY (MRI_Session_ID) )'  ] 
END = ')''';

DCTL_CMD = [ DC_INIT DC_COLS  END ] 

