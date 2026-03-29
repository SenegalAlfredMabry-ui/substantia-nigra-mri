%%
%%This is the 4th of 7 necessary steps in producing scaled neuromelanin
%%intensity data **for the substantia nigra**
%%Written by EBR 2024


clear all
close all
%corrects by z-scoring within-slice in the pons mask only

pathprompt='Please enter local path to TSE folder. Including leading and trailing slash: ';
localpath=getenv('BLACK_GOLIATH_TSE');
if isempty(localpath)
    localpath='';
end
localpath=convertCharsToStrings(localpath);
if strlength(localpath)==0
    localpath=convertCharsToStrings('${MOTIP_ROOT:-/path/to/MOTIP2018}/MRI_data/TSE/SN/');
end
pathitems=split(localpath,'/'); %split the path so you can derive all the others
addpath(localpath)

%build a path to the age data file
agefolders=pathitems;
agefolders(end-3)='';
agefolders(end-2)='';
agefolders(end-1)='';
agepath=join(agefolders,'/');
addpath(agepath)
load([char(agepath),'age_and_subject_number_20240704.mat'])

%build a path to the ntools folder
ntoolsfolders=pathitems;
ntoolsfolders(end-2)='scripts';
ntoolsfolders(end-1)='tse';
ntoolspath=join(ntoolsfolders,'/');

%build a path to original TSe data
tsefolders=pathitems;
tsefolders(end-1)='';
tsepath=join(tsefolders,'/');

cutoff1=40;
cutoff2=70;

young=[];
middle=[];
old=[];
for i=1:length(age_and_subject_number(:,1))
    subject=age_and_subject_number(i,1);
    age=age_and_subject_number(i,2);
    if age<cutoff1
        if isfile(strcat(tsepath,'zeropad_tse_',num2str(subject),'.nii'))
            young=[young subject];
        end
    elseif age>=cutoff1 && age<cutoff2
        if isfile(strcat(tsepath,'zeropad_tse_',num2str(subject),'.nii'))
            middle=[middle subject];
        end
    elseif age>=cutoff2
        if isfile(strcat(tsepath,'zeropad_tse_',num2str(subject),'.nii'))
            old=[old subject];
        end
    end
end

override = strtrim(getenv('BLACK_GOLIATH_SUBJECTS'));
if strlength(override) > 0
    tokens = strsplit(override);
    cleanup_subjects = [];
    for idx = 1:numel(tokens)
        token = tokens{idx};
        match = regexp(token, '(\d+)$', 'tokens', 'once');
        if ~isempty(match)
            cleanup_subjects(end+1) = str2double(match{1}); %#ok<AGROW>
        end
    end
    cleanup_subjects = unique(cleanup_subjects);
    cleanup_subjects(isnan(cleanup_subjects)) = [];
    if isempty(cleanup_subjects)
        cleanup_subjects = [530];
    else
        fprintf('Black Goliath restricted to subjects: %s\n', strjoin(string(cleanup_subjects)));
    end
else
    cleanup_subjects=[530];
end
everybody=cleanup_subjects;

%%

for i=1:length(everybody)
    subject=everybody(i);
    addpath(strcat(ntoolspath,'NIfTI_20140122'))

    %if isfile(strcat(tsepath,'corrected_tse_SN_',num2str(subject),'.nii'))
   %disp(['The corrected tse for SN already exists for ',num2str(subject),', skipping...'])
    %else
    
    subjnm=niftiread(strcat(tsepath,'zeropad_tse_',num2str(subject),'.nii'));
    nm=double(subjnm); %actual neuromelanin data for this person
    nminfo=niftiinfo(strcat(tsepath,'zeropad_tse_',num2str(subject),'.nii'));
    try
        subjmask=niftiread(strcat(localpath,'DC_mask_inTSE_',num2str(subject),'.nii'));
        mask=double(subjmask); %can't do these transformations on int16
        mask=rdivide(mask,128); %brainstem and 4th ventricle mask, loaded
    
        DC=nm.*mask;
        DC(DC==0)=NaN; %pull out neuromelanin data within the brainstem only
    
%   figure() %use this if you need to check on progress
%   image(pons(:,:,5)) %can also be moved around the script
        overallmean=nanmean(nanmean(nanmean(DC)));
        % checker=[];
        % rechecker=[];
        for j=1:size(DC,3)
            slice=DC(:,:,j); %loop through brainstem/4v by slice

            if length(find(~isnan(slice)))>0

                %checker=[checker nanmean(nanmean(slice))]; %see what the average value is
                slice=slice-overallmean; %NORMALIZE THE SLICE
                DC(:,:,j)=slice; %re-insert the normalized slice
                %rechecker=[rechecker nanmean(nanmean(slice))]; %see what the average value is now            end
            end
        end
        DC(isnan(DC))=0; %return NaN values to 0 because NaNs aren't handled in NIFTI
        subjnm=DC; %return corrected DC to its file

        
        niftiwrite(int16(subjnm),strcat(localpath,'corrected_tse_DC_',num2str(subject),'.nii'),nminfo)
        disp(['Saved output for ',num2str(subject)])

        catch error
        warning(strcat('Failed on subject ',num2str(subject)))
        fprintf(2,'Error!!',error.message)
        %rethrow(error)
    end
    end
%end %put this back when you want your else loop again
