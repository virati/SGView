%%
%Written by: Vineet Tiruvadi 2013-2016
%SG Viewer
%Main Gui for rapid viewing and spectral figures of BrainRadio and hdEEG
%data for MaybergLab.

function varargout = SG_view(varargin)
% SG_VIEW MATLAB code for SG_view.fig
%      SG_VIEW, by itself, creates a new SG_VIEW or raises the existing
%      singleton*.DBS901-Chronic-5day_perday.fig
%
%      H = SG_VIEW returns the handle to a new SG_VIEW or the handle to
%      the existing singleton*.
%
%      SG_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SG_VIEW.M with the given input arguments.
%
%      SG_VIEW('Property','Value',...) creates a new SG_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SG_view_OpeningFcn gets called.  Anhttp://proxy.library.emory.edu/login?url=
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SG_view_OpeningFcn via varargin.
% 
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SG_view

% Last Modified by GUIDE v2.5 06-Jun-2017 18:02:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SG_view_OpeningFcn, ...
                   'gui_OutputFcn',  @SG_view_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SG_view is made visible.
function SG_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SG_view (see VARARGIN)

% Choose default command line output for SG_view
handles.output = hObject;
handles.svInts.num = 1;

handles.interval_colors = {'k','b','m','c','r','g','y','k'};

%Hardload settings
handles.LFP_HardLoad = '';
handles.EEG_HardLoad = '';

%General directory settings
handles.LFP_dir = '~/local_data/BR';
handles.EEG_dir = '~/local_data/hdEEG';
handles.LFP_curr_dir = handles.LFP_dir;

handles.welchPSD = [];
handles.welchBounds = [];
handles.welchColor = {};
handles.pmap = colormap('lines');
colormap('default');

%Figure variables
handles.tdFig = 0;

%What should every modality be sampled to? This INTERPOLATES, need to make
%sure at every step that this is ok
handles.DATA.common_Fs = 422;

handles.BLP_pow = [];

handles.brLFP.Fs = str2num(get(handles.edtLFP_Fs,'String'));
handles.hdEEG.Fs = str2num(get(handles.edtEEG_Fs,'String'));

%% NEW STRUCTURES

%Assume 260 channels, set up the UI channel list
for ii = 1:260
    handles.DATA.UI_TS_List{ii} = [num2str(ii) ': Empty'];
end
handles.DATA.TS_List = zeros(260,1);

%variable for address of current active channel
%Start at 0, so throws exception if not properly handled after
%initialization
handles.UI.active_chann = 0;

%%
%Analysis tracking
    handles.SG.SG_Labels = {};
    handles.SG.Done = zeros(260,1);
    handles.ANALYSIS.ds_factor = 1;
    
    %Bigger one for the bar plots across files/patients
    handles.ANALYSIS.Aggr.idx = 0;
%%

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SG_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function LoadEEG(hObject,eventdata,handles,eeg_fname,eeg_dir)
    handles.RAW.EEG = [];
    
    disp('Loading EEG...');
    
    %Linux
    handles.DATA.EEG.locs = readlocs('~/Dropbox/projects/Research/MDD-DBS/Ephys/hdEEG/GSN257.sfp');
    %Windows
    %handles.DATA.EEG.locs = readlocs('C:\Users\Vineet\Dropbox\projects\Research\MDD-DBS\Ephys\hdEEG\GSN257.sfp');
    
    handles.RAW.EEG.rawdata = load([eeg_dir eeg_fname{1}]);
    handles.RAW.EEG.Fs = handles.RAW.EEG.rawdata.EEGSamplingRate;
    
    %Display information in UI
    set(handles.edtEEGfname,'String',eeg_fname{1});
    set(handles.edtEEG_Fs,'String',num2str(handles.RAW.EEG.Fs));
    
    %Find the fieldnames
    list_fields = fieldnames(handles.RAW.EEG.rawdata);
    handles.DATA.EEG.ts = handles.RAW.EEG.rawdata.(list_fields{1});
    ts_size = size(handles.DATA.EEG.ts);
    avg_idx = ts_size(1)+1;
    
    %This will be expanded and incorporate a configuration file/viewer so
    %can quickly switch between configs/montages
    
    %channel_list = {[32],[37],[25],[18],[241],[244],[238],[234],[89],[88],[130],[142],[136],[135],[148],[157],[9],[45],[186],[132]};
    
    %Differential channels
    % Subset indeed
    handles.DATA.EEG.SubSet =  {[32,37],[25,18],[241,244],[238,234],[89,88],[130,142],[136,135],[148,157],[9,45],[186,132]};
    %handles.DATA.EEG.SubSetLabels = {};
    
    %Single channels
    %handles.DATA.EEG.SubSet = {[32],[37],[25],[18],[241],[244],[238],[234],[89],[88],[130],[142],[136],[135],[148],[157],[9],[45],[186],[132]};
    
    %Avg Reference
    %handles.DATA.EEG.SubSet = {[32,avg_idx],[37,avg_idx],[25,avg_idx],[18,avg_idx],[241,avg_idx],[244,avg_idx],[238,avg_idx],[234,avg_idx],[89,avg_idx],[88,avg_idx],[130,avg_idx],[142,avg_idx],[136,avg_idx],[135,avg_idx],[148,avg_idx],[157,avg_idx],[9,avg_idx],[45,avg_idx],[186,avg_idx],[132,avg_idx]};
        
    %Just load all
%     for ii = 1:257
%         handles.DATA.EEG.SubSet{ii} = [ii];
%     end
    %Compute the channel average reference
    handles.DATA.EEG.ts(end+1,:) = mean(handles.DATA.EEG.ts,1);
    
    %Load in the EEG channels from the Subset variable
    for ii = 1:length(handles.DATA.EEG.SubSet)
        eeg_chann = cell2mat(handles.DATA.EEG.SubSet(ii));
        if numel(eeg_chann) == 1
            handles.DATA.chann{ii+2}.ts = detrend(handles.DATA.EEG.ts(eeg_chann,:));
            handles.DATA.chann{ii+2}.label = ['EEG Channel '  num2str(eeg_chann)];
        else
            handles.DATA.chann{ii+2}.ts = detrend(handles.DATA.EEG.ts(eeg_chann(1),:) - handles.DATA.EEG.ts(eeg_chann(2),:));
            handles.DATA.chann{ii+2}.label = ['EEG Diff Channel ' num2str(eeg_chann(1)) ' : ' num2str(eeg_chann(2))];
        end
        handles.DATA.chann{ii+2}.Fs = handles.RAW.EEG.Fs;
        handles.DATA.chann{ii+2}.num = ii+2;
        
        %Setup all EEG channels to be shifted when choosing windows
        handles.DATA.chann{ii+2}.doShift = 1;
    end
    
    add_channs = length(handles.DATA.EEG.SubSet)-1
    handles.DATA.TS_List(3:3+add_channs) = 1;
    
  
    %Now do UI stuff, like populating the channel popdown for timedomain
    for ii = 1:length(handles.DATA.EEG.SubSet)
        handles.DATA.UI_TS_List{ii+2} = [num2str(ii+2) ':' handles.DATA.chann{ii+2}.label];
    end
    
    set(handles.popChannTs,'String',handles.DATA.UI_TS_List);

    disp('...EEG File Loaded');
    
    guidata(hObject,handles);
    

%Loading of data function, this is for the brainradio LFP, 
%SINGLE FILES ONLY/EXPERIMENTS
    
function LoadLFP(hObject, eventdata, handles, data_loc)
    
    handles.RAW.LFP = [];

    disp('Loading LFP...');
        
    %Load in data from data_loc(ation) structure argument
    raw_data = dlmread([data_loc.data_path data_loc.data_file]);
    handles.RAW.LFP.rawdata = raw_data(:,[1,3]);
    handles.RAW.LFP.file = data_loc;
    %!!CHANGEBACK TO 422
    handles.RAW.LFP.rawFs = str2num(get(handles.edtLFP_Fs,'String'));
    
    %Reflect changes in UI
    set(handles.editFname,'String',data_loc.data_file);
    set(handles.cmdAvgSpectr,'Enable','off');
    handles.LFP_curr_dir = data_loc.data_path;
    
    %Interpolation and channel timeseries structures
    %handles.CHANN.ts will ALWAYS sampled at 1000Hz; though, here, the sizes
    %may not be consistent...
    for ii = 1:2
        handles.DATA.chann{ii}.ts = resample(handles.RAW.LFP.rawdata(:,ii),handles.DATA.common_Fs,handles.RAW.LFP.rawFs);
        handles.DATA.chann{ii}.Fs = handles.RAW.LFP.rawFs * handles.DATA.common_Fs / handles.RAW.LFP.rawFs;
        
        handles.DATA.chann{ii}.doShift = 0;
        handles.DATA.chann{ii}.num = ii;
    end    
    
    handles.DATA.TS_List(1:2) = [1 1];
  
    handles.DATA.chann{1}.label = ['LFP Channel Left'];
    handles.DATA.chann{2}.label = ['LFP Channel Right'];
    
    handles.DATA.UI_TS_List{1} = [num2str(1) ':' handles.DATA.chann{1}.label];
    handles.DATA.UI_TS_List{2} = [num2str(2) ':' handles.DATA.chann{2}.label];
    
    
    
    set(handles.popChannTs,'String',handles.DATA.UI_TS_List);
    
    disp('...Raw LFP loaded');
    
    handles.UI.active_chann = 1;
    
    guidata(hObject,handles);



% --- Outputs from this function are returned to the command line.
function varargout = SG_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%This function is for LOADING THE CURRENT DATA
%Includes BR and (aspirationally) hdEEG data

% --- Executes on button press in cmdLoad.
function cmdLoad_Callback(hObject, eventdata, handles)
% hObject    handle to cmdLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB30816
% handles    structure with handles and user data (see GUIDATA)
    [data_file data_path] = uigetfile({'*.txt'},'Pick BR Data File',handles.LFP_curr_dir,'MultiSelect','On')
    
    
    data_info.data_file = data_file;
    data_info.data_path = data_path;
    
    LoadLFP(hObject,eventdata,handles,data_info);


%This button ALWAYS redoes analysis on all channels
%And resets analysis sets (display sets)
% --- Executes on button press in cmdLFP_TF.
% function cmdLFP_TF_Callback(hObject, eventdata, handles)
% % hObject    handle to cmdLFP_TF (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
%     disp('Generating [Spectrograms] for all channels...');
%     
%     ds_fact = str2num(get(handles.edtDownsample,'String'));
%     %Preprocessing steps here
%     
%     %Clear current downsampled signals and spectrograms
%     handles.brLFP.ds = [];
%     handles.SG = [];    
%     
%     handles.brLFP.ds(:,1) = decimate(handles.brLFP.rawdata{1}(1,:),ds_fact);
%     handles.brLFP.ds(:,2) = decimate(handles.brLFP.rawdata{1}(2,:),ds_fact);
%     handles.brLFP.aFs = handles.brLFP.Fs/ds_fact;
%     
%     handles.anls.data = {};
%     
%     if (get(handles.cbFilter,'Value'))
%         disp('Filtering...');
%         lpf_50 = lowpass_gen(handles.brLFP.aFs);
%         handles.anls.data{end}(1,:) = filtfilthd(lpf_50,handles.anls.data{end}(1,:));
%         handles.anls.data{end}(2,:) = filtfilthd(lpf_50,handles.anls.data{end}(2,:));
%     end
%     
%     sgWin = str2num(get(handles.edSgWin,'String'));
%     sgDispl = str2num(get(handles.edSgDispl,'String'));
%     sgNfft = str2num(get(handles.edSgNFFT,'String'));
%     
%     %Does the spectrogram analysis
%     for c_num = 1:2
%         [handles.CHANN.SG{c_num}.S handles.CHANN.SG{c_num}.F handles.CHANN.SG{c_num}.T] = spectrogram(handles.brLFP.ds(:,c_num),blackmanharris(sgWin),sgDispl,sgNfft,handles.brLFP.aFs);
%     end
% 
%     handles.SG.Method = 'PWelch STFT';
%     
%     %Which channel to plot here
%     colormap('jet');
%     achann = handles.UI.active_chann;
%     imagesc(handles.CHANN.SG{achann}.T, handles.CHANN.SG{achann}.F, 10*log10(abs(handles.CHANN.SG{achann}.S)),[-50 0]);set(gca,'YDir','normal');colorbar();
%     title(['Channel ' num2str(achann)]);
%     
%     disp('...Done [Spectrograms]');
%     
%     handles.run_mean = [];
%     handles.run_length = 0;
%     handles.run_mean_color = {};
%     
%     guidata(hObject,handles);
%        

% --- Executes during object creation, after setting all properties.
function editFname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edSgWin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSgWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edSgDispl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSgDispl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edSgNFFT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edSgNFFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edtIStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtIEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edtI2Start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtI2Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edtI2End_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtI2End (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%MAIN 1-D PLOTTING FUNCTION BELOW
% --- Executes on button press in cmdSelect.
%OLDOLDOLD
% function cmdSelect_Callback(hObject, eventdata, handles)
% % hObject    handle to cmdSelect (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% sig_cont = select_window(hObject,eventdata,handles);
% 
% PSD_estimate(handles,handles.UI.active_chann,sig_cont);
% 
% 
% colormap('jet');
% 
% axes(handles.axes1);
% iv_rect = getrect();
% 
% use_welch = 1;
% 
% clear l_hpsd r_hpsd;
% [itv,t_vect] = iv_to_ts(handles,iv_rect);
% 
% if use_welch
%     %Uses the Welch Estimate
%     [l_hpsd] = intrv_spectra(hObject,eventdata,handles,itv,'left');
%     [r_hpsd] = intrv_spectra(hObject,eventdata,handles,itv,'right');
% end
% 
% rectangle('Position',[iv_rect(1), 0, iv_rect(3), 100],'LineWidth',4,'EdgeColor',handles.pmap(handles.svInts.num+1,:));
% 
% PSD_plotting(hObject,eventdata,handles,l_hpsd,'Left');



function PSD_plotting(hObject,eventdata,handles,PSD,fig_title)
handles.svInts.num = handles.svInts.num + 1;

figure(handles.ag_spect);
%suptitle(fig_title);

% %Display spectrogram for each interval
% subplot(3,2,1);
% imagesc(PSD.sgT,PSD.sgF,10*log10(abs(PSD.sg)));
% set(gca,'YDir','normal');

%pwelch estimate
subplot(3,2,2:4);
hold on;
for ii = 1:2
    b(:,1,ii) = abs(10*log10(PSD.P(:,ii)) - 10*log10(PSD.PC(:,1,ii)));
    b(:,2,ii) = abs(10*log10(PSD.P(:,ii)) - 10*log10(PSD.PC(:,2,ii)));
end

%Need to do running matrix of all PSDs being plotted and all bs

handles.runPSD = [handles.runPSD, PSD.P];
handles.runBounds = cat(3,handles.runBounds,b);
handles.runColor{handles.svInts.num} = handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1};


pmap = handles.pmap;

boundedline(PSD.welchF,10*log10(handles.welchPSD),handles.welchBounds,'cmap',pmap,'alpha');
%boundedline(PSD.welchF,10*log10(PSD.welch),b,handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1},'alpha');

title('Welch Estimate');
xlabel('Freq (Hz)');ylabel('Power (dB)');
xlim([0 25]);

%Save to handles.spectra{associated index}
handles.PSD_estimates{handles.svInts.num}.LEFT = PSD.welch;
handles.PSD_estimates{handles.svInts.num}.Freqs = PSD.welchF;

% subplot(3,1,2);
% boundedline(PSD.sgF,10*log10(PSD.sg_m),1.96 * (10*log10(PSD.sg_s))/sqrt(PSD.sg_n),handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1},'alpha');hold on;
% 
% xlim([0 25]);
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% title('Spectrogram Average Estimate');


% subplot(3,1,3);
% b_mts(:,1) = abs(10*log10(LEFT.l_mts) - 10*log10(LEFT.l_mts_ci(1,:)'));
% b_mts(:,2) = abs(10*log10(LEFT.l_mts) - 10*log10(LEFT.l_mts_ci(1,:)'));
% 
% boundedline(LEFT.l_mts_f,10*log10(LEFT.l_mts),b_mts,handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1},'alpha');
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% xlim([0 50]);s
% title('Multitaper Estimate');


plot_ratio = 0;
subplot(3,2,5:6);
if plot_ratio
    %subplot(3,2,5:6);

    if handles.svInts.num >= 2
        plot(handles.PSD_estimates{1}.Freqs,10*log10(handles.PSD_estimates{2}.LEFT ./ handles.PSD_estimates{1}.LEFT));hold on;
    end
    xlim([0 25]);
    title('Log Power Ratio')
else
    %Plot band limited power
    alphalim = [8,14];
    thetalim = [4,8];
    
    alpha_mask = find(PSD.welchF > alphalim(1) & PSD.welchF < alphalim(2));
    theta_mask = find(PSD.welchF > thetalim(1) & PSD.welchF < thetalim(2));
    
    alpha_pow = mean(PSD.welch(alpha_mask));
    theta_pow = mean(PSD.welch(theta_mask));
    
    handles.Alpha_pow = [handles.Alpha_pow,alpha_pow];
    handles.Theta_pow = [handles.Theta_pow,theta_pow];
    bar_h = bar([handles.Alpha_pow;handles.Theta_pow]);
    for kk = 1:length(bar_h)
        set(bar_h(kk),'FaceColor',pmap(kk,:));
    end
    set(gca,'XTickLabel',{'Alpha','Theta'});
    
end

set(findall(gcf,'-property','FontSize'),'FontSize',24);
guidata(hObject,handles);        

function [nbounds,tvect] = iv_to_ts(handles,iv_rect,lim_time)
iv_start = iv_rect(1);
if lim_time == 0
    iv_end = iv_rect(1)+iv_rect(3);
else
    iv_end = iv_rect(1)+lim_time;
end

if iv_start < 0
    iv_start = 0;
end
if iv_end > (length(handles.DATA.chann{handles.UI.active_chann}.ts) / handles.DATA.chann{handles.UI.active_chann}.Fs)
    iv_end = length(handles.DATA.chann{handles.UI.active_chann}.ts) / handles.DATA.chann{handles.UI.active_chann}.Fs;
end

nbounds = [iv_start iv_end] * handles.DATA.chann{handles.UI.active_chann}.Fs + 1; %+1 to account for no 0 index
tvect = linspace(iv_start,iv_end,(iv_end - iv_start) *  handles.DATA.chann{handles.UI.active_chann}.Fs);


% --- Executes on button press in cmdSvInt.
function cmdSvInt_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSvInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.svInts.num = handles.svInts.num + 1;
handles.svInts.Int{handles.svInts.num}.file = get(handles.editFname,'String');
handles.svInts.Int{handles.svInts.num}.tInt = handles.curr_Int.t_Int;
handles.svInts.Int{handles.svInts.num}.label = get(handles.edtIntLabel,'String');

handles.svInts.Int{handles.svInts.num}.l_psd = handles.curr_Int.l_hpsd;
handles.svInts.Int{handles.svInts.num}.r_psd = handles.curr_Int.r_hpsd;

base_list = get(handles.lsbIvs,'String');
base_list = strvcat(char(base_list),handles.svInts.Int{handles.svInts.num}.label);

set(handles.lsbIvs,'String',cellstr(base_list));

guidata(hObject,handles);

% --- Executes on selection change in lsbIvs.
function lsbIvs_Callback(hObject, eventdata, handles)
% hObject    handle to lsbIvs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lsbIvs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lsbIvs


% --- Executes during object creation, after setting all properties.
function lsbIvs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lsbIvs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdLoadInt.
function cmdLoadInt_Callback(hObject, eventdata, handles)
% hObject    handle to cmdLoadInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ANALYSIS.active_figure = figure();

%Refresh SG main axes
%Plot the active chann
SG_Plot_Chann(hObject,handles,handles.UI.active_chann);


%Reset the interval count
handles.ANALYSIS.active_interval = 0;
handles.ANALYSIS.intv = {};
handles.ANALYSIS.Aggr.idx = handles.ANALYSIS.Aggr.idx + 1;

%
%handles.ag_spect = figure;
% handles.svInts.num = 0;
% handles.sel = 0;
% handles.freq_means = [];
% 
% handles.welchPSD = [];
% handles.welchBounds = [];
% handles.Alpha_pow = [];
% handles.Theta_pow = [];

guidata(hObject,handles);



function edtIntLabel_Callback(hObject, eventdata, handles)
% hObject    handle to edtIntLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtIntabel as text
%        str2double(get(hObject,'String')) returns contents of edtIntLabel as a double


% --- Executes during object creation, after setting all properties.
function edtIntLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIntLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbFilter.
function cbFilter_Callback(hObject, eventdata, handles)
% hObject    handle to cbFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbFilter


% --- Executes on button press in cmdTime.
%Displays the time series for the rectangle chosen
function cmdTime_Callback(hObject, eventdata, handles)
% hObject    handle to cmdTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
iv_rect = getrect();

iv_start = iv_rect(1);
iv_end = iv_rect(1)+iv_rect(3);

if iv_start < 0
    iv_start = 0;
end

if iv_end > (length(handles.brLFP.ds(:,1)) / handles.brLFP.aFs)
    iv_end = length(handles.brLFP.ds(:,1)) / handles.brLFP.aFs;
end

iv = [iv_start iv_end] * handles.brLFP.aFs + 1;

t_vect = linspace(iv_start,iv_end,(iv_end - iv_start) *  handles.brLFP.aFs + 1);


t_sig = handles.brLFP.ds(iv(1):iv(2),:);


figure;
subplot(3,1,1);
plot(t_vect,t_sig(:,1)); hold on;
plot(t_vect,t_sig(:,2),'r');



set(findall(gcf,'-property','FontSize'),'FontSize',24);
%set(findall(gcf,'-property','FontName'), 'Helvetica');



subplot(3,1,3);
[acor, lag] = xcorr(t_sig(:,1),t_sig(:,2),'coeff');
plot(lag,acor);

write_file = 1;

if write_file
    dlmwrite('/tmp/left_time_sig.txt',[t_vect;t_sig(:,1)'],'delimiter',',');
    dlmwrite('/tmp/right_time_sig.txt',[t_vect;t_sig(:,2)'],'delimiter',',');
end


% --- Executes on button press in cmdSvBL.
function cmdSvBL_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSvBL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.baselineInt.left = handles.curr_Int.left;
handles.baselineInt.right = handles.curr_Int.right;
guidata(hObject,handles);

function edtIntLength_Callback(hObject, eventdata, handles)
% hObject    handle to edtIntLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtIntLength as text
%        str2double(get(hObject,'String')) returns contents of edtIntLength as a double


% --- Executes during object creation, after setting all properties.
function edtIntLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtIntLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtDownsample_Callback(hObject, eventdata, handles)
% hObject    handle to edtDownsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDownsample as text
%        str2double(get(hObject,'String')) returns contents of edtDownsample as a double


% --- Executes during object creation, after setting all properties.
function edtDownsample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDownsample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in cmdAvgSpectr.
function cmdAvgSpectr_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAvgSpectr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.svInts.num = handles.svInts.num + 1;
%Periodogram
for ii = 1:length(handles.anls.multi_data)
    %precondition signal
    sig_n = detrend(handles.anls.multi_data{ii}(1,:),'constant');
    %Take just the last part of sig_n, minus the initial transient
    sig_n = sig_n(end - 5*handles.brLFP.aFs:end);
    %[Lpxx(ii,:),Lf] = periodogram(handles.anls.multi_data{ii}(1,:),blackmanharris(length(handles.anls.multi_data{ii}(1,:))),1024*2,handles.brLFP.aFs);
    [Lpxx(ii,:),Lf] = pwelch(sig_n,blackmanharris(5 * handles.brLFP.aFs),[],2056,handles.brLFP.aFs);
    disp(num2str(size(sig_n)));
end

if ~isfield(handles,'avgSpect')
    handles.avgSpect = figure();
end

avgLp = mean(Lpxx,1);
%stdLp = std(10*log10(Lpxx),[],1);
stdLp = std(Lpxx,[],1);
rngLp = range(Lpxx,1);

figure(handles.avgSpect);

subplot(4,1,1);
plot(Lf,10*log10(avgLp),handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1});hold on;
xlim([0 50]);

subplot(4,1,2);
boundedline(Lf,10*log10(avgLp),2*(10*log10(stdLp))/sqrt(length(handles.anls.multi_data)),handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1},'alpha')%,handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1},'alpha');hold on;
xlim([0 50]);

%Subplot for theta, alpha, beta band limited
theta_mask = (Lf > 4 ) & (Lf < 8);
theta_band_spectr = avgLp .* theta_mask';
theta_mean = mean(theta_band_spectr(theta_mask ~= 0));

alpha_mask = (Lf > 8 ) & (Lf < 14);
alpha_band_spectr = avgLp .* alpha_mask';
alpha_mean = mean(alpha_band_spectr(alpha_mask ~= 0));

handles.run_length = handles.run_length + 1;
handles.run_mean(:,handles.run_length) = [theta_mean,alpha_mean]';
%handles.run_mean = [handles.run_mean';theta_mean,alpha_mean]';
handles.run_mean_color = [handles.run_mean_color,handles.interval_colors{mod(handles.svInts.num,length(handles.interval_colors))+1}];

set(findall(gcf,'-property','FontSize'),'FontSize',24);

guidata(hObject,handles);

% --- Executes on button press in cmdReload.
function cmdReload_Callback(hObject, eventdata, handles)
% hObject    handle to cmdReload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
      LoadLFP(hObject,eventdata,handles);
    guidata(hObject,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over cmdSelect.
function cmdSelect_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to cmdSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on cmdSelect and none of its controls.
function cmdSelect_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to cmdSelect (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) https://docs.google.com/presentation/d/1F5ppNZjxxtwawCghksNKbVTURH4TyB8eCgRIxSMwlfs/edit?usp=sharingof the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in SaveSGCmd.
function SaveSGCmd_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSGCmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Fig2 = figure;
copyobj(handles.axes1,Fig2);colormap('jet');

data_file = get(handles.editFname,'String');

%hgsave(Fig2, ['/tmp/' handles.data_file '_spectrogram.fig']);
pause(2);
print(Fig2,'-dpng',['/tmp/' data_file '_Chann' num2str(handles.UI.active_chann) '_spectrogram.png']);


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in cmdClearSG.
function cmdClearSG_Callback(hObject, eventdata, handles)
% hObject    handle to cmdClearSG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
cla;


% --- Executes on button press in cmdMeta.
function cmdMeta_Callback(hObject, eventdata, handles)
% hObject    handle to cmdMeta (see GCBO)%
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Textprompts for ease
pt = input('What patient is this?')
exp = input('What experiment is this?')
timept = input('What timepoint is this?')
condit = input('What condition is this?')
trial = input('What trial is this?')

base_dir = '/tmp/Results/SysID/';
out_dir = [base_dir '/' exp '/' timept '/' pt '/']
mkdir(out_dir)
output_file = [out_dir condit '_T' trial '.xml']


% --- Executes on button press in cmdPathClip.
function cmdPathClip_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPathClip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clipboard('copy',[handles.RAW.LFP.file.data_path handles.RAW.LFP.file.data_file]);


function SG_Plot_Chann(hObject,handles,achann)

axes(handles.axes1);cla;
colormap('jet');
%achann = handles.UI.active_chann;
colormap('jet');
imagesc(handles.TF.chann{achann}.T, handles.TF.chann{achann}.F, 10*log10(abs(handles.TF.chann{achann}.S)));
%caxis [-20,20] on the z-scored input signal seems like a great
%VISUALIZATION technique; keep this for the main window spectrogram/view
%But might need to be tweaked for the EEG data
set(gca,'YDir','normal');colorbar();ylim([0 handles.TF.Fs / 2]);caxis([-20,20]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['Channel ' num2str(achann) ' - ' handles.DATA.chann{achann}.label]);


% --- Executes during object creation, after setting all properties.
function popChann_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popChann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popTFMethod.
function popTFMethod_Callback(hObject, eventdata, handles)
% hObject    handle to popTFMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popTFMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popTFMethod
handles.TF_Method.idx = get(hObject,'Value');
temp_string = get(hObject,'String');
handles.TF_Method.name = temp_string{handles.TF_Method.idx};

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popTFMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popTFMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdEEGLoad.
function cmdEEGLoad_Callback(hObject, eventdata, handles)
% hObject    handle to cmdEEGLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [handles.EEG_file handles.EEG_path] = uigetfile({'*.mat'},'Pick EEG Data File',handles.EEG_dir,'MultiSelect','Off')
    handles.EEG_file = cellstr(handles.EEG_file);
    handles.EEG_dir = handles.EEG_path;
    
    eeg_fname = handles.EEG_file;eeg_dir = handles.EEG_path;
    
    guidata(hObject,handles);
    
    LoadEEG(hObject,eventdata,handles,eeg_fname,eeg_dir);
        
function edtEEGfname_Callback(hObject, eventdata, handles)
% hObject    handle to edtEEGfname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtEEGfname as text
%        str2double(get(hObject,'String')) returns contents of edtEEGfname as a double


% --- Executes during object creation, after setting all properties.
function edtEEGfname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtEEGfname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdPAC.
function cmdPAC_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPAC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Interval return function
axes(handles.axes1);
iv_rect = getrect();

if get(handles.chkLimIntv,'Value')
    lim_time = str2num(get(handles.edtIntLength,'String'));
else
    lim_time = 0;
end
[itv,t_vect] = iv_to_ts(handles,iv_rect,lim_time);

disp('Starting PAC');

PAC.sig = handles.brLFP.ds(itv(1):itv(2),:);

disp('Running PAC...');
%Do PAC Calculation
%CFC Section
a_range = 1:1:100; 
p_range = 1:1:25;
bw = 2; pbr = 0.02;
n = 36; Fb = 2; Fc = 1;
    

% %GPU version
% %GPU VERSION
% s1c = reshape(repmat(PAC.sig(:,1)',100,1),10,10,[]);
% s2c = reshape(repmat(PAC.sig(:,2)',100,1),10,10,[]);
% 
% S1 = gpuArray(PAC.sig(:,1)');
% S2 = gpuArray(PAC.sig(:,2)');


%Non-GPU version

CFC.MI{1}.MIs = GLMcomodulogram(PAC.sig(:,1)',PAC.sig(:,2)',a_range,p_range,handles.brLFP.aFs,bw,pbr,'No');
CFC.MI{2}.MIs = GLMcomodulogram(PAC.sig(:,2)',PAC.sig(:,1)',a_range,p_range,handles.brLFP.aFs,bw,pbr,'No');

CFC.MI{3}.MIs = GLMcomodulogram(PAC.sig(:,1)',PAC.sig(:,1)',a_range,p_range,handles.brLFP.aFs,bw,pbr,'No');
CFC.MI{4}.MIs = GLMcomodulogram(PAC.sig(:,2)',PAC.sig(:,2)',a_range,p_range,handles.brLFP.aFs,bw,pbr,'No');


handles.CFC = CFC;

guidata(hObject,handles);




figure;
colormap('jet');
subplot(2,2,1);
imagesc(p_range,a_range,handles.CFC.MI{1}.MIs',[0 0.8]);set(gca,'ydir','Normal');colorbar();
subplot(2,2,2);
imagesc(p_range,a_range,handles.CFC.MI{2}.MIs',[0 0.8]);set(gca,'ydir','Normal');colorbar();
subplot(2,2,3);
imagesc(p_range,a_range,handles.CFC.MI{3}.MIs',[0 0.8]);set(gca,'ydir','Normal');colorbar();
subplot(2,2,4);
imagesc(p_range,a_range,handles.CFC.MI{4}.MIs',[0 0.8]);set(gca,'ydir','Normal');colorbar();

disp('DONE PAC');





% --- Executes on selection change in popPSDEst.
function popPSDEst_Callback(hObject, eventdata, handles)
% hObject    handle to popPSDEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popPSDEst contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popPSDEst
handles.PSD_Method.idx = get(hObject,'Value');
temp_string = get(hObject,'String');
handles.PSD_Method.name = temp_string{handles.PSD_Method.idx};

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popPSDEst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popPSDEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 



function edtLFP_Fs_Callback(hObject, eventdata, handles)
% hObject    handle to edtLFP_Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtLFP_Fs as text
%        str2double(get(hObject,'String')) returns contents of edtLFP_Fs as a double
handles.brLFP.Fs = str2num(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edtLFP_Fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtLFP_Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtEEG_Fs_Callback(hObject, eventdata, handles)
% hObject    handle to edtEEG_Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtEEG_Fs as text
%        str2double(get(hObject,'String')) returns contents of edtEEG_Fs as a double
handles.hdEEG.Fs = str2num(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edtEEG_Fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtEEG_Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popPACc1.
function popPACc1_Callback(hObject, eventdata, handles)
% hObject    handle to popPACc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popPACc1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popPACc1


% --- Executes during object creation, after setting all properties.
function popPACc1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popPACc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popPACc2.
function popPACc2_Callback(hObject, eventdata, handles)
% hObject    handle to popPACc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popPACc2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popPACc2


% --- Executes during object creation, after setting all properties.
function popPACc2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popPACc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkLimIntv.
function chkLimIntv_Callback(hObject, eventdata, handles)
% hObject    handle to chkLimIntv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkLimIntv


% --- Executes on button press in cmdTimeFile.
function cmdTimeFile_Callback(hObject, eventdata, handles)
% hObject    handle to cmdTimeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%time_file = dlmread('/tmp/Andrea_imagined_movement/DBS906_emoLFPpilot_11092015.txt');
%time_file = dlmread('/tmp/Andrea_imagined_movement/DBS905_motorLFPpilot_11102015.txt');
meta_file = [handles.RAW.LFP.file.data_path handles.RAW.LFP.file.data_file(1:end-3) '.xml']

meta_data = dlmread(meta_file);



timings = time_file(:,2);
stim_art = time_file(:,3);

disp('Timings Loaded...');

% --- Executes on button press in boxResamp1K.
stim_art = time_file(:,3);


disp('Timings Loaded...');

% --- Executes on button press in boxResamp1K.
function boxResamp1K_Callback(hObject, eventdata, handles)
% hObject    handle to boxResamp1K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boxResamp1K


% --- Executes on selection change in popChannTs.
function popChannTs_Callback(hObject, eventdata, handles)
% hObject    handle to popChannTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popChannTs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popChannTs

%When this value changes, update the current active channel
handles.UI.active_chann = get(hObject,'Value');

disp(['Changing Active Channel to ' num2str(handles.UI.active_chann)]);


guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popChannTs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popChannTs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cmdTsShow.
function cmdTsShow_Callback(hObject, eventdata, handles)
% hObject    handle to cmdTsShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curr_chann = get(handles.popChannTs,'Value');
figure(100);
subplot(3,1,1);
plot(handles.DATA.chann{curr_chann}.ts);hold on;
title(['Channel ' num2str(curr_chann) ' :: Time Domain']);legend();
subplot(3,1,2);
plot(handles.RAW.LFP.rawdata(:,curr_chann));hold on;
title(['Raw Channel ' num2str(curr_chann) ' at ' num2str(handles.RAW.LFP.rawFs) ' :: T Domain']);
%subplot(3,1,3);
%spectrogram(handles.RAW.LFP.rawdata(:,curr_chann),blackmanharris(512),500,2^10,handles.RAW.LFP.rawFs);
%title(['Spectrogram of RAW Channel ' num2str(curr_chann) ' :: TF Domain']);legend();
guidata(hObject,handles);


% --- Executes on button press in cmdTFPlot.
function cmdTFPlot_Callback(hObject, eventdata, handles)
% hObject    handle to cmdTFPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = find(handles.DATA.TS_List == 1)';

disp('Computing T-F for Loaded Channels...');

handles.SG.Done = [];
handles.SG.SG_Labels = {};

%Make a local container for the raw data
chann_ts{1} = [];
chann_fs = handles.DATA.chann{min(a)}.Fs;
ds_fact = 1;

for ii = a
    dummy_chann = handles.DATA.chann{ii};
    chann_ts{ii}.ts = decimate(dummy_chann.ts,ds_fact);
end
chann_fs = chann_fs ./ ds_fact;

%Without downsampling ahead of time
for ii = a
    [handles.TF.chann{ii}.S, handles.TF.chann{ii}.F, handles.TF.chann{ii}.T] = spectrogram(chann_ts{ii}.ts,blackmanharris(512),500,2^10,422)%chann_fs); 
end
%Update the decimated Fs
handles.TF.Fs = chann_fs;

%Set UI active channel to lowest available channel
handles.UI.active_chann = min(a);

%Plot the active chann
SG_Plot_Chann(hObject,handles,handles.UI.active_chann);

%Update the popdown text
for ii = a
    SG_Labels{ii} = ['Channel ' num2str(ii)];
    handles.SG.SG_Labels = SG_Labels;
    handles.SG.Done = [handles.SG.Done ii];
end

set(handles.popChann,'String',SG_Labels);

disp('...Done computing T-F for Channel Set');
guidata(hObject,handles);




% --- Executes on button press in cmdCoher.
function cmdCoher_Callback(hObject, eventdata, handles)
% hObject    handle to cmdCoher (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in popChann.
function popChann_Callback(hObject, eventdata, handles)
% hObject    handle to popChann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popChann contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popChann
handles.UI.active_chann = get(hObject,'Value');
%Is value gotten in the current list of SG'ed channels?
guidata(hObject, handles);

if handles.SG.Done(handles.UI.active_chann)
    SG_Plot_Chann(hObject,handles,handles.UI.active_chann);
else
    disp('That Channel Was Not Computed');
end


%Function that finds the stim periods in ALL CHANNELS, and aligns them
%Do this on raw data, sampled at highest Fs

% --- Executes on button press in cmdStimAlign.
function cmdStimAlign_Callback(hObject, eventdata, handles)
% hObject    handle to cmdStimAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Filter all channels at 130 Hz
loaded_channs = find(handles.DATA.TS_List == 1)';

d = fdesign.bandpass('N,F3dB1,F3dB2',20,129,131,1000);
Hd = design(d,'butter');
fvtool(Hd);

for ii = loaded_channs
    %filted_sig{ii}.stim_sig = eegfilt(handles.DATA.chann{ii}.ts,handles.DATA.chann{ii}.Fs,129,130);
    a{ii}.ts = filtfilthd(Hd,handles.DATA.chann{ii}.ts);
end

disp('...Done Aligning');

function intv = select_window(hObject,eventdata,handles)

axes(handles.axes1);
iv_rect = getrect();
%Bounds should be the boundary TIMES
[intv.nbounds,intv.tvect] = iv_to_ts(handles,iv_rect,0); %0 is for the lim time ability, to keep fixed window length
intv.fromChannel = handles.UI.active_chann;


% --- Executes on button press in cmdChannPCA.
function cmdChannPCA_Callback(hObject, eventdata, handles)
% hObject    handle to cmdChannPCA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Do PCA on all channels
%Choose your epoch
curr_epoch = select_window(hObject,eventdata,handles)

epoch_data = extractEpoch(handles,curr_epoch);
ChannBand(handles,epoch_data);

%Do a banding out

%Just do raw PCA
[coeff,score,latent] = pca(epoch_data,'Algorithm','eig');
%or do PCA on banded data

% 
% figure;
% subplot(3,2,1);
% plot(score);
% subplot(3,2,2);
% scatter(epoch_data(:,1),epoch_data(:,2),'.');
% 
% subplot(3,2,3);
% hist(epoch_data(:,1));
% subplot(3,2,4);
% hist(epoch_data(:,2));
% 
% subplot(3,2,5);
% [s,f,t] = spectrogram(score(:,1),blackmanharris(512),500,2^10,1000);
% imagesc(t,f,10*log10(abs(s)));colormap('jet');set(gca,'YDir','normal');
% subplot(3,2,6);
% [s,f,t] = spectrogram(score(:,2),blackmanharris(512),500,2^10,1000);
% imagesc(t,f,10*log10(abs(s)));colormap('jet');set(gca,'YDir','normal');
% 
% disp('Done with time-series PCA');

function ChannBand(handles,epoch_data)
%Break into bands
band_epoch_data = extractBands(handles,epoch_data);

%find timepoints
bed_size = size(band_epoch_data);

%Reshape so it's timex(channelxbands)
PCA_input = reshape(permute(band_epoch_data,[2,3,1]),[],bed_size(1))'; %YES! this reshapes this matrix properly, as you would expect
PCA_input = abs(PCA_input)

[coeff,score,latent] = pca(PCA_input,'Algorithm','eig');

score = real(score)

figure;
subplot(2,2,1);
%plot(abs(PCA_input));
imagesc(abs(coeff))
subplot(2,2,2);
%plot(abs(score));
plot(latent)

subplot(2,2,3);
scatter3(abs(PCA_input(:,1)),abs(PCA_input(:,2)),abs(PCA_input(:,3)),'.');
subplot(2,2,4);
scatter3(score(:,1),score(:,2),score(:,3),'.');

%subplot(3,2,3);
%hist(PCA_input(:,1));
%subplot(3,2,4);
%hist(PCA_input(:,2));

%this stuff is a bit meaningless when banded out; why take spectrogram of
%banded power; ehh, there are reasons, but not right now
%subplot(3,2,5);
%[s,f,t] = spectrogram(score(:,1),blackmanharris(512),500,2^10,422);
%imagesc(t,f,10*log10(abs(s)));colormap('jet');set(gca,'YDir','normal');
%subplot(3,2,6);
%[s,f,t] = spectrogram(score(:,2),blackmanharris(512),500,2^10,422);
%imagesc(t,f,10*log10(abs(s)));colormap('jet');set(gca,'YDir','normal');

disp('Done with time-series PCA');

% figure;
% subplot(3,1,1);
% imagesc(abs(coeff));
% subplot(3,1,2);
% plot(score);

% %Do covariance matrix
% 
% %Plotting, bullshit right now
% figure;
% subplot(4,2,1);plot(epoch_data);
% subplot(4,2,3);imagesc(coeff);
% subplot(4,2,5);plot(latent);
% topovect = linspace(1,260,260);
% subplot(4,2,7);topoplot(topovect,handles.DATA.EEG.locs);
% 
% %Plot the head and EEG channels
% %figure;
% %topoplot([],handles.DATA.EEG.locs,'style','blank','electrodes','labelpoint');
% %disp('PCA Analysis Done...');
% 
% %Blindly plot the first channel theta power vs third channel theta power
% figure;title('Theta Power Space');
% %scatter3(PCA_input(:,5*(1-1) + 2),PCA_input(:,5*(3-1) + 2),PCA_input(:,5*(5-1) + 2));hold on;
% line(PCA_input(:,5*(1-1) + 2),PCA_input(:,5*(3-1) + 4),PCA_input(:,5*(5-1) + 4));
% xlabel('Channel 1 Theta');
% ylabel('Channel 3 Theta');


function banded_data = extractSGBands(handles,epoch_data)
%Break into traditional bands
osc_bands.F = {[1,4],[4,8],[8,14],[15,30],[30,50],[1,20]};
osc_bands.name = {'Delta','Theta','Alpha','Beta','Broad Gamma','Norm'};

%Dimensionality of the data
%Row is going to be observation, column is channel
data_dim = size(epoch_data);

for ii = 1:data_dim(2)
    %For each channel in the above
    %Do multitaper for spectrogram
        
    %Do pwelch spectrogram
    [S,F,T] = spectrogram(epoch_data(:,ii),blackmanharris(256),250,2^10,1000); %!!! Fs hardcoded here, change this
        
    for jj = 1:length(osc_bands.name)
        %For each band
        %banded_data;
        bandlims = osc_bands.F{jj};
        %Log transform HERE
        %banded_data(:,jj,ii) = mean(10*log10(abs(S(F > bandlims(1) & F < bandlims(2),:))));        
        %OR Save log transform for viz
        banded_data(:,jj,ii) = sum(abs(S(F > bandlims(1) & F < bandlims(2),:)));      
    end
end

function hilb_data = extractBands(handles,epoch_data)
%Break into traditional bands
osc_bands.F = {[1,4],[4,8],[8,14],[15,30],[30,50],[1,20]};
osc_bands.name = {'Delta','Theta','Alpha','Beta','Broad Gamma','Norm'};

%Dimensionality of the data
%Row is going to be observation, column is channel
data_dim = size(epoch_data);

for ii = 1:length(osc_bands.name)
    f(ii) = fdesign.bandpass('N,Fc1,Fc2',100,osc_bands.F{ii}(1),osc_bands.F{ii}(2),1000);
    Fd(ii) = design(f(ii),'butter');    
end

for ii = 1:data_dim(2)
    %For each channel in the above
    %Do multitaper for spectrogram
    
    %Do FILTERING approach
    for jj = 1:length(osc_bands.name)
        %Make the filter
        fsig(:,jj,ii) = filtfilthd(Fd(jj),epoch_data(:,ii));        
        %banded_data(:,jj,ii) = fsig(:,jj,ii).^2;
        hilb_data(:,jj,ii) = hilbert(fsig(:,jj,ii));
    end
end


%Do any band normalization here
disp('Done Banding...');

function epochts = extractEpoch(handles,curr_epoch)
%Go through all loaded channels and extract the needed epoch
for ii = 1:length(handles.DATA.chann)
    [aligned_epoch,achann] = align_epoch(handles,curr_epoch, handles.DATA.chann{ii});
    
    achann.ds = decimate(achann.ts,handles.ANALYSIS.ds_factor);
    %Now take the the interval chunk we wanted
    achann.intv{1}.ds = achann.ds(round(aligned_epoch.ivn));
    achann.intv{1}.tvect = linspace(aligned_epoch.tvect(1),aligned_epoch.tvect(end),length(achann.intv{1}.ds));

    epochts(:,ii) = achann.intv{1}.ds;
end

function [out_epoch,achann] = align_epoch(handles,in_epoch,achann)
    %Grab active channel and decimate it
    
    achann.dsFs = achann.Fs ./ handles.ANALYSIS.ds_factor;
    
    out_epoch = in_epoch;
    out_epoch.ivn = (in_epoch.tvect(1)*achann.dsFs) : (in_epoch.tvect(end)*achann.dsFs);
    if in_epoch.fromChannel <= 2
        %window is in the reference of the LFPs
        if achann.num > 2
            out_epoch.ivn = (in_epoch.tvect(1)*achann.dsFs - handles.DATA.lag_val) : (in_epoch.tvect(end)*achann.dsFs - handles.DATA.lag_val);
        end
    else
        %Window is in the reference of the EEGs
        if achann.num <= 2
            out_epoch.ivn = (in_epoch.tvect(1)*achann.dsFs + handles.DATA.lag_val) : (in_epoch.tvect(end)*achann.dsFs + handles.DATA.lag_val);
        end
    end



% --- Executes on button press in cmdOscModel.
function cmdOscModel_Callback(hObject, eventdata, handles)
% hObject    handle to cmdOscModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Do PCA on all channels
%Choose your epoch
curr_epoch = select_window(hObject,eventdata,handles)

channel_Ts = extractEpoch(handles,curr_epoch);
band_Ts = extractBands(handles,channel_Ts);

%F(round(length(band_Ts(:,1,1))./10)) = struct('cdata',[],'colormap',[]);

iband = 2;

%Compute PLV for all pairwise channels
figure();
nchann = length(band_Ts(1,1,:));

PLV = zeros(length(band_Ts(:,1,1)),nchann,nchann);

for ii = 1:length(band_Ts(1,1,:))
    for jj = 1:ii-1
        subplot(4,1,1);plot(real(band_Ts(:,iband,ii)));hold on;plot(real(band_Ts(:,iband,jj)),'r');
        phasediff = phase(band_Ts(:,iband,ii)) - phase(band_Ts(:,iband,jj));
        PLV(:,ii,jj) = sin(phasediff);
        subplot(4,1,2);plot(phasediff);hold on;title('Phase Difference In Theta');
        subplot(4,1,3);plot(diff(PLV(:,ii,jj)));hold on;title('Derivative of Phase Difference');
        
    end
end


PLV2 = sin(phase(band_Ts(:,iband,1)) - phase(band_Ts(:,iband,2)));
pval = angle(band_Ts(:,iband,:));
subplot(4,1,4);plot(PLV2);ylim([-2,2]);title('Sin Phase Difference');


% figure(81);
% plot(channel_Ts);

%figure();
%plot(PLV);title('Left and Right PLV');
%ylim([-2,2]);

% Animation for phase offset
% for kk = 1:50:length(band_Ts(:,1,1));
%     figure(82);
%     clf;
%     for jj = iband;
%         subplot(3,1,1);
%         %plot(channel_Ts);
%         line([kk,kk],[0,2]);
%         
%         subplot(2,1,1);
%         %plot(imag(PLV));
%         plot(PLV);
%         %plot(angle(band_Ts(:,jj,1)) - angle(band_Ts(:,jj,2)));
%         %plot([diff(sin(phase(band_Ts(:,jj,1)) - phase(band_Ts(:,jj,2))));0],'r');
%         line([kk,kk],[-pi,pi]);
%         
%         subplot(2,1,2);
% 
%         for ss = 1:length(band_Ts(1,1,:))
%             polar(pval(kk,ss),1,'o');hold on;
%         end
%     end
%     pause(0.001);
%     %F(kk) = getframe(gcf);
% end

disp('Done with osc model');

function [delta,theta,alpha,beta,bbgamma] = breakBands(ts)
disp('Band Analysis...');
%bp.delta = fdesign.bandpass('',)


function edtDSFactor_Callback(hObject, eventdata, handles)
% hObject    handle to edtDSFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtDSFactor as text
%        str2double(get(hObject,'String')) returns contents of edtDSFactor as a double
handles.ANALYSIS.ds_factor = str2double(get(hObject,'String'));

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edtDSFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtDSFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%MAIN 1-D PLOTTING FUNCTION BELOW
% --- Executes on button press in cmdSelect.
function cmdSelect_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sig_cont = select_window(hObject,eventdata,handles);
solo_chann = 0;

%Increment the current interval
handles.ANALYSIS.active_interval = handles.ANALYSIS.active_interval + 1;

%Get channel data for the interval of interest
if solo_chann
    %if you want a single channel (the active one)
    achann = get_intv(handles,handles.UI.active_chann,sig_cont);
else
    %if you want all the loaded channels
    achann = get_intv(handles,find(handles.DATA.TS_List==1)',sig_cont);
end

handles.ANALYSIS.intv{handles.ANALYSIS.active_interval}.achann = calc_Oscil(handles,achann);

%Plot oscillatory data
band_bars = plot_Oscil(handles);

handles.ANALYSIS.Aggr.blp{handles.ANALYSIS.Aggr.idx} = band_bars;

%Paint the window onto the SGview
axes(handles.axes1);
hold on;
c = lines;
for ii = 1:length(handles.ANALYSIS.intv)
    line([handles.ANALYSIS.intv{ii}.achann{1}.tvect(1),handles.ANALYSIS.intv{ii}.achann{1}.tvect(end)],[0,0],'LineWidth',10,'Color',c(ii,:));
end


guidata(hObject,handles);

%Calculate the PSD and also oscillatory bands; these will go hand in hand
%for these analyses
function achann = calc_Oscil(handles,achann)

    for ii = 1:length(achann)
        
        %Now calculate the PSD using pwelch
        %[P,F,B] = pwelch(achann.intv{1}.ds,512,500,2^10,achann.dsFs,'ConfidenceLevel',0.95);
        %Or with multitaper
        [achann{ii}.P,achann{ii}.F,achann{ii}.B] = pmtm(achann{ii}.ds,10,2^10,achann{ii}.dsFs,'ConfidenceLevel',0.95);

        %Extract the slope of the PSD
        achann{ii}.soi = achann{ii}.F > 1 & achann{ii}.F < 20;
        achann{ii}.slop = (10*log10(achann{ii}.P)) \ (10 * log10(achann{ii}.F));
        achann{ii}.slop20 = (10*log10(achann{ii}.P(achann{ii}.soi))) \ (10 * log10(achann{ii}.F(achann{ii}.soi)));
        
        %Now calculate the oscillatory power in each band
        %achann{cidx}.intv{1}.bands = 
        achann{ii}.bandData = extractSGBands(handles,achann{ii}.ds);
        achann{ii}.oscData = extractBands(handles,achann{ii}.ds);
    end

function band_bars = plot_Oscil(handles,achann)
    %Start plotting
    figure(handles.ANALYSIS.active_figure);
    clf;
    numints = length(handles.ANALYSIS.intv);
    
    cmap = [];
    
    for jj = 1:numints
        achann = handles.ANALYSIS.intv{jj}.achann;
        for ii = 1:length(achann)
            subplot(3,length(achann),ii);
            colormap('lines');
            h = plot(achann{ii}.tvect,achann{ii}.ds);hold on;
            xlabel('Time (s)');ylabel('Voltage (mV)');
            title(['Raw Downsampled Signal - Channel ' num2str(ii)]);

            %set up color scheme properly
            %c = get(h,'Color');
            cmap = colormap;

            subplot(3,length(achann),2+ii);
            bplot(achann{ii}.F,achann{ii}.P,achann{ii}.B,cmap(jj,:));hold on;
            title('PSD Plots');
            %plot(achann{ii}.F(achann{ii}.soi),achann{ii}.slop20 * achann{ii}.F(achann{ii}.soi),zeros(size(achann{ii}.F,1),2),[1 0 0]);
        end
    end
    figure(handles.ANALYSIS.active_figure);
    achann = handles.ANALYSIS.intv{1}.achann;
    
    for ii = 1:length(achann)
        for jj = 1:numints
            %Normalize the bands
            handles.ANALYSIS.osc_bands(:,jj,ii) = 10*log10(mean(handles.ANALYSIS.intv{jj}.achann{ii}.bandData ./ mean(handles.ANALYSIS.intv{jj}.achann{ii}.bandData(:,6),1)));
            %Don't normalize the bands
            handles.ANALYSIS.osc_bands(:,jj,ii) = mean(abs(handles.ANALYSIS.intv{jj}.achann{ii}.bandData));

            handles.ANALYSIS.psd.P(:,:,jj,ii) = handles.ANALYSIS.intv{jj}.achann{ii}.P;
            handles.ANALYSIS.psd.F(:,:,jj,ii) = handles.ANALYSIS.intv{jj}.achann{ii}.F;
            handles.ANALYSIS.psd.B(:,:,jj,ii) = handles.ANALYSIS.intv{jj}.achann{ii}.B;
        end
        %subplot(1,length(achann),ii)
        %bplot(handles.ANALYSIS.psd.F,handles.ANALYSIS.psd.P,handles.ANALYSIS.psd.B,colormap('jet'));
        %title('PSD Plots');
        
        subplot(3,length(achann),4+ii)
        b = bar(handles.ANALYSIS.osc_bands(:,:,ii));
        for kk = 1:length(b)
            b(kk).FaceColor = cmap(kk,:);
        end
        set(gca,'XTickLabels',{'Delta','Theta','Alpha','Beta','Gamma','NORM'});
        title('Banded Power');
    end
    %Store banded data and colors into larger array
    band_bars = handles.ANALYSIS.osc_bands;
    

function achann = get_intv(handles,onchanns,epoch)
    %Check if we want to downsample the data
    ds_fact = handles.ANALYSIS.ds_factor;
    cidx = 0;

    for activec = onchanns
        cidx = cidx + 1;

        disp(['Calculating PSD for chann ' num2str(activec)]);
        %Grab active channel and decimate it
        achann{cidx}.ts = handles.DATA.chann{activec}.ts;
        achann{cidx}.dsFs = handles.DATA.chann{activec}.Fs ./ ds_fact;
        set(handles.txtDSFs,'String',['Fs: ' num2str(achann{cidx}.dsFs)]);

        %Decimate the signal
        achann{cidx}.ds = decimate(achann{cidx}.ts,ds_fact);

        %IF WE'RE DEALING WITH AN EEG SIGNAL, we need to shift it
        epoch.ivn = (epoch.tvect(1)*achann{cidx}.dsFs) : (epoch.tvect(end)*achann{cidx}.dsFs);
        if epoch.fromChannel <= 2
            %window is in the reference of the LFPs
            if activec > 2
                epoch.ivn = (epoch.tvect(1)*achann{cidx}.dsFs - handles.DATA.lag_val) : (epoch.tvect(end)*achann{cidx}.dsFs - handles.DATA.lag_val);
            end
        else
            %Window is in the reference of the EEGs
            if activec <= 2
                epoch.ivn = (epoch.tvect(1)*achann{cidx}.dsFs + handles.DATA.lag_val) : (epoch.tvect(end)*achann{cidx}.dsFs + handles.DATA.lag_val);
            end
        end
        %Ignoring doShift in this implementation

        %Now take the the interval chunk we wanted
        achann{cidx}.ds = achann{cidx}.ds(round(epoch.ivn));
        achann{cidx}.tvect = linspace(epoch.tvect(1),epoch.tvect(end),length(achann{cidx}.ds));
    end

%figure;
%pwelch(achann.intv{1}.ds,512,500,2^10,achann.dsFs,'ConfidenceLevel',0.95);

%subplot(2,1,2);
%plot(b,10*log10(abs(a)));hold on;

function bplot(F,P,B,c)
%B comes directly from the computation functions
b(:,1) = abs(10*log10(abs(squeeze(P))) - 10*log10(abs(B(:,1))));
b(:,2) = abs(10*log10(abs(squeeze(P))) - 10*log10(abs(B(:,2))));

h = boundedline(F,10*log10(abs(P)),b,'alpha','cmap',c);
xlabel('Frequency (Hz)');ylabel('Power (dB)');
%set(h,'Colormap',c);
xlim([0 100]);


% --- Executes on button press in cmdAChannTF.
function cmdAChannTF_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAChannTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a = handles.UI.active_chann;

disp('Computing T-F for Active Channel...');

%Make a local container for the raw data
chann_ts = [];
chann_fs = handles.DATA.chann{min(a)}.Fs;
ds_fact = 1; %Purely for display reasons; need not be put on front-UI

for ii = a
    dummy_chann = handles.DATA.chann{ii};
    chann_ts{ii}.ts = decimate(dummy_chann.ts,ds_fact);
end
chann_fs = chann_fs ./ ds_fact;

%Without downsampling ahead of time??
%z-score the input signals
for ii = a
    chann_ts{ii}.ts = zscore(chann_ts{ii}.ts);
    [handles.TF.chann{ii}.S, handles.TF.chann{ii}.F, handles.TF.chann{ii}.T] = spectrogram(chann_ts{ii}.ts,blackmanharris(512),500,2^10,chann_fs); 
end
%Update the decimated Fs
handles.TF.Fs = chann_fs;

%Plot the active chann
SG_Plot_Chann(hObject,handles,handles.UI.active_chann);

%Update that the SG for the active channel has been done
handles.SG.Done(a) = 1;
handles.SG.SGLabels{a} = ['Channel ' num2str(a)];

set(handles.popChann,'String',handles.SG.SGLabels);
set(handles.popChann,'Value',a);

disp(['...Done computing T-F for Channel ' num2str(handles.UI.active_chann)]);
guidata(hObject,handles);


function update_UI_SG_Pop(handles)

%check the handles SG tracking variables and populate the list
%appropriately


% --- Executes on button press in cmdAutoAlignStims.
function cmdAutoAlignStims_Callback(hObject, eventdata, handles)
% hObject    handle to cmdAutoAlignStims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Automated version here; move to other button when ready
disp('Aligning channels');
d = fdesign.bandpass('N,F3dB1,F3dB2',20,129,131,1000);
Hd = design(d,'butter');

%Compare the LFP signal with an EEG signal
%I added the intrinsic differential signal for common mode rejection of
%noise
%This was more robust for SINGLE channels from the EEG
%And seemed to stay the same for DIFFERENTIAL channels from the EEG
c1 = zscore(filtfilthd(Hd,handles.DATA.chann{1}.ts));
c2 = zscore(filtfilthd(Hd,handles.DATA.chann{3}.ts - handles.DATA.chann{4}.ts));

%Try a cross correlation approach
[cc,lags] = xcorr(c1,c2);

[~,max_cc_idx] = max(cc);
%We know that the lag should be limited, it's not going to be > 
max_cc_loc = lags(max_cc_idx)

% Other way of doing it, trying to find the first point of large overlap;
% this is shit

% %find max region
% %first z-score
ccz = zscore(cc);

%Try plotting the aligned signals;
figure;
subplot(3,1,1);
c1t = 1:length(handles.DATA.chann{1}.ts);
c2t = 1:length(handles.DATA.chann{3}.ts);
c2ts = c2t + max_cc_loc;

plot(c1t,zscore(handles.DATA.chann{1}.ts));hold on;
plot(c2ts,zscore(handles.DATA.chann{3}.ts));

%Plot filtered stim artifact one
subplot(3,1,2);
plot(c1);hold on;plot(c2);

subplot(3,1,3);
plot(lags,ccz);

%Store lag value for EEG channel set
handles.DATA.lag_val = max_cc_loc;

disp('Done Aligning');
guidata(hObject,handles);


% --- Executes on button press in cmdDiffChann.
function cmdDiffChann_Callback(hObject, eventdata, handles)
% hObject    handle to cmdDiffChann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdCmpLoad.
function cmdCmpLoad_Callback(hObject, eventdata, handles)
% hObject    handle to cmdCmpLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdSGSett.
function cmdSGSett_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSGSett (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chann_sig = extract_window(handles,sig_cont)

%get channel data
chann_sig = handles.DATA.chann{sig_cont.fromChannel}.ts(sig_cont.nbounds(1):sig_cont.nbounds(2));


% --- Executes on button press in cmdMedFilt.
function cmdMedFilt_Callback(hObject, eventdata, handles)
% hObject    handle to cmdMedFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sig_cont = select_window(hObject,eventdata,handles);
sig_oi = extract_window(handles,sig_cont);

figure;
subplot(2,2,1);
plot(sig_oi);title('Unfiltered interval of interest on channel');

%Filter the data
for ii = 1:length(sig_cont.fromChannel)
    %filted(:,ii) = medfilt1(handles.DATA.chann{ii}.ts,50);
    hampeled(:,ii) = hampel(1:length(sig_oi),sig_oi,1,3,'Adaptive',0.1);
end

subplot(2,2,2);
plot(hampeled);title('Filtered signal');

subplot(2,2,3);
[s,f,t] = spectrogram(sig_oi,blackmanharris(512),500,2^10,1000);
imagesc(t,f,10*log10(abs(s)));colormap('jet');set(gca,'YDir','normal');

subplot(2,2,4);
[s,f,t] = spectrogram(hampeled,blackmanharris(512),500,2^10,1000);
imagesc(t,f,10*log10(abs(s)));colormap('jet');set(gca,'YDir','normal');

% --- Executes on button press in cmdBreak.
function cmdBreak_Callback(hObject, eventdata, handles)
% hObject    handle to cmdBreak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cmdPlotAggrBands.
function cmdPlotAggrBands_Callback(hObject, eventdata, handles)
% hObject    handle to cmdPlotAggrBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure;
cmap = colormap('colorcube');
%Plot the left signals first (channel 1)
for ii = 1:2 %channel number
    subplot(1,2,ii);
    for jj = 2 %band number
        for ll = 1:length(handles.ANALYSIS.Aggr.blp)
            dummysize = size(handles.ANALYSIS.Aggr.blp{ll});
            intvs = dummysize(2);
            for kk = 1:intvs %interval number
                %plot all patients next to each other
                osc(ll,:,:,:) = handles.ANALYSIS.Aggr.blp{ll};
                b = bar(squeeze(osc(:,jj,:,ii))');
                xlabel('Interval Number');ylabel('Power');
                leg_end{kk} = ['Patient #' num2str(kk)];
            end
            for pp = 1:length(b)
                b(pp).FaceColor = cmap(pp,:);
            end
        end
    end
    title(['Channel ' num2str(ii) ' Theta power']);
    legend(leg_end);
end
suptitle('TurnOn Powers');


% --- Executes on button press in cmdLowPassExtr.
function cmdLowPassExtr_Callback(hObject, eventdata, handles)
% hObject    handle to cmdLowPassExtr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Just take channel 1
chann = 1;

raw_sig = handles.DATA.chann{chann}.ts;

%Band it below
lph = fdesign.lowpass('N,F3dB',20,20,handles.DATA.chann{chann}.Fs);
lpf = design(lph);

filted_sig = filtfilthd(lpf,raw_sig);
disp('Filtered...');


% --- Executes on button press in cmdSelectInfo.
function cmdSelectInfo_Callback(hObject, eventdata, handles)
% hObject    handle to cmdSelectInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sig_cont = select_window(hObject,eventdata,handles);
disp('Interval Info:');

if sig_cont.tvect(end) - sig_cont.tvect(1) > 15
    str_ok = 'OK!';
else
    str_ok = 'SMALL!!';
end
disp(['Length of this interval is: ' num2str(sig_cont.tvect(end) - sig_cont.tvect(1)) ' which is ' str_ok]);
disp([num2str(sig_cont.tvect(1)) ' to ' num2str(sig_cont.tvect(end))]);
