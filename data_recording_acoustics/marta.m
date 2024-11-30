function marta(varargin)
%MARTA  - experiment sequencing tool
%
%	usage:  marta(trials, info, ...)
%
% TRIALS, INFO specify the experiment structure; see PARSEEXPFILE
%
% Recognized 'NAME',VALUE parameters (defaults shown within {}):
%   AUDIO   - override audio recording parameters; e.g. struct('DEVIDX',1,'NCHAN',2,'SRATE',22050)
%	LAYOUT  - load window configuration from specified variable
%	RECORD  - record data {1} or display only (0)
%	SAVELAY - save window layout to specified variable
%   VOXTHR  - VOX-activation threshold
%
% examples:
%  [trials,info] = ParseExpFile('test')
%  marta(trials, info)   % record by default
%  marta(trials, info, 'RECORD', 0)   % view stimuli w/o recording

% INFO specifies experiment information (except for CSS all written to log):
%  SESSION     - session name
%  PREFIX      - prefixed to each output file
%  LOG         - log filename (diary of session)
%  TIMINGS     - MAT filename for recording trial timestamps
%  PARTICIPANT - participant(s) info (array defines dyadic expt)
%  ACQHW       - hardware interface module(s) for recording and their parameters
%  CSS         - prefixed to each trial HTML if present
%
% TRIALS defines trials sequence
%  TYPE        - PAUSE | RECORD | BLKREC
%  FNAME       - trial filename 
%  PROMPT      - experimenter prompt (defaults to FNAME)
%  DUR         - explicit reset of trial duration, if any
%  STIM        - one or more trial stimuli
%   PART       - partner seeing stimulus:  1 | 2 | [] both (default)
%   DELAY      - delay to display of stimulus (msecs from start of trial)
%   HTML       - delivered to browser prefixed by CSS
%   EXTRA      - may include HANDLER (executable function plus arg), additional flags
%
% to view available AUDIO input devices and obtain device ID  use
%  d = audiodevinfo; d = d.input; [{d.ID}',{d.Name}']

% mkt 03/08 v0.1
% mkt 08/23 v2.0 rewritten; not bw compatible

vers = 'mkt 08/23 v2.0';	% current version

% PAUSE trials stop active recording and autocycling and persist until cleared
% RECORD trials start and stop active recording and initiate ISI under autocycling
% BLKREC trials start recording if necessary but do not stop it, but display timings
%  function as with RECORD trials

% trial recording timestamps are recorded if the TIMINGS tag is included in INFO
% and saved as a variable of that name to a mat file of that name
% timestamps: FNAME, RECSTART, START, END, COMPLETED, STIM1..STIMn
% RECSTART is "now" following return of START recording from all ACQHW modules
% START is "now" at onset of trial
% END is "now" at end of trial (precedes call to STOP recording for all ACQHW modules)
% COMPLETED is true for normal completion, false if trial aborted
% STIM1,..STIMn give timstamps of onsets of displayed stimuli for each trial
% dt = @(a,b) seconds(time(between(datetime(a,'convertFrom','datenum'),datetime(b,'convertFrom','datenum'))))
persistent timings

% temporary HTML filename for rendering media
persistent htmlFname

% internal state
persistent state	

if nargin < 1, help marta; return; end

%% branch by action
if isstruct(varargin{1})
	action = 'INIT';
elseif ischar(varargin{1})
	action = varargin{1};
else
	varargin(1:2) = [];		% delete src, evt
	action = varargin{1};
end

switch action

% ABOUT  - show current marta and handler version info
	case 'ABOUT'
		fprintf('MARTA  - experiment sequencing tool (%s)\n', vers)
		if ~isempty(state)
			for mi = 1 : length(state.ACQHW), state.ACQHW(mi).FUN('ABOUT'); end
		end

% CLOSE	 - shut down	
	case 'CLOSE'
		if ~strcmp('PAUSED',state.MODE), return; end
		if strcmp('Cancel',uiconfirm(state.CTL.FH,'Exit MARTA?','Verify...','Icon','question')), return; end
		ShutDown
		
% INIT  - initialize
	case 'INIT'
		if verLessThan('matlab','9.14')
			error('This version of MARTA requires Matlab R2023a or later')
		end
		htmlFname = 'tempH3gX2t9.html';		% temporary HTML filename for rendering media
		Initialize(varargin{1}, varargin{2}, varargin{3:end})
		
% LAYOUT  - specify window configuration		
	case 'LAYOUT'
		layout = varargin{2};
		for li = 1 : length(layout)
			k = find(strcmp(layout(li).NAME,{state.WINDOWS.NAME}),1);
			if isempty(k), continue; end
			state.WINDOWS(k).HANDLE.Position = layout(li).POSITION;
		end

% LIST  - list update, either through click on list element	or menu entry	
	case 'LIST'
		if nargin < 4	% click on list; value set directly
			SetTrial(get(gcbo,'value'));	
		else			% menu item - pass increment (+/-1)
			SetTrial([],varargin{2});
		end	

% SAVELAY  - save window configuration		
	case 'SAVELAY'
		layout(length(state.WINDOWS),1) = struct('NAME',[],'POSITION',[]);
		for wi = 1 : length(state.WINDOWS)
			w = state.WINDOWS(wi);
			layout(wi) = struct('NAME',w.NAME,'POSITION',w.HANDLE.Position);
		end
		vs = GetName('Save layout as...','layout');
		if ~isempty(vs)
			assignin('base',vs,layout);
			fprintf('current window layout saved as %s in base workspace\n',vs)
		end

% TOGGLE  - toggle sequencing (called by control button press)
	case 'TOGGLE'
		switch state.MODE
			case 'PAUSED', StartTrial	% start execution of current trial -> ACTIVE
			case 'ACTIVE', AbortTrial	% abort execution of current trial -> PAUSED
			case 'ISI',    StopISI		% abort ISI delay                  -> PAUSED
		end			
		
% VOX  - display VOX-activated stimulus (called from acq_AUDIO)
	case 'VOX'
		if ~isempty(state.VOXSTIM)
			UpdateDisplay(state.VOXSTIM{1},state.VOXSTIM{2});
			state.VOXSTIM = [];
		end
		
% ERROR trap		
	otherwise
		error('MARTA: unrecognized parameter (%s)', action)

end  % switch

%% ========== LOCAL FUNCTIONS ============================================================

%% ----- ABORTTRIAL  - interrupt active trial

function AbortTrial

state.CTL.AUTO.Value = 0;
t = timerfind('tag','STIM'); stop(t); delete(t)
state.TRIALTIMER.StopFcn = ''; stop(state.TRIALTIMER)
state.COMPLETED = false;

EndTrial

end % AbortTrial

%% ----- BROWSERCHANGEDSIZE  - browser changed size; update uihtml child

function BrowserChangedSize(fh, ~, bh)

p = fh.Position;
bh.Position = [1 1 p(3:4)];

end % BrowserChangedSize

%% ----- ENDTRIAL  - end-of-trial processing; called on TRIALTIMER timeout or by AbortTrial

function EndTrial(~,~)

% blank screen(s)
set(state.BROWSER,'HTMLSource','<BODY bgcolor="darkgrey" />')

% stop recording unless BLKREC recording in progress
if state.RECORD, StopRecording; end

% stop in-progress stimulus handler
if ~isempty(state.STIMFUN), 
	state.STIMFUN('STOP',state.STIMARG); 
	state.STIMFUN = [];
	state.STIMARG = [];
end

% increment current trial unless on last trial
nTrials = length(state.TRIALS);
oldIdx = state.CURIDX;
if state.CURIDX < nTrials, 
	SetTrial(state.CURIDX+1); 
elseif state.CURIDX > nTrials,	% "extra" trial
	SetTrial(nTrials+1);		% this resets the filename query 
end

% PAUSE and RECORD trials stop all ongoing (BLKREC) recording
if state.RECORD && state.RECORDING
	switch state.CURTRIAL.TYPE
		case {'PAUSE','RECORD'}, 
			StopRecording
		otherwise, % ignore
	end
end

% if autocycling set up for next trial or ISI delay
if state.CTL.AUTO.Value 
	if state.ISI
		state.MODE = 'ISI';
		set(state.ISITIMER, 'StartDelay',state.ISI, 'executionMode','singleShot', 'timerFcn',@ISITimeOut);
		start(state.ISITIMER);
	else
		StartTrial
	end

% else pause and wait for input
else
	set(state.CTL.CTLB,'String','Start','backgroundColor',[.8 .9 .8])
	state.CTL.LIST.Enable = 'on';
	set(state.CTL.LISTMEN,'Enable','on')
	state.CTL.PROGRESS.Visible = 'off';
	state.MODE = 'PAUSED';
	if strcmp('BLKSTRT',state.CURTRIAL.TYPE), state.CTL.AUTO.Value = 1; end
end

end % EndTrial

%% ----- GETNAME  - get name interactively

function name = GetName(promptStr, defStr)

pos = state.CTL.FH.Position;
width = 300;
height = 100;
figPos = [pos(1)+(pos(3)-width)/2, pos(2)+(pos(4)-height)/2, width, height];

fh = uifigure('name', promptStr, ...
	'Position', figPos, ...
	'WindowStyle', 'alwaysontop', ...
	'UserData', 0);

% name field
eh = uicontrol(fh, ...
	'Position', [20 60 width-40 20], ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', defStr);

% OK, cancel buttons
uicontrol(fh, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback',{@Resume,1});
uicontrol(fh, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback',{@Resume,0});

% wait for input
uiwait(fh);
if ishandle(fh) && fh.UserData,
	name = strtok(get(eh, 'string'));
else,
	name = [];
end;
if ishandle(fh), delete(fh); end

	function Resume(~,~,choice)
		fh.UserData = choice; uiresume(fh)
	end

end % GetName

%% ----- INITIALIZE  - initialization sequence

function Initialize(trials, info, varargin)

% defaults
dur = 3;		% default duration (secs)
ISI = 1;		% default ISI (secs)
record = 1;		% record by default
layout = [];	% use default layout
audioParams=[];	% use info-specified audio params
VOXthr = [];	% override stimulus threshold
html = ['<HEAD><STYLE>#txt {text-align:center;margin-top:100px;font-style:italic;font-family:Arial;font-size:64;font-weight:bold;color:#000088}</STYLE>', ...
		'</HEAD><BODY bgcolor="darkgrey"><div id="txt">Don''t w<FONT color="#C80000">o</FONT>rry, be h<FONT color="#00C800">a</FONT>ppy...</div></BODY>'];

record = ~isempty(info.ACQHW);		% recording requires a HW module

% ensure not already running
h = findobj(findall(0),'name','MARTA','type','figure');
if ~isempty(h), figure(h); return; end

% parse additional arguments
for ai = 2 : 2 : length(varargin)
	switch upper(varargin{ai-1})
		case 'AUDIO', audioParams = varargin{ai};
		case 'LAYOUT', layout = varargin{ai};
		case 'RECORD', record = (varargin{ai} == 1);
		case 'VOXTHR', VOXthr = varargin{ai};
		otherwise, error('unrecognized MARTA parameter (%s)',varargin{ai-1});
	end
end

% disable backtracing (enable for debugging)
warning('off','backtrace')	

% init control window
controller = InitController(dur, ISI, trials);
windows = struct('NAME','CONTROLLER','HANDLE',controller.FH);

% init logging
if record, 
	diary(info.LOG)		% start logging
	if exist(sprintf('%s.mat',info.TIMINGS),'file') == 2 && strcmp('Cancel', ...
			uiconfirm(controller.FH,sprintf('%s.mat exists',info.TIMINGS),'Verify...', ...
					'Icon','warning','Options',{'Overwrite','Cancel'}))
		delete(controller.FH)
		return
	end
end
ls = repmat('=',1,60);
fprintf('\n\n%s\n  %s initiated %s\n', ls, info.SESSION, upper(char(datetime('now','format','ddMMMyy HH:mm:ss'))))
if record,
	fprintf('%s\n',repmat('-',1,60)) 
	InitLog(rmfield(info,'CSS'))
end
fprintf('%s\n\n',ls)

% init HW modules
ACQHW = [];
if record
	if isempty(info.ACQHW)
		record = 0;
	else
		for mi = 1 : length(info.ACQHW)
			n = info.ACQHW(mi).NAME;
			f = str2func(n);
			p = info.ACQHW(mi).PARAMS;
			k = find(strcmp('AUDIO',{info.ACQHW(mi).PARAMS.NAME}));
			if ~isempty(audioParams) && ~isempty(k)		% replace info AUDIO with command line override
				info.ACQHW(mi).PARAMS(k).ARGS = audioParams;
			end
			try
				cfg = f('INIT',info.ACQHW(mi).PARAMS);
			catch
				fprintf('%s\n',lasterr)
				diary off
				delete(controller.FH)
				delete(findobj('TAG','MARTA'))
				return
			end
			hw = struct('NAME',n,'FUN',f,'CFG',cfg);
			if isempty(ACQHW), ACQHW = hw; else, ACQHW(end+1) = hw; end
			if ~isempty(hw.CFG) && isfield(hw.CFG,'WINDOW')
				w = hw.CFG.WINDOW;
				for wi = 1 : length(w)
					if isempty(windows), windows = w(wi); else, windows(end+1) = w(wi); end
				end
			end
		end
	end
end	

% init browser window(s)
fh = uifigure('tag','MARTA','name','Participant 1','closeRequestFcn',{@marta,'CLOSE'});
p = fh.Position; p = [1 p(2) 940 p(4)]; fh.Position = p;
bh = uihtml(fh,'position',[1,1,p(3:4)]);
set(fh,'AutoResizeChildren','off','SizeChangedFcn',{@BrowserChangedSize,bh});
bh.HTMLSource = html;
windows(end+1) = struct('NAME','BROWSER','HANDLE',fh);
if length(info.PARTICIPANT) > 1
	p = [p(1)+50 p(2)-50 p(3:4)];
	fh = uifigure('tag','MARTA','position',[p(1)+50 p(2)-50 p(3:4)],'name','Participant 2');
	bh(2) = uihtml(fh,'position',[1,1,p(3:4)]);
	set(fh,'AutoResizeChildren','off','SizeChangedFcn',{@BrowserChangedSize,bh(2)});
	bh(2).HTMLSource = html;
	windows(end+1) = struct('NAME','BROWSER2','HANDLE',fh);
end

% preallocate trial recording timings (double number of trials)
if record
	n = max(cellfun(@length,{trials.STIM}));
	varNames = ["FNAME","RECSTART","START","STOP","COMPLETED"];
	varTypes = ["string","double","double","double","logical"];
	for k = 1 : n
		varNames = [varNames , sprintf("STIM%d",k)];
		varTypes = [varTypes , "double"];
	end
	timings = table('Size',[length(trials)*2 , length(varNames)],'VariableNames',varNames,'VariableTypes',varTypes);
end

% set initial state
state = struct( ...
				'ACQHW', ACQHW, ...				% acquisition HW modules
				'BROWSER', bh, ...				% browser(s) (uihtml object)
				'COMPLETED', [], ...			% trial completion status (0 aborted, 1 normal)
				'CTL', controller, ...			% controller object
				'CURIDX', [], ...				% current trial index
				'CURTRIAL', [], ...				% current trial
				'DUR', dur, ...					% active duration (secs)
				'FNAME', [], ...				% trial filename
				'INFO', info, ...				% experient info
				'ISI', ISI, ...					% active ISI (secs)
				'ISITIMER', timer, ...			% ISI timer object
				'LATENCY', 100, ...				% latency added to requested recording duration to ensure completion (ms)
				'MAXDUR', 300, ...				% maximum supported BLKREC duration (secs)
				'MODE', 'PAUSED', ...			% display status:  PAUSED, ACTIVE, ISI
				'TRIALTIMER', timer, ...		% trial timer object
				'RECORD', record, ...			% data are recorded if nonzero
				'RECORDING', 0, ...				% recording status: 0 off; 1 on; 2 BLKREC (stopped only by PAUSE trial)
				'STIMARG', [], ...				% current stimulus handler argument
				'STIMFUN', [], ...				% current stimulus handler
				'TIDX', 1, ...					% completed trials index
				'TRIALS', trials, ...			% experiment trials
				'VOXSTIM', [], ...				% VOX-activated stimulus
				'VOXTHR', VOXthr, ...			% overrides stim-specified threshold
				'WINDOWS', windows);

% update layout if necessary
if ~isempty(layout)
	for li = 1 : length(layout)
		k = find(strcmp(layout(li).NAME,{windows.NAME}),1);
		if isempty(k), continue; end
		windows(k).HANDLE.Position = layout(li).POSITION;
	end
end

% set initial trial
SetTrial(1)

% focus control window
figure(controller.FH)

end % Initialize

%% ----- INITCONTROLLER  - initialize control window
function controller = InitController(dur, ISI, trials)

trialList = {trials.FNAME};
trialList{end+1} = '<< Select Interactively >>';

color = [0.9255 0.9137 0.8471];
width = 400;
height = 400;
pos = get(0,'monitorPositions');
k = find(pos(:,1)==1);
pos = pos(k,:);
pos = [(pos(3)-width)-10, 45, width, height];
fh = uifigure('name','MARTA', ...
			'position', pos, ...
			'closeRequestFcn', {@marta,'CLOSE'}, ...
			'color', color, ...		
 			'windowStyle','alwaysontop', ...
			'resize', 'off');

% menu items
mh = uimenu(fh, 'label', '&File');
uimenu(mh, 'label', 'About MARTA...', ...
			'callback', {@marta,'ABOUT'});
uimenu(mh, 'label', 'E&xit', ...
			'separator', 'on', ...
			'accelerator', 'X', ...
			'callback', {@marta,'CLOSE'});

mh = uimenu(fh, 'label', '&Control');
uimenu(mh, 'label', 'Save &Layout', ...
			'accelerator', 'L', ...
			'callback', {@marta,'SAVELAY'});
menuH = uimenu(mh, 'label', '&Backup', ...
			'separator', 'on', ...
			'accelerator', 'B', ...
			'callback', {@marta,'LIST',-1});
menuH(2) = uimenu(mh, 'label', '&Next Trial', ...
			'accelerator', 'N', ...
			'callback', {@marta,'LIST',1});

% border
ph = uipanel(fh, ...
			'units', 'pixels', ...
			'position', [13 106 376 290], ...
			'fontName', 'Helvetica', ...
			'fontSize', 10, ...
			'backgroundColor', color, ...
			'title', ' Current Trial ');

% trial list
listH = uilistbox(fh, ...
			'position',[13 187 376 192], ...
			'items', trialList, ...	
			'itemsData', [1:length(trialList)], ...
			'value', 1, ...
			'backgroundColor', 'w', ...
			'valueChangedFcn', {@marta,'LIST'});

% prompt
promptH = uicontrol(fh, ...
			'style','edit', ...
			'units','pixels', ...
			'position',[19.8 139 362 45], ...
			'fontName', 'Helvetica', ...
			'fontSize', 16, ...
			'horizontalAlignment', 'left', ...
			'enable', 'inactive', ...
			'backgroundColor', [1 1 .9], ...
			'hitTest', 'off', ...
			'string', trials(1).PROMPT);
			
% progress bar
uicontrol(fh, ...
			'style','text', ...
			'units','pixels', ...
			'position',[17 110.5 60 21], ...
			'fontName', 'Helvetica', ...
			'fontSize', 10, ...
			'horizontalAlignment', 'right', ...
			'backgroundColor', color, ...
			'hitTest', 'off', ...
			'string', 'Progress:');
pH = axes('parent',ph,'units','pixels', 'position',[67 7 301 21], 'hitTest','off');
set(pH, 'color','w', 'box','on', 'xtick',[], 'ytick',[], 'xlim',[0 1],'ylim',[0 1]);
progressH = patch([0 0 1 1],[0 1 1 0],'r','edgeColor','w','visible','off','parent',pH);

% acquisition duration
uicontrol(fh, ...
			'style','text', ...
			'units','pixels', ...
			'position',[17 67 200 24], ...
			'fontName', 'Helvetica', ...
			'fontSize', 12, ...
			'horizontalAlignment', 'right', ...
			'backgroundColor', color, ...
			'string', 'Acquisition Duration (secs):');
durH = uicontrol(fh, ...
			'style','edit', ...
			'units','pixels', ...
			'position',[221 67 60 30], ...
			'fontName', 'Helvetica', ...
			'fontSize', 12, ...
			'horizontalAlignment', 'center', ...
			'backgroundColor', 'w', ...
			'string', sprintf('%.2g',dur));

% ISI
uicontrol(fh, ...
			'style','text', ...
			'units','pixels', ...
			'position',[292 67 25 24], ...
			'fontName', 'Helvetica', ...
			'fontSize', 12, ...
			'horizontalAlignment', 'right', ...
			'backgroundColor', color, ...
			'string', 'ISI:');
isiH = uicontrol(fh, ...
			'style','edit', ...
			'units','pixels', ...
			'position',[321 67 60 30], ...
			'fontName', 'Helvetica', ...
			'fontSize', 12, ...
			'horizontalAlignment', 'center', ...
			'backgroundColor', 'w', ...
			'string', sprintf('%.2g',ISI));

% autocycling checkbox
autoH = uicontrol(fh, ...
			'style','checkbox', ...
			'units','pixels', ...
			'position',[49 25 120 30], ...
			'fontName', 'Helvetica', ...
			'fontSize', 12, ...
			'backgroundColor', color, ...
			'string', 'Auto Cycling', ...
			'value', 0);

% control button
ctlH = uicontrol(fh, ...
			'style','pushbutton', ...
			'units','pixels', ...
			'position',[195 16 165 42], ...
			'fontName', 'Helvetica', ...
			'fontSize', 18, ...
			'fontWeight', 'bold', ...
			'backgroundColor', [.8 .9 .8], ...
			'string', 'Start', ...
			'interruptible', 'off', ...
			'callback', {@marta,'TOGGLE'});

% controller
controller = struct('AUTO', autoH, ...			% auto cycling checkbox handle
					'CTLB', ctlH, ...			% sequencing control button
					'DUR', durH, ...			% duration field
					'FH', fh, ...				% controller figure handle
					'ISI', isiH, ...			% ISI field
					'LIST', listH, ...			% trial listbox
					'LISTMEN', menuH, ...		% menu list increment items
					'PROGRESS', progressH, ...	% progress bar
					'PROMPT', promptH);			% prompt field

end % InitController

%% ----- INITLOG  - write info to log

function InitLog(info, indent)

if nargin < 2, indent = 4; end
fl = 25;

fn = fieldnames(info);
for fi = 1 : length(fn)
	n = fn{fi};
	v = info.(n);
	is = repmat(' ',1,indent); 
	ds = repmat('.',1,fl-length(n)-indent);
	if isstruct(v)
		for vi = 1 : length(v)
			fprintf('%s%s\n', is, n)
			InitLog(v(vi),indent+2)
		end		
	elseif ischar(v)
		fprintf('%s%s%s %s\n', is, n, ds, v)
	elseif isscalar(v)
		fprintf('%s%s%s %g\n', is, n, ds, v)
	else
		fs = sprintf('%s%s%s [', is, n, ds);
		for k = 1 : length(v)-1, fs = [fs,'%g ']; end
		fs = [fs,'%g]\n'];
		fprintf(fs, v)
	end
end

end % InitLog

%% ----- ISITIMEOUT  - ISI delay timeout handler

function ISITimeOut(~,~)

if state.CTL.AUTO.Value
	StartTrial
else
	set(state.CTL.CTLB,'String','Start','backgroundColor',[.8 .9 .8])
	state.CTL.LIST.Enable = 'on';
	set(state.CTL.LISTMEN,'Enable','on')
	state.CTL.PROGRESS.Visible = 'off';
	state.MODE = 'PAUSED';
end

end % ISITimeOut

%% ----- SETTRIAL  - update current trial (prepare for execution)

function SetTrial(newTrial, increment)

if nargin > 1	% set by menu
	newTrial = state.CTL.LIST.Value + increment;
	if newTrial < 1, return; end
	if newTrial > length(state.TRIALS)+1 
		newTrial = length(state.TRIALS)+1; 	% << Select Interactively >>
		state.CTL.AUTO.Value = 0;
	end
end

% set up extra trial
if newTrial <= length(state.TRIALS)
	trial = state.TRIALS(newTrial);
else
	trial = state.TRIALS(end);
	trial.FNAME = GetName('Filename for extra trial...',sprintf('%s_EXTRA',state.INFO.PARTICIPANT.ID));
	if isempty(trial.FNAME)		% cancel
		newTrial = length(state.TRIALS);
	else
		trial.TYPE = 'RECORD';
		trial.PROMPT = [];
		trial.DUR = [];
		trial.STIM = trial.STIM(1);
		trial.STIM.DELAY = 0;
		trial.STIM.EXTRA = [];
		trial.STIM.HTML = '<BODY bgcolor="white" />';
	end
end

% append version number to filename; use FNAME.wav to test for existing reps
d = dir(sprintf('%s*.wav', trial.FNAME));
if isempty(d), r = 1; else, r = length(d)+1; end	% increment version number
state.FNAME = sprintf('%s_%02d', trial.FNAME, r);

% update state with current (newly selected) trial
state.CTL.LIST.Value = newTrial;
scroll(state.CTL.LIST, newTrial)
state.CURIDX = newTrial;
state.CURTRIAL = trial;	

% update prompt
if isempty(trial.PROMPT), trial.PROMPT = trial.FNAME; end
state.CTL.PROMPT.String = trial.PROMPT;

% update trial duration, ISI
curDur = str2double(state.CTL.DUR.String);
if isempty(curDur)
	curDur = state.DUR;
	state.CTL.DUR.String = sprintf('%.2g',curDur);
elseif curDur ~= state.DUR
	state.DUR = curDur;
end
newDur = trial.DUR;
if ~isempty(newDur) && newDur ~= curDur
	state.DUR = newDur;
	switch trial.TYPE
		case {'BLKSTRT','BLKEND'}, % don't update
		otherwise, state.CTL.DUR.String = sprintf('%.2g',newDur); 
	end
end
ISI = str2double(state.CTL.ISI.String);
if isempty(ISI), 
	state.CTL.ISI.String = sprintf('%.2g',state.ISI);
elseif ISI ~= state.ISI
	state.ISI = ISI;
end

end % SetTrial

%% ----- SHUTDOWN  - shutdown processing

function ShutDown

t = timerfind('tag','STIM'); stop(t); delete(t) 
state.ISITIMER.StopFcn = ''; stop(state.ISITIMER); delete(state.ISITIMER)
state.TRIALTIMER.StopFcn = ''; stop(state.TRIALTIMER); delete(state.TRIALTIMER)
for mi = 1 : length(state.ACQHW), try, state.ACQHW(mi).FUN('CLOSE'); catch, lasterr, end; end
delete(findobj(findall(0),'tag','MARTA'))
try, closereq; catch, close all force ; end
if state.RECORD && isfield(state.INFO,'TIMINGS')
	timings(ismissing(timings.FNAME),:) = [];		% clean up prealloc
	tName = state.INFO.TIMINGS;
	save(state.INFO.TIMINGS,'timings')
	assignin('base',tName,timings)
	evalin('base',sprintf('save %s %s ; clear %s', tName, tName, tName))
	fprintf('\n  Trial timings written to %s.mat\n', tName)
end

ls = repmat('=',1,60);
fprintf('\n%s\n  %s terminated %s\n%s\n', ls, state.INFO.SESSION, upper(char(datetime('now','format','ddMMMyy HH:mm:ss'))),ls)
diary off		% stop logging

if exist(htmlFname,'file'), delete(htmlFname); end

warning('on','backtrace')

end % ShutDown

%% ----- STARTRECORDING  - start recording

function StartRecording(newMode,dur)

if state.RECORDING == 0
	ACQHW = state.ACQHW;
	for mi = 1 : length(ACQHW), ACQHW(mi).FUN('START',state.FNAME,dur); end
	t = now;		% start of recording
	timings.RECSTART(state.TIDX) = t;
	if newMode == 2
		fprintf('BLKREC recording initiated %s\n', char(datetime(t,'convertfrom','datenum','format','HH:mm:ss.SSS')))
	end
	state.RECORDING = newMode;
else
	t = now;		% start of trial
end

timings.START(state.TIDX) = t;			
timings.FNAME(state.TIDX) = state.FNAME;
state.COMPLETED = true;
fprintf('%40s  [%s]  ', state.FNAME, char(datetime(t,'convertfrom','datenum','format','HH:mm:ss.SSS')))

end % StartRecording

%% ----- STARTTRIAL  - start execution of current trial

function StartTrial

trial = state.CURTRIAL;
state.VOXSTIM = [];

% handle possibly updated trial DURation
dur = str2double(state.CTL.DUR.String);
if isempty(dur)
	dur = state.DUR;
	state.CTL.DUR.String = sprintf('%.2g',dur);
elseif dur ~= state.DUR
	state.DUR = dur;
end
if ~isempty(trial.DUR) && trial.DUR ~= dur
	dur = trial.DUR;
	state.DUR = dur;
	switch trial.TYPE
		case {'BLKSTRT','BLKEND'},	% don't update
		otherwise, state.CTL.DUR.String = sprintf('%.2g',dur); 
	end
end

% handle PAUSE trial (turns off autocycling)
if strcmpi('PAUSE',trial.TYPE)
	fprintf('%s\n',trial.PROMPT)
	UpdateDisplay(trial.STIM,[])
	state.CTL.AUTO.Value = 0;

% handle RECORD or BLKREC trial
else

% begin recording
	if state.RECORD
		switch trial.TYPE
			case 'BLKSTRT'
				StartRecording(2,state.MAXDUR)
			case 'RECORD'
				StartRecording(1,dur)
			case 'BLKREC'
				if state.RECORDING
					StartRecording(state.RECORDING,dur)
				else
					StartRecording(1,dur)	% this handles a BLKREC trial run individually
				end
			otherwise 
				StartRecording(state.RECORDING,dur)
		end
	end
	
% display initial stimulus if no delay
	if trial.STIM(1).DELAY==0, UpdateDisplay(trial.STIM(1),1); end
	
% start trial timer:  updates progress bar 20x / sec
% include latency to allow recording to complete full duration
	nUpdates = ceil((state.DUR + state.LATENCY/1000) * 20);
	set(state.TRIALTIMER, 'timerFcn',{@UpdateProgress,state.CTL.PROGRESS,nUpdates}, ...
		'busyMode','drop', 'executionMode','fixedRate', 'period',1/20, ...
		'tasksToExecute',nUpdates, 'stopFcn',@EndTrial)
	state.CTL.PROGRESS.Visible = 'on';
	start(state.TRIALTIMER)
		
% start delayed stimulus timers
% negative delay encodes VOX threshold
	for si = 1 : length(trial.STIM)
		stim = trial.STIM(si);
		if stim.DELAY == 0		% immediate
			continue
		elseif stim.DELAY > 0	% delayed
			t = timer('StartDelay',stim.DELAY/1000,'tag','STIM','timerFcn',{@StimTimeOut,stim,si});
			start(t);
		else					% init VOX-activated
			k = find(strcmp('acq_AUDIO',{state.ACQHW.NAME}));
			thresh = state.VOXTHR;
			if isempty(thresh), thresh = -stim.DELAY; end
			state.ACQHW(k).FUN('VOX',thresh);
			state.VOXSTIM = {stim,si};
		end
	end
end

state.CTL.LIST.Enable = 'off';
set(state.CTL.LISTMEN,'Enable','off')
set(state.CTL.CTLB,'string','Abort','backgroundColor',[.9 .8 .8]);
state.MODE = 'ACTIVE';

end % StartTrial

%% STIMTIMEOUT  - process stimulus after its onset delay (relative to start of trial)

function StimTimeOut(t, ~, stim, si)

UpdateDisplay(stim, si)
stop(t); delete(t)

end % StimTimeOut

%% ----- STOPISI  - abort executing ISI

function StopISI

state.ISITIMER.StopFcn = '';
stop(state.ISITIMER);
state.CTL.AUTO.Value = 0;
set(state.CTL.CTLB,'String','Start','backgroundColor',[.8 .9 .8])
state.MODE = 'PAUSED';

end % StopISI

%% ----- STOPRECORDING  - stop recording

function StopRecording

t = now;
timings.STOP(state.TIDX) = t;
timings.COMPLETED(state.TIDX) = state.COMPLETED;
dt = datetime(table2array(timings(state.TIDX,3:4)),'convertfrom','datenum');
state.TIDX = state.TIDX + 1;	% increment completed trial counter
if state.RECORDING
	dt = seconds(time(between(dt(1),dt(2))));
	fprintf('[%s] %6.3f s\n', char(datetime(t,'convertfrom','datenum','format','HH:mm:ss.SSS')),dt)
	if state.RECORDING == 1 || strcmp('BLKEND',state.CURTRIAL.TYPE)
		ACQHW = state.ACQHW;
		for mi = length(ACQHW) : -1 : 1, ACQHW(mi).FUN('STOP'); end
		state.RECORDING = 0;			% set not recording
		if strcmp('BLKEND',state.CURTRIAL.TYPE), state.CTL.AUTO.Value = 0; end
	end
end

end % StopRecording

%% ----- UPDATEDISPLAY  - update browser(s) with current stimulus HTML

function UpdateDisplay(stim, si)

if state.RECORD && ~isempty(si), timings(state.TIDX,si+5) = table(now); end	% STIMn onset

% format HTML
if strcmp('BLKPAD',stim.HTML)		% BLKREC start, end trials
	html = sprintf('<BODY bgcolor="black">%s</BODY>',stim.HTML);
elseif isempty(state.INFO.CSS)
	html = sprintf('<BODY bgcolor="white">%s</BODY>',stim.HTML);
else
	html = sprintf('<HEAD><STYLE TYPE="text/css">%s</STYLE></HEAD><BODY>%s</BODY>',state.INFO.CSS,stim.HTML);

% anything more complicated than text seems inaccessible to uihtml unless loaded from a file
	if ~isempty(stim.EXTRA)
		if isfield(stim.EXTRA,'WRITEHTML')
			fid=fopen(htmlFname,'wt'); fprintf(fid,'%s',html); fclose(fid);
			html = htmlFname;
		end
		if isfield(stim.EXTRA,'HANDLER') && ~isempty(stim.EXTRA.HANDLER)
			state.STIMFUN = str2func(stim.EXTRA.HANDLER.FUN);
			state.STIMARG = state.STIMFUN('START',stim.EXTRA.HANDLER.ARG);	% returns handle
		else
			state.STIMFUN = [];
			state.STIMARG = [];
		end
	end
end

% pump to browser(s)
if isempty(stim.PART)
	set(state.BROWSER,'HTMLSource',html)
else
	state.BROWSER(stim.PART).HTMLSource = html;
end

end % UpdateDisplay

%% ----- UPDATEPROGRESS  - update progress bar every 50 ms through trial

function UpdateProgress(t, ~, h, nReps)

h.XData = [0 0 1 1] * t.TasksExecuted/nReps;

end % UpdateProgress

end % marta
