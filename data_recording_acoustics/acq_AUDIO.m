function r = acq_AUDIO(action, varargin)
%ACQ_AUDIO  - MARTA generic audio recording handler
%
% records mono or stereo audio using audiorecorder mechanism
% note that input device must support two channel recording for stereo
%
% to view available AUDIO input devices and obtain device ID  use
%  d = audiodevinfo; d = d.input; [{d.ID}',{d.Name}']
%
% expects two optional XML-specified parameters passed through INFO; e.g.
%   <ACQHW>
%     <MODULE name="acq_AUDIO">
%       <AUDIO devidx="2" nchan="1" srate="22050" record="0" />
%       <DISPLAY color="0 0 1" />
%     </MODULE>
%   </ACQHW>
%
% where all parameter attributes are optional:
%	DEVIDX defaults to "-1"
%	NCHAN to "1"
%	SRATE to "22050"
%	RECORD to "1" (record files; 0 displays audio for monitoring but does not recode)
%	COLOR = "0 0 0"
%
% AUDIO parameters may also be specified as a MARTA command line argument; e.g.
% marta(trials,info,'AUDIO',struct('DEVIDX',2,'NCHAN',2))

% mkt 10/18
% mkt 08/23 rewrite for marta 2.0

vers = 'mkt 08/23 v2.0';	% current version

%% persistent parameters
persistent sRate		% sampling rate
persistent nChan		% # channels
persistent idx			% sample index
persistent offset		% offset into animated line
persistent ah			% audio plot axis
persistent lh			% audio plot line handle(s)
persistent rh			% audiorecorder handle
persistent fName		% saved trial filename
persistent saveFile		% true if recording files, false if monitoring only
persistent maxPts		% max # plotted samples
persistent VOX			% tests for VOX-activation if nonempty
persistent win			% VOX envelope window length

%% branch by action
switch action

%% ABOUT  - show current marta and handler version info
	case 'ABOUT'
		fprintf('acq_AUDIO  - MARTA generic audio recording handler (%s)\n', vers)

%% CLOSE  - shutdown handler
	case 'CLOSE'
		stop(rh)

%% INIT  - initialize handler (varargin: params)
	case 'INIT'
		params = varargin{1};
		
% parse expected audio parameters
		audio = struct('DEVIDX',-1,'NCHAN',1,'SRATE',22050,'RECORD',1);
		k = find(strcmp('AUDIO',{params.NAME}),1);
		if ~isempty(k)
			audioIn = params(k).ARGS;
			fn = upper(fieldnames(audioIn));
			for k = 1 : length(fn), audio.(fn{k}) = audioIn.(fn{k}); end
		end
		k = find(strcmp('DISPLAY',{params.NAME}),1);
		if isempty(k), col = 'k'; else, col = params(k).ARGS.COLOR; end
		                
% initialize audio
		if audio.DEVIDX >= 0
			d = audiodevinfo;
			k = find(audio.DEVIDX == cell2mat({d.input.ID}));
			if isempty(k), error('audio input device %d not available', audio.DEVIDX), end
			rh = audiorecorder(audio.SRATE,16,audio.NCHAN,audio.DEVIDX);	% 16 bits
			fprintf('using %s (%d) for input audio', d.input(k).Name, audio.DEVIDX)
		else
			rh = audiorecorder(audio.SRATE,16,audio.NCHAN);
			fprintf('using default audio input (DeviceID = -1)')
		end
		if audio.RECORD, fprintf('\n'); else, fprintf('  (monitoring only)\n'); end
		set(rh, 'TimerPeriod',0.02, 'TimerFcn',@AudioCB);		% 20 ms period
		sRate = audio.SRATE;
		nChan = audio.NCHAN;
		saveFile = audio.RECORD;
		win = round(25*sRate/1000);		% VOX envelope window (samps)
		
% init audio monitoring window
		maxPts = 1e6;	% max # plotted samples
		fh = figure('name','Audio', ...
						'tag','MARTA', ...
						'color','w', ...
						'menubar','none', ...
						'closeRequestFcn',{@marta,'CLOSE'}, ...
						'numberTitle','off');
		m = get(0,'monitorPositions'); k = find(m(:,1)==1); m = m(k,:);
		p = fh.Position; fh.Position = [m(3)-p(3) p(2)+50 p(3) p(4)-50];

		if nChan == 1,
			ah = axes('position',[.01 .03 .98 .94],'xtick',[],'ytick',[],'box','on');
			lh = animatedline('color',col,'MaximumNumPoints',maxPts);
		else,
			ah = axes('position',[.05 .54 .93 .44],'xtick',[],'ytick',[],'box','on');
			lh = animatedline('color',col,'MaximumNumPoints',maxPts);
			ah(2) = axes('position',[.05 .1 .93 .44],'xtick',[],'ytick',[],'box','on');
			lh(2) = animatedline('color',col,'MaximumNumPoints',maxPts);
		end;

% return configuration
		r = struct('WINDOW',struct('NAME','AUDIO','HANDLE',fh));
		
%% START  - init recording (varargin: fName, dur)
	case 'START'
        fName = varargin{1};	% saved trial filename
        dur = varargin{2};		% recording duration (secs)   

% init audio plotting
        clearpoints(lh(1));
		if nChan > 1, clearpoints(lh(2)); end;
        xl = [1 min([floor(sRate*dur),maxPts])];
        set(ah,'xlim',xl,'ylim',[-1 1]);
        drawnow;
        
% init audio recording
		idx = 0;				% sample index
		offset = 0;				% offset into animated line
		record(rh, dur);

%% STOP  - terminate recording (may be truncated)
	case 'STOP'
		stop(rh);
		s = getaudiodata(rh);
		s = s * 2;				% map to +/-1
		buf = s(idx+1:end,:);
		if ~isempty(buf)
			x = [1:size(buf,1)]' + idx;
			addpoints(lh(1), x, buf(:,1));			% update plot (mono/L)
			if nChan > 1, addpoints(lh(2), x, buf(:,2)); end;			% stereo (R)
		end

% write acquired samples
		if saveFile, audiowrite([fName,'.wav'],s,sRate); end

%% VOX  - set awaiting voice activation
	case 'VOX'
		VOX = varargin{1};		% threshold (% of RMS envelope)

%% error trap
	otherwise
		error('unrecognized action (%s) in ACQ_AUDIO', action);
end	% switch

%% ========== LOCAL FUNCTIONS ============================================================

%% ----- AUDIOCB  - audio monitoring timer callback
function AudioCB(h,~)

buf = getaudiodata(h);			% retrieve audio samples
ns = size(buf,1) - idx;			% # new samples
if ns < 1, return; end;			% empty buffer
buf(1:idx,:) = [];				% trim to new samples
buf = buf*2;					% map to +/-1
x = [1:ns]' + idx - offset;
if x(end) > maxPts
	offset = offset + maxPts;
	clearpoints(lh(1));
	if nChan > 1, clearpoints(lh(2)); end;
	k = find(x<=offset);
	buf(k,:) = [];
	x(k) = [];
	x = x - offset + 1;
end;
idx = idx + ns;
addpoints(lh(1), x, buf(:,1));						% update plot (mono/L)
if nChan > 1, addpoints(lh(2), x, buf(:,2)); end;	% stereo (R)

if ~isempty(VOX)
	if any(envelope(buf,win,'peak') > VOX)
		marta('VOX')
		VOX = [];
	end
end

end	% AudioCB

end % acq_audio

