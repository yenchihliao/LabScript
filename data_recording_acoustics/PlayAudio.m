function h = PlayAudio(action, varargin)
%PLAYAUDIO  - load & play audio file (MARTA stimulus handler)

switch action
	case 'START', 
		fName = varargin{1};
		[p,f,e] = fileparts(fName);
		if isempty(e), fName = fullfile(p,[f,'.wav']); end;
		try
			[s,sr] = audioread(fName);
		catch
			fprintf('Warning:  %s not found\n', fName);
			h = [];
			return
		end
		h = audioplayer(s,sr);
		play(h)
	case 'STOP'
		h = varargin{1};
		if ishandle(h), stop(h); end
end
