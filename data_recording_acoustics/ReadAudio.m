function [s,sr,hdr] = ReadAudio(fn)
%READAUDIO  - read audio format data
%
%	usage:  [s,sr,hdr] = ReadAudio(fn)
%
% returns signal S from AUDIO format file FN.audio as [nSamps x nChan] 'single' samples
% optionally returns sampling rate SR (Hz) and file header HDR
%
% see also AUDIO2WAV

% v2 AUDIO files have a 16 byte header:
%   bytes 1:4 are signature char 'AUDIOBIN' (char)
%   bytes 5:8 gives audio file version number (single)
%   bytes 9:12 gives sampling rate in Hz (single)
%	bytes 13:16 gives number of channels (single)
% the remainder encode audio data with channels interleaved (if more than one)
% in little-endian single format
%
% the current version is v2

% mkt 08/23

if nargin < 1, help ReadAudio; return; end

AUDIO_SIGNATURE = 'AUDB';	% audio file signature

[p,f,e] = fileparts(fn);
if isempty(e), fn = fullfile(p,[f,'.audio']); end

% open file
fid = fopen(fn, 'rb', 'l');
assert(fid > 0, 'unable to open %s', fn);

% attempt to read header
h = single(fread(fid,4,'single'));
sig = char(typecast(h(1),'uint8'));
if strcmp(sig,AUDIO_SIGNATURE)		% v2+
	sr = h(3);
	hdr = struct('VERSION',h(2), 'SRATE',sr, 'NCHAN',h(4));
else
	error('error parsing header for %s', fn)
end

% load data
s = fread(fid,inf,'single')';
fclose(fid);
if hdr.NCHAN > 1, s = reshape(s,hdr.NCHAN,[])'; end
