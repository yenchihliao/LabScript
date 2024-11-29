function [trials,info] = ParseExpFile(fName)
%PARSEEXPFILE  - parse marta experiment file in XML format
%
%	usage:  [trials,info] = ParseExpFile(fName)
%
% this prcedure reads marta experiment specification file FNAME (default .xml extension) and
% returns array-or-structs TRIALS and experiment INFO used by MARTA to sequence an experiment
%
% for a description of the XML format expected see MARTAEXPFILESPEC.pdf
%
% example:
%  [trials,info] = ParseExpFile('test')
%  marta(trials, info)

% mkt 07/23 rewrite

if nargin < 1, help ParseExpFile; return; end
[p,f,e] = fileparts(fName);
if isempty(e), fName = fullfile(p,[f,'.xml']); end

% parse the XML document
tree = xmlread(fName);

% construct a nested structure by recursion over child nodes
% cf. https://www.mathworks.com/help/matlab/ref/xmlread.html
s = ParseChildNodes(tree);

%% check for expected general structure
if length(s)>1 || ~strcmpi('SESSION',s.NAME) || isempty(s.ATTRIBUTES), error('missing valid SESSION element'); end
nodes = {s.CHILDREN.NAME};
k = find(strcmpi('INFO',nodes));
if isempty(k), error('missing valid INFO element'); end
INFO = s.CHILDREN(k).CHILDREN;
k = find(strcmpi('DEFS',nodes));
if isempty(k), DEFS = []; else, DEFS = s.CHILDREN(k).CHILDREN; end
k = find(strcmpi('ORDER',nodes));
if isempty(k), error('missing valid ORDER element'); end
ORDER = s.CHILDREN(k).CHILDREN;

%% construct INFO variable
info.SESSION = s.ATTRIBUTES.VALUE;
nodes = {INFO.NAME};

k = find(strcmpi('PREFIX',nodes));
if isempty(k), info.PREFIX = info.SESSION; else, info.PREFIX = INFO(k).ATTRIBUTES.VALUE; end

k = find(strcmpi('LOG',nodes));
if isempty(k), info.LOG = info.SESSION; else, info.LOG = INFO(k).ATTRIBUTES.VALUE; end
[~,f,e] = fileparts(info.LOG);
if isempty(e), info.LOG = fullfile(sprintf('%s.log',f)); end

k = find(strcmpi('TIMINGS',nodes));
if ~isempty(k), info.TIMINGS = INFO(k).ATTRIBUTES.VALUE; end

idx = find(strcmpi('PARTICIPANT',nodes));
if ~isempty(idx)
	info.PARTICIPANT = [];
	for k = 1 : length(idx)
		a = MapAttributes(INFO(idx(k)));
		if isempty(info.PARTICIPANT), info.PARTICIPANT = a; else, info.PARTICIPANT(end+1) = a; end
	end
end
k = find(strcmpi('ACQHW',nodes));
if isempty(k)
	info.ACQHW = []; 
else 
	ACQHW = INFO(k).CHILDREN;
	idx = find(strcmpi('MODULE',{ACQHW.NAME}));
	if isempty(idx)
		modules = [];
	else
		for mi = 1 : length(idx)
			a = MapAttributes(ACQHW(idx(mi)));
			modules(mi) = struct('NAME',a.NAME,'PARAMS',[]);
			params = [];
			p = ACQHW(idx(mi)).CHILDREN;
			for xi = 1 : length(p)
				switch upper(p(xi).NAME)
					case {'#TEXT','REM'}  % ignore
					otherwise
						pp = struct('NAME',p(xi).NAME,'ARGS',MapAttributes(p(xi)));
						if isempty(params), params = pp; else, params(end+1) = pp; end
				end
			end
			modules(mi).PARAMS = params;
		end
	end
	info.ACQHW = modules;
end
k = find(strcmpi('CSS',nodes));
if ~isempty(k)
	CSS = INFO(k).CHILDREN;
	if length(CSS) > 1		% if more than one offspring (as with embedded REM) use longest
		[~,k] = max(cellfun(@length,{CSS.CONTENT}));
		CSS = CSS(k);
	end
	info.CSS = replace(CSS.CONTENT,newline,'');	% kill CR
end

% define blocks if necessary
if ~isempty(DEFS)
	DB = DEFS(strcmpi('DEFBLOCK',{DEFS.NAME}));
	blocks(1,length(DB)) = struct('NAME',[],'CODE',[],'NREPS',[],'RAND',[],'DUR',[],'BLKREC',[],'BLKDLY',[],'TRIALS',[]);
	for bi = 1 : length(DB)
		trials = [];
		blocks(bi).NREPS = 1;
		blocks(bi).RAND = 0;
		blocks(bi).BLKREC = 0;
		blocks(bi).BLKDLY = 0;
		attr = DB(bi).ATTRIBUTES;
		for ai = 1 : length(attr)
			switch upper(attr(ai).NAME)
				case {'CODE','NAME'}, blocks(bi).(upper(attr(ai).NAME)) = attr(ai).VALUE;
				otherwise, blocks(bi).(upper(attr(ai).NAME)) = str2double(attr(ai).VALUE); 
			end
		end
		if isempty(blocks(bi).CODE), blocks(bi).CODE = blocks(bi).NAME; end
		temps = DB(bi).CHILDREN(find(strcmpi('TEMPLATE',{DB(bi).CHILDREN.NAME})));
		tokens = DB(bi).CHILDREN(find(strcmpi('TOKEN',{DB(bi).CHILDREN.NAME})));
		
% expand tokens
		for ti = 1 : length(tokens)
			tid = 1;	% default
			dur = [];	% unspecified
			args = {tokens(ti).CHILDREN.CONTENT};% a0 %originally >1 leading to no space betwn words, changed to >2 by Fangchi 20241125	
			if length(regexp(args{1},'\s','split')) > 2, error('illegal whitespace in token %d block %d content (a0)', ti, bi); end
			for xi = 1 : length(tokens(ti).ATTRIBUTES)
				switch upper(tokens(ti).ATTRIBUTES(xi).NAME)
					case 'TID', tid = str2double(tokens(ti).ATTRIBUTES(xi).VALUE);
					case 'DUR', dur = str2double(tokens(ti).ATTRIBUTES(xi).VALUE);
					otherwise, args{end+1} = tokens(ti).ATTRIBUTES(xi).VALUE;
				end
			end
			t = temps(tid);
			if ~isempty(t.ATTRIBUTES), dur = str2double(t.ATTRIBUTES.VALUE); end
			s = t.CHILDREN(strcmpi('STIMULUS',{t.CHILDREN.NAME}));

% expand stimuli
			stim = [];
			for si = 1 : length(s)
				sa = MapAttributes(s(si));
				delay = 0; 
				part = [];		% display to both partners if relevant
				if ~isempty(sa)
					if isfield(sa,'DELAY')
						delay = sa.DELAY; 	% two element vector specifies range for randomized delay
						sa = rmfield(sa,'DELAY');
					end
					if isfield(sa,'PART')
						part = sa.PART;
						sa = rmfield(sa,'PART');
					end
					if isfield(sa,'VOX')	% VOX threshold mapped to negative DELAY
						if isempty(info.ACQHW) || ~strcmp('acq_AUDIO',{info.ACQHW.NAME}), error('VOX requires acq_AUDIO plugin'); end
						delay = -sa.VOX;
						sa = rmfield(sa,'VOX');
					end
				end
				k = find(strcmpi('HTML',{s(si).CHILDREN.NAME}),1);
				if isempty(k), error('missing required HTML element (block %d stimulus %d)', bi, si); end
				HTML = s(si).CHILDREN(k).CHILDREN.CONTENT;
				for ai = 1 : length(args)	% substitute argument placeholders
					HTML = replace(HTML,sprintf('@%d',ai-1),args{ai});
				end
				handler = [];
				k = find(strcmpi('HANDLER',{s(si).CHILDREN.NAME}));
				if ~isempty(k) 
					a = MapAttributes(s(si).CHILDREN(k));
					arg = a.ARG;
					for ai = 1 : length(args)	% substitute argument placeholders
						arg = replace(arg,sprintf('@%d',ai-1),args{ai});
					end
					handler = struct('FUN',a.FUNCTION,'ARG',arg);
				end
				if isempty(sa) && isempty(handler)
					extra = [];
				elseif isempty(sa)
					extra = struct('HANDLER',handler);
				else
					extra = sa;
					extra.HANDLER = handler;
				end
				ss = struct('PART',part,'DELAY',delay,'HTML',HTML,'EXTRA',extra);
				if isempty(stim), stim = ss; else, stim(end+1) = ss; end
			end
			if ~all(diff(sort(cell2mat({stim.DELAY})))), error('identical delays found for stimuli in block %s token %d',blocks(bi).NAME,ti); end
			fn = sprintf('%s_%s_%s_BXX_RXX',info.PREFIX,args{1},blocks(bi).CODE);
			if blocks(bi).BLKREC, rt = 'BLKREC'; else, rt = 'RECORD'; end
			tt = struct('TYPE',rt,'FNAME',fn,'PROMPT',[],'DUR',dur,'STIM',stim);
			if isempty(trials), trials = tt; else, trials(end+1) = tt; end
		end
		if blocks(bi).BLKREC,	% prefix and append BLKREC trials
			fn = sprintf('%s_BLKREC_%s_BXX',info.PREFIX,blocks(bi).CODE);
			dur = blocks(bi).BLKDLY/1000 + .5;
			s = struct('PART',[],'DELAY',blocks(bi).BLKDLY,'HTML','BLKPAD','EXTRA',[]);
			t0 = struct('TYPE','BLKSTRT','FNAME',fn,'PROMPT',[],'DUR',dur,'STIM',s);
			t1 = struct('TYPE','BLKEND','FNAME',[fn,'_END'],'PROMPT',[],'DUR',dur,'STIM',s);
			blocks(bi).TRIALS = [t0,trials,t1];
		else
			blocks(bi).TRIALS = trials;
		end
	end
end

% expand trials
trials = [];
nn = {ORDER.NAME};
if ~isempty(DEFS)
	BLOCK = string({blocks.NAME}');
	COUNT = zeros(length(BLOCK),1);
	bc = table(BLOCK,COUNT);	% tracks block repetitions
end
for ni = 1 : length(nn)
	switch upper(nn{ni})
		case 'PAUSE'
			a = MapAttributes(ORDER(ni));
			if isfield(a,'DUR'), dur = a.DUR; else, dur = []; end
			s = struct('PART',[],'DELAY',[],'HTML',GetCDATA(ORDER(ni)),'EXTRA',[]);
			tt = struct('TYPE','PAUSE','FNAME','PAUSE','PROMPT',a.PROMPT,'DUR',dur,'STIM',s);
			if isempty(trials), trials = tt; else, trials(end+1) = tt; end
		case 'TRIAL'
			a = MapAttributes(ORDER(ni));
			idx = find(strcmpi('STIMULUS',{ORDER(ni).CHILDREN.NAME}));
			if isempty(idx), error('explicit TRIAL missing STIMULUS'); end
			ss = ORDER(ni).CHILDREN(idx);
			stim = [];
			for si = 1 : length(ss)
				sa = MapAttributes(ss(si));
				if isfield(sa,'DELAY'), delay = sa.DELAY; else, delay = 0; end
				if isfield(sa,'PART'), part = sa.PART; else, part = []; end
				k = find(strcmpi('HTML',{ss(si).CHILDREN.NAME}),1);
				if isempty(k), error('explicit TRIAL missing STIMULUS HTML'); end
				HTML = ss(si).CHILDREN(k).CHILDREN.CONTENT;
				k = find(strcmpi('HANDLER',{ss(si).CHILDREN.NAME}),1);
				if isempty(k) 
					extra = []; 
				else
					xa = MapAttributes(ss(si).CHILDREN(k));
					extra = struct('HANDLER', struct('FUN',xa.FUNCTION,'ARG',xa.ARG));
				end
				s = struct('PART',part,'DELAY',delay,'HTML',HTML,'EXTRA',extra);
				if isempty(stim), stim = s; else, stim(end+1) = s; end
			end
			tt = struct('TYPE','RECORD','FNAME',a.FNAME,'PROMPT',a.PROMPT,'DUR',a.DUR,'STIM',stim);
			if isempty(trials), trials = tt; else, trials(end+1) = tt; end
		case 'BLOCK'
			a = MapAttributes(ORDER(ni));
			if isempty(DEFS), error('BLOCK %s has not been defined',a.NAME); end
			bi = find(strcmp(a.NAME,bc.BLOCK));
			n = bc.COUNT(bi) + 1;
			bc.COUNT(bi) = n;
			tt = ExpandTrials(blocks(bi),n);
			if isempty(tt(1).DUR), tt(1).DUR = blocks(bi).DUR; end
			if isempty(trials), trials = tt; else, trials = [trials,tt]; end
		otherwise  % ignore
	end
end

% add trial number
for ti = 1 : length(trials)
	trials(ti).FNAME = sprintf('%s_%04d', trials(ti).FNAME, ti);
end


%% ----- local functions

%% EXPANDTRIALS - expand block trials
function trials = ExpandTrials(b, n)	% block, block repetition count
t = b.TRIALS;
bn = sprintf('B%02d',n);
for ti = 1 : length(t), t(ti).FNAME = replace(t(ti).FNAME,'BXX',bn); end
if b.BLKREC
	t0 = t(1);
	t1 = t(end);
	t([1 end]) = [];
end
trials = [];
for ri = 1 : b.NREPS
	tt = t;
	rn = sprintf('R%02d',ri);
	for ti = 1 : length(t)
		tt(ti).FNAME = replace(t(ti).FNAME,'RXX',rn); 
		for si = 1 : length(tt(ti).STIM)
			if ~isscalar(tt(ti).STIM(si).DELAY)	% two element vector specifies range for randomized delay
				tt(ti).STIM(si).DELAY = randi(tt(ti).STIM(si).DELAY,1);
			end
		end
	end
	if b.RAND, tt = tt(randperm(length(t))); end
	if isempty(trials), trials = tt; else, trials = [trials , tt]; end
end
if b.BLKREC
	trials = [t0 , trials , t1];
end

%% GETCDATA - get CDATA content
function s = GetCDATA(node)
k = find(strncmp('#cdata',{node.CHILDREN.NAME},6));
if isempty(k), error('malformed CDATA section'); end
s = node.CHILDREN(k).CONTENT;

%% MAKESTRUCTFROMNODE - create structure from node info
function nodeStruct = MakeStructFromNode(theNode)
nodeStruct = struct(                        ...
				'NAME', char(theNode.getNodeName),       ...
				'ATTRIBUTES', ParseAttributes(theNode),  ...
				'CONTENT', '',                           ...
				'CHILDREN', ParseChildNodes(theNode));
if any(strcmp(methods(theNode), 'getData'))
	nodeStruct.CONTENT = char(theNode.getData); 
else
	nodeStruct.CONTENT = '';
end

%% MAPATTRIBUTES - map node attributes to struct
function a = MapAttributes(node)
attr = node.ATTRIBUTES;
if isempty(attr), a = []; return; end
for ai = 1 : length(attr)
	v = str2num(attr(ai).VALUE);
	if isempty(v), v = attr(ai).VALUE; end
	a.(upper(attr(ai).NAME)) = v;
end

%% PARSEATTRIBUTES - create attributes structure from node info
function attributes = ParseAttributes(theNode)
attributes = [];
if theNode.hasAttributes
	theAttributes = theNode.getAttributes;
	numAttributes = theAttributes.getLength;
	allocCell = cell(1, numAttributes);
	attributes = struct('NAME',allocCell,'VALUE',allocCell);
	for count = 1 : numAttributes
		attrib = theAttributes.item(count-1);
		attributes(count).NAME = char(attrib.getName);
		attributes(count).VALUE = char(attrib.getValue);
	end
end

%% PARSECHILDNODES - recurse over node children
function children = ParseChildNodes(theNode)
children = [];
if theNode.hasChildNodes
	childNodes = theNode.getChildNodes;
	numChildNodes = childNodes.getLength;
	allocCell = cell(1, numChildNodes);
	children = struct('NAME',allocCell,'ATTRIBUTES',allocCell,'CONTENT',allocCell,'CHILDREN',allocCell);
	for count = 1 : numChildNodes
		theChild = childNodes.item(count-1);
		children(count) = MakeStructFromNode(theChild);
	end
end

