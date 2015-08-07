function [gs ap ir] = cigar2gappedsequence(seqs, cigars, varargin)
%CIGAR2GAPPEDSEQUENCE convert unaligned sequences to aligned sequences using CIGAR strings.
%
%   GS = CIGAR2GAPPEDSEQUENCE(S, C) converts the unaligned sequences in
%   cell array S into gapped sequences using the information stored in the
%   cell array C containing CIGAR strings. GC, S and C are all cell arrays of
%   the same size.
%
%   CIGAR2GAPPEDSEQUENCE(S, C, 'GapsInRef', TF) specifies whether positions
%   corresponding to gaps in the reference sequence should be included or
%   not. Default is false.
%
%   CIGAR2GAPPEDSEQUENCE(S, C, 'SoftClipping', TF) specifies whether
%   positions corresponding to soft clippling ends should be included or
%   not. Default is false.
%
%   CIGAR2GAPPEDSEQUENCE(S, C, 'Quality', TF) specifies whether the input
%   sequences are the base qualities. Default is false. When TF is set to
%   true the character symbol used for representing deletions (D) and
%   skipped positions (N) is char(0).
%
%   [GS AP] = CIGAR2GAPPEDSEQUENCE(S, C) returns the anchor position AP in
%   the gapped sequences GS. An anchor position is the position in a gapped
%   alignment that should be aligned to its respective start position in 
%   the hypothetical reference when combining the gapped sequences into a
%   multiple alignment. Note: when using the default inputs for GapsInRef
%   and SoftClipping no soft-clippling (S), insert (I) or padding (P)
%   symbols occur before the anchor position, therefore, the AP output is 
%   a column vector of ones.
%
%   [GS AP IR] = CIGAR2GAPPEDSEQUENCE(S, C) returns the number of gaps that
%   need to be inserted in the hypothetical ungapped reference for every
%   sequence. When GapsInRef is set to false, IR contains vectors with only
%   zeros.

%   Copyright 2010 The MathWorks, Inc.


mfilename = 'cigar2align';

%=== Check inputs
bioinfochecknargin(nargin, 2, mfilename)

if ~iscellstr(seqs) || isempty(seqs)
	error(message('bioinfo:cigar2align:InvalidSequenceInput'));
end

if ~iscellstr(cigars)|| isempty(cigars)
	error(message('bioinfo:cigar2align:InvalidCigarInput'));
end

if numel(seqs) ~= numel(cigars)
	error(message('bioinfo:cigar2align:InvalidInputSize'));
end

%=== Parse PVP
[GapsInRef, SoftClipping, Quality, Accelerator] = parse_inputs(varargin{:});

if Accelerator
    if nargout > 2
        [gs ap ir] = bioinfoprivate.cigar2gappedsequencemex(seqs,cigars,SoftClipping,GapsInRef,Quality);
    elseif nargout > 1
        [gs ap] = bioinfoprivate.cigar2gappedsequencemex(seqs,cigars,SoftClipping,GapsInRef,Quality);
    else
        gs = bioinfoprivate.cigar2gappedsequencemex(seqs,cigars,SoftClipping,GapsInRef,Quality);
    end
else
    %=== Initialize variables
    if Quality
        gapSymbol = char(0);
        skipSymbol = char(0);
    else
        gapSymbol  = '-';
        skipSymbol = '.';
    end
    sV('DIMNPSH') = uint8([1 2 3 4 5 6 7]);
    sH = [false true true false false true];
    sI = [true false true false false false];
    sK = [false false false true false false];
    sG = [false true false false true false];
    if GapsInRef
        if SoftClipping
            sL = [false true true false false true];
            sJ = [true true true true true true];
            sS = [false false false false false true];
        else
            sL = [false true true false false false];
            sJ = [true true true true true false];
        end
    else
        if SoftClipping
            sL = [false false true false false true];
            sJ = [true false true true false true];
            sS = [false false false false false true];
        else
            sL = [false false true false false false];
            sJ = [true false true true false false];
        end
    end
    
    n = numel(cigars);
    
    gs = cell(n,1);
    ap = ones(n,1);
    ir = cell(n,1);
    
    validCigar = false; % keep track of whether a valid cigar string has been found
    T = regexprep(cigars,'[^MDIPSHN]','');
    
    for i = 1:n
        t = sV(T{i});
        lent = length(t);
        %=== determine where and how many gaps in the current read (D and P)
        if ~lent % uninformative cigar
            gs{i} = '*';
            ap(i) = 1;
            ir{i} = 0;
        else
            validCigar = true;
            try
                [v lenv] = sscanf(cigars{i}, '%d%*c');
                
                if lenv~=lent
                    error('bioinfo:cigar2align:InvalidCigar','Different length')
                end
                %=== remove hard clipped positions H
                hard = t==7;  %'H'
                v(hard) = [];
                t(hard) = [];
                
                %=== idx -> indices to the cigar symbols in t for every
                %    position in the gapped sequence (size of idx is still
                %    congruent to the CIGAR string, i.e. contains gaps in
                %    reference and soft clipping )
                acv = accumarray(cumsum(v),1);
                idx = cumsum(acv)-acv+1;
                
                %=== tt -> a cigar symbol for every position in the gapped
                %    sequence (size of tt is still congruent to the CIGAR
                %    string, i.e. contains gaps in reference and soft clipping)
                tt = t(idx);
                
                %=== lidx -> logical vector indicating the positions in the
                %    gapped sequence that need to be filled from characters in
                %    the sequence, the logical values are set depending on
                %    GapsInRef and SoftClipping, but the size is still
                %    congruent to the full CIGAR string
                lidx = sL(tt);
                
                %=== jidx -> logical vector indicating the positions that will
                %    remain in the gapped sequence according with GapsInRef and
                %    SoftClipping. For example if GapsInRef and SoftClipping
                %    are both set to true, lidx will be all true.
                jidx = sJ(tt);
                
                %=== hidx -> logical vector indicating -ALL- the candidate
                %    positions that could take a symbol from the sequence
                %    (independently of GapsInRef and SoftClipping)
                hidx = sH(tt);
                
                %=== Initialize to all gapSymbol
                gst = gapSymbol(uint8(jidx(jidx)));
                
                %=== Copy symbols from sequence
                gst(lidx(jidx)) = seqs{i}(lidx(hidx));
                
                %=== Insert skipSymbols which now are marked as the gapSymbol
                %    since we use that for init
                kidx = sK(tt);
                gst(kidx(jidx)) = skipSymbol;
                
                %=== Soft clippling bases are changed to lower case
                if SoftClipping && ~Quality
                    sidx = sS(tt);
                    gst(sidx(jidx)) = lower(gst(sidx(jidx)));
                end
                
                %=== Store gapped sequence
                gs{i} = gst;
                
                if GapsInRef || SoftClipping
                    
                    %=== Find anchor
                    iidx = sI(tt);
                    ap(i) = find(iidx(jidx),1);
                    
                    %=== gidx -> logical vector indicating the positions that
                    %    insert gaps in the reference
                    gidx = sG(tt);
                    q = cumsum(gidx(jidx));
                    
                    % Gaps that need to be inserted in the reference positions that
                    % exist in the output gapped sequence
                    ir{i} = diff(q(~gidx(jidx)));
                    
                else
                    ir{i} = 0;
                end
                
            catch ME
                error('bioinfo:cigar2align:InvalidCigar',...
                    'Invalid CIGAR string found.');
            end
        end
    end
    
    if ~validCigar
        error(message('bioinfo:cigar2align:AllUninformativeCigar'));
    end

end % if else Accelerator

%--------------------------------------------------------------------------

function [GapsInRef, SoftClipping, Quality, Accelerator] = parse_inputs(varargin)
% Parse input PV pairs.

mfilename = 'cigar2align';

%=== defaults
SoftClipping = false;
GapsInRef = false;
Quality = false;
Accelerator = true;

%=== check for the right number of inputs
if rem(nargin,2) == 1
	error(message('bioinfo:cigar2align:IncorrectNumberOfArguments', mfilename));
end

%=== allowed parameters
okargs = {'GapsInRef', 'SoftClipping', 'Quality', 'Accelerator'};

%=== parse inputs
for j = 1:2:nargin
	pname = varargin{j};
	pval = varargin{j+1};
	
	if ~ischar(pname)
		error(message('bioinfo:cigar2align:InvalidParameter'));
	end
	
	k = find(strncmpi(pname, okargs, numel(pname)));
	if isempty(k)
		error(message('bioinfo:cigar2align:UnknownParameterName', pname));
	elseif length(k)>1
		error(message('bioinfo:cigar2align:AmbiguousParameterName', pname));
	else
        switch(k)
            case 1 % GapsInRef
                GapsInRef = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            case 2 % SoftClipping
                SoftClipping = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            case 3 % Quality
                Quality = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
            case 4 % Accelerator
                Accelerator = bioinfoprivate.opttf(pval, okargs{k}, mfilename);
        end
	end
end
