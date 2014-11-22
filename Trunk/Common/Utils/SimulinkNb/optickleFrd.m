function frdOut = optickleFrd(opt,f,varargin)
    % creates a frd object from optickle object, for use in SimulinkNB
    % Syntax: optickleFrd(opt,f) 
    % or optickleFrd(opt,f,sigAC) if sigAC is precomputed
    %
    % If not called from a FlexTf block, you must add a cell array of drive
    % names and probe names. Example:
    % optickleFrd(opt,f,{'W PM'},{'W REFL I'});
    % or optickleFrd(opt,f,sigAC,{'W PM'},{'W REFL I'});

    parseBlock = 0;
    % see if we have a sigAC
    if nargin<3 || iscell(varargin{1})
        % compute from optickle
        [~,~,sigAC] = tickle(opt,[],f);
        if ~isempty(varargin)
            drives = varargin{1};
            probes = varargin{2};
        else 
            parseBlock = 1;
        end
    else
        sigAC = varargin{1};
        if length(varargin)>1
            drives = varargin{2};
            probes = varargin{3};
        else
            parseBlock = 1;
        end
    end
    
    % get the drives and probes from the simulink block
    if parseBlock
        % from http://www.mathworks.com/matlabcentral/answers/12810
        h_Inports = find_system(gcbh,'FindAll','On','SearchDepth',1,'BlockType','Inport');
        drives = get(h_Inports,'Name');
        
        h_Outports = find_system(gcbh,'FindAll','On','SearchDepth',1,'BlockType','Outport');
        probes = get(h_Outports,'Name');
    end
    
    if isempty(drives) || isempty(probes)
        error('optickleFrd:missingdriveprobes','No probe or drive names provided, nor found in simulink block');
    end
    
    driveIndex = cellfun(@(drive) getDriveNum(opt,drive),drives);
    probeIndex = cellfun(@(probe) getProbeNum(opt,probe),probes);
    
    frdOut = frd(sigAC(probeIndex,driveIndex,:),f,'Units','Hz');
end
