classdef Session
    % Data for a given session.
    % The idea here is that by using a class, one can easily manipulate the
    % data files (add, resort spikes, etc.) using the structures functions.
    % "Day" structure is a container for a days worth of sessions.
    
    properties
        % Information about build files
        navDat % Navigation data file
        neuralDat % BR file extensions
        numchannels = 36;
        % Data
        shape
        time % vector of time data, sectioned by trials
        posdat % vector of position data, sectioned by trials
        dirdat % vector of direction data, sectioned by trials
        eyedat % vector of eye data, sectioned by trials
        LFPdat % vector of LFP data, sectioned by trials
        spikes
        banpos
        baneat
    end
    
    methods
        % Constructor
        function S = Session(varargsin)
            % Creates a new session object.
            % S = Session() creates an empty object
            %
            % UNDER CONSTRUCTION!!
            % S = Session(navDat,neuralDat) creates as session from the
            % specified data files. Best called from the "Day" class, as
            % this contains the requisite folders. Alignes data and etc.
            %
            % S = Session(datFile) creates a session object from an
            % existing "FieldTrip" type file
            %
            % S = Session(datFile,navDat,neuralDat) creats a session object
            % from an existing "Filed Trip" type file, but with specified
            % navigation and neural data tags.
            
            switch nargin
                % No arguments specified, builds an empty session.
                case 0
                    S.navDat = '';
                    S.neuralDat = '';
                    S.shape = '';
                    S.time = [];
                    S.posdat = [];
                    S.dirdat = [];
                    S.eyedat = [];
                    S.LFPdat = [];
                    S.spikes = [];
                    S.banpos = [];
                    S.baneat = [];
                    % Only one argument specified. Assumes that this argument
                    % is an existing "Field Trip" style format, builds session
                    % accordingly.
                case 1
                    data = varargsin;
                    try
                        load(data)
                    catch
                    end
                    % Assign each of the existing filds in the "FieldTrip" like
                    % structure to the session.
                    S.time = data.time;
                    S.posdat = data.posdat;
                    % Autodetect shape.
                    highest = 0;
                    for trl = 1:length(data.posdat);
                        highest = max(max(abs(S.posdat{trl})));
                    end
                    if highest>10
                        S.shape = 'circle';
                    else
                        S.shape = 'square';
                    end
                    S.dirdat = data.dirdat;
                    S.eyedat = data.eyedat;
                    S.LFPdat = data.neural;
                    S.banpos = data.banpos;
                    S.baneat = data.baneat;
                    S.spikes = data.spikes;
                    % Three arguments specified. Recursivly builds a
                case 3
                    S = Session(varargsin{1});
                    S.navDat = varargsin{2};
                    S.neuralDat = varargsin{3};
                case 2
                    S.navDat  = varargsin{1};
                    S.neuralDat = varargsin{2};
            end
        end
        
        function  [S2] = getcontinuous(S)
            % Returns data as continuouos variables, rather than in the
            % "trial sorted" way mandated by fieldtrip.
            %
            % [tme,pos,dir,eye,lfp,spk] = Session.getcontinous returns
            % specific data vectors
            %
            % S = Session.getcontinous returns data as a new Session object
            % with one trial. This is the default if no output is
            % specified.
            
            tme = [];
            pos = [];
            dir = [];
            eye = [];
            lfp = [];
            spk = [];
            S2 = S;
            if iscell(S.time)
                for trl = 1:length(S.time)
                    tme = [tme S.time{trl}];
                    pos = [pos S.posdat{trl}];
                    dir = [dir S.dirdat{trl}];
                    eye = [eye S.eyedat{trl}];
                    lfp = [lfp S.LFPdat{trl}];
                    spk = [spk S.spikes{trl}];
                end
                S2.time = tme;
                S2.posdat = pos;
                S2.dirdat = dir;
                S2.eyedat = eye;
                S2.LFPdat = lfp;
                S2.spikes = spk;
            end
        end
        
        function R = makeratemap(S,type,channel,unit,size,rez)
            % Creates a Gaussian Kernal rate map for the specified session.
            %
            % Session.makeratemap('pos',size) creates a rate map for the VR
            % position data w/ the specified filter STD "size" and
            % resolution "rez"
            %
            % Session.makeratemap('eye',size) creates a rate map for the
            % eye data with the specified filter STD "size" and resolution
            % "rez".
            switch type
                case 'pos'
                    switch S.shape
                        case 'circle'
                            range = [-15 15 -15 15];
                        case 'square'
                            range = [-10 10 -10 10];
                    end
                    postrain = S.getcontinuous.posdat;
                case 'eye'
                    range = [-600 600 -500 500];
                    postrain = S.getcontinuous.eyedat;
            end
            spiketrain = S.getcontinuous.spikes(channel,:) == unit;
            filt_sigma = size;
            
            % Set up space
            x = linspace(range(1),range(2),rez);
            y = linspace(range(3),range(4),rez);
            [X,Y] = meshgrid(x,y);
            posX = postrain(1,:);
            posY = postrain(2,:);
            
            % Define Gaussian Kernal
            mu = [0,0];
            sigma = filt_sigma*eye(2);
            gk = mvnpdf([X(:) Y(:)],mu,sigma);
            gk = reshape(gk,length(x),length(y));
            
            % Calculate histograms
            H_spike = hist3([posX(spiketrain)',posY(spiketrain)'],{x,y});
            H_space = hist3([posX(:), posY(:)],{x,y});
            
            % Smooth histograms
            F_spike = conv2(H_spike,gk,'same');
            F_space = conv2(H_space,gk,'same');
            
            %Generate "rate map"
            R = F_spike./F_space;
            
            % Plot
            pcolor(X,Y,R);
            colorbar
            shading interp
            axis('square');
        end
        
        function mergeunits(S,channel,spks)
            % combine units sorted by klusta into one unit w/in this
            % session.
            % 
            % Session.mergeunits(channel, [spks]) modifies "Session" while
            % s.t. units in [spks] vector are combined on channel "channel"
            newID = spks(1);
            for s = 2:length(spks);
                S.spikes(channel,S.spikes(channel,:)==spks(s)) = newID;
            end
        end
        
    end
end


