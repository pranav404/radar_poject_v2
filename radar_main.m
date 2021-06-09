clc;
clearvars;
close all;

fc     = 2e9;                 % Radar carrier frequency (Hz)
c      = 3e8;                 % Propagation speed (m/s)
lambda = c/fc;                % Radar wavelength (m)

maxrng = 70e3;               % Maximum range (m)
minrng = 2e3;                 % Minimum range (m)
bw     = 1e6;
fs     = 3*bw;
prf    = 1/range2time(maxrng,c);
dcycle = 0.1;

wav = phased.LinearFMWaveform('SampleRate', fs, ...
    'DurationSpecification', 'Duty cycle', 'DutyCycle', dcycle, ...
    'PRF', prf, 'SweepBandwidth', bw);
rngres = bw2range(bw,c)



arraysz   = 50;                                                  
ant       = phased.URA('Size',arraysz,'ElementSpacing',lambda/2);   
ant.Element.BackBaffled = true;

arraystv  = phased.SteeringVector('SensorArray',ant,'PropagationSpeed',c);
radiator  = phased.Radiator('OperatingFrequency',fc, ...
    'PropagationSpeed', c, 'Sensor',ant, 'WeightsInputPort', true);
collector = phased.Collector('OperatingFrequency',fc, ...
    'PropagationSpeed', c, 'Sensor',ant);

beamw = rad2deg(lambda/(arraysz*lambda/2))


pd      = 0.9;                     % Probability of detection
pfa     = 1e-6;                    % Probability of false alarm
snr_min = albersheim(pd, pfa, 1);
ampgain = 20;
tgtrcs  = 1;
ant_snrgain = pow2db(arraysz^2);

ppower  = radareqpow(lambda,maxrng,snr_min,wav.PulseWidth,...
    'RCS',tgtrcs,'Gain',ampgain+ant_snrgain);

tx = phased.Transmitter('PeakPower',ppower,'Gain',ampgain,'InUseOutputPort',true);
rx = phased.ReceiverPreamp('Gain',ampgain,'NoiseFigure',5,'EnableInputPort',true);


mfcoeff = getMatchedFilter(wav);
mf      = phased.MatchedFilter('Coefficients',mfcoeff,'GainOutputPort', true);

% time varying gain
tgrid   = unigrid(0,1/fs,1/prf,'[)');
rgates  = c*tgrid/2;
rngloss = 2*fspl(rgates,lambda);
refloss = 2*fspl(maxrng,lambda);
tvg     = phased.TimeVaryingGain('RangeLoss',rngloss,'ReferenceLoss',refloss);

% monopulse
monfeed = phased.MonopulseFeed('SensorArray',ant,'PropagationSpeed',c,...
    'OperatingFrequency',fc,'SquintAngle',3);
monest  = getMonopulseEstimator(monfeed);


tracker = trackerGNN('FilterInitializationFcn',@initMPARGNN,...
    'ConfirmationThreshold',[2 3], 'DeletionThreshold',3,...
    'HasDetectableTrackIDsInput',true,'AssignmentThreshold',100,...
    'MaxNumTracks',5,'MaxNumSensors',1);


mfradar.Tx      = tx;
mfradar.Rx      = rx;
mfradar.TxAnt   = radiator;
mfradar.RxAnt   = collector;
mfradar.Wav     = wav;
mfradar.RxFeed  = monfeed;
mfradar.MF      = mf;
mfradar.TVG     = tvg;
mfradar.DOA     = monest;
mfradar.STV     = arraystv;
mfradar.Tracker = tracker;
mfradar.IsTrackerInitialized = false;
t = 0:0.1:20;

tgtpos = [39875;4225; 1000]; 
yveldat = 344;
tgtvel = [200; 0.5*yveldat; 100]; 
X1 = tgtpos(1,1) + tgtvel(1,1)*t + (5/2)*t.^2;
%X2 = tgtpos(1,2) + tgtvel(1,2)*t + (5/2)*t.^2;
%X3 = tgtpos(1,3) + tgtvel(1,3)*t + (5/2)*t.^2;
Y1 = tgtpos(2,1) + tgtvel(2,1)*t + (5/2)*t.^2;
%Y2 = tgtpos(2,2) + tgtvel(2,2)*t + (5/2)*t.^2;
%Y3 = tgtpos(2,3) + tgtvel(2,3)*t + (5/2)*t.^2;
Z1 = tgtpos(3,1) + tgtvel(3,1)*t + (-9.8/2)*t.^2;
%Z2 = 0*t;
%Z3 = 0*t;

wpts = [t.' X1.' Y1.' Z1.'];
ntgt = size(tgtpos,2);                                      
tgtmotion = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts);   
target = phased.RadarTarget('MeanRCS',tgtrcs*ones(1,ntgt),'OperatingFrequency',fc);    


channel = phased.FreeSpace('SampleRate',fs,'TwoWayPropagation',true,'OperatingFrequency',fc);


env.Target       = target;
env.TargetMotion = tgtmotion;
env.Channel      = channel;



scanregion   = [-30, 30, 0, 20]; 
azscanspan   = diff(scanregion(1:2));
numazscan    = ceil(azscanspan/beamw);
azscanangles = linspace(scanregion(1),scanregion(2),numazscan);
elscanspan   = diff(scanregion(3:4));
numelscan    = ceil(elscanspan/beamw);
elscanangles = linspace(scanregion(3),scanregion(4),numelscan);
[elscangrid,azscangrid] = meshgrid(elscanangles,azscanangles);  
scanangles   = [azscangrid(:) elscangrid(:)].';


sceneplot = helperMPARTaskPlot('initialize',scanangles,azscanangles,maxrng,beamw,tgtpos,X1,Y1,Z1);


searchq = struct('JobType','Search','BeamDirection',num2cell(scanangles,1),...
    'Priority',1000,'WaveformIndex',1);
current_search_idx = 1;


disp(searchq(current_search_idx))





trackq(10) = struct('JobType',[],'BeamDirection',[],'Priority',3000,'WaveformIndex',[],...
    'Time',[],'Range',[],'TrackID',[]);
num_trackq_items = 0;
disp(trackq(1))







jobq.SearchQueue  = searchq;
jobq.SearchIndex  = current_search_idx;
jobq.TrackQueue   = trackq;
jobq.NumTrackJobs = num_trackq_items;





rng(2018);
current_time = 0;
Npulses      = 10;
numdwells    = 200;
dwelltime    = 0.01;

jobload.num_search_job = zeros(1,numdwells);
jobload.num_track_job  = zeros(1,numdwells);

xef =  0;



for dwell_idx = 1:200
   %{
    if xef == 600
        yveldat = -1*yveldat;
        tgtvel =[-400 400; 0.5*yveldat 0.5*yveldat; 0 0]; 
        tgtmotion = phased.Platform('InitialPosition',tgtpos,'Velocity',tgtvel); 
        xef = 0;
    end
    %}
    %xef = xef + 10;
    searchidx = jobq.SearchIndex;
    searchqidx = mod(searchidx-1,numel(jobq.SearchQueue))+1;
    readyidx = find([jobq.TrackQueue(1:jobq.NumTrackJobs).Time]<=current_time);
    [~,maxpidx] = max([jobq.TrackQueue(readyidx).Priority]);
    taskqidx = readyidx(maxpidx);
    if ~isempty(taskqidx) && jobq.TrackQueue(taskqidx).Priority >= jobq.SearchQueue(searchqidx).Priority
        current_job = jobq.TrackQueue(taskqidx);
        for m = taskqidx+1:jobq.NumTrackJobs
            jobq.TrackQueue(m-1) = jobq.TrackQueue(m);
        end
        jobq.NumTrackJobs = jobq.NumTrackJobs - 1;
        jobq.SearchQueue(searchqidx).Priority = jobq.SearchQueue(searchqidx).Priority+100;
    else
        current_job = jobq.SearchQueue(searchqidx);
        searchidx = searchqidx+1;
    end
    jobq.SearchIndex = searchidx;
    % Simulate the received I/Q signal
    radarpos = [0;0;0];
    radarvel = [0;0;0];
    Npulse = 10;
    for i = 1:Npulse
        x = wav();
        [tgtposition,tgt_velocity] = tgtmotion(1/prf);
        %tgtposition
        [~,tgtangle] = rangeangle(tgtposition);
        [xtran,inuseflag] = tx(x);
        w = arraystv(fc,current_job.BeamDirection);
        xtran = radiator(xtran,tgtangle,conj(w));
        xpropagated = channel(xtran,radarpos,tgtposition,radarvel,tgt_velocity);
        xpropagated = target(xpropagated);
        xrec = collector(xpropagated,tgtangle);
        [xrs,xrdaz,xrdel] = monfeed(xrec,current_job.BeamDirection);
        if i == 1
            xrsint   = mfradar.Rx(xrs,~(inuseflag>0));
            xrdazint = mfradar.Rx(xrdaz,~(inuseflag>0));
            xrdelint = mfradar.Rx(xrdel,~(inuseflag>0));
        else
            xrsint   = xrsint+mfradar.Rx(xrs,~(inuseflag>0));
            xrdazint = xrdazint+mfradar.Rx(xrdaz,~(inuseflag>0));
            xrdelint = xrdelint+mfradar.Rx(xrdel,~(inuseflag>0));
        end
            
    end
    xsum = xrsint; xdaz = xrdazint; xdel = xrdelint;
    % Signal processor to extract detection
    nbw         = mfradar.Rx.SampleRate/(mfradar.Wav.SampleRate/mfradar.Wav.SweepBandwidth);
    npower      = noisepow(nbw,mfradar.Rx.NoiseFigure,mfradar.Rx.ReferenceTemperature);
    pfa         = 1e-6;
    threshold   = npower * db2pow(npwgnthresh(pfa,1,'noncoherent'));
    arraysz     = mfradar.TxAnt.Sensor.Size(1);
    ant_snrgain = pow2db(arraysz^2);
    mfcoeff     = getMatchedFilter(mfradar.Wav);
    mfgain      = pow2db(norm(mfcoeff)^2);
    threshold   = threshold * db2pow(mfgain+2*ant_snrgain); 
    threshold   = sqrt(threshold);    
    tgrid   = unigrid(0,1/mfradar.Wav.SampleRate,1/mfradar.Wav.PRF,'[)');
    rgates  = mfradar.TxAnt.PropagationSpeed*tgrid/2;
    xrsmf = mfradar.TVG(mfradar.MF(xrsint));
    if any(abs(xrsmf)>threshold)
        [~,tgtidx] = findpeaks(abs(xrsmf),'MinPeakHeight',threshold,...
            'Sortstr','Descend','NPeaks',1);
        rng_est = rgates(tgtidx-(numel(mfcoeff)-1));
        ang_est = mfradar.DOA(xrsint(tgtidx-1),xrdazint(tgtidx-1),xrdelint(tgtidx-1),current_job.BeamDirection);
        % Form the detection object.  
        measNoise = diag([0.1, 0.1, 150].^2);           % Measurement noise matrix
        detection = objectDetection(current_time,...
            [ang_est(1);ang_est(2);rng_est], 'MeasurementNoise', measNoise,...
            'MeasurementParameters',struct('Frame','spherical', 'HasVelocity', false));  
    else
        detection = objectDetection.empty;
    end
    fprintf('\n%f sec:\t%s\t[%f %f]',current_time,current_job.JobType,current_job.BeamDirection(1),...
            current_job.BeamDirection(2));
    
    % Radar manager to perform data processing and update track queue
    trackq           = jobq.TrackQueue;
    num_trackq_items = jobq.NumTrackJobs;

    % Execute current job
    switch current_job.JobType
        case 'Search'
            % For search job, if there is a detection, establish tentative
            % track and schedule a confirmation job
            if ~isempty(detection)
                ang_est = detection.Measurement(1:2);
                rng_est = detection.Measurement(3);
                if ~mfradar.IsTrackerInitialized
                    [~,~,allTracks] = mfradar.Tracker(detection,current_time,uint32([]));
                    mfradar.IsTrackerInitialized = true;
                else
                    [~,~,allTracks] = mfradar.Tracker(detection,current_time,uint32([]));
                end
                num_trackq_items = num_trackq_items+1;
                trackq(num_trackq_items) = struct('JobType','Confirm','Priority',2000,...
                    'BeamDirection',ang_est,'WaveformIndex',1,'Time',current_time+dwelltime,...
                    'Range',rng_est,'TrackID',allTracks(~[allTracks.IsConfirmed]).TrackID);
                if current_time < 0.3 || strcmp(current_job.JobType,'Track')
                    fprintf('\tTarget detected at %f m',rng_est);
                end
            else
                allTracks = [];
            end

        case 'Confirm'
            % For confirm job, if the detection is confirmed, establish a track
            % and create a track job corresponding to the revisit time
            if ~isempty(detection)
                trackid = current_job.TrackID;
                [~,~,allTracks] = mfradar.Tracker(detection,current_time,trackid);

                rng_est = detection.Measurement(3);
                if rng_est >= 50e3
                    updateinterval = 0.5;
                else
                    updateinterval = 0.1;
                end
                revisit_time = current_time+updateinterval;
                predictedTrack = predictTracksToTime(mfradar.Tracker,trackid,revisit_time);
                xpred = predictedTrack.State([1 3 5]);
                [phipred,thetapred,rpred] = cart2sph(xpred(1),xpred(2),xpred(3));
                num_trackq_items = num_trackq_items+1;
                trackq(num_trackq_items) = struct('JobType','Track','Priority',3000,...
                    'BeamDirection',rad2deg([phipred;thetapred]),'WaveformIndex',1,'Time',revisit_time,...
                    'Range',rpred,'TrackID',trackid);
                if current_time < 0.3 || strcmp(current_job.JobType,'Track')
                    fprintf('\tCreated track %d at %f m',trackid,rng_est);
                end
            else
                allTracks = [];
            end

        case 'Track'
            % For track job, if there is a detection, update the track and
            % schedule a track job corresponding to the revisit time. If there
            % is no detection, predict and schedule a track job sooner so the
            % target is not lost.
            if ~isempty(detection)
                trackid = current_job.TrackID;
                [~,~,allTracks] = mfradar.Tracker(detection,current_time,trackid);

                rng_est = detection.Measurement(3);
                if rng_est >= 50e3
                    updateinterval = 0.5;
                else
                updateinterval = 0.1;
                end

                revisit_time = current_time+updateinterval;
                predictedTrack = predictTracksToTime(mfradar.Tracker,trackid,revisit_time);
                xpred = predictedTrack.State([1 3 5]);
                [phipred,thetapred,rpred] = cart2sph(xpred(1),xpred(2),xpred(3));
                num_trackq_items = num_trackq_items+1;
                trackq(num_trackq_items) = struct('JobType','Track','Priority',3000,...
                    'BeamDirection',rad2deg([phipred;thetapred]),'WaveformIndex',1,'Time',revisit_time,...
                    'Range',rpred,'TrackID',trackid);

                if current_time < 0.3 || strcmp(current_job.JobType,'Track')
                    fprintf('\tTrack %d at %f m',trackid,rng_est);
                end
            else
                trackid = current_job.TrackID;
                [~,~,allTracks] = mfradar.Tracker(detection,current_time,trackid);

                updateinterval = 0.1;  % revisit sooner
                revisit_time = current_time+updateinterval;
                predictedTrack = predictTracksToTime(mfradar.Tracker,trackid,revisit_time);
                if(~isempty(predictedTrack))
                    xpred = predictedTrack.State([1 3 5]);
                    [phipred,thetapred,rpred] = cart2sph(xpred(1),xpred(2),xpred(3));
                    num_trackq_items = num_trackq_items+1;
                    trackq(num_trackq_items) = struct('JobType','Track','Priority',3000,...
                    'BeamDirection',rad2deg([phipred;thetapred]),'WaveformIndex',1,'Time',revisit_time,...
                    'Range',rpred,'TrackID',trackid);
                    current_job.JobType = 'Search';
                    fprintf('\tNo detection, track %d predicted, prediction details theta = %f,phi = %f',...
                        current_job.TrackID,rad2deg(phipred),rad2deg(thetapred));
                else
                    num_trackq_items = num_trackq_items+1;
                    trackq(num_trackq_items) = struct('JobType',[],'Priority',[],...
                    'BeamDirection',[],'WaveformIndex',[],'Time',[],...
                    'Range',[],'TrackID',[]);
                    current_job.JobType = 'Search';
                    fprintf('TrackID %d has been dropped. Continuing Search.',current_job.TrackID);
                end
            end

    end
    
    jobq.TrackQueue   = trackq;
    jobq.NumTrackJobs = num_trackq_items;
    % Update time
    tgtpos = tgtmotion(dwelltime-Npulses/mfradar.Wav.PRF);
    helperMPARTaskPlot('update',sceneplot,current_job,maxrng,beamw,tgtpos,allTracks,detection.Measurement);
    current_time = current_time+dwelltime;
end