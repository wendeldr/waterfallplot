function rwave = ecgLms(ecg,sampling,b_butter_ecg4mwi,a_butter_ecg4mwi,dbg) 

mwitholdfract = .25;  % normally 0.30 , low # = very tolerant, easy detection
mwitholdff    = .80;   % usually .80

ecg = ecg - mean(ecg); % neha - baseline wander removal?

origecg=ecg; sl = sampling/1000 ;
mwiwidthmsec = 175  ; mwiwidthpts = round( mwiwidthmsec * sl ) ;  % usually 175/...
refractmsec  = 250 ;   refractpts = round(  refractmsec * sl );   % usually 250/...
% refract should be greater than mwiwidth

% filter, differentiate and square ecg
filtecg = filtfilt(b_butter_ecg4mwi,a_butter_ecg4mwi,origecg); % zero phase
difffiltecg = diff2(filtecg) ; % '84 paper, 5pt diff-ation% perhaps i will use this signal later
sqdifffiltecg = difffiltecg.^2 ; 
ll = numel(sqdifffiltecg);

% create Moving-Window-Integration							   
mwisignal=zeros(1,ll);
for cnt=1:mwiwidthpts
    mwisignal(cnt) = sum(sqdifffiltecg(1:cnt)); 
end

l = mwiwidthpts ;
wholesum = cumsum(sqdifffiltecg); 
if max(wholesum) > realmax/100
    % just crash it for now; if this ever actually happened, then make a fix/work-around
    disp('vvvvvvvvvvvvvvvvv');
    length_wholesum = length(wholesum);
    max_val_wholesum = max(wholesum);
    error('in ecgLms.m, variable "wholesum" (cumsum) comes close to exceeding max allowed Matlab value - quitting - create a work-around now . . .');
end

mwisignal(l+1:ll) = wholesum(l+1:ll)  -  wholesum(1:ll-l);  

% the begin of mwisignal always starts near zero, this is sometimes bad
mwisignal(1:mwiwidthpts) = mwisignal(mwiwidthpts+1)*ones(1,mwiwidthpts) ;                      
% "mwisignal" is what is used in the detection algorithm 


%######################################
% initiate the algorithm with the firts 4 seconds, 10 seconds

pt1=1; pt2 = min([ 6*sampling   numel(mwisignal)  ]);  pt3 = min([ 8*sampling   numel(mwisignal)  ]);

x = sort( mwisignal(pt1:pt2)  );  lx = numel(x);  % they go from small to large values  
mvimaxval = mean(x( round(.90*lx):round(.95*lx) ))  ;
mwithold  = mwitholdfract * mvimaxval  ;

rpos = ecglaux(    origecg(pt1:pt3) , mwisignal(pt1:pt3), sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff);
rneg = ecglaux( -1*origecg(pt1:pt3) , mwisignal(pt1:pt3), sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff) ;


if numel(rpos)<2
    rpos = [ pt1 pt3 ]; 
end

if numel(rneg)<2
    rneg = [ pt1 pt3 ]; 
end

if numel(rpos)<3
    ecg = -1*origecg;
elseif numel(rneg)<3   
    ecg = origecg;
elseif  abs(median(origecg(rneg)))  >  abs(median(origecg(rpos)))  
    ecg = -1*origecg;
end

rwave = ecglaux(ecg,mwisignal,sampling,mwithold,mvimaxval,mwiwidthpts,refractpts,mwitholdfract,mwitholdff);

end

