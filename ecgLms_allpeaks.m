function [rwave, rwave_oppdirUD] = ecgLms_allpeaks(ecg , sampling, b_butter_ecg4mwi , a_butter_ecg4mwi )

mwitholdfract = .25;  % normally 0.30 , low # = very tolerant, easy detection
mwitholdff    = .80;  % usually .80

ecg = ecg - mean(ecg); % clear Mfactor ;
origecg=ecg; sl = sampling/1000 ;
mwiwidthmsec = 175  ; mwiwidthpts = round( mwiwidthmsec * sl ) ;  % usually 175/...
refractmsec  = 250 ;   refractpts = round(  refractmsec * sl );   % usually 250/...
% refract should be greater than mwiwidth

filtecg = filtfilt( b_butter_ecg4mwi , a_butter_ecg4mwi , origecg ); % zero phase
difffiltecg = diff2(filtecg) ; % '84 paper, 5pt diff-ation% perhaps i will use this signal later
sqdifffiltecg = difffiltecg.^2 ; 
ll = numel(sqdifffiltecg) ;

% create Moving-Window-Integration
% why doesn't this work, to pre-make : 
mwisignal=zeros(1,ll);
for cnt=1:mwiwidthpts
    mwisignal(cnt) = sum(sqdifffiltecg(1:cnt)); 
end
l = mwiwidthpts ; 
wholesum = cumsum(sqdifffiltecg); 
mwisignal(l+1:ll) = wholesum(l+1:ll)  -  wholesum(1:ll-l);  

% the begin of mwisignal always starts near zero, this is sometimes bad
mwisignal(1:mwiwidthpts) = mwisignal(mwiwidthpts+1)*ones(1,mwiwidthpts) ; 
                       
                       
% "mwisignal" is what is used in the detection algorithm 

clear wholesum difffiltecg

%######################################
% initiate the algorithm with the firts 4 seconds, 10 seconds

pt1=1; 
pt2 = min([ 6*sampling   numel(mwisignal)  ]);  
pt3 = min([ 8*sampling   numel(mwisignal)  ]);

x = sort( mwisignal(pt1:pt2)  );  
lx = numel(x);  % they go from small to large values  
mvimaxval = mean(x( round(.90*lx):round(.95*lx) ))  ;
mwithold  = mwitholdfract * mvimaxval  ;

% do this tester with the usual, still need to learn if the "right" way is UP or DOWN
rpos = ecglaux(    origecg(pt1:pt3) , mwisignal(pt1:pt3), sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff) ;
rneg = ecglaux( -1*origecg(pt1:pt3) , mwisignal(pt1:pt3), sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff) ;

if numel(rpos)<2
    rpos=[ pt1 pt3 ]; 
end
if numel(rneg)<2
    rneg=[ pt1 pt3 ]; 
end

if numel(rpos)<3  
    ecg = -1*origecg;
elseif numel(rneg)<3  
    ecg = origecg;
elseif  abs(median(origecg(rneg)))  >  abs(median(origecg(rpos)))   
    ecg = -1*origecg;
end

[ rwave , rwave_oppdirUD , mwitholdATEND, mvimaxvalATEND ]  = ...
          ecglaux_allpeaks( ecg    , mwisignal, sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff);

% "rwave" will contain both the good "rwave"s
% _and_ rwave_oppdirUD 
% because in the end, I just have 1 "rwave" vector, with an accompanying listcycok vector,
% but I cannot make the listcycok vector right here and now

MinAllowMsec = 50 ;
MinAllowPts  = round( MinAllowMsec * sampling / 1000 ) ;

toss_i=[];
if numel(rwave)>0
    if numel(rwave_oppdirUD)>0   
        for perOppDir = 1:numel(rwave_oppdirUD)
            d_npts = min(   abs(rwave_oppdirUD(perOppDir)-rwave)  ) ;        
            if d_npts<=MinAllowPts
                toss_i = [ toss_i perOppDir ];
            end

        end % for perOppDir = 1:numel(rwave_oppdirUD),

        if numel(toss_i)==numel(rwave_oppdirUD)
            rwave_oppdirUD=[]; 
        elseif numel(toss_i)>0     
            rwave_oppdirUD(toss_i)=[]; 
        end

    end % if numel(rwave_oppdirUD)>0,
end % if numel(rwave)>0,

rwave = unique([ round(rwave)   round(rwave_oppdirUD)  ]) ; % unique does sort

