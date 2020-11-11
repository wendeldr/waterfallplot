function rwave = ecgLmsDog1(ecg , sampling, b_butter_ecg4mwi , a_butter_ecg4mwi , MoreSensitive4GP ); 


mwitholdfract = .25;  % normally 0.30 , low # = very tolerant, easy detection
mwitholdff    = .80;   % usually .80


if exist('MoreSensitive4GP','var')==1,
if numel( MoreSensitive4GP )==1,
if        MoreSensitive4GP == 'y' , 
mwitholdfract = .125 ; % special for Guinea Pig - more sensitive to picka peak as an R, less of a Twave to worry about
mwitholdff    = .65 ; % less history going into each new value of withold
end ; end ; end ; 


ecg = ecg - mean(ecg); 


origecg=ecg; sl = sampling/1000 ;
mwiwidthmsec = 115  ; mwiwidthpts = round( mwiwidthmsec * sl ) ;  % usually HUMAN = 175
refractmsec  = 135  ;  refractpts = round(  refractmsec * sl );   % usually HUMAN = 250
% refract should be greater than mwiwidth



filtecg = filtfilt( b_butter_ecg4mwi , a_butter_ecg4mwi , origecg ); % zero phase

difffiltecg = diff2(filtecg) ; % '84 paper, 5pt diff-ation% perhaps i will use this signal later
sqdifffiltecg = difffiltecg.^2 ; ll = numel(sqdifffiltecg) ;

% create Moving-Window-Integration

							  % why doesn't this work, to pre-make : 
mwisignal=zeros(1,ll);
for cnt=1:mwiwidthpts, mwisignal(cnt) = sum(  sqdifffiltecg(  1: cnt  )   ) ; end;
l = mwiwidthpts ; %%%%%%%%%%%%% duh : prev bug : l = numel(mwisignal);
wholesum = cumsum(sqdifffiltecg); 
% 3/29/14 : 
% no need to clear (see other for note) : clear realmax ;
if max(wholesum) > realmax/100,
    % just crash it for now; if this ever actually happened, then make a fix/work-around
    disp('vvvvvvvvvvvvvvvvv');
    length_wholesum = length(wholesum) ,
    max_val_wholesum = max(wholesum) ,
    error('in ecgLmsDog1.m, variable "wholesum" (cumsum) comes close to exceeding max allowed Matlab value - quitting - create a work-around now . . .');
end ;

					   mwisignal(l+1:ll) = wholesum(l+1:ll)  -  wholesum(1:ll-l);  

% the begin of mwisignal always starts near zero, this is sometimes bad
                       mwisignal(1:mwiwidthpts) = mwisignal(mwiwidthpts+1)*ones(1,mwiwidthpts) ; 
                       
                       
% "mwisignal" is what is used in the detection algorithm 

% % % % % % % % % % % % backto=gcf; figure; newfig=gcf;
% % % % % % % % % % % % subplot(411); hold on ; plot(ecg,'b'); axis('tight');
% % % % % % % % % % % % subplot(412); hold on ; plot(filtecg); axis('tight');
% % % % % % % % % % % % subplot(413); hold on ; plot(difffiltecg); axis('tight');
% % % % % % % % % % % % subplot(414); hold on ; plot(mwisignal); axis('tight');
% % % % % % % % % % % % figure(backto);

% no need to clear, save time : clear wholesum difffiltecg



%######################################
% initiate the algorithm with the firts 4 seconds, 10 seconds

pt1=1; pt2 = min([ 6*sampling   numel(mwisignal)  ]);  pt3 = min([ 8*sampling   numel(mwisignal)  ]);

x = sort( mwisignal(pt1:pt2)  );  lx = numel(x);  % they go from small to large values  
mvimaxval = mean(x( round(.90*lx):round(.95*lx) ))  ;
mwithold  = mwitholdfract * mvimaxval  ;

rpos = ecglauxDog1(    origecg(pt1:pt3) , mwisignal(pt1:pt3), sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff);
rneg = ecglauxDog1( -1*origecg(pt1:pt3) , mwisignal(pt1:pt3), sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff) ;

%%%%%%%%%%%%%%%%%%%%%%%% save temp.mat rpos rneg ; error('crash');


if numel(rpos)<2, rpos=[ pt1 pt3 ]; end;
if numel(rneg)<2, rneg=[ pt1 pt3 ]; end;

% a PVC peak will still be taller than a normal peak
% my method was that 1 PVC would be picked, plus many mini peaks near 0 line, so median val will be low

% what's my assumption, that there are fewer PVC's, or that
%drpos = diff(rpos);  toss_i = find(drpos>1800*sampling/1000); rpos(toss_i)=[]; if numel(rpos)<2, rpos=[ pt1 pt3 ]; end;
%drneg = diff(rneg);  toss_i = find(drneg>1800*sampling/1000); rneg(toss_i)=[]; if numel(rneg)<2, rneg=[ pt1 pt3 ]; end;

%[ abs(median(origecg(rpos)))    abs(median(origecg(rneg)))        mean(origecg(rpos))   mean(origecg(rneg))  ]
%[ origecg(rpos) ]
%[ origecg(rneg)  ]

%yn = [ abs(median(origecg(rneg)))  >  abs(median(origecg(rpos)))  ]
%yn = [ abs(median(origecg(rpos)))  <  abs(median(origecg(rneg)))  ]
%dadiff = abs(median(origecg(rpos))) - abs(median(origecg(rneg)))

% this should be better, for something so important :

    if numel(rpos)<3,   ecg = -1*origecg;
elseif numel(rneg)<3,   ecg = origecg;
elseif  abs(median(origecg(rneg)))  >  abs(median(origecg(rpos)))  ,   ecg = -1*origecg;
end;



% can force up or down here: 
% ecg = origecg ; disp('Force R up, ecgLms.m')
% ecg = -1*origecg ; disp('Force R down, ecgLms.m')

rwave = ecglauxDog1(  ecg , mwisignal, sampling, mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff);


% % % % % % figure(newfig);
% % % % % % subplot(411); hold on ; plot( rwave ,       ecg(rwave) , 'r+' ) ; 
% % % % % % subplot(414); hold on ; plot( rwave , mwisignal(rwave) , 'r+' ) ; 
% % % % % % figure(backto);

