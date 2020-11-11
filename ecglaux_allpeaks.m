function [rwave , rwave_oppdirUD , mwitholdATEND, mvimaxvalATEND ] = ecglaux_allpeaks(ecg, mwisignal, sampling, ... 
		mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff )

% clear Mfactor ; 
mwitholdForOppDir = mwithold * .75 ; 

sl = sampling/1000 ;
% now i will examine the mwisignal for threshold detection
examwindow = round(200*sl) ; % should be <50% of a typical RRint, so as not to capture 2 R peaks in it
sub1 = round(275*sl) ; 	% subtract val 1, "225" in '86 paper
sub2 = round(150*sl) ; 	% subtract val 2, "125" in '86 paper
lookmorepts = round(0*sl) ;
ifno = round(25*sl); 		% 25 msec ahead if not picked

bufindA = round(6*sampling/120);
bufindB = round(2*sampling/120);
bufindC = round(4*sampling/120);
divbufindC = bufindC+1; 

% pre-make with zeros :
Rpickval = zeros(1, round(100*numel(ecg)/sampling/60)   );  % overestimate for hr=100
Rpickind = Rpickval ; 
Rpickind(1)=-9999*refractpts;
prevslopeup = 0 ;  

RpickvalTotDown = zeros(1, round(100*numel(ecg)/sampling/60/5)   );  % overestimate for hr=100, div5
RpickindTotDown = RpickvalTotDown ; 
RpickindTotDown(1)=-9999*refractpts;

prevRpickindEitherUPorDOWN=-9999*refractpts;


limitofwhile = numel(mwisignal)-examwindow-1.25*mwiwidthpts-1  - lookmorepts   -  bufindA ; 
mwiwidthpts1pt25 = 1.25*mwiwidthpts ; 
sl10 = 10*sl; sl20 = 20*sl ; 

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
perpt = sub1+1+bufindA ; Rcount = 0 ; RcountDOWN=0;

while perpt <  limitofwhile   % examwindow = min([ examwindow  l-perpt-1]);
	taken='N'; takenDOWN='N';
	[val,ind] = max(mwisignal(perpt:perpt+examwindow)) ; absind = ind+perpt-1;
	% i will catch the upstroke of the MWI, at the extreme right of this window
    if val>mwithold       &  ( mwisignal(round(absind+sl10))>val  |  mwisignal(round(absind+sl20))>val ) % i caught the upstroke in the right window, 
		% find the max of this bump
        [val2,ind2] = max(mwisignal(absind : round(absind+mwiwidthpts1pt25)  )); Aind2 = ind2+absind-1;
		% this is the max of the MWI bump	% do further testing, or just take the max ECG's vals in a fixed window from here?
		% usual : exclude VERY big peaks, are prob PVCs : if val2<3.0*mvimaxval ,
        % for all-peaks : 
        if val2<10000*3.0*mvimaxval
			% positive = take
			pt1 = Aind2-sub1   ; 	
            pt2 = Aind2+lookmorepts ;		           
			[RpickvalUP, RpickindUP] = max(ecg(pt1:pt2)); 
            RpickindUP = RpickindUP + (pt1) - 1 ;	
			nowslopeup  =  abs(  ecg(RpickindUP-0)-ecg(RpickindUP-bufindC)  ) / (divbufindC)  ;
            % need to decide whether to keep slope crit or not. default has been yes in the past.

            if ecg(RpickindUP-bufindA)<ecg(RpickindUP-0) & ecg(RpickindUP+bufindA)<ecg(RpickindUP+0)   ... 
               &  (RpickindUP > prevRpickindEitherUPorDOWN+refractpts) &   (    nowslopeup   >   0.33*prevslopeup     ) % all but last new 4/14/00
            % positive = take
            taken='Y'; 	
            end
            
            if taken=='N',
		    % new for all_peaks - if no "right" peak was chosen on this interval, try a "wrong" peak - on (-1)*ecg
            % retain the same pt1 and pt2
			[RpickvalDOWN, RpickindDOWN] = max(  (-1)*ecg(pt1:pt2)   ); 
            RpickindDOWN= RpickindDOWN + (pt1) - 1 ;	
			nowslopeDOWN  =  abs(  (-1)*ecg(RpickindDOWN-0)  -  (-1)*ecg(RpickindDOWN-bufindC)  ) / (divbufindC)  ;
            % need to decide whether to keep slope crit or not. default has been yes in the past.
                if (-1)*ecg(RpickindDOWN-bufindA)  <  (-1)*ecg(RpickindDOWN-0) ...
                 & (-1)*ecg(RpickindDOWN+bufindA)  <  (-1)*ecg(RpickindDOWN+0)   ...  %  & ecg(RpickindUP-bufindB)<ecg(RpickindUP-0)  & ecg(RpickindUP+bufindB)<ecg(RpickindUP+0)
                 &         (RpickindDOWN > prevRpickindEitherUPorDOWN+refractpts) 
                % toss last one : 	&   (    nowslopeupDOWN   >   0.33*prevslope??up??     ), % all but last new 4/14/00
                % positive = take
                takenDOWN='Y'; 	
                end
            end % if taken=='N', 

         end % if val2<...

     end % if val>...

	if taken=='Y'
		Rcount = Rcount + 1 ;  
        Rpickval(Rcount)          = RpickvalUP ;  
        Rpickind(Rcount)          = RpickindUP ;
        prevRpickindEitherUPorDOWN    = RpickindUP ;
		perpt = Rpickind(Rcount) + refractpts ; 
        mwithold = mwitholdfract*val2 + mwitholdff*(mwithold-mwitholdfract*val2) ; 
		mvimaxval= val2 + mwitholdff*(mvimaxval-val2) ;
		prevslopeup = nowslopeup +  mwitholdff*(prevslopeup-nowslopeup) ;   
    elseif takenDOWN=='Y'         
		RcountDOWN = RcountDOWN + 1 ;  
        RpickvalTotDown(RcountDOWN)  = RpickvalDOWN ;  
        RpickindTotDown(RcountDOWN)  = RpickindDOWN ;
        prevRpickindEitherUPorDOWN       = RpickindDOWN ;
		perpt = RpickindTotDown(RcountDOWN) + refractpts ;                 
        % do _not_ update trigger params based on a "wrong" detection                                
    else
		perpt = perpt + ifno ;
	end  %if val

end

rwave = Rpickind(1:Rcount);   mwitholdATEND = mwithold;  mvimaxvalATEND = mvimaxval; 

rwave_oppdirUD=[];
if RcountDOWN>0
    rwave_oppdirUD = RpickindTotDown(1:RcountDOWN) ; 
end 

% ^ this can, will be often, null
% bring this idea back from the dead - now they are specifically "rwave_oppdirUD" :
% % % % % % do one more check thru the set, infer if any R peaks were missed
moreOppDirR=[];  
for perrwave=2:numel(rwave)-1 
    % don't do if there already is an OppDir R in this cyc, obv.
    hit_i = [];
    if numel(rwave_oppdirUD)>0
        hit_i = find( rwave_oppdirUD>=rwave(perrwave) & rwave_oppdirUD<=rwave(perrwave+1)-1 ) ; 
    end 
    if numel(hit_i)>0
        continue
    end
        
    SubtractN = min([ 9   perrwave-1  ]);
    rrmedian = median(1+diff(rwave( perrwave-SubtractN : perrwave  ))); 
    rrdist=numel( rwave(perrwave) : rwave(perrwave+1)-1 );
	if rrdist>1.80*rrmedian
        %%% DOWN =  MIN of the up-"ecg" = NEGATIVE VALUE
        [val,ind] = min(ecg( rwave(perrwave)+round(.20*rrmedian)  : rwave(perrwave+1)-1-round(.20*rrmedian)  )); 
		ind  = ind +    rwave(perrwave)+round(.20*rrmedian)   - 1 ;
	    m = mean([ ecg(rwave(perrwave))  ecg(rwave(perrwave+1))  ]); 	    
        if abs(val)>.35*m  
            moreOppDirR = [  moreOppDirR   ind  ] ; 
        end % take it
    end % if rrdist >

end % for

if numel(moreOppDirR)>0
rwave_oppdirUD = unique([  rwave_oppdirUD   moreOppDirR   ]);  % 'unique' does 'sort' ; 'sort' doesn't do 'unique'
end


