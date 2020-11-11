function [rwave , mwitholdATEND, mvimaxvalATEND ] = ecglaux(ecg, mwisignal, sampling, ... 
		mwithold, mvimaxval, mwiwidthpts, refractpts, mwitholdfract, mwitholdff )

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
Rpickind = Rpickval ; Rpickind(1)=-9999*refractpts;
prevslopeup = 0 ;  

limitofwhile = numel(mwisignal)-examwindow-1.25*mwiwidthpts-1  - lookmorepts   -  bufindA ; 
mwiwidthpts1pt25 = 1.25*mwiwidthpts ; 
sl10 = 10*sl; sl20 = 20*sl ; 

perpt = sub1+1+bufindA ; Rcount = 0 ; 

while perpt <  limitofwhile
	taken='N';
	[val,ind] = max(mwisignal(perpt:perpt+examwindow)) ; absind = ind+perpt-1;
    % i will catch the upstroke of the MWI, at the extreme right of this window
    if val>mwithold  &  ( mwisignal(round(absind+sl10))>val  |  mwisignal(round(absind+sl20))>val ) 
       % i caught the upstroke in the right window, 
       % find the max of this bump
       [val2,ind2] = max(mwisignal(absind : round(absind+mwiwidthpts1pt25)  ));
        Aind2 = ind2+absind-1;
       % this is the max of the MWI bump	% do further testing, or just take the max ECG's vals in a fixed window from here?
        if val2<3.0*mvimaxval
			pt1 = Aind2-sub1   ; 	
            pt2 = Aind2+lookmorepts ;		
			[ RpickvalUP, RpickindUP] = max(ecg(pt1:pt2)); 
            RpickindUP = RpickindUP + (pt1) - 1 ;	
			nowslopeup  =  abs(  ecg(RpickindUP-0)-ecg(RpickindUP-bufindC)  ) / (divbufindC)  ;
            % need to decide whether to keep slope crit or not. default has been yes in the past.
            if ecg(RpickindUP-bufindA)<ecg(RpickindUP-0) & ecg(RpickindUP+bufindA)<ecg(RpickindUP+0)   ...  %  & ecg(RpickindUP-bufindB)<ecg(RpickindUP-0)  & ecg(RpickindUP+bufindB)<ecg(RpickindUP+0)
                        &  (RpickindUP > Rpickind(max([Rcount 1]))+refractpts)   ... 
                        &   (    nowslopeup   >   0.33*prevslopeup     ) % all but last new 4/14/00
                taken='Y'; 	
            end					
		end 

    end 
	
    if taken=='Y'
		Rcount = Rcount+1;  Rpickval(Rcount)=RpickvalUP;  Rpickind(Rcount)=RpickindUP;
		perpt = Rpickind(Rcount) + refractpts ; 
		mwithold = mwitholdfract*val2 + mwitholdff*(mwithold-mwitholdfract*val2) ; 
		mvimaxval= val2 + mwitholdff*(mvimaxval-val2) ;
		prevslopeup = nowslopeup +  mwitholdff*(prevslopeup-nowslopeup) ;
    else
		perpt = perpt + ifno ;
    end  %if val

end
rwave = Rpickind(1:Rcount);   
mwitholdATEND = mwithold;  
mvimaxvalATEND = mvimaxval; 

return ;
