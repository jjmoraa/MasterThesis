function [amids,apmids,cpmids,ctmids,cmmids,ocp,oct,ocm,dT,dQ,totalpwr,moment,rsize,rmids]=AD_performance(chordVector,relwnds,aind,apind,aoas,Fs,cls,cds,cpvector,wndspeed,ltsr,R,span)

%ask this to lackner
aind(end)=0;
%aind(end)=aind(end-1);
apind(end)=0;
%apind(end)=apind(end-1);
%Fs(end)=1;
rotspeed=ltsr.*wndspeed/R;
%So, due to the "flawed" evaluation of the final section of the blade, I'm
%going to assume temporaly no tip losses on this forced condiiton
for i=1:(length(chordVector)-1) %If we get the xfoil working (length(blade.ispan)-1)
    amids(i)=(aind(i)+aind(i+1))/2;
    apmids(i)=(apind(i)+apind(i+1))/2;
    rmids(i)=(span(i)+span(i+1))/2;
    rsize(i)=span(i+1)-span(i);
    Fmids(i)=(Fs(i+1)+Fs(i))/2;
    clsmid(i)=(cls(i+1)+cls(i))/2;
    cdsmid(i)=(cds(i+1)+cds(i))/2;
    chordmid(i)=(chordVector(i+1)+chordVector(i))/2;
    relwndsmid(i)=(relwnds(i+1)+relwnds(i))/2;
    %rotspeedmid(i)=(rotspeed(i+1)+rotspeed(i))/2;
end

%Actuator disk
    dT=Fmids.*1.225*(wndspeed^2)*4.*amids.*(1-amids)*pi.*rmids.*rsize;
    dQ=Fmids.*4.*apmids.*(1-amids).*1.225*wndspeed.*rotspeed*pi().*(rmids.^3).*rsize;
    
%Blade element momentum %%%%%%%its urel!!!!!! fixxxx
    dFn=(3/2)*1.225*(wndspeed^2).*(clsmid.*cos(relwndsmid)+cdsmid.*sin(relwndsmid)).*chordmid.*rsize;
    dQp=Fmids.*(3/2)*1.225*(wndspeed^2).*((clsmid).*sin(relwndsmid)-(cdsmid).*cos(relwndsmid)).*chordmid.*rmids.*rsize;

cpmids=4*amids.*((1-amids).^2);
ctmids=4*amids.*(1-amids);
cmmids=(8/3)*amids.*(1-amids);

annulusAreas=pi()*((rmids+rsize/2).^2-(rmids-rsize/2).^2);
ocp=sum(cpmids.*annulusAreas)/(pi()*R^2);
oct=sum(ctmids.*annulusAreas)/(pi()*R^2);
ocm=sum(cmmids.*annulusAreas)/(pi()*R^2);
dpower=dQ.*rotspeed;

moment=0.5*1.225*(wndspeed^2)*pi()*(R^3)*ocm;
%totalpwr=sum(dpower(2:end))*0.96; %JJM: Is this correct? I'm assuming some generator efficiency;
totalpwr=sum(dpower(2:end));
% also with more discretization dpower will converge to a more accurate value
