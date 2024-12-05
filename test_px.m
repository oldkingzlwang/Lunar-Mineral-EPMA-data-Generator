function StrctFrm=test_px(data)
% Modified from Jesse B. Walters's program: MinPlot (Mineralogia, vol. 53, 2022, 51-66)

[m,n]=size(data);
for k = 1:n
    data(isnan(data(:, k)), k) = 0;
end
MC=zeros(size(data));
SiO2=60.083;TiO2=79.865;Al2O3=101.961;
Cr2O3=151.989;FeO=71.8442;MnO=70.937;
MgO=40.304;CaO=56.0774;Na2O=61.979;
K2O=94.195;P2O5=141.943;
W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O,P2O5];
for i=1:size(W,2)
    if i==3 || i==4 || i==9 || i==10 || i==11
        MC(:,i)=data(:,i)./W(:,i).*2;
    else
        MC(:,i)=data(:,i)./W(:,i);
    end
end
cat=4.0; %cations per formula unit
Opfu=6.0; %oxygens per formula unit
MCnormfact=cat./sum(MC,2); %normalization factor
MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations
O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5); %for FeO
O2(:,6)=MCnorm(:,6); %for MnO
O2(:,7)=MCnorm(:,7); %for MgO
O2(:,8)=MCnorm(:,8); %for CaO
O2(:,9)=MCnorm(:,9)./2; %for Na2O
O2(:,10)=MCnorm(:,10)./2; %for K2O
O2(:,11)=MCnorm(:,11).*(5/2); %for P2O5

O2total=sum(O2,2); %O2 totals

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr
APFU(:,7)=MCnorm(:,6); %for Mn
APFU(:,8)=MCnorm(:,7); %for Mg
APFU(:,9)=MCnorm(:,8); %for Ca
APFU(:,10)=MCnorm(:,9); %for Na
APFU(:,11)=MCnorm(:,10); %for K
APFU(:,12)=MCnorm(:,11); %for P

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 6
%if so, then there is no Fe3+
%if totalO2 < 6, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(6-totalO2) then the amount
%of Fe3+ = 2*(6-totalO2), if false then, all Fe is Fe3+
for c=1:m
    if (Opfu-O2total(c,1)) > 0
        if MCnorm(c,5) > (Opfu-O2total(c,1))*2
            APFU(c,4)=(Opfu-O2total(c,1))*2; 
        else
            APFU(c,4)=MCnorm(c,5);
        end
    else
        APFU(c,4)=0;
    end
end

APFU(:,6)=MCnorm(:,5)-APFU(:,4); %the APFU of Fe2+ = totalFe-Fe3+

APFU(:,13)=sum(APFU,2); %calculations the total, which should be 4
% Oxygen deficiency 
APFU(:,14)=Opfu-O2total; %must be greater than zero

XMg=APFU(:,8)./(APFU(:,8)+APFU(:,6)); 

%% structural formula calculation

%T SITE
%Si+P
for c=1:m
    if APFU(c,1)+APFU(c,12)<2.000
        StrctFrm(c,1)=APFU(c,1)+APFU(c,12); %If Si+P < 2, then Si+P(T) = the measured Si+P content
    else
        StrctFrm(c,1)=2; %If Si+P is in excess, then Si+P(T) = 2
    end
end

%Al(T)
for c=1:m
    if 2-StrctFrm(c,1)>0 %Is 2-(Si+P) > 0? If y, then some Al goes into T
        if 2-StrctFrm(c,1)>APFU(c,3) %For low Al cpx, 2-(Si+P) may be > Al
            StrctFrm(c,2)=APFU(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=2-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to M1
        end
    else
        StrctFrm(c,2)=0; %if Si=2, then no Al goes into T
    end
end

%Fe3+(T)
for c=1:m
    if 2-StrctFrm(c,1)-StrctFrm(c,2)>0 %Is 2-(Si+Al) > 0? If y, then some Fe3+ goes into T
        if 2-StrctFrm(c,1)-StrctFrm(c,2)>APFU(c,4) %For low Fe3+ cpx, 2-(Si+Al) may be > Fe3+
            StrctFrm(c,3)=APFU(c,4); %All Fe3+ goes into T
        else
            StrctFrm(c,3)=2-StrctFrm(c,1)-StrctFrm(c,2); %if there isn't enough space in T for all Fe3+, the rest will go to M1
        end
    else
        StrctFrm(c,3)=0; %if Si+Al=2, then no Fe3+ goes into T
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3); %Si + P + Al + Fe3+ in T

%M1 SITE

%Al(M1)
StrctFrm(:,5)=APFU(:,3)-StrctFrm(:,2); %Al(M1) = Total Al - Al(T)

%Ti (M1)
StrctFrm(:,6)=APFU(:,2);

%Fe3+ (M1)
StrctFrm(:,8)=APFU(:,4)-StrctFrm(:,3); %Fe3+(M1) = Total Fe3+ - Fe3+(T)

%Cr3+ (M1)
StrctFrm(:,7)=APFU(:,5);

%Mn (M1)
StrctFrm(:,9)=APFU(:,7);

%Mg (M1)
for c=1:m
    if XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5))<APFU(c,8) %if XMg*(1-Sum(Al to Mn) in M1) is < Mg
        StrctFrm(c,10)=XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5)); %Mg(M1)=XMg*(1-Sum(Al to Mn) in M1) and some Mg goes into M2
    else
        StrctFrm(c,10)=APFU(c,8); %if not, all Mg goes into M1
    end
end

%Fe2+ (M1)
for c=1:m
    if 1-(StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,5))>0 %Is 1-Sum(Al to Mg) in M1 >0?, if yes then:
        if 1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5)>APFU(c,6) %Is 1-Sum(Al to Mg) in M1 > Fe2+?, if yes then:
            StrctFrm(c,11)=APFU(c,6); %all Fe2+ goes into M1
        else
            StrctFrm(c,11)=1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5); %if there isn't enough space in M1 for all Fe2+, some goes into M2
        end
    else
        StrctFrm(c,11)=0; %If M1 is already filled, then no Fe2+ goes into M1
    end
end

%Sum of M1 site
StrctFrm(:,12)=sum(StrctFrm(:,5:11),2);

%M2 SITE
%Mg (M2)
for c=1:m
    if APFU(c,8)-StrctFrm(c,10)>0 %Is Mgtotal-Mg(M1) > 0?
        StrctFrm(c,13)=APFU(c,8)-StrctFrm(c,10); %if yes, then some Mg goes into M2
    else
        StrctFrm(c,13)=0; %if no, then all Mg is in M1
    end
end

%Fe2+ (M2)
for c=1:m
    if APFU(c,6)-StrctFrm(c,11)>0 %Is Fe2+total-Fe2+(M1) > 0?
        StrctFrm(c,14)=APFU(c,6)-StrctFrm(c,11); %if yes, then some Fe2+ goes into M2
    else
        StrctFrm(c,14)=0; %if no, then all Fe2+ is in M1
    end
end

%Ca (M2)
StrctFrm(:,15)=APFU(:,9);

%Na (M2)
StrctFrm(:,16)=APFU(:,10);

%K (M2)
StrctFrm(:,17)=APFU(:,11);

%Sum of M2 site
StrctFrm(:,18)=StrctFrm(:,13)+StrctFrm(:,14)+StrctFrm(:,15)+StrctFrm(:,16)+StrctFrm(:,17);

%           1    2     3       4    5      6      7        8        9
% StrctFrm: Si+P Al(T) Fe3+(T) sumT Al(M1) Ti(M1) Cr3+(M1) Fe3+(M1) Mn(M1)
% 10     11       12    13     14       15     16     17    18
% Mg(M1) Fe2+(M1) sumM1 Mg(M2) Fe2+(M2) Ca(M2) Na(M2) K(M2) sumM2
end