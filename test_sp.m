function APFU=test_sp(data)
% Modified from Jesse B. Walters's program: MinPlot (Mineralogia, vol. 53, 2022, 51-66)

[m,n]=size(data);
for k = 1:n
    data(isnan(data(:, k)), k) = 0;
end
MC=zeros(size(data));
SiO2=60.083;TiO2=79.865;Al2O3=101.961;
Cr2O3=151.989;FeO=71.8442;MnO=70.937;
MgO=40.304;CaO=56.0774;
W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO];
for i=1:size(W,2)
    if i==3 || i==4
        MC(:,i)=data(:,i)./W(:,i).*2;
    else
        MC(:,i)=data(:,i)./W(:,i);
    end
end
cat=3.0; %cations per formula unit
Opfu=4.0; %oxygens per formula unit

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

O2total=sum(O2,2); %O2 totals

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr
APFU(:,7)=MCnorm(:,6); %for Mn
APFU(:,8)=MCnorm(:,7); %for Mg
APFU(:,9)=MCnorm(:,8); %for Ca

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 4
%if so, then there is no Fe3+
%if totalO2 < 4, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(4-totalO2) then the amount
%of Fe3+ = 2*(4-totalO2), if false then, all Fe is Fe3+

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

APFU(:,10)=sum(APFU,2); %calculations the total, which should be 4
% Oxygen deficiency 
APFU(:,11)=Opfu-O2total; %must be greater than zero

%       1  2  3  4    5  6    7    8   9  10 11
% APFU: Si Ti Al Fe3+ Cr Fe2+ Mn   Mg  Ca sum O
end