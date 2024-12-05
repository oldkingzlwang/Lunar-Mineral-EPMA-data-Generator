function StrctFrm1=MWU_test(filtered_samples,test_data,mineral)
% Perform the Mann-Whitney U test

%% Perform the Mann-Whitney U test for pyroxene
if strcmp(mineral,'px')==1
    StrctFrm1=test_px(filtered_samples);
    StrctFrm2=test_px(test_data);

    [p1, hyp1, ~] = ranksum(StrctFrm1(:,6)./(StrctFrm1(:,2)+StrctFrm1(:,5)), ...
        StrctFrm2(:,6)./(StrctFrm2(:,2)+StrctFrm2(:,5))); % Test 1: Ti/Al
    fprintf('P-value: %.4f\n', p1);
    if hyp1 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p2, hyp2, ~] = ranksum((StrctFrm1(:,10)+StrctFrm1(:,13))./(StrctFrm1(:,11)+StrctFrm1(:,14)), ...
        (StrctFrm2(:,10)+StrctFrm2(:,13))./(StrctFrm2(:,11)+StrctFrm2(:,14))); % Test 2: Mg/Fe
    fprintf('P-value: %.4f\n', p2);
    if hyp2 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p3, hyp3, ~] = ranksum(StrctFrm1(:,2)+StrctFrm1(:,3)-StrctFrm1(:,16), ...
        StrctFrm2(:,2)+StrctFrm2(:,3)-StrctFrm2(:,16)); % Test 3: Al6+Fe3-Na
    fprintf('P-value: %.4f\n', p3);
    if hyp3 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p4, hyp4, ~] = ranksum(StrctFrm1(:,16)+StrctFrm1(:,17)+StrctFrm1(:,15), ...
        StrctFrm2(:,16)+StrctFrm2(:,17)+StrctFrm2(:,15)); % Test 4: Na+K+Ca
    fprintf('P-value: %.4f\n', p4);
    if hyp4 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

%% Perform the Mann-Whitney U test for feldspar
elseif strcmp(mineral,'pl')==1 || strcmp(mineral,'kfs')==1
    StrctFrm1=test_fsp(filtered_samples);
    StrctFrm2=test_fsp(test_data);

    [p1, hyp1, ~] = ranksum(StrctFrm1(:,8)./(StrctFrm1(:,8)+StrctFrm1(:,9)+StrctFrm1(:,10)), ...
        StrctFrm2(:,8)./(StrctFrm2(:,8)+StrctFrm2(:,9)+StrctFrm2(:,10))); % Test 1: An#
    fprintf('P-value: %.4f\n', p1);
    if hyp1 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p2, hyp2, ~] = ranksum(StrctFrm1(:,12), StrctFrm2(:,12)); % Test 2: total
    fprintf('P-value: %.4f\n', p2);
    if hyp2 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

%% Perform the Mann-Whitney U test for olivine
elseif strcmp(mineral,'ol')==1
    StrctFrm1=test_ol(filtered_samples);
    StrctFrm2=test_ol(test_data);

    [p1, hyp1, ~] = ranksum(StrctFrm1(:,11)./StrctFrm1(:,10), ...
        StrctFrm2(:,11)./StrctFrm2(:,10)); % Test 1: Mg/Fe
    fprintf('P-value: %.4f\n', p1);
    if hyp1 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

%% Perform the Mann-Whitney U test for ilmenite
elseif strcmp(mineral,'ilm')==1
    StrctFrm1=test_ilm(filtered_samples);
    StrctFrm2=test_ilm(test_data);

    [p1, hyp1, ~] = ranksum(StrctFrm1(:,4)./(StrctFrm1(:,4)+StrctFrm1(:,2)), ...
        StrctFrm2(:,4)./(StrctFrm2(:,4)+StrctFrm2(:,2))); % Test 1: Fe3/(Fe3+Ti)
    fprintf('P-value: %.4f\n', p1);
    if hyp1 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p2, hyp2, ~] = ranksum(StrctFrm1(:,6)./(StrctFrm1(:,6)+StrctFrm1(:,7)+StrctFrm1(:,8)), ...
        StrctFrm2(:,6)./(StrctFrm2(:,6)+StrctFrm2(:,7)+StrctFrm2(:,8))); % Test 2: Fe2/(Fe2+Mn+Mg)
    fprintf('P-value: %.4f\n', p2);
    if hyp2 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p3, hyp3, ~] = ranksum(StrctFrm1(:,2)./(StrctFrm1(:,2)+StrctFrm1(:,3)+StrctFrm1(:,5)), ...
        StrctFrm2(:,2)./(StrctFrm2(:,2)+StrctFrm2(:,3)+StrctFrm2(:,5))); % Test 3: Ti/(Ti+Al+Cr)
    fprintf('P-value: %.4f\n', p3);
    if hyp3 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

%% Perform the Mann-Whitney U test for spinel
elseif strcmp(mineral,'crsp')==1 || strcmp(mineral,'mgsp')==1
    StrctFrm1=test_sp(filtered_samples);
    StrctFrm2=test_sp(test_data);

    [p1, hyp1, ~] = ranksum(StrctFrm1(:,4)./(StrctFrm1(:,4)+StrctFrm1(:,2)), ...
        StrctFrm2(:,4)./(StrctFrm2(:,4)+StrctFrm2(:,2))); % Test 1: Fe3/(Fe3+Ti)
    fprintf('P-value: %.4f\n', p1);
    if hyp1 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p2, hyp2, ~] = ranksum(StrctFrm1(:,6)./(StrctFrm1(:,6)+StrctFrm1(:,7)+StrctFrm1(:,8)), ...
        StrctFrm2(:,6)./(StrctFrm2(:,6)+StrctFrm2(:,7)+StrctFrm2(:,8))); % Test 2: Fe2/(Fe2+Mn+Mg)
    fprintf('P-value: %.4f\n', p2);
    if hyp2 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

    [p3, hyp3, ~] = ranksum(StrctFrm1(:,2)./(StrctFrm1(:,2)+StrctFrm1(:,3)+StrctFrm1(:,5)), ...
        StrctFrm2(:,2)./(StrctFrm2(:,2)+StrctFrm2(:,3)+StrctFrm2(:,5))); % Test 3: Ti/(Ti+Al+Cr)
    fprintf('P-value: %.4f\n', p3);
    if hyp3 == 1
        disp('There is a significant difference between the two datasets.');
    else
        disp('There is no significant difference between the two datasets.');
    end

else
    error('The set of mineral type is incorrect!');
end