clear;clc;

mode=1; % Specify the calculation mode. 

mineral='px'; % Specify the mineral type. 
% Choose from: px ol pl ilm crsp mgsp kfs

numValidSamples = 100; % Specify the number of generated EPMA data you want
% If mode = 0 or 1, then the number will change with the number of test_data

% Import the basic database
if strcmp(mineral,'px')==1
    prmp = input('Which is the type of the pyroxene? If CaO<10 wt%, input L; otherwise, input H.',"s");
    if strcmp(prmp,'L')==1
        data=readmatrix('Highlands.xlsx','Sheet','LPx');
    elseif strcmp(prmp,'H')==1
        data=readmatrix('Highlands.xlsx','Sheet','HPx');
    else
        error('Input character must be one of H or L!');
    end
elseif strcmp(mineral,'pl')==1
    data=readmatrix('Highlands.xlsx','Sheet','Pl');
elseif strcmp(mineral,'ol')==1
    data=readmatrix('Highlands.xlsx','Sheet','Ol');
elseif strcmp(mineral,'ilm')==1
    data=readmatrix('Highlands.xlsx','Sheet','Ilm');
elseif strcmp(mineral,'crsp')==1
    data=readmatrix('Highlands.xlsx','Sheet','CrSp');
elseif strcmp(mineral,'mgsp')==1
    data=readmatrix('Highlands.xlsx','Sheet','MgSp');
elseif strcmp(mineral,'kfs')==1
    data=readmatrix('Highlands.xlsx','Sheet','Kfs');
else
    error('The set of mineral type is incorrect!');
end

% If modeling ilmenite and spinel data, only SiO2, TiO2, Al2O3, Cr2O3,
% FeO, MnO, MgO, CaO are needed.
if strcmp(mineral,'ilm')==1 || strcmp(mineral,'crsp')==1 || strcmp(mineral,'mgsp')==1
    data=data(:,5:12);
    data(data == 0) = NaN;
else
    data=data(:,5:15);
    data(data == 0) = NaN;
end

% Remove outliers from the EPMA database.
data=rmoutliers(data,'median');

if mode==0
    % Split the EPMA database into train data (train_data) and test data (test_data)
    % Statistics and Machine Learning Toolbox is needed.
    cv = cvpartition(size(data, 1), 'Holdout', 0.3);
    % The number of train_data versus test data is 7:3 in this default set.
    train_indices = training(cv);
    test_indices = test(cv);
    train_data = data(train_indices, :);
    test_data = data(test_indices, :);

    % Calculate the Spearman correlation coefficients (rho) of the train_data
    [rho, P] = corr(train_data, 'Type', 'Spearman', 'Rows', 'pairwise');

    % Calculate the mean and standard deviation values of the test_data
    test_mean=mean(test_data,'omitnan');
    test_std=std(test_data,'omitmissing');

elseif mode==1
    % Import the test_data from Test.xlsx. Make sure you have put your EPMA
    % data into the Test.xlsx datafile previously.

    % Calculate the Spearman correlation coefficients (rho) of the train_data
    % Using the whole basic database as the train_data
    [rho, P] = corr(data, 'Type', 'Spearman', 'Rows', 'pairwise');
   
    % Import the basic database
    if strcmp(mineral,'px')==1
        if strcmp(prmp,'L')==1
            test_data=readmatrix('Test.xlsx','Sheet','LPx');
        elseif strcmp(prmp,'H')==1
            test_data=readmatrix('Test.xlsx','Sheet','HPx');
        else
            error('Input character must be one of H or L!');
        end
    elseif strcmp(mineral,'pl')==1
        test_data=readmatrix('Test.xlsx','Sheet','Pl');
    elseif strcmp(mineral,'ol')==1
        test_data=readmatrix('Test.xlsx','Sheet','Ol');
    elseif strcmp(mineral,'ilm')==1
        test_data=readmatrix('Test.xlsx','Sheet','Ilm');
    elseif strcmp(mineral,'crsp')==1
        test_data=readmatrix('Test.xlsx','Sheet','CrSp');
    elseif strcmp(mineral,'mgsp')==1
        test_data=readmatrix('Test.xlsx','Sheet','MgSp');
    elseif strcmp(mineral,'kfs')==1
        test_data=readmatrix('Test.xlsx','Sheet','Kfs');
    else
        error('The set of mineral type is incorrect!');
    end

    if strcmp(mineral,'ilm')==1 || strcmp(mineral,'crsp')==1 || strcmp(mineral,'mgsp')==1
        test_data=test_data(:,5:12);
        test_data(test_data == 0) = NaN;
    else
        test_data=test_data(:,5:15);
        test_data(test_data == 0) = NaN;
    end

    % Calculate the mean and standard deviation values of the test_data
    test_mean=mean(test_data,'omitnan');
    test_std=std(test_data,'omitmissing');

elseif mode==2
    % Calculate the Spearman correlation coefficients (rho) of the train_data
    % Using the whole basic database as the train_data
    [rho, P] = corr(data, 'Type', 'Spearman', 'Rows', 'pairwise');
    
    % Input mean and std values of unknown dataset manually
    % SiO2 TiO2 Al2O3 Cr2O3 FeO MnO MgO CaO Na2O K2O P2O5
    test_mean=[53.25 0.5 0.75 0.5 20.4 0 23.2 1.25 0 0 0];
    test_std=[0.65 0.1 0.05 0.1 0.5 0.1 0.5 0.35 0.1 0.1 0.1];
end

cov_matrix = diag(test_std) * rho * diag(test_std);

% This step is to ensure the cov_matrix is positive semi-definite
% Function mvnrnd needs a positive semi-definite cov_matrix input
eig_values = eig(cov_matrix);  % Calculate eigenvalues
if all(eig_values >= 0)
    disp('The cov_matrix is OK.')
else
    warning('The covariance matrix is not positive semi-definite');
    [V, D] = eig(cov_matrix);  % Get eigenvalues and eigenvectors
    D(D < 0) = 0;              % Set negative eigenvalues to zero
    cov_matrix = V * D * V';   % Reconstruct the matrix
    warning('The covariance matrix has been reconstructed');
end

% If mode = 0 or 1, then numValidSamples will change with the number of test_data
% You can comment out these three lines of code if you want to use the 
% previously customized numValidSamples value.
if mode==0 || mode==1
    numValidSamples = size(test_data,1);
end

% Initialize variables
xValidGenerated = [];
numAttempts = 0;

% Regard test_data is an unknown dataset, generate random data using the 
% mvnrnd function based on the mean and std values of test_data
while size(xValidGenerated, 1) < numValidSamples
    num_samples = numValidSamples * (numAttempts+1);  % Adjust factor as needed
    samples = mvnrnd(test_mean, cov_matrix, num_samples);

    % Delete the data that do not satisfy the EPMA standards
    % Condition 1: Non-negative values in all dimensions
    nonNegativeSamples = all(samples >= 0, 2);

    % Condition 2: Sum within [98.5, 101.5]
    rowSums = sum(samples, 2, 'omitnan');
    sumWithinRange = (rowSums >= 98.5) & (rowSums <= 101.5);

    % Combine both conditions
    validSamples = nonNegativeSamples & sumWithinRange;
    xValidGenerated = samples(validSamples, :);
    numAttempts = numAttempts + 1;

    % Optionally, prevent infinite loop
    if numAttempts > 1000  % Adjust as needed
        warning('Could not generate enough valid samples after multiple attempts.');
        break;
    end
end

% Select the required number of samples
if size(xValidGenerated, 1) >= numValidSamples
    filtered_samples = xValidGenerated(1:numValidSamples, :);
else
    error('Not enough valid samples generated. Try increasing numNewSamples or numAttempts.');
end

% If it is in mode 0 or 1, then compare the generated data 
% (filtered_samples) with the test_data both qualitatively and
% quantitatively.
if mode==0 || mode==1
    % Quantitatively comparison: the Mann-Whitney U test
    Ionic_occup=MWU_test(filtered_samples,test_data,mineral);

    % Qualitatively comparison: visual comparison using diagrams
    combined_data = [filtered_samples; test_data];
    group_labels = [ones(size(filtered_samples, 1), 1); 2 * ones(size(test_data, 1), 1)];
    colors = 'rb';  markers = '.';
    figure;
    [h, ax, bigax] = gplotmatrix(combined_data, [], group_labels, colors, markers, [], 'on', '', '');
    axgd=["SiO2","TiO2","Al2O3","Cr2O3","FeO","MnO","MgO","CaO","Na2O","K2O","P2O5"];
    numVars = size(combined_data, 2);
    for i = 1:numVars
        xlabel(ax(numVars, i), sprintf(axgd(i)));
        ylabel(ax(i, 1), sprintf(axgd(i)));
    end
end

% Output the modeling results in .xlsx file
rowSums = sum(filtered_samples, 2, 'omitnan');
if strcmp(mineral,'px')==1
    if mode==2
        Ionic_occup=test_px(filtered_samples);
    end
    results=array2table(round([filtered_samples,rowSums,Ionic_occup],3),'VariableNames',...
        {'SiO2','TiO2','Al2O3','Cr2O3','FeO','MnO','MgO','CaO','Na2O','K2O','P2O5','Total',...
        'Si+P', 'Al(T)', 'Fe3(T)', 'sumT', 'Al(M1)', 'Ti(M1)', 'Cr(M1)', 'Fe3(M1)', 'Mn(M1)',...
        'Mg(M1)', 'Fe2+(M1)', 'sumM1', 'Mg(M2)', 'Fe2+(M2)', 'Ca(M2)', 'Na(M2)', 'K(M2)', 'sumM2'});
    writetable(results,['Output_', prmp, mineral,'.xlsx']);
elseif strcmp(mineral,'pl')==1 || strcmp(mineral,'kfs')==1
    if mode==2
        Ionic_occup=test_pl(filtered_samples);
    end
    results=array2table(round([filtered_samples,rowSums,Ionic_occup],3),'VariableNames',...
        {'SiO2','TiO2','Al2O3','Cr2O3','FeO','MnO','MgO','CaO','Na2O','K2O','P2O5','Total',...
        'Si', 'Ti', 'Al', 'Cr', 'Fe', 'Mn', 'Mg', 'Ca', 'Na', 'K', 'P', 'sum'});
    writetable(results,['Output_', mineral,'.xlsx']);
elseif strcmp(mineral,'ol')==1
    if mode==2
        Ionic_occup=test_ol(filtered_samples);
    end
    results=array2table(round([filtered_samples,rowSums,Ionic_occup],3),'VariableNames',...
        {'SiO2','TiO2','Al2O3','Cr2O3','FeO','MnO','MgO','CaO','Na2O','K2O','P2O5','Total',...
        'Si+P', 'Al(T)', 'Fe3(T)', 'sumT', 'Al', 'Ti', 'Cr', 'Fe3', 'Mn',...
        'Fe2', 'Mg', 'Ca', 'Na', 'K', 'sum'});
    writetable(results,['Output_', mineral,'.xlsx']);
elseif strcmp(mineral,'ilm')==1 || strcmp(mineral,'crsp')==1 || strcmp(mineral,'mgsp')==1
    if mode==2
        if strcmp(mineral,'ilm')==1
            Ionic_occup=test_ilm(filtered_samples);
        else
            Ionic_occup=test_sp(filtered_samples);
        end
    end
    results=array2table(round([filtered_samples,rowSums,Ionic_occup],3),'VariableNames',...
        {'SiO2','TiO2','Al2O3','Cr2O3','FeO','MnO','MgO','CaO','Total',...
        'Si', 'Ti', 'Al', 'Fe3', 'Cr', 'Fe2', 'Mn', 'Mg', 'Ca', 'sum', 'O'});
    writetable(results,['Output_', mineral,'.xlsx']);
end