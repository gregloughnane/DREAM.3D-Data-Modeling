function [mbcfull, mbclow, mbchigh, ds_meanerrors] = ComputeGrainStatistics(no_bins_shape, no_bins_quantiles, qlow, qhigh)
% This function computes the modified bhattacharyya coefficient for grain
% size distributions from output .csv files created using PipelineRunner
% no_bins_shape is the number of bins for full distribution histograms of b/a, c/a
% no_bins_quantiles is the number of bins for partial distribution histograms
% qlow is the low quantile of interest for MBC computation
% qhigh is the high "quantile" of interest for MBC computation
% Notes:
% The number of bins for grain size distribution comparisons is hard-coded
% The columns from which to retrieve data from CSV files are hard-coded
% This code is running with version 6.0-Beta-7afa1de-Win64
 
format compact
[RefFile, DSFiles ] = GetCSVs( ); 
% Ref = Reference Structure, DS = Downsampled
 
%% Get Data From ORIGINAL Structure
orig_data = csvread(RefFile,3,0);
[numrows,numcols] = size(orig_data);
cnt=0;
for j = 1:numrows
% Note (cnt,:) makes the resulting vectors to be column vectors
    if orig_data(j,14) == 0; %sort using bounding box algorithm
        cnt = cnt+1;
        orig_grainids(cnt,:) = orig_data(j,1);
        orig_esd(cnt,:) = orig_data(j,18);
        orig_ba(cnt,:) = orig_data(j,2);
        orig_ca(cnt,:) = orig_data(j,3);
        orig_nnn(cnt,:) = orig_data(j,24);
    end
end
orig_esd_sorted = orig_esd(orig_esd~=0); % Remove any zero volume grains
orig_ba_sorted = orig_ba(orig_ba~=0);
orig_ca_sorted = orig_ca(orig_ca~=0); 
orig_nnn_sorted = orig_nnn(orig_nnn~=0); 
 
%% Define binning for histograms 
mid_bins_esd = [0.25:0.5:13.75]; %%%  Define bins based on grain sizes
bins_esd = [mid_bins_esd inf];
bins_shape = [0:1/no_bins_shape:1]; %%%  b/a, c/a
bins_nnn = [0:max(orig_nnn_sorted)]; %%%  Number of Nearest Neighbors
 
%%% Compute Mean & Standard Deviation for each reference data set    
orig_mean_esd = mean(orig_esd_sorted); orig_std_esd = std(orig_esd_sorted);
orig_mean_ba = mean(orig_ba_sorted); orig_std_ba = std(orig_ba_sorted);
orig_mean_ca = mean(orig_ca_sorted); orig_std_ca = std(orig_ca_sorted);
orig_mean_nnn = mean(orig_nnn_sorted); orig_std_nnn = std(orig_nnn_sorted);
 
%%% Determine quantiles and bins for quantile analysis
tp_orig_low_esd = quantile(orig_esd_sorted, qlow);
tp_orig_high_esd = quantile(orig_esd_sorted, qhigh);
tp_orig_low_ba = quantile(orig_ba_sorted, qlow);
tp_orig_high_ba = quantile(orig_ba_sorted, qhigh);
tp_orig_low_ca = quantile(orig_ca_sorted, qlow);
tp_orig_high_ca = quantile(orig_ca_sorted, qhigh);
tp_orig_low_nnn = quantile(orig_nnn_sorted, qlow);
tp_orig_high_nnn = quantile(orig_nnn_sorted, qhigh);
 
LowBins_esd = linspace(0,tp_orig_low_esd,no_bins_quantiles);
LowBins_ba = linspace(0,tp_orig_low_ba,no_bins_quantiles);
LowBins_ca = linspace(0,tp_orig_low_ca,no_bins_quantiles);
LowBins_nnn = linspace(0,tp_orig_low_nnn,no_bins_quantiles);
 
HighBins_esd = linspace(tp_orig_high_esd,max(mid_bins_esd),no_bins_quantiles);
HighBins_ba = linspace(tp_orig_high_ba,1,no_bins_quantiles);
HighBins_ca = linspace(tp_orig_high_ca,1,no_bins_quantiles);
HighBins_nnn = linspace(tp_orig_high_nnn,max(orig_nnn_sorted),no_bins_quantiles);
 
%%% Partition original data set into low and high quantile data sets
for w = 1:length(orig_esd_sorted)
    if orig_esd_sorted(w) < tp_orig_low_esd;
        orig_esd_low(w) = orig_esd_sorted(w); % Partition for low
    end
    if orig_esd_sorted(w) > tp_orig_high_esd;
        orig_esd_high(w) = orig_esd_sorted(w); % Partition for high
    end
end
for w = 1:length(orig_ba_sorted);
    if orig_ba_sorted(w) < tp_orig_low_ba;
        orig_ba_low(w) = orig_ba_sorted(w); % Partition for low
    end
    if orig_ba_sorted(w) > tp_orig_high_ba;
        orig_ba_high(w) = orig_ba_sorted(w); % Partition for high
    end
end
for w = 1:length(orig_ca_sorted);
    if orig_ca_sorted(w) < tp_orig_low_ca;
        orig_ca_low(w) = orig_ca_sorted(w); % Partition for low 
    end
    if orig_ca_sorted(w) > tp_orig_high_ca;
        orig_ca_high(w) = orig_ca_sorted(w); % Partition for high 
    end
end
for w = 1:length(orig_nnn_sorted);
    if orig_nnn_sorted(w) < tp_orig_low_nnn;
        orig_nnn_low(w) = orig_nnn_sorted(w);; % Partition for low 
    end
    if orig_nnn_sorted(w) > tp_orig_high_nnn;
        orig_nnn_high(w) = orig_nnn_sorted(w); % Partition for high  
    end
end
orig_esd_low_sorted = orig_esd_low(orig_esd_low~=0); % Remove any zeros
orig_esd_high_sorted = orig_esd_high(orig_esd_high~=0); % Remove any zeros
orig_ba_low_sorted = orig_ba_low(orig_ba_low~=0); % Remove any zeros 
orig_ba_high_sorted = orig_ba_high(orig_ba_high~=0); % Remove any zeros 
orig_ca_low_sorted = orig_ca_low(orig_ca_low~=0); % Remove any zeros 
orig_ca_high_sorted = orig_ca_high(orig_ca_high~=0); % Remove any zeros 
orig_nnn_low_sorted = orig_nnn_low(orig_nnn_low~=0); % Remove any zeros 
orig_nnn_high_sorted = orig_nnn_high(orig_nnn_high~=0); % Remove any zeros 
 
%% Construct original structure histograms
%%% FULL
esd_o = histc(orig_esd_sorted,bins_esd);
norm_esd_o = esd_o./(length(orig_esd_sorted)); % normalize
ba_o = histc(orig_ba_sorted,bins_shape);
norm_ba_o = ba_o./(length(orig_ba_sorted)); % normalize
ca_o = histc(orig_ca_sorted,bins_shape);
norm_ca_o = ca_o./(length(orig_ca_sorted)); % normalize 
nnn_o = histc(orig_nnn_sorted,bins_nnn);
norm_nnn_o = nnn_o./(length(orig_nnn_sorted)); % normalize
%%% LOW
esd_o_low = histc(orig_esd_low_sorted,LowBins_esd);
norm_esd_o_low = esd_o_low./(length(orig_esd_low_sorted)); % normalize
ba_o_low = histc(orig_ba_low_sorted,LowBins_ba);
norm_ba_o_low = ba_o_low./(length(orig_ba_low_sorted)); % normalize 
ca_o_low = histc(orig_ca_low_sorted,LowBins_ca);
norm_ca_o_low = ca_o_low./(length(orig_ca_low_sorted)); % normalize 
nnn_o_low = histc(orig_nnn_low_sorted,LowBins_nnn);
norm_nnn_o_low = nnn_o_low./(length(orig_nnn_low_sorted)); % normalize 
%%% HIGH
esd_o_high = histc(orig_esd_high_sorted,HighBins_esd);
norm_esd_o_high = esd_o_high./(length(orig_esd_high_sorted)); % normalize
ba_o_high = histc(orig_ba_high_sorted,HighBins_ba);
norm_ba_o_high = ba_o_high./(length(orig_ba_high_sorted)); % normalize 
ca_o_high = histc(orig_ca_high_sorted,HighBins_ca);
norm_ca_o_high = ca_o_high./(length(orig_ca_high_sorted)); % normalize 
nnn_o_high = histc(orig_nnn_high_sorted,HighBins_nnn);
norm_nnn_o_high = nnn_o_high./(length(orig_nnn_high_sorted)); % normalize 
 
[numrows,numcols] = size(DSFiles);
for k = 1:numrows;
    %% Get Data From Structures with Errors
    ds_file = DSFiles(k,:);
    ds_data = csvread(ds_file,3,0);
    [numrows,numcols] = size(ds_data);
    cnt=0;
    for j = 1:numrows;
    % Note (cnt,:) makes the resulting vectors to be column vectors
        if ds_data(j,10) == 0; %sort using bounding box algorithm
            cnt = cnt+1;
            ds_grainids(cnt,:) = ds_data(j,1);
            ds_esd(cnt,:) = ds_data(j,14);
            ds_ba(cnt,:) = ds_data(j,2);
            ds_ca(cnt,:) = ds_data(j,3);
            ds_nnn(cnt,:) = ds_data(j,16);
        end
    end
    ds_esd_sorted = ds_esd(ds_esd~=0); % Remove any zeros
    ds_ba_sorted = ds_ba(ds_ba~=0); 
    ds_ca_sorted = ds_ca(ds_ca~=0); 
    ds_nnn_sorted = ds_nnn(ds_nnn~=0); 
 
    %% Construct down-sampled structure histograms
    %%% Construct down-sampled structure histogram (FULL)
    esd_ds = histc(ds_esd_sorted,bins_esd); 
    norm_esd_ds = esd_ds./(length(ds_esd_sorted)); % normalize
    ba_ds = histc(ds_ba_sorted,bins_shape); 
    norm_ba_ds = ba_ds./(length(ds_ba_sorted)); % normalize
    ca_ds = histc(ds_ca_sorted,bins_shape); 
    norm_ca_ds = ca_ds./(length(ds_ca_sorted)); % normalize
    nnn_ds = histc(ds_nnn_sorted,bins_nnn); 
    norm_nnn_ds = nnn_ds./(length(ds_nnn_sorted)); % normalize
    
    %%% Partition original data set into low and high quantile data sets
    for w = 1:length(ds_esd_sorted);
        if ds_esd_sorted(w) < tp_orig_low_esd;
            ds_esd_low_sorted(w) = ds_esd_sorted(w); % Partition for low
        end
        if ds_esd_sorted(w) > tp_orig_high_esd;
            ds_esd_high_sorted(w) = ds_esd_sorted(w); % Partition for high
        end
    end
    for w = 1:length(ds_ba_sorted);
        if ds_ba_sorted(w) < tp_orig_low_ba;
            ds_ba_low_sorted(w) = ds_ba_sorted(w);
        end
        if ds_ba_sorted(w) > tp_orig_high_ba;
            ds_ba_high_sorted(w) = ds_ba_sorted(w);
        end
    end    
    for w = 1:length(ds_ca_sorted);
        if ds_ca_sorted(w) < tp_orig_low_ca;
            ds_ca_low_sorted(w) = ds_ca_sorted(w);
        end
        if ds_ca_sorted(w) > tp_orig_high_ca;
            ds_ca_high_sorted(w) = ds_ca_sorted(w);
        end
    end        
    for w = 1:length(ds_nnn_sorted);
        if ds_nnn_sorted(w) < tp_orig_low_nnn;
            ds_nnn_low_sorted(w) = ds_nnn_sorted(w);
        end
        if ds_nnn_sorted(w) > tp_orig_high_nnn;
            ds_nnn_high_sorted(w) = ds_nnn_sorted(w);
        end
    end       
    ds_esd_low_sorted = ds_esd_low_sorted(ds_esd_low_sorted~=0); % Remove any zeros
    ds_esd_high_sorted = ds_esd_high_sorted(ds_esd_high_sorted~=0); % Remove any zeros 
    ds_ba_low_sorted = ds_ba_low_sorted(ds_ba_low_sorted~=0); % Remove any zeros
    ds_ba_high_sorted = ds_ba_high_sorted(ds_ba_high_sorted~=0); % Remove any zeros
    ds_ca_low_sorted = ds_ca_low_sorted(ds_ca_low_sorted~=0); % Remove any zeros 
    ds_ca_high_sorted = ds_ca_high_sorted(ds_ca_high_sorted~=0); % Remove any zeros
    ds_nnn_low_sorted = ds_nnn_low_sorted(ds_nnn_low_sorted~=0); % Remove any zeros 
    ds_nnn_high_sorted = ds_nnn_high_sorted(ds_nnn_high_sorted~=0); % Remove any zeros 
    %%% LOW
    esd_ds_low = histc(ds_esd_low_sorted,LowBins_esd);
    norm_esd_ds_low = esd_ds_low./(length(ds_esd_low_sorted)); % normalize
    ba_ds_low = histc(ds_ba_low_sorted,LowBins_ba);
    norm_ba_ds_low = ba_ds_low./(length(ds_ba_low_sorted)); % normalize
    ca_ds_low = histc(ds_ca_low_sorted,LowBins_ca);
    norm_ca_ds_low = ca_ds_low./(length(ds_ca_low_sorted)); % normalize 
    nnn_ds_low = histc(ds_nnn_low_sorted,LowBins_nnn);
    norm_nnn_ds_low = nnn_ds_low./(length(ds_nnn_low_sorted)); % normalize 
    %%% HIGH
    esd_ds_high = histc(ds_esd_high_sorted,HighBins_esd);
    norm_esd_ds_high = esd_ds_high./(length(ds_esd_high_sorted)); % normalize
    ba_ds_high = histc(ds_ba_high_sorted,HighBins_ba);
    norm_ba_ds_high = ba_ds_high./(length(ds_ba_high_sorted)); % normalize
    ca_ds_high = histc(ds_ca_high_sorted,HighBins_ca);
    norm_ca_ds_high = ca_ds_high./(length(ds_ca_high_sorted)); % normalize
    nnn_ds_high = histc(ds_nnn_high_sorted,HighBins_nnn);
    norm_nnn_ds_high = nnn_ds_high./(length(ds_nnn_high_sorted)); % normalize
    
    %% Compute MBCs
    mbc_esd(k,:) = real(sqrt(1-sum(sqrt(norm_esd_o.*norm_esd_ds))));
    mbc_esd_low(k,:) = real(sqrt(1-sum(sqrt(norm_esd_o_low.*norm_esd_ds_low))));
    mbc_esd_high(k,:) = real(sqrt(1-sum(sqrt(norm_esd_o_high.*norm_esd_ds_high))));
    mbc_ba(k,:) = real(sqrt(1-sum(sqrt(norm_ba_o.*norm_ba_ds))));
    mbc_ba_low(k,:) = real(sqrt(1-sum(sqrt(norm_ba_o_low.*norm_ba_ds_low))));
    mbc_ba_high(k,:) = real(sqrt(1-sum(sqrt(norm_ba_o_high.*norm_ba_ds_high))));
    mbc_ca(k,:) = real(sqrt(1-sum(sqrt(norm_ca_o.*norm_ca_ds))));
    mbc_ca_low(k,:) = real(sqrt(1-sum(sqrt(norm_ca_o_low.*norm_ca_ds_low))));
    mbc_ca_high(k,:) = real(sqrt(1-sum(sqrt(norm_ca_o_high.*norm_ca_ds_high))));
    mbc_nnn(k,:) = real(sqrt(1-sum(sqrt(norm_nnn_o.*norm_nnn_ds))));
    mbc_nnn_low(k,:) = real(sqrt(1-sum(sqrt(norm_nnn_o_low.*norm_nnn_ds_low))));
    mbc_nnn_high(k,:) = real(sqrt(1-sum(sqrt(norm_nnn_o_high.*norm_nnn_ds_high))));
    
    %% Compute Mean & error relative to reference for each down-sampled data set
    ds_mean_esd(k,:) = mean(ds_esd_sorted); ds_std_esd(k,:) = std(ds_esd_sorted);
    ds_mean_error_esd(k,:) = 100*(ds_mean_esd(k,:) - orig_mean_esd)/orig_mean_esd;
    ds_mean_ba(k,:) = mean(ds_ba_sorted); ds_std_ba(k,:) = std(ds_ba_sorted);
    ds_mean_error_ba(k,:) = 100*(ds_mean_ba(k,:) - orig_mean_ba)/orig_mean_ba;
    ds_mean_ca(k,:) = mean(ds_ca_sorted); ds_std_ca(k,:) = std(ds_ca_sorted);
    ds_mean_error_ca(k,:) = 100*(ds_mean_ca(k,:) - orig_mean_ca)/orig_mean_ca;
    ds_mean_nnn(k,:) = mean(ds_nnn_sorted); ds_std_nnn(k,:) = std(ds_nnn_sorted);
    ds_mean_error_nnn(k,:) = 100*(ds_mean_nnn(k,:) - orig_mean_nnn)/orig_mean_nnn;
    
    %% Store dist data for each down-sampled set
    esd_ds_stack(:,k) = norm_esd_ds;
    ba_ds_stack(:,k) = norm_ba_ds;
    ca_ds_stack(:,k) = norm_ca_ds;
    nnn_ds_stack(:,k) = norm_nnn_ds;
    
    clear ds_esd ds_esd_low ds_esd_high ds_ba ds_ba_low ds_ba_high ds_ca ds_ca_low ds_ca_high
end
%% Output
% Avg & CI for each full MBC
mbc_esd_avg = mean(mbc_esd); 
mbc_esd_moe = std(mbc_esd)/sqrt(length(mbc_esd)) * tinv(0.975,length(mbc_esd)-1);
mbc_ba_avg = mean(mbc_ba); 
mbc_ba_moe = std(mbc_ba)/sqrt(length(mbc_ba)) * tinv(0.975,length(mbc_ba)-1);
mbc_ca_avg = mean(mbc_ca); 
mbc_ca_moe = std(mbc_ca)/sqrt(length(mbc_ca)) * tinv(0.975,length(mbc_ca)-1);
mbc_nnn_avg = mean(mbc_nnn); 
mbc_nnn_moe = std(mbc_nnn)/sqrt(length(mbc_nnn)) * tinv(0.975,length(mbc_nnn)-1);
mbcfull = [mbc_esd_avg mbc_esd_moe
            mbc_ba_avg  mbc_ba_moe
            mbc_ca_avg  mbc_ca_moe
            mbc_nnn_avg  mbc_nnn_moe];
% Avg & CI for each low MBC        
mbc_esd_avg_low = mean(mbc_esd_low); 
mbc_esd_moe_low = std(mbc_esd_low)/sqrt(length(mbc_esd_low)) * tinv(0.975,length(mbc_esd_low)-1);
mbc_ba_avg_low = mean(mbc_ba_low); 
mbc_ba_moe_low = std(mbc_ba_low)/sqrt(length(mbc_ba_low)) * tinv(0.975,length(mbc_ba_low)-1);
mbc_ca_avg_low = mean(mbc_ca_low); 
mbc_ca_moe_low = std(mbc_ca_low)/sqrt(length(mbc_ca_low)) * tinv(0.975,length(mbc_ca_low)-1);
mbc_nnn_avg_low = mean(mbc_nnn_low); 
mbc_nnn_moe_low = std(mbc_nnn_low)/sqrt(length(mbc_nnn_low)) * tinv(0.975,length(mbc_nnn_low)-1);
mbclow = [mbc_esd_avg_low mbc_esd_moe_low
            mbc_ba_avg_low  mbc_ba_moe_low
            mbc_ca_avg_low  mbc_ca_moe_low
            mbc_nnn_avg_low  mbc_nnn_moe_low];
% Avg & CI for each high MBC                
mbc_esd_avg_high = mean(mbc_esd_high); 
mbc_esd_moe_high = std(mbc_esd_high)/sqrt(length(mbc_esd_high)) * tinv(0.975,length(mbc_esd_high)-1);
mbc_ba_avg_high = mean(mbc_ba_high); 
mbc_ba_moe_high = std(mbc_ba_high)/sqrt(length(mbc_ba_high)) * tinv(0.975,length(mbc_ba_high)-1);
mbc_ca_avg_high = mean(mbc_ca_high); 
mbc_ca_moe_high = std(mbc_ca_high)/sqrt(length(mbc_ca_high)) * tinv(0.975,length(mbc_ca_high)-1);
mbc_nnn_avg_high = mean(mbc_nnn_high); 
mbc_nnn_moe_high = std(mbc_nnn_high)/sqrt(length(mbc_nnn_high)) * tinv(0.975,length(mbc_nnn_high)-1);
mbchigh = [mbc_esd_avg_high mbc_esd_moe_high
            mbc_ba_avg_high  mbc_ba_moe_high
            mbc_ca_avg_high  mbc_ca_moe_high
            mbc_nnn_avg_high  mbc_nnn_moe_high];           
% Avg & CI for each percentage error in the mean
ds_mean_error_esd_avg = mean(ds_mean_error_esd); 
ds_mean_error_esd_moe = std(ds_mean_error_esd)/sqrt(length(ds_mean_error_esd)) * tinv(0.975,length(ds_mean_error_esd)-1);
ds_mean_error_ba_avg = mean(ds_mean_error_ba); 
ds_mean_error_ba_moe = std(ds_mean_error_ba)/sqrt(length(ds_mean_error_ba)) * tinv(0.975,length(ds_mean_error_ba)-1);
ds_mean_error_ca_avg = mean(ds_mean_error_ca); 
ds_mean_error_ca_moe = std(ds_mean_error_ca)/sqrt(length(ds_mean_error_ca)) * tinv(0.975,length(ds_mean_error_ca)-1);
ds_mean_error_nnn_avg = mean(ds_mean_error_nnn); 
ds_mean_error_nnn_moe = std(ds_mean_error_nnn)/sqrt(length(ds_mean_error_nnn)) * tinv(0.975,length(ds_mean_error_nnn)-1);
ds_meanerrors = [ds_mean_error_esd_avg ds_mean_error_esd_moe
                 ds_mean_error_ba_avg  ds_mean_error_ba_moe
                 ds_mean_error_ca_avg  ds_mean_error_ca_moe
                 ds_mean_error_nnn_avg ds_mean_error_nnn_moe];
% Avg & CI for each distribution      d
ds_esd_dist_avg = mean(esd_ds_stack,2);
ds_esd_dist_moe = std(esd_ds_stack,0,2)/sqrt(min(size(esd_ds_stack))) * tinv(0.975,min(size(esd_ds_stack))-1);
ds_ba_dist_avg = mean(ba_ds_stack,2);
ds_ba_dist_moe = std(ba_ds_stack,0,2)/sqrt(min(size(ba_ds_stack))) * tinv(0.975,min(size(ba_ds_stack))-1);
ds_ca_dist_avg = mean(ca_ds_stack,2); 
ds_ca_dist_moe = std(ca_ds_stack,0,2)/sqrt(min(size(ca_ds_stack))) * tinv(0.975,min(size(ca_ds_stack))-1);
ds_nnn_dist_avg = mean(nnn_ds_stack,2); 
ds_nnn_dist_moe = std(nnn_ds_stack,0,2)/sqrt(min(size(nnn_ds_stack))) * tinv(0.975,min(size(nnn_ds_stack))-1);              
 
% Compute percent difference for each bin relative to reference
ds_esd_dist_bin_error = 100*abs((norm_esd_o - ds_esd_dist_avg));
ds_ba_dist_bin_error = 100*abs((norm_ba_o - ds_ba_dist_avg));
ds_ca_dist_bin_error = 100*abs((norm_ca_o - ds_ca_dist_avg));
ds_nnn_dist_bin_error = 100*abs((norm_nnn_o - ds_nnn_dist_avg));
 
ResultsToWrite = [mbcfull
                  1e+100*ones(1,2)  
                  mbclow
                  1e+100*ones(1,2) 
                  mbchigh
                  1e+100*ones(1,2) 
                  ds_meanerrors];
DistDataToWrite = [bins_esd' norm_esd_o ds_esd_dist_avg ds_esd_dist_moe ds_esd_dist_bin_error
                   1e+100*ones(1,5) 
                   bins_shape' norm_ba_o ds_ba_dist_avg  ds_ba_dist_moe ds_ba_dist_bin_error
                   1e+100*ones(1,5)
                   bins_shape' norm_ca_o ds_ca_dist_avg  ds_ca_dist_moe ds_ca_dist_bin_error
                   1e+100*ones(1,5)
                   bins_nnn' norm_nnn_o ds_nnn_dist_avg ds_nnn_dist_moe ds_nnn_dist_bin_error];
 
DescribeResults = input('What do you want to call these results? : ','s')
%Write MBC and Percent Error Results to .csv
filename = strcat(sprintf('%s',DescribeResults),sprintf('%s', '_MBCandPercErrMeanData'),sprintf('%s', '.csv'));
csvwrite(filename,ResultsToWrite,1,1);
%Write distribution data to .csv
filename = strcat(sprintf('%s',DescribeResults),sprintf('%s', '_DistData'),sprintf('%s', '.csv'));
csvwrite(filename,DistDataToWrite,1,1);
