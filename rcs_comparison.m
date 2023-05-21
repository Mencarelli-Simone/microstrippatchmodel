clear all 
close all 
clc
%%
t_rcs = readtable("rcs_parametric_sweep.txt");

nsize = 101;
preamble = 2;
tables = (height(t_rcs) + 2) / (nsize + preamble);
epsilon = [3, 3, 4, 4, 5, 5, 6];
patchwidth = [ 7, 8, 7, 8, 7, 8, 7];

%% array extraction
frequency = table2array(t_rcs(1:nsize,1));
rcs = zeros(nsize, tables);

for i = 1:tables
    rcs(:,i) = table2array(t_rcs(((i-1) * (nsize + preamble) + 1):(i * (nsize + preamble) - preamble),2));
end


%% example plot
figure(1)
hold on
% for i = 1:tables
%     plot(frequency,rcs(:,i))
% end

%% interpolation
% desired resolution
delta_f_min = 0.000005; %GHz
Ni = ceil(length(frequency) * (frequency(2)-frequency(1))/delta_f_min);
% new axis
f_i = linspace(frequency(1), frequency(length(frequency)), Ni);
s_i = zeros(length(f_i), tables);
leg = [];
for i = 1:tables
    s_i(:,i) = interp1(frequency, rcs(:,i), f_i, 'cubic' );
    leg = [leg, 'er=' + string(epsilon(i))+ ', L=' +  string(patchwidth(i))];
end
plot(f_i,s_i)
legend(leg)
lgd = legend;
lgd.NumColumns = 2;
%% derivatives
delta = f_i(2) - f_i(1);
ds = diff(s_i) / delta;
ff = f_i + delta/2;
ff = ff(2:length(ff));
figure(2)
plot(ff, abs(ds))
hold on
% dds = diff(ds) / delta;
% fff = ff + delta/2;
% fff = fff(2:length(fff));
% plot(fff, dds)
%plot(ff, abs(ds)-ds)

%% Peakfinder
threshold = 0.0001;
indexes = {};
for i = 1:tables
    indexes = [indexes; find(abs(ds(:,i)) < threshold)];
    if isempty(find(abs(ds(:,i)) < threshold))
        indexes = [indexes; 0];
    end
end

% %% points visualization
% figure(1)
% for i = 1:tables
%     for j = indexes{i}
%         if j~=0
%             scatter(ff(j), (s_i(j,i) + s_i(j+1,i)) / 2,'+k') % point
%             text(ff(j) + 0.001, s_i(j,i) + 0.001, string(ff(j)))
%         end
%     end
% end
% legend off
% legend(leg)

%% line fit
slope = mean(ds , 1);
% offset fit
offset=[];
fitlines = {};
for i = 1:tables
    of = @(q) mean((s_i(:,i)' - (f_i * slope(i) + q)).^2);
    x_o = -10;
    offset = [offset,fminsearch(of, x_o)]; 
end
%% lines
figure(3)
hold on
f_i_m = ones(tables, 1) .* f_i;
lreg = slope' .* f_i_m  + offset';
for  i = 1:tables
    plot(f_i_m(i,:), s_i(:,i)' - lreg(i,:))
    %plot(f_i_m(i,:), lreg(i,:))
end

%% the peak should be the minimum of the new signal
for i = 1:tables
    [c, ii] = min(s_i(:,i)' - lreg(i,:));
    if i == 6
        [c, ii] = min(s_i(1:200000,i)' - lreg(i,1:200000));% refined to exclde secondary peaks
    end
    plot(f_i_m(i,ii), s_i(ii,i)' - lreg(i,ii), 'ok')
    text(f_i_m(i,ii) + 0.001, s_i(ii,i)' - lreg(i,ii) + 0.001, string(i)+'= '+string(f_i_m(i,ii)))
end

legend(leg)





