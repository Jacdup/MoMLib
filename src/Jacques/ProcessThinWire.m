numfiles = 9;
% FF_MBF{1} = importdata('FF_MBF_391_1.mat');
for k = 1:12
    filename = sprintf('FF_MBF_103_%d.mat',k);
    FF_MBF{k} = importdata(filename);
end

% Import FEKO RWG
for k = 1:12
%     filename = sprintf('FF_RWG_714_%d.mat',k);
    filename = sprintf('C:\\Users\\19083688\\Desktop\\Masters\\MoMLib\\Results\\Thin Wire vs MBF\\FEKO\\XY\\FEKO_RWG_%d.txt', k);
    FF_RWG{k} = feko_farfield_extract(filename);
    FF_RWG{k} = FF_RWG{k}(1:180,:);
%     FF_RWG{k} = importdata(filename);
end

% Import Thin wire
for k = 1:numfiles
    numbers = {'0.0001','0.0004', '0.0008', '0.001', '0.004', '0.008','0.01', '0.04', '0.08', '0.1','0.14', '0.18'};
    filename = sprintf('C:\\Users\\19083688\\Desktop\\Masters\\MoMLib\\Results\\Thin Wire vs MBF\\FEKO\\XY\\FEKO_%d_%s.txt', k, char(numbers{k}));
    FF_FEKO{k} = feko_farfield_extract(filename);
    FF_FEKO{k} = FF_FEKO{k}(1:180,:);
end

% Calculate RCS & Errors
for k = 9:12
%     avgMBF(k)  = sum(FF_MBF{k}(:,3),1)/180;
%      avgRWG(k)  = sum(FF_RWG{k}(:,3),1)/180;
%       avgFEKO(k)  = sum(FF_FEKO{k}(:,3),1)/180;
    
%     rcsMBF(k) = (4*pi*((abs(FF_MBF{k}(1,3))^2)/(1)));
%     rcsThin(k) = (4*pi*((abs(FF_FEKO{k}(46,3))^2)/(1)));
    rcsRWG(k)  = (4*pi*((abs(FF_RWG{k}(91,3))^2)/(1)));
%      [errAbsFF(k), errRelFF(k)] = pNormError(rcsMBF, rcsRWG, 2);
%      [errAbsFFThin(k), errRelFFThin(k)] = pNormError(rcsThin, rcsRWG, inf);
    
%     [errAbsFF(k), errRelFF(k)] = pNormError(sqrt(FF_MBF{k}(:,3).^2 + FF_MBF{k}(:,5).^2), sqrt(FF_RWG{k}(:,3).^2 + FF_RWG{k}(:,5).^2), 2);
%      [errAbsFFThin(k), errRelFFThin(k)] = pNormError(sqrt(FF_FEKO{k}(:,3).^2 + FF_FEKO{k}(:,5).^2), sqrt(FF_RWG{k}(:,3).^2 + FF_RWG{k}(:,5).^2), 2);
end

% Plot farfields
for k = 1:12
        figure
     plot(FF_MBF{k}(:,1),sqrt(FF_MBF{k}(:,3).^2 + FF_MBF{k}(:,5).^2));
     hold on
     plot(FF_RWG{k}(:,1)*(pi/180),sqrt(FF_RWG{k}(:,3).^2 + FF_RWG{k}(:,5).^2));
%        hold on
%      plot(FF_FEKO{k}(:,1)*(pi/180),sqrt(FF_FEKO{k}(:,3).^2 + FF_FEKO{k}(:,5).^2));
end
 temp =[ 0.0001,0.0004, 0.0008, 0.001, 0.004, 0.008, 0.01, 0.04, 0.08, 0.1, 0.14, 0.2];
% Plot errors
    figure
    semilogx(temp(1:9), errRelFF);
    hold on
    semilogx(temp(1:9), errRelFFThin);
    
    % Plot RCS
    figure
        semilogx(temp(1:12), rcsMBF(1:12));
    hold on
     semilogx(temp(1:9), rcsThin(1:9));
 semilogx(temp(1:12), rcsRWG(1:12));
    legend("MBF","Thin Wire", "RWG");
%         loglog(temp(1:12), rcsMBF(1:12));
%     hold on
%     loglog(temp(1:9), rcsThin(1:9));
%     loglog(temp(1:12), rcsRWG(1:12));