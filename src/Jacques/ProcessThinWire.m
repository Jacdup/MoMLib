numfiles = 9;
% FF_MBF{1} = importdata('FF_MBF_391_1.mat');
for k = 1:12
    filename = sprintf('FF_MBF_103_%d.mat',k);
    FF_MBF{k} = importdata(filename);
end

% Import FEKO RWG
for k = 1:18
%     filename = sprintf('FF_RWG_714_%d.mat',k);
    filename = sprintf('C:\\Users\\19083688\\Desktop\\Masters\\MoMLib\\Results\\Thin Wire vs MBF\\FEKO\\XY\\Theta45\\FEKO_RWG_%d.ffe', k);
    fileID = fopen(filename,'r');
    % formatSpec = '%f %f %f %f %f %f %*f %*f %*f %*f %*f %*s';
    % formatSpec = '%f %f %f %f %f %f %*f %*f %*f %*f %*f %*s';
    formatSpec = '%f %f %f %f %f %f %f %f %f';
    sizeA = [9 inf];
    FF_RWG{k} = textscan(fileID,formatSpec, 'HeaderLines',16 );
    FF_RWG{k} = cell2mat(FF_RWG{k});
    fclose(fileID);
    rcsRWG(k,:) = FF_RWG{k}(1:180,8);
%     FF_RWG{k} = feko_farfield_extract(filename);
%     FF_RWG{k} = FF_RWG{k}(1:180,:);

%     FF_RWG{k} = importdata(filename);
end

% Import Thin wire
for k = 1:numfiles
    numbers = {'0.0001','0.0004', '0.0008', '0.001', '0.004', '0.008','0.01', '0.04', '0.08', '0.1','0.14', '0.18'};
%     filename = sprintf('C:\\Users\\19083688\\Desktop\\Masters\\MoMLib\\Results\\Thin Wire vs MBF\\FEKO\\XY\\FEKO_%d_%s.txt', k, char(numbers{k}));
filename = sprintf('C:\\Users\\19083688\\Desktop\\Masters\\MoMLib\\Results\\Thin Wire vs MBF\\FEKO\\XY\\Theta45\\FEKO_Thin_%d.ffe', k);   
 fileID = fopen(filename,'r');
formatSpec = '%f %f %f %f %f %f %f %f %f';
    sizeA = [9 inf];
    FF_Thin{k} = textscan(fileID,formatSpec, 'HeaderLines',16 );
    FF_Thin{k} = cell2mat(FF_Thin{k});
    fclose(fileID);

% FF_FEKO{k} = feko_farfield_extract(filename);
%     FF_FEKO{k} = FF_FEKO{k}(1:180,:);
end

% Calculate RCS & Errors
for k = 1:18
    
    rcsRWG(k,:) = FF_RWG{k}(1:180,8);

%     rcsThin(k,:) = FF_Thin{k}(1:180,8);
    
%     rcsMBF(k) = (4*pi*((abs(sum(FF_MBF{k}(:,5)))^2)/(180)));
%     rcsMBF(k) = (4*pi*((abs(FF_MBF{k}(46,5))^2)/(1)));
%     rcsThin(k) = (4*pi*((abs(FF_FEKO{k}(22,3))^2)/(1)));% 45 Degrees sidelobe
%     rcsThin(k) = (4*pi*((abs(FF_FEKO{k}(46,3))^2)/(1))); % Mainlobe
%     rcsThin(k) = (4*pi*((abs(sum(FF_FEKO{k}(:,3)))^2)/(180)));
%     rcsRWG(k)  = (4*pi*((abs(FF_RWG{k}(66,3))^2)/(1))); % 45 Degrees sidelobe
%      rcsRWG(k)  = (4*pi*((abs(FF_RWG{k}(91,3))^2)/(1))); % 180 Degrees mainlobe
%      rcsRWG(k)  = (4*pi*((abs(sum(FF_RWG{k}(:,3)))^2)/(180))); 
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
%  temp =[ 0.0001,0.0004, 0.0008, 0.001, 0.004, 0.008, 0.01, 0.04, 0.08, 0.1, 0.14, 0.2];
 temp =[ 0.0001,0.0004, 0.0008, 0.001, 0.004, 0.008, 0.01, 0.04, 0.08, 0.1, 0.14, 0.18, 0.2, 0.24, 0.28, 0.3, 0.34, 0.38];
% Plot errors
    figure
    semilogx(temp(1:9), errRelFF);
    hold on
    semilogx(temp(1:9), errRelFFThin);
    
    % Plot RCS as surf
    lim = 1:18;
    figure
    surf(1:2:360,temp(lim),abs(rcsMBF(lim,:) -rcsRWG(lim,:)));
        figure
    surf(1:180,temp(lim),rcsRWG(lim,:));
    surf(1:180,temp(1:9),rcsThin(1:9,:));
    
    % Plot single RCS
    figure
    semilogx(temp(lim), rcsMBF(lim,70));
    hold on
%      loglog(temp(4:9), rcsThin(4:9));
    semilogx(temp(lim), rcsRWG(lim,70));
     semilogx(temp(1:9), rcsThin(1:9,70));
    legend("MBF", "RWG","Thin Wire");
    ylabel('$\sigma/\lambda^2$');
    xlabel('$r/\lambda$');
   
%     hold on
%     loglog(temp(1:9), rcsThin(1:9));
%     loglog(temp(1:12), rcsRWG(1:12));