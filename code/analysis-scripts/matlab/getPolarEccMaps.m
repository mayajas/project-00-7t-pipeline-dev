function [Srf,Idx] = getPolarEccMaps(Srf,Thrsh,Type)

%% Remove rubbish
if strcmpi(Srf.Values{1}, 'R^2') || strcmpi(Srf.Values{1}, 'nR^2')
    r = Srf.Data(1,:) <= Thrsh(1) | isnan(Srf.Data(1,:));
    Srf.Data(1,r) = NaN; % Set rubbish to NaN
else 
    r = isnan(Srf.Data(1,:));
end


%% Calculate colours for the map 
if strcmpi(Type, 'Polar') 
    % Polar map
    Pha = atan2(Srf.Data(3,:), Srf.Data(2,:)) / pi * 180;
    Data = Pha;
    if Srf.Hemisphere(1) == 'l'
        % Left hemisphere
        Pha = mod(ceil(Pha + 270), 360) + 1;
    elseif Srf.Hemisphere(1) == 'b'
        % Bilateral maps
        Pha(1:Srf.Nvert_Lhem) = mod(ceil(Pha(1:Srf.Nvert_Lhem) + 270), 360) + 1;
        Pha(Srf.Nvert_Lhem+1:end) = -Pha(Srf.Nvert_Lhem+1:end);
        Pha(Srf.Nvert_Lhem+1:end) = mod(ceil(Pha(Srf.Nvert_Lhem+1:end) + 90), 360) + 1;
    else
        % Right hemisphere
        Pha = -Pha;
        Pha = mod(ceil(Pha + 90), 360) + 1;
    end
    Pha(Pha == 0) = 360;
    Pha(r) = 360;
    
    Srf.Data(end+1,:) = Pha;
    Srf.Values{end+1} = Type;
    Idx = size(Srf.Data,1);
%     % Colourmap
%     cstr = ['colormap(' def_cmap_angle(2:end) '(360));'];
%     Cmap = eval(cstr);        
%     if def_cmap_angle(1) == '-'
%         Cmap = flipud(Cmap);
%     end
%     
%     % Determine colours
%     Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
%     if isnan(PathColour) 
%         PathColour = [1 1 1];
%     end
    
elseif strcmpi(Type, 'Eccentricity')
    % Eccentricity map
    Rho = sqrt(Srf.Data(2,:).^2 + Srf.Data(3,:).^2);
    Data = Rho;
    Rho(r) = 0;
 
%     % If logarithmic maps
%     if def_logmaps
%     	Rho = log2(Rho);
%     end	
        
    % Set all below minimum to minimum
    Rho(Rho < Thrsh(2)) = Thrsh(2);
    % Adjust minimum
    AdjThr = Thrsh(3) - Thrsh(2);  
    Rho = Rho - Thrsh(2);
    Rho(Rho < 0) = 0;
    % Set all above maximum to maximum
    Rho = Rho / AdjThr;
    Rho(Rho > 1) = 1;   
    % Convert to integers
    Pha = round(Rho * 360);

    Srf.Data(end+1,:) = Pha;
    Srf.Values{end+1} = Type;
    Idx = size(Srf.Data,1);
    
    
%     % Colourmap
%     cstr = ['[colormap(' def_cmap_eccen(2:end) '(360)); CurvGrey];'];
%     Cmap = eval(cstr);        
%     if def_cmap_eccen(1) == '-'
%         Cmap = flipud(Cmap);
%     end
%     
%     % Determine colours
%     Pha = mod(Pha, 360);
%     Pha(Pha==0) = 360;
%     Pha(r|isnan(Pha)) = 360;
%     Colours = Cmap(Pha,:).*Alpha + CurvGrey(Curv,:).*(1-Alpha); % Colours transparently overlaid onto curvature
%     if isnan(PathColour) 
%         PathColour = [0 0 0];
%     end
end