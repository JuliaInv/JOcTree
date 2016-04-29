function plotOcTreeMesh(f, i, j, k, bsz, n, h, x0, cdata)

% Julia provides indices as int64 but Matlab requires double
f     = double(f);
i     = double(i);
j     = double(j);
k     = double(k);
bsz   = double(bsz);
n     = double(n);
cdata = double(cdata);

nx  = n(1);
ny  = n(2);
nz  = n(3);
hx  = h(1);
hy  = h(2);
hz  = h(3);
z0  = x0(3);
y0  = x0(2);
x0  = x0(1);

clim = sort(cdata(isfinite(cdata)));
clim = clim([1 end]);
c    = mean(clim);
d    = diff(clim);
if d < eps(single(c))
  if abs(c) < realmin('single')
    clim = [-1 1];
  else
    clim = c + eps(single(c)) * [-1 1];
  end
end
cmap = jet(320);
cmap = cmap(33:288, :);
clear c d

% variables accessed by draw and callback functions
vertices  = [];
faces     = [];
slices    = [];
iSlice    = [];
nSlices   = [];
edgeColor = [0.5 0.5 0.5];
hPatch    = [];
ratio     = [];
tfun      = [];
iPlane    = 0;

% plot
figure(f)
clf
axis equal
caxis(clim);
colormap(cmap);
box on
set(gca, ...
  'Layer', 'top');
colorbar;
% We have to deactivate interactive plotting before changing the scroll
% callback.
plotedit off
zoom off
pan off
rotate3d off
datacursormode off
brush off
set(gcf, ...
  'KeyPressFcn', @keyFcn, ...
  'WindowScrollWheelFcn', @scrollFcn, ...
  'Toolbar', 'figure');

% checkbox show mesh
uicontrol( ...
  'Style', 'checkbox', ...
  'Units', 'normalized', ...
  'Position', [0.02 0.10 0.24 0.04], ...
  'Value', true, ...
  'HorizontalAlignment', 'left', ...
  'String', 'Show mesh', ...
  'TooltipString', 'Check to show mesh', ...
  'Callback', @meshFcn, ...
  'Tag', 'Mesh');

% checkbox aspect ratio
uicontrol( ...
  'Style', 'checkbox', ...
  'Units', 'normalized', ...
  'Position', [0.02 0.06 0.24 0.04], ...
  'Value', true, ...
  'HorizontalAlignment', 'left', ...
  'String', 'Axis equal', ...
  'TooltipString', 'Check for axis equal', ...
  'Callback', @axisFcn, ...
  'Tag', 'AxisEqual');

% checkbox vertical axis direction up
uicontrol( ...
  'Style', 'checkbox', ...
  'Units', 'normalized', ...
  'Position', [0.02 0.02 0.24 0.04], ...
  'Value', false, ...
  'HorizontalAlignment', 'left', ...
  'String', 'Reverse y-axis', ...
  'TooltipString', 'Check to reverse y-axis', ...
  'Callback', @yAxisFcn, ...
  'Tag', 'YDir');

% popup for x|y|z plane
h = uicontrol( ...
  'Style', 'popup', ...
  'Units', 'normalized', ...
  'Position', [0.86 0.02 0.12 0.04], ...
  'HorizontalAlignment', 'left', ...
  'String', {'x','y','z'}, ...
  'TooltipString', 'Select slice direction', ...
  'Value', 3, ...
  'Callback', @planeFcn, ...
  'Tag','Plane');

% finally, trigger planeFcn to draw
feval(get(h, 'Callback'), h);

% nested callback functions
  function draw

    cla
    hPatch = patch( ...
      'Vertices', vertices, ...
      'Faces', faces(slices{iSlice}, :), ...
      'FaceVertexCData', cdata(slices{iSlice}), ...
      'FaceColor', 'flat', ...
      'EdgeColor', edgeColor, ...
      'ButtonDownFcn', @resetFcn);
    title(tfun(iSlice));
    drawnow
    
  end

  function keyFcn(hObject, hEvent)
    
    if hObject ~= gcf
      figure(hObject);
    end

    switch hEvent.Key
        case 'uparrow'
        iSlice = min(iSlice + 1, nSlices);
      case 'downarrow'
        iSlice = max(1, iSlice - 1);
      otherwise
        return
    end

    draw;

  end

  function scrollFcn(hObject, hEvent)
    
    if hObject ~= gcf
      figure(hObject);
    end

    inc = hEvent.VerticalScrollCount;
    iSlice = min(max(1, iSlice + inc), nSlices);
    draw;

  end

  function meshFcn(hObject, ~, ~)
    
    if get(hObject, 'Value') == get(hObject, 'Max')
      edgeColor = [0.5 0.5 0.5];
    else
      edgeColor = 'none';
    end
    set(hPatch, 'EdgeColor', edgeColor);
    
  end

  function axisFcn(hObject, ~, ~)
    
    if get(hObject, 'Value') == get(hObject, 'Max')
      set(gca, 'DataAspectRatio', [1 1 1]);
    else
      set(gca, 'DataAspectRatio', [ratio 1]);
    end
    
    
  end

  function resetFcn(hObject, ~)
    
    % hObj is handle to patch. 'SelectionType' is defined for figure that
    % is two levels up from patch.
    type = get(get(get(hObject, 'Parent'), 'Parent'), 'SelectionType');
    
    if strcmp(type, 'open')
      
      iSlice = floor((nSlices + 1) / 2);
      draw;
      
    end
    
  end

  function planeFcn(hObject, ~)
    
    jPlane = get(hObject, 'Value');
    if iPlane == jPlane
      % same plane, nothing to do
      return
    end
    
    % get current axis limits
    alim = axis;
    switch iPlane
      case 0 % initialization
        xmin = x0;
        xmax = x0 + nx * hx;
        ymin = y0;
        ymax = y0 + ny * hy;
        zmin = z0;
        zmax = z0 + nz * hz;
      case 1
        ymin = alim(1);
        ymax = alim(2);
        zmin = alim(3);
        zmax = alim(4);
        dy   = ymax - ymin;
        dz   = zmax - zmin;
        dx   = max(dy, dz);
        xc   = x0 + (iSlice - 0.5) * hx;
        xmin = xc - 0.5 * dx;
        xmax = xc + 0.5 * dx;
      case 2
        xmin = alim(1);
        xmax = alim(2);
        zmin = alim(3);
        zmax = alim(4);
        dx   = xmax - xmin;
        dz   = zmax - zmin;
        dy   = max(dx, dz);
        yc   = y0 + (iSlice - 0.5) * hy;
        ymin = yc - 0.5 * dy;
        ymax = yc + 0.5 * dy;
      case 3
        xmin = alim(1);
        xmax = alim(2);
        ymin = alim(3);
        ymax = alim(4);
        dx   = xmax - xmin;
        dy   = ymax - ymin;
        dz   = max(dx, dy);
        zc   = z0 + (iSlice - 0.5) * hz;
        zmin = zc - 0.5 * dz;
        zmax = zc + 0.5 * dz;
    end
    
    iPlane = jPlane;
    
    hYDir = findobj(gcf, 'Type', 'uicontrol', 'Tag', 'YDir');

    switch iPlane

      case 1

        vertices = [
          j     k
          j+bsz k
          j+bsz k+bsz
          j     k+bsz];
        [vertices, ~, faces] = unique(vertices, 'rows');
        faces = reshape(faces, [], 4);

        vertices(:,1) = y0 + (vertices(:,1) - 1) * hy;
        vertices(:,2) = z0 + (vertices(:,2) - 1) * hz;

        xlim([ymin ymax]);
        ylim([zmin zmax]);
        xlabel('y in m');
        ylabel('z in m');
        tfun  = @(i)(sprintf('x = %gm', x0 + (i - 0.5) * hx));
        % update checkbox for vertical axis
        set(hYDir, ...
          'String', 'Reverse z-axis', ...
          'TooltipString', 'Check to reverse z-axis');
        ratio = [hy hz];

        m = i;
        nSlices = nx;
        iSlice  = round((0.5 * (xmin + xmax) - x0) / hx + 0.5);
        iSlice  = min(max(1, iSlice), nSlices);

      case 2

        vertices = [
          i     k
          i+bsz k
          i+bsz k+bsz
          i     k+bsz];
        [vertices, ~, faces] = unique(vertices, 'rows');
        faces = reshape(faces, [], 4);

        vertices(:,1) = x0 + (vertices(:,1) - 1) * hx;
        vertices(:,2) = z0 + (vertices(:,2) - 1) * hz;

        xlim([xmin xmax]);
        ylim([zmin zmax]);
        xlabel('x in m');
        ylabel('z in m');
        tfun  = @(j)(sprintf('y = %gm', y0 + (j - 0.5) * hy));
        % update checkbox for vertical axis
        set(hYDir, ...
          'String', 'Reverse z-axis', ...
          'TooltipString', 'Check to reverse z-axis');
        ratio = [hx hz];

        m = j;
        nSlices = ny;
        iSlice  = round((0.5 * (ymin + ymax) - y0) / hy + 0.5);
        iSlice  = min(max(1, iSlice), nSlices);

      case 3

        vertices = [
          i     j
          i+bsz j
          i+bsz j+bsz
          i     j+bsz];
        [vertices, ~, faces] = unique(vertices, 'rows');
        faces = reshape(faces, [], 4);

        vertices(:,1) = x0 + (vertices(:,1) - 1) * hx;
        vertices(:,2) = y0 + (vertices(:,2) - 1) * hy;

        xlim([xmin xmax]);
        ylim([ymin ymax]);
        xlabel('x in m');
        ylabel('y in m');
        tfun  = @(k)(sprintf('z = %gm', z0 + (k - 0.5) * hz));
        % update checkbox for vertical axis
        set(hYDir, ...
          'String', 'Reverse y-axis', ...
          'TooltipString', 'Check to reverse y-axis');
        ratio = [hx hy];

        m = k;
        nSlices = nz;
        iSlice  = round((0.5 * (zmin + zmax) - z0) / hz + 0.5);
        iSlice  = min(max(1, iSlice), nSlices);

    end
    
    bval = reshape(unique(bsz), 1, []);
    nbsz = length(bval);

    slices = cell(nSlices, 1);

    for ibsz = 1:nbsz

      b = bval(ibsz);
      M = find(bsz == b);
      N = repmat(m(M), 1, b) + repmat(0:b-1, length(M), 1);
      M = repmat(M, 1, b);

      s = accumarray(N(:), M(:), [nSlices, 1], @(x)({x}));
      slices = cellfun(@(x,y)([x;y]), slices, s, 'UniformOutput', false);

    end
    
    draw;
    
  end

  function yAxisFcn(hObject, ~, ~)

    if get(hObject, 'Value') == get(hObject, 'Max')
      set(gca, 'YDir', 'reverse');
    else
      set(gca, 'YDir', 'normal');
    end

  end

end