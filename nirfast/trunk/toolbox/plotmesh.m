function plotmesh(mesh, plotflag)

% plotmesh(mesh, plotflag)
%
% Allows fast and easy viewing of mesh
% 
% mesh is the input mesh (variable or filename)
% plotflag is optional, if it is 1, the source/detectors will show




% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end


%****************************************
% Plotting
%****************************************


if nargin == 1
   plotflag = 0;
end
if plotflag == 1
    figure;
    hold on;
    ind = find(mesh.bndvtx==1);
    if mesh.dimension == 2
        if isfield(mesh,'source') == 1
            plot(mesh.source.coord(:,1),...
                mesh.source.coord(:,2),'ro','LineWidth',2,'MarkerSize',8);


        end
        if isfield(mesh,'meas') == 1
            plot(mesh.meas.coord(:,1),...
                mesh.meas.coord(:,2),'bx','LineWidth',2,'MarkerSize',8);
        end
        plot(mesh.nodes(ind,1),mesh.nodes(ind,2),'c.');
        axis equal;
        legend('Source','Detector');
    elseif mesh.dimension == 3
        if isfield(mesh,'source') == 1
            plot3(mesh.source.coord(:,1),...
                mesh.source.coord(:,2),...
                mesh.source.coord(:,3),'ro',...
                'LineWidth',2,'MarkerSize',8);
        end
        if isfield(mesh,'meas') == 1
            plot3(mesh.meas.coord(:,1),...
                mesh.meas.coord(:,2),...
                mesh.meas.coord(:,3),'bx',...
                'LineWidth',2,'MarkerSize',8);
        end
        plot3(mesh.nodes(ind,1),...
            mesh.nodes(ind,2),...
            mesh.nodes(ind,3),'c.');
        axis equal;
        legend('Source','Detector');
    end
end

figure;
set(gca,'FontSize',28)

% STANDARD
if strcmp(mesh.type,'stnd') == 1
  subplot(1,2,1);
  plotim(mesh,mesh.mua);
  title('\mu_a','FontSize',20);
  colorbar('horiz');
  
  subplot(1,2,2);
  plotim(mesh,mesh.mus);
  title('\mu_s''','FontSize',20);
  colorbar('horiz');
  
% STANDARD SPN
elseif strcmp(mesh.type,'stnd_spn') == 1
  subplot(2,2,1);
  plotim(mesh,mesh.mua);
  title('\mu_a','FontSize',20);
  colorbar('horiz');
  
  subplot(2,2,2);
  plotim(mesh,mesh.mus);
  title('\mu_s''','FontSize',20);
  colorbar('horiz');
  
  subplot(2,2,3);
  plotim(mesh,mesh.g);
  title('g','FontSize',20);
  colorbar('horiz');
  
% STANDARD BEM
elseif strcmp(mesh.type,'stnd_bem') == 1
  subplot(1,2,1);
  val = zeros(size(mesh.nodes,1),1);
  for i=1:size(mesh.mua,1)
      if i==1
          ind = find(mesh.region(:,2)==0);
          ind = unique(mesh.elements(ind,:));
        val(ind) = mesh.mua(i);
      else
          ind = find(mesh.region(:,2)==i);
          ind = unique(mesh.elements(ind,:));
        val(ind) = mesh.mua(i);
      end
  end
  plotim(mesh,val);
  title('\mu_a','FontSize',20);
  colorbar('horiz');
  
  subplot(1,2,2);
  val = zeros(size(mesh.nodes,1),1);
  for i=1:size(mesh.mus,1)
      if i==1
          ind = find(mesh.region(:,2)==0);
          ind = unique(mesh.elements(ind,:));
        val(ind) = mesh.mus(i);
      else
          ind = find(mesh.region(:,2)==i);
          ind = unique(mesh.elements(ind,:));
        val(ind) = mesh.mus(i);
      end
  end
  plotim(mesh,val);
  title('\mu_s''','FontSize',20);
  colorbar('horiz');

% FLUORESCENCE
elseif strcmp(mesh.type,'fluor') == 1
  subplot(3,2,1);
  plotim(mesh,mesh.muax);
  title('\mu_{ax}','FontSize',10);
  colorbar;

  subplot(3,2,2);
  plotim(mesh,mesh.musx);
  title('\mu_{sx}''','FontSize',10);
  colorbar;

  subplot(3,2,3);
  plotim(mesh,mesh.muam);
  title('\mu_{am}','FontSize',10);
  colorbar;

  subplot(3,2,4);
  plotim(mesh,mesh.musm);
  title('\mu_{sm}''','FontSize',10);
  colorbar;
  
  subplot(3,2,5);
  if isfield(mesh,'etamuaf') == 1
      plotim(mesh,mesh.etamuaf);
  else
      plotim(mesh,mesh.muaf.*mesh.eta);
  end
  title('\eta\mu_{fl}','FontSize',10);
  colorbar;

  subplot(3,2,6);
  plotim(mesh,mesh.tau);
  title('\tau','FontSize',10);
  colorbar;

% SPECTRAL
elseif strcmp(mesh.type,'spec') == 1
  [nc,junk]=size(mesh.chromscattlist);
  
  if isfield(mesh,'etamuaf')
    n = ceil((nc-2)/2)+2;
  else
    n = ceil((nc-2)/2)+1;
  end
  k = 0;
  for i = 1 : nc-2
    k = k + 1;
    subplot(n,2,k);
    plotim(mesh,mesh.conc(:,i));
    t = char(mesh.chromscattlist(i,1));
    title(t,'FontSize',10);
    colorbar;
  end
  subplot(n,2,k+1);
  plotim(mesh,mesh.sa);
  title('Scatter Amplitude','FontSize',10);
  colorbar;
  subplot(n,2,k+2);
  plotim(mesh,mesh.sp);
  title('Scatter Power','FontSize',10);
  colorbar;
  if isfield(mesh,'etamuaf')
    subplot(n,2,k+3);
    plotim(mesh,mesh.etamuaf);
    title('etamuaf','FontSize',10);
    colorbar;
  end
else
    errordlg('Mesh type not supported','NIRFAST Error');
    error('Mesh type not supported');
end

function plotim(mesh,val)
if mesh.dimension == 3 && strcmp(mesh.type,'stnd_bem')
    h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val,'FaceAlpha',0.5);
else
    h = trisurf(mesh.elements,...
	    mesh.nodes(:,1),...
	    mesh.nodes(:,2),...
	    mesh.nodes(:,3),...
	    val);
end
shading interp;
view(2);
axis equal; 
axis off;
colormap hot;
