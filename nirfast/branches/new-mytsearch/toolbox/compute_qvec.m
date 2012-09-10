function [qvec mesh] = compute_qvec(mesh,s_ind,i)

ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
nsource = length(source);
nnodes=size(mesh.nodes,1);

if nargin >= 3 && ~isempty(s_ind) && ~isempty(i)
    if mesh.source.fwhm(s_ind) == 0
        if ~isfield(mesh.source,'qvec') ...
                || length(mesh.source.qvec) < i || ...
                isempty(mesh.source.qvec{i}) || ...
                size(mesh.source.qvec{i},1) ~= nnodes
            qvec = ...
                gen_source_point(mesh,mesh.source.coord(s_ind,1:3));
            mesh.source.qvec{i} = qvec;
        else
            qvec = mesh.source.qvec{i};
        end
    else
        error(' Report a bug')
    end
else
    qvec = spalloc(nnodes,nsource,nsource*100);
    for i = 1 : nsource
        s_ind = mesh.source.num == source(i);
        if mesh.source.fwhm(s_ind) == 0
            if isfield(mesh,'force_qvec_compute') && ...
                mesh.force_qvec_compute ~= 0
                qvec(:,i) = ...
                    gen_source_point(mesh,mesh.source.coord(s_ind,1:3));
                mesh.source.qvec{i} = qvec(:,i);
            elseif ~isfield(mesh.source,'qvec') ...
                    || length(mesh.source.qvec) < i || ...
                    isempty(mesh.source.qvec{i}) || ...
                    size(mesh.source.qvec{i},1) ~= nnodes;
                qvec(:,i) = ...
                    gen_source_point(mesh,mesh.source.coord(s_ind,1:3));
                mesh.source.qvec{i} = qvec(:,i);
            end
    %     else
    %         qvec(:,i) = gen_source(mesh.nodes,...
    %             sort(mesh.elements')',...
    %             mesh.dimension,...
    %             mesh.source.coord(s_ind,:),...
    %             mesh.source.fwhm(s_ind));
        end
    end
end