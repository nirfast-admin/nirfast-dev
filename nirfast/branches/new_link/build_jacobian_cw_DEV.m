function [J] = build_jacobian_cw(mesh,data)

% J = build_jacobian(mesh,data)
%
% Builds Jacobian, both complex and in terms of log amplitude and
% phase using direct field, data.phi, and adjoint field,
% data.aphi. For structure of Jacobian see any of Dartmouth
% publications.
%
% mesh is the mesh
% data is the data
% J is the resulting Jacobian



[ncol,junk] = size(mesh.nodes);
[nrow] = length(find(mesh.link(:,3)~=0));
[nsd, msd] = size(mesh.link); 

J.complex = zeros(nrow,ncol);
J.complete = zeros(nrow,ncol);

% create a fake imaginary part here as mex files assume complex
% numbers
fake_i = ones(ncol,1).*1e-20;

k = 1;
 for i = 1 : nsd
     if mesh.link(i,3) == 1 
      if mesh.dimension == 2
	
	% Calculate the absorption part here
	J.complex(k,:) = ...
	    -IntFG(mesh.nodes(:,1:2),...
		   sort(mesh.elements')',...
		   mesh.element_area,...
		   complex(full(data.phi(:,mesh.link(i,1))),fake_i),...
		   conj(complex(full(data.aphi(:,mesh.link(i,2))),fake_i)));
	% Extract log amplitude
	J.complete(k,:) = ...
	    real(J.complex(k,:)./data.complex(k));
      elseif mesh.dimension == 3
	
	% Calculate the absorption part here
	J.complex(k,:) = ...
	    -IntFG_tet4(mesh.nodes,...
			sort(mesh.elements')',...
			mesh.element_area,...
			complex(full(data.phi(:,mesh.link(i,1))),fake_i),...
			conj(complex(full(data.aphi(:,mesh.link(i,2))),fake_i)));
	% Extract log amplitude
	J.complete(k,:) = ...
	    real(J.complex(k,:)./data.complex(k));
      end
      k = k + 1;
    end
  end
end
