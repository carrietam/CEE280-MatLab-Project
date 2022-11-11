classdef AJCT_Analysis < handle
% Analysis class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        all_ele
        all_nodes
        nnodes
        coord
        concen
        nele
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = AJCT_Analysis(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Azz, Iyy, Izz, J, E, v, webdir, w, truss)
            self.nnodes = nnodes;
            self.coord = coord;
            self.concen = concen;
            self.nele = nele;
            self.CreateNodes(nnodes,coord);
            self.CreateElements(A, Ayy, Azz, E, ends, Izz, Iyy, J, nele, v, webdir, w);
        end
        
        %% Run the analysis
        function RunAnalysis(self)
            self.CreateStiffnessMatrix();
            self.CreateLoadVectors();
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Create Nodes
        % Create a vector of node objects, size nnodes
        function CreateNodes(self,nnodes,coord)
            % Loop over number of nodes
            for i = 1:nnodes
                % Create node & add to all_nodes vector
                node_i = AJCT_Node(coord(i,1:3),i);
                all_nodes(i) = node_i;
            end
            self.all_nodes = all_nodes;
        end
        
        %% Create Elements
        % Create a vector of element objects, size nele
        function CreateElements(self,A, Ayy, Azz, E, ends, Izz, Iyy, J, nele, v, webdir, w)
            % Loop over number of elements
            for i = 1:nele
                % Create element & add to all_ele vector
                element_i = AJCT_Element(A(i), Ayy(i), Azz(i), E(i),...
                    [self.all_nodes(ends(1)); self.all_nodes(ends(2))], Izz(i), Iyy(i), J(i), v(i),...
                    webdir(i), w(i,1:3));
                all_ele(i) = element_i;
            end
            self.all_ele = all_ele;
        end
        
        %% Create Stiffness Matrix
        % Assemble system stiffness matrix
        function CreateStiffnessMatrix(self)
            nDOF = self.nnodes*6; % Calculate the number of degrees of freedom

            for i = 1:self.nele
                % Get necessary element properties
                ele_dof = self.all_ele(i).GetDOF(); % relevant DOFs
                K_ele = self.all_ele(i).GetKGloEleMatrix(); % global element...
                % stiffness matrix
                
                % Create necessary vectors to use sparse fnct
                % e1 repeats ele_dof over and over (e.g. 
                % [1,2,3,1,2,3,1,2,3...])
                e1 = ele_dof;
                for j = 1:length(ele_dof)-1
                    e1 = [e1,ele_dof];
                end
                
                % e2 repeats every element in ele_dof (e.g.
                % [1,1,1,2,2,2,3,3,3...])
                e2 = repelem(ele_dof,length(ele_dof));
                
                % Convert element global matrix as row vector
                K_ele = reshape(K_ele,1,[]);
                
                % Create mapped element stiffness matrix using sparse fnct
                K_mapped = sparse(e1,e2,K_ele,nDOF,nDOF);
                
                % Add mapped element matrix to existing global system matrix
                % Try to add to system matrix. If system matrix has not
                % been made yet, system matrix = 1st element matrix
                try
                    K_g = K_g + K_mapped;
                catch
                    K_g = K_mapped;
                end
            end
            disp('This is the global stiffness matrix:')
            disp(K_g)
        end
        
        function CreateLoadVectors(self)
           nDOF = self.nnodes*6; % Calculate the number of degrees of freedom
           concen_g = zeros(nDOF,1); % Initialize concen_g
           
           % Create concentrated load vector
           for i = 1:self.nnodes
                % Get necessary element properties
                node_DOF = self.all_nodes(i).GetNodeDOF();
                
                % Map out element concen vector onto global concen vector
                concen_g(node_DOF:node_DOF+5,1) = [self.concen(i,1:6)]';
           end
           
           % Create fixed end force vector
           nDOF = self.nnodes*6; % Calculate the number of degrees of freedom
           FEF_sys = zeros(nDOF,1); % Blank vector of FEF
           
           % Loop over number of elements
           for j = 1:self.nele
               % Get necessary element properties
               ele_dof = self.all_ele(j).GetDOF(); % relevant DOFs
               ele_FEF = self.all_ele(j).GetFEF(); % global FEF vector
               
               % Loop over element degrees of freedom
               for k = 1:length(ele_dof)
                   dof = ele_dof(k); % Find current relevant DOF
                   
                   % Map out element FEF vector
                   FEF_sys(dof-1+k) = FEF_sys(dof-1+k)+ele_FEF(k);
               end
           end
           disp('This is the concentrated load vector:')
           disp(concen_g)
           disp('This is the fixed end force vector:')
           disp(FEF_sys)
        end
    end
end
