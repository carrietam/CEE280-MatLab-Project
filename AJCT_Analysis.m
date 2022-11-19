classdef AJCT_Analysis < handle
% Analysis class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        all_ele
        all_nodes
        
        % Note that properties under this comment can probably be
        % eliminated, but code needs to be modified to do this
        nnodes
        coord
        concen
        nele
        freeDOF
        suppDOF
        psDOF
        K_g
        KFF
        KFN
        KFS
        PF
        PS
        PP
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
            self.ClassifyDOF(fixity);
        end
        
        %% Run the analysis
        function RunAnalysis(self)
            self.CreateStiffnessMatrix();
            self.ComputeStiffnessSubMatrices();
            self.CheckKffMatrix();
            self.CreateLoadVectors();
            self.ComputeDisplacementReactions();
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Create Nodes
        % Create a vector of node objects, size nnodes
        function CreateNodes(self,nnodes,coord)
            self.all_nodes = AJCT_Node.empty;
            % Loop over number of nodes
            for i = 1:nnodes
                % Create node & add to all_nodes vector
                node_i = AJCT_Node(coord(i,1:3),i);
                self.all_nodes(i) = node_i;
            end

        end
        
        %% Create Elements
        % Create a vector of element objects, size nele
        function CreateElements(self,A, Ayy, Azz, E, ends, Izz, Iyy, J, nele, v, webdir, w)
            self.all_ele = AJCT_Element.empty;
            % Loop over number of elements
            for i = 1:nele
                % Create element & add to all_ele vector
                element_i = AJCT_Element(A(i), Ayy(i), Azz(i), E(i),...
                    [self.all_nodes(ends(i,1)); self.all_nodes(ends(i,2))], Izz(i), Iyy(i), J(i), v(i),...
                    webdir(i,1:3), w(i,1:3));
                self.all_ele(i) = element_i;
            end
        end
        
        %% Create Stiffness Matrix
        % Assemble system stiffness matrix
        function CreateStiffnessMatrix(self)
            nDOF = self.nnodes*6; % Calculate the number of degrees of freedom

            for i = 1:self.nele
                % Get necessary element properties
                disp('this is i')
                disp(i)
                disp('this is self ele')
                disp(self.all_ele)
                ele_dof = self.all_ele(i).GetDOF() % relevant DOFs
                K_ele = self.all_ele(i).GetKGloEleMatrix() % global element...
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
                K_mapped = sparse(e1,e2,K_ele,nDOF,nDOF)
                
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
            disp(full(K_g))
            self.K_g = K_g;
        end
        
        %% Create load vectors
        function CreateLoadVectors(self)
           nDOF = self.nnodes*6; % Calculate the number of degrees of freedom
           concen_g = zeros(nDOF,1); % Initialize concen_g
           
           %% Create concentrated load vectors
           for i = 1:self.nnodes
               % Get necessary element properties
               node_DOF = self.all_nodes(i).GetNodeDOF();
                
               % Map out element concen vector onto global concen vector
               concen_g(node_DOF,1) = [self.concen(i,1:6)]';
               % need a line to continuously add these
           end
           
           %% Classify concentrated load vectors
           % Initialize classified concen vectors
           concen_f = zeros(length(self.freeDOF),1); % free DOFs
           concen_s = zeros(length(self.suppDOF),1); % support DOFs
           concen_p = zeros(length(self.psDOF),1); % prescribed DOFs

           % Create concentrated load vectors for classified DOFs
           concen_f = concen_g(self.freeDOF);
           concen_s = concen_g(self.suppDOF);
           concen_p = concen_g(self.psDOF);

           %% Create fixed end force vector
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
                   FEF_sys(dof) = FEF_sys(dof)+ele_FEF(k);
               end
           end

           %% Classify fixed end force vectors
           % Initialize classified concen vectors
           FEF_f = zeros(length(self.freeDOF),1); % free DOFs
           FEF_s = zeros(length(self.suppDOF),1); % support DOFs
           FEF_p = zeros(length(self.psDOF),1); % prescribed DOFs

           % Create fixed end force vectors for classified DOFs
           

           % Loop over free element degrees of freedom
           for i = 1:nfdof
               FEF_f(i,1) = [FEF_sys(self.freeDOF(i))];
           end

           % Loop over support element degrees of freedom
           for i = 1:nsdof
               FEF_s(i,1) = [FEF_sys(self.suppDOF(i))];
           end

           % Loop over prescribed element degrees of freedom
           for i = 1:npdof
               FEF_p(i,1) = [FEF_sys(self.psDOF(i))];
           end

           %% Find total force vectors based on classified DOFs
           self.PF = concen_f + FEF_f; % free DOFs
           self.PS = concen_s + FEF_s; % support DOFs
           self.PP = concen_p + FEF_p; % prescribed DOFs

           %% Display results
%            disp('This is the concentrated load vector:')
%            disp(concen_g)
%            disp('This is the global fixed end force vector:')
%            disp(FEF_sys)
%            disp(size(FEF_sys))
% 
%            disp('This is the free fixed end force vector:')
%            disp(FEF_f)
%            disp(size(FEF_f))
% 
%            disp('This is the support fixed end force vector:')
%            disp(FEF_s)
%            disp(size(FEF_s))
% 
%            disp('This is the prescribed fixed end force vector:')
%            disp(FEF_p)
%            disp(size(FEF_p))
        end

        %% Classify Degrees of Freedom
        function ClassifyDOF(self, fixity)
            % Transform the fixity matrix so that different columns 
            % correspond to different nodes
            fix_t = fixity';

            % Classify DOFs as free, supports, or prescribed
            self.freeDOF = find(isnan(fix_t)); % Free DOFs
            self.suppDOF = find(fix_t == 0); % Support DOFs
            self.psDOF = find(fix_t ~=0 & ~isnan(fix_t)); % Prescribed DOFs
        end

        %% Compute Stiffness Submatrices
        function ComputeStiffnessSubMatrices(self)
            % Retrieve necessary classifications of degrees of freedom
            freeDOF = self.freeDOF;
            suppDOF = self.suppDOF;
            psDOF = self.psDOF;

            % Create submatrices from global K matrix
            self.KFF = self.K_g(freeDOF,freeDOF);
            self.KFS = self.K_g(suppDOF,suppDOF);
            self.KFN = self.K_g(psDOF,psDOF);
        end

        %% Check KFF Matrix to see if structure is stably supported
        function CheckKffMatrix(self)
            % Calculate and display condition number
            cond_kff = condest(self.KFF);
            disp('The condition number for the KFF matrix is:')
            disp(cond_kff)

            % Display the number of sig figs that are to be lost
            lost_sigfigs = log10(cond_kff);
            disp('The estimated number of significant figures to be lost is:')
            disp(lost_sigfigs)

            if log10(cond_kff)<10^12
                AFLAG = 1;
            else
                AFLAG = 0;
                disp('Unstable Structure. Unsuccessful Analysis')
                quit
            end
        end

        %% Compute global displacements & reactions at each DOF
        function ComputeDisplacementReactions(self)
            % Compute displacements at free degrees of freedom
            nDOF = self.nnodes*6;
            DN = zeros(length(self.psDOF));

            try
                DF = self.KFF\(self.PF-self.KFN*DN);
            catch
%                 disp('KFF');
%                 disp(self.KFF);
%                 disp(size(self.KFF));
%                 disp(full(self.KFF));
%                 disp('PF');
%                 disp(self.PF);
%                 disp(size(self.PF));
                DF = self.KFF\self.PF;
            end

            % Initalize global displacement vector
            disp_glo_vector = zeros(nDOF,1);
            
            % Loop over free degrees of freedom and put them in a global
            % displacement vector
            for i = 1:length(self.freeDOF)
                % Create a vector of displacements
                disp_glo_vector(self.freeDOF(i)) = DF(i);
            end

            % Loop over prescribed degrees of freedom and put them in the
            % global displacement vector
            for i = 1:length(self.psDOF)
                disp_glo_vector(self.psDOF(i)) = DN(i);
            end

            % Format displacement vector to be outputted back to MASTAN2
            DEFL = reshape(disp_glo_vector,self.nnodes,6);

%             REACT(i,1:6)
        end
    end
end
