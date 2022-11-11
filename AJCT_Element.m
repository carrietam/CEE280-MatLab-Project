classdef AJCT_Element < handle
% Element class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        A % Total Area
        Ayy % Shear Area in y axis
        Azz % Shear Area in z axis
        d % Vector containing distance between Nodes
        DOF % Degrees of freedom for element i
        E % Modulus of Elasticity
        element_nodes % vector of 2 Node objects
        FEF_G % Global element fixed end force vector
        gam % Small gamma matrix
        gamt % Small gamma transpose matrix
        GAM % Big gamma matrix
        GAMT % Big gamma transpose matrix
        Iyy % Moment of inertia in the y axis
        Izz % Moment of inertia in the z axis
        J % Torsional constant
        Keg % Global element stiffness matrix
        L % Length
        v % Poisson's ratio
        webdir % Unit web vector for element i
        w % 3x1 vector of uniform loads in 3 coord. dir.s

    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        % Pass element properties to object
        function self = AJCT_Element(A, Ayy, Azz, E, element_nodes, Izz, Iyy, J, v, webdir, w)

            % Initialize given properties
            self.A = A;
            self.Ayy = Ayy;
            self.Azz = Azz;
            self.E = E;
            self.element_nodes = element_nodes;
            self.Iyy = Iyy;
            self.Izz = Izz;
            self.J = J;
            self.v = v;
            self.webdir = webdir;
            self.w = w;

            % Compute any other necessary properties
            self.ComputeLength(); % Length
            self.ComputeTransformationMatrix(); % Gamma Matrix
            self.ComputeElasticStiffnessMatrix(); % Stiffness Matrix
            self.RetrieveDOF(); % Degrees of Freedom
            self.ComputeFixedEndForces(); % Fixed End Forces
        end
        
        %% Make degrees of freedom public
        function eDOF = GetDOF(self)
            eDOF = self.DOF;
        end
        
        %% Make global stiffness matrix public
        function Kglo = GetKGloEleMatrix(self)
            Kglo = self.Keg;
        end
        
        %% Make global fixed end force vector public
        function fefg = GetFEF(self)
            fefg = self.FEF_G;
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Compute the element's length
        function ComputeLength(self)
            % Pass two Node objects to function
            node1 = self.element_nodes(1,1);
            node2 = self.element_nodes(2,1);

            % Get coordinates of Node objects
            ncoord1 = node1.GetNodeCoord();
            ncoord2 = node2.GetNodeCoord();

            % Calculate length of object
            d = ncoord2 - ncoord1; % calc distance between nodes
            d = d';
            self.d = d;
            length = sqrt((d(1,1))^2+(d(2,1))^2+(d(3,1))^2); % 3D distance formula
            self.L = length; % assign length property to object
        end
        
        %% Compute the element's geometric transformation matrix
        function ComputeTransformationMatrix(self)
            % Create small gamma matrix with placeholder 0's
            gam = zeros(3,3);

            % Calculate actual small gamma matrix values
            nor = self.d/self.L; % Normalize vector containing distance btwn nodes
            gam(1,:) = transpose(nor); % Set first row of small gamma to...
            % normalized distance vector
            gam(2,:) = self.webdir; % Set second row of small...
            % gamma to cross product of third and first row
            gam(3,:) = (cross(gam(1,:),gam(2,:))); % Set third row of small...
            % gamma to cross product of normalized distance vector & unit 
            % y-axis vector
            self.gam = gam;

            % Create small gamma transpose matrix
            gamt = transpose(gam);
            self.gamt = gamt;

            % Create big gamma matrix
            blanks = zeros(size(gam)); % Create matrix of zeros the size...
            % of small gamma matrix
            GAM = blkdiag(gam,gam,gam,gam);
            self.GAM = GAM;

            % Create big gamma transpose matrix
            GAMT = blkdiag(gamt,gamt,gamt,gamt);
            self.GAMT = GAMT;
        end
        
        %% Compute the element's elastic stiffness matrix in local and global coordinates
        function ComputeElasticStiffnessMatrix(self)
            % Compute any variables used in stiffness matrix
            G = self.E/(2*(1+self.v));
            phiY = (12*self.E*self.Izz)/(G*self.Ayy*self.L^2);
            phiZ = (12*self.E*self.Iyy)/(G*self.Azz*self.L^2);
            
            % Calculate stiffness matrix
            Ke = zeros(12);
            Ke(1,1) = self.E*self.A/self.L;
            Ke(7,1) = -self.E*self.A/self.L;
            Ke(2,2) = (12*self.E*self.Izz)/(self.L^3*(1+phiY));
            Ke(6,2) = (6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(8,2) = -(12*self.E*self.Izz)/(self.L^3*(1+phiY));
            Ke(12,2) = (6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(3,3) = (12*self.E*self.Iyy)/(self.L^3*(1+phiZ));
            Ke(5,3) = -(6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(9,3) = -(12*self.E*self.Iyy)/(self.L^3*(1+phiZ));
            Ke(11,3) = -(6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(4,4) = G*self.J/self.L;
            Ke(10,4) = -G*self.J/self.L;
            Ke(3,5) = -(6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(5,5) = ((4+phiZ)*self.E*self.Iyy)/(self.L*(1+phiZ));
            Ke(9,5) = (6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(11,5) = ((2-phiZ)*self.E*self.Iyy)/(self.L*(1+phiZ));
            Ke(2,6) = (6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(6,6) = ((4+phiY)*self.E*self.Izz)/(self.L*(1+phiY));
            Ke(8,6) = -(6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(12,6) = ((2-phiY)*self.E*self.Izz)/(self.L*(1+phiY));
            Ke(1,7) = -self.E*self.A/self.L;
            Ke(7,7) = self.E*self.A/self.L;
            Ke(2,8) = -(12*self.E*self.Izz)/(self.L^3*(1+phiY));
            Ke(6,8) = -(6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(8,8) = (12*self.E*self.Izz)/(self.L^3*(1+phiY));
            Ke(12,8) = -(6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(3,9) = -(12*self.E*self.Iyy)/(self.L^3*(1+phiZ));
            Ke(5,9) = (6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(9,9) = (12*self.E*self.Iyy)/(self.L^3*(1+phiZ));
            Ke(11,9) = (6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(4,10) = -G*self.J/self.L;
            Ke(10,10) = G*self.J/self.L;
            Ke(3,11) = -(6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(11,11) = ((4+phiZ)*self.E*self.Iyy)/(self.L*(1+phiZ));
            Ke(9,11) = (6*self.E*self.Iyy)/(self.L^2*(1+phiZ));
            Ke(5,11) = ((2-phiZ)*self.E*self.Iyy)/(self.L*(1+phiZ));
            Ke(2,12) = (6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(12,12) = ((4+phiY)*self.E*self.Izz)/(self.L*(1+phiY));
            Ke(8,12) = -(6*self.E*self.Izz)/(self.L^2*(1+phiY));
            Ke(6,12) = ((2-phiY)*self.E*self.Izz)/(self.L*(1+phiY));

            % Convert local stiffness matrix to global stiffness matrix
            self.ComputeTransformationMatrix();
            Keg = self.GAMT*Ke*self.GAM;
            self.Keg = Keg;

            % Output length results
            output_l = ['The length of the element is: ',num2str(self.L)];
            disp(output_l)

            % Output gamma matrix results
            output_g = ['The small gamma matrix for this element is:'];
            output_G = ['The big gamma matrix for this element is:'];
            disp(output_g)
            disp(self.gam)
            disp(output_G)
            disp(self.GAM)

            % Output stiffness matrix results
            output_k = ['The local stiffness matrix of the element is:'];
            output_K = ['The global stiffness matrix of the element is:'];
            disp(output_k)
            disp(Ke)
            disp(output_K)
            disp(Keg)
        end
        
        %% Retrieve the degrees of freedom associated with the element nodes
        function RetrieveDOF(self)
            % Pass two Node objects to function
            node1 = self.element_nodes(1,1);
            node2 = self.element_nodes(2,1);
            
            % Get degrees of freedom of nodes
            dof = [GetNodeDOF(node1);GetNodeDOF(node2)];
            self.DOF = dof;
        end
        
        %% Compute fixed end forces of the element
        function ComputeFixedEndForces(self)
            % Pull in any necessary properties
            L = self.L;
            GAMT = self.GAMT;
            
            % Pull uniform loads in each of the coordinate directions...
            % from w vector
            self.w = self.w'
            wx = self.w(1,1);
            wy = self.w(2,1);
            wz = self.w(3,1);
            
            % Compute FEF / FEM in local coordinates
            FEF_LOC = zeros(12,1);
            
            FEF_LOC(1,1) = wy*L/2;
            FEF_LOC(2,1) = wx*L/2; % axial?
            FEF_LOC(3,1) = wz*L/2;
            FEF_LOC(4,1) = wz*L^2/12;
            FEF_LOC(5,1) = 0; % no moment?
            FEF_LOC(6,1) = wy*L/2;
            FEF_LOC(7,1) = wy*L^2/12;
            FEF_LOC(8,1) = wx*L/2; % axial?
            FEF_LOC(9,1) = wz*L/2;
            FEF_LOC(10,1) = wz*L^2/12;
            FEF_LOC(11,1) = 0; % no moment?
            FEF_LOC(12,1) = wy*L^2/12;
                            
            % Compute FEF / FEM in global coordinates
            FEF_G = GAMT*FEF_LOC;
            self.FEF_G = FEF_G;
                
        end
    end
end
