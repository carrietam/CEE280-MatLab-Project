classdef AJCT_Element < handle
% Element class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        DOF % Degrees of freedom for element i
        FEF_G % Global element fixed end force vector
        FEF_LOC % Local element fixed end force vector
        GAM % Big gamma matrix
        Keg % Global element stiffness matrix
        Ke_loc % Local element stiffness matrix

    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        % Pass element properties to object
        function self = AJCT_Element(A, Ayy, Azz, E, element_nodes, Izz, Iyy, J, v, webdir, w)

            % Compute any other necessary properties
            L = ComputeLength(self,element_nodes); % Length
            ComputeTransformationMatrix(self,webdir,element_nodes,L); % Gamma Matrix
            ComputeElasticStiffnessMatrix(self, E, v, A, Ayy, Azz, Iyy, Izz, J, L); % Stiffness Matrix
            RetrieveDOF(self,element_nodes); % Degrees of Freedom
            ComputeFixedEndForces(self,L,w); % Fixed End Forces
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

        %% Compute local forces
        function local_forces = ComputeForces(self, dispglobal)
            % Calculate local forces using system displacement vector,
            % element local stiffness matrix, and element local FEF vector
            local_forces = self.Ke_loc*self.GAM*dispglobal + self.FEF_LOC;
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Compute the element's length
        function L = ComputeLength(self, element_nodes)

            % Get coordinates of Node objects
            ncoord1 = element_nodes(1,1).GetNodeCoord();
            ncoord2 = element_nodes(2,1).GetNodeCoord();

            % Calculate length of object
            dist = ncoord2 - ncoord1; % calc distance between nodes, % 1x3 matrix of distance between nodes
            L = sqrt((dist(1,1))^2+(dist(1,2))^2+(dist(1,3))^2); % 3D distance formula
        end
        
        %% Compute the element's geometric transformation matrix
        function ComputeTransformationMatrix(self,webdir,element_nodes,L)
            % Create small gamma matrix with placeholder 0's
            small_gam = zeros(3,3);

            % Get Node Coordinates
            ncoord1 = element_nodes(1,1).GetNodeCoord();
            ncoord2 = element_nodes(2,1).GetNodeCoord();

            % Calculate length of object
            dist = ncoord2 - ncoord1; % calc distance between nodes

            % Calculate actual small gamma matrix values
            small_gam(1,:) = dist/L; % Set first row of small gamma to...
            % normalized distance vector
            small_gam(2,:) = webdir; % Set second row of small...
            % gamma to cross product of third and first row
            small_gam(3,:) = (cross(small_gam(1,:),small_gam(2,:))); % Set third row of small...
            % gamma to cross product of normalized distance vector & unit 
           
            % Create big gamma matrix
            self.GAM =  blkdiag(small_gam,small_gam,small_gam,small_gam);
        end
        
        %% Compute the element's elastic stiffness matrix in local and
        % global coordinates
        function ComputeElasticStiffnessMatrix(self, E, v, A, Ayy, Azz,...
                Iyy, Izz, J, L )
            % Initiate sparse matrix
            Ke_loc = sparse(12,12); 
            % Note: if you want less degrees of freedom,
            % this will have to be manually changed later
     
            % Compute any variables used in stiffness matrix
            G = E/(2*(1+v));
            phiY = (12*E*Izz)/(G*Ayy*L^2);
            phiZ = (12*E*Iyy)/(G*Azz*L^2);
            
            % Calculate terms of stiffness matrix
            k_axial = E*A/L;
            k_flex_z1 = (12*E*Izz)/(L^3*(1+phiY));
            k_flex_z2 = (6*E*Izz)/(L^2*(1+phiY));
            k_flex_z3 = ((4+phiY)*E*Izz)/(L*(1+phiY));
            k_flex_z4 = ((2-phiY)*E*Izz)/(L*(1+phiY));
            k_flex_y1 = (12*E*Iyy)/(L^3*(1+phiZ));
            k_flex_y2 = (6*E*Iyy)/(L^2*(1+phiZ));
            k_flex_y3 = ((4+phiZ)*E*Iyy)/(L*(1+phiZ));
            k_flex_y4 = ((2-phiZ)*E*Iyy)/(L*(1+phiZ));
            k_torsion = G*J/L;

            % Create local axial stiffness matrix
            k_mat_axial = k_axial*[1,-1;-1,1];

            % Create local flexural z stiffness matrix
            k_mat_flexz = [k_flex_z1,k_flex_z2,-k_flex_z1,k_flex_z2;...
                k_flex_z2,k_flex_z3,-k_flex_z2,k_flex_z4;...
                -k_flex_z1,-k_flex_z2,k_flex_z1,-k_flex_z2;...
                k_flex_z2,k_flex_z4,-k_flex_z2,k_flex_z3];

            % Create local flexural y stiffness matrix
            k_mat_flexy = [k_flex_y1,-k_flex_y2,-k_flex_y1,-k_flex_y2;...
                -k_flex_y2,k_flex_y3,k_flex_y2,k_flex_y4;...
                -k_flex_y1,k_flex_y2,k_flex_y1,k_flex_y2;...
                -k_flex_y2,k_flex_y4,k_flex_y2,k_flex_y3];

            % Create local torsional stiffness matrix
            k_mat_torsion = k_torsion*[1,-1;-1,1];

            % Map local submatrices to local element stiffness matrix
            self.Ke_loc([1,7],[1,7]) = k_mat_axial;
            self.Ke_loc([2,6,8,12],[2,6,8,12]) = k_mat_flexz;
            self.Ke_loc([3,5,9,11],[3,5,9,11]) = k_mat_flexy;
            self.Ke_loc([4,10],[4,10]) = k_mat_torsion;
 
            % Convert local stiffness matrix to global stiffness matrix
            self.Keg = (self.GAM)'*self.Ke_loc*self.GAM;
        end

        %% Retrieve the degrees of freedom associated with the 
        % element nodes
        function RetrieveDOF(self,element_nodes)
            % Pass two Node objects to function
            node1 = element_nodes(1,1);
            node2 = element_nodes(2,1);
            
            % Get degrees of freedom of nodes
            self.DOF = [GetNodeDOF(node1);GetNodeDOF(node2)];
        end
        
        %% Compute fixed end forces of the element
        function ComputeFixedEndForces(self,L,w)
            % Pull uniform loads in each of the coordinate directions...
            % from w vector
            w = w';
            wx = w(1,1);
            wy = w(2,1);
            wz = w(3,1);
            
            % Compute fixed end forces of element in local coordinates
            Fx = -wx*L/2;
            Fy = -wy*L/2;
            Fz = -wz*L/2;
            Mx = 0;
            My = -wy*L^2/12;
            Mz = -wz*L^2/12;

            % Arrange the local fixed end forces/moments in a vector
            self.FEF_LOC = [Fx; Fy; Fz; Mx; -Mz; My;...
                            Fx; Fy; Fz; Mx;Mz;-My];

            % Translate local FEF/FEM to global coordinates
            self.FEF_G = (self.GAM')*(self.FEF_LOC);
        end
    end
end
