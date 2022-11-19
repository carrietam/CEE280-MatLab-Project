classdef AJCT_Node < handle
% Node class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        node_coord % 1x3 vector containing the x, y, and z coordinates of the node
        node_dof % degrees of freedom for the node
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %  node_coord:  1x3 vector containing the x, y, and z coordinates of the node
        function self = AJCT_Node(node_coord, node_num)
            self.node_coord = node_coord;
            self.assignDOF(node_num);
        end
        
        %% Get Node Coordinates
        %  Return "node_coord"
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end
        
        %% Get Node Degrees of freedom
        function pubdf = GetNodeDOF(self)
            pubdf = self.node_dof;
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Assign degrees of freedom for each node
        % Returns degrees of freedom
        function assignDOF(self,node_num)
               p = (node_num-1)*6;
               self.node_dof= (1:6)'+ [p;p;p;p;p;p];
        end
    end
end
