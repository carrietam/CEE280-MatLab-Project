classdef AJCT_Node < handle
% Node class for a 3-dimensional framed structure
    
    % Private properties go here
    properties (Access = private)
        node_coord % 3x1 vector containing the x, y, and z coordinates of the node
        node_num % assigned node number from mastan
        node_dof % degrees of freedom for the node
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        %  node_coord:  3x1 vector containing the x, y, and z coordinates of the node
        function self = AJCT_Node(node_coord, node_num)
            self.node_coord = node_coord;
            self.node_num = node_num;
            self.assignDOF();
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
        function assignDOF(self)
            for i=1:6
                df(i,1) = [(self.node_num-1)*6+i];
                self.node_dof = df;
            end
        end
    end
end
