%Cecilia Doyle
%Structural Systems 1
%Final Project

%these coordinates are based on a dome of radius 1
data = readtable('dome_geometry.xlsx');
data = table2array(data);

%nodes are row 3-36,[x y z]
node_geometry = data(3:63,1:3);
%get x, y, and z coordinates
x = node_geometry(:,1);
z = node_geometry(:,2);
y = node_geometry(:,3);
%flip y and z to correct visualization (given coordinates were flipped)
node_geometry = [x,y,z];
%scale dome to radius of 20 feet
scaled = [];
S = [20 0 0 0
    0 20 0 0
    0 0 20 0
    0 0 0 1];
for i = 1:length(node_geometry)
    xyz = [x(i),y(i),z(i),1];
    new = xyz * S;
    scaled = [scaled; new];
end
%replace coordinates with coordinates of scaled dome
node_geometry = scaled(:,1:3); %cords = [x' y' z']
x = node_geometry(:,1);
y = node_geometry(:,2);
z = node_geometry(:,3);
%make bottom flat, set all base z values to lowest value
z(47:61) = z(61);
node_geometry(:,3) = z;

%members are row 67 - 231 of data
members = data(67:231,1:2);
members = members + 1; %lists starts with node 0, shift all #s by 1 so starts at 1

%get lengths of members
lengths = get_member_lengths(members, node_geometry);


%%%%%GRAPH ORIGINAL SHAPE%%%%%
%plot member 1 to start 3d graph
%get nodes and their coordinates
nodei = members(1,1);
nodej = members(1,2);
xi = node_geometry(nodei, 1);
yi = node_geometry(nodei, 2);
zi = node_geometry(nodei, 3);
xj = node_geometry(nodej, 1);
yj = node_geometry(nodej, 2);
zj = node_geometry(nodej, 3);

figure()
plot3([xi,xj],[yi,yj],[zi,zj], '-b');
title('Original Geodesic Dome');
xlabel('x');
ylabel('y');
zlabel('z');

%plot rest of members
for i = 1:length(members)
    nodei = members(i,1);
    nodej = members(i,2);
    
    xi = node_geometry(nodei, 1);
    yi = node_geometry(nodei, 2);
    zi = node_geometry(nodei, 3);
    xj = node_geometry(nodej, 1);
    yj = node_geometry(nodej, 2);
    zj = node_geometry(nodej, 3);
    
    hold on
    plot3([xi,xj],[yi,yj],[zi,zj], '-b');
    
end

%plot again for combined graph
%plot member 1 to start 3d graph
%get nodes and their coordinates
nodei = members(1,1);
nodej = members(1,2);
xi = node_geometry(nodei, 1);
yi = node_geometry(nodei, 2);
zi = node_geometry(nodei, 3);
xj = node_geometry(nodej, 1);
yj = node_geometry(nodej, 2);
zj = node_geometry(nodej, 3);

figure()
plot3([xi,xj],[yi,yj],[zi,zj], '-b');
title('Combined Geodesic Dome');
xlabel('x');
ylabel('y');
zlabel('z');

%plot rest of members
for i = 1:length(members)
    nodei = members(i,1);
    nodej = members(i,2);
    
    xi = node_geometry(nodei, 1);
    yi = node_geometry(nodei, 2);
    zi = node_geometry(nodei, 3);
    xj = node_geometry(nodej, 1);
    yj = node_geometry(nodej, 2);
    zj = node_geometry(nodej, 3);
    
    hold on
    plot3([xi,xj],[yi,yj],[zi,zj], '-b');
    
end


%properties
A = 20; %in^2
E = 1600; %ksi, E of wood

%num elements and nodes
num_elem = length(members);
num_nodes = length(node_geometry);
%external forces
external_forces = zeros(num_nodes*3,1); %sets ext forces to 0, has x,y,z directions
highest = find(z == max(z)); %top node
external_forces(highest*3) = -200; %kips


%set all supports as 0
supports = zeros(num_nodes*3,1);
%base nodes are pinned - all three dof are fixed
fixed = find(z == min(z));
for i = 1:length(fixed)
    %get node
    j = fixed(i);
    %set that node's dof to 1
    supports((3*j)-2,:) = 1;
    supports((3*j)-1,:) = 1;
    supports((3*j),:) = 1;
end


%build global stiffness matrices
components = [];

for i = 1:length(members)
    %get the 2 nodes that the member is connecting
    iNode = members(i,1);
    jNode = members(i,2);
    %get the coordinates of each node
    iNodePosition=node_geometry(iNode,:);
    jNodePosition=node_geometry(jNode,:);

    xj = jNodePosition(1);
    xi = iNodePosition(1);
    yj = jNodePosition(2);
    yi = iNodePosition(2);
    zj = jNodePosition(3);
    zi = iNodePosition(3);
    %length of current member
    IJ = lengths(i);
    %matrix components (eq 5.11)
    l11 = (xj - xi) / IJ;
    l12 = (yj - yi) / IJ;
    l13 = (zj - zi) / IJ;
    %add to list
    add = [l11,l12,l13];

    components = [components; add];

end


%Global Stiffness Matrices
k_global = cell(num_elem,1);
for i = 1:num_elem
    mem_components = components(i,:);
    L = lengths(i);
    
    k_global{i} = get_global_k(E, A, L, mem_components);
end


%Combined Stiffness Matrix
K = zeros(num_nodes*3);
for i = 1:num_elem
    node_i = members(i,1);
    node_j = members(i,2);
    
    dof_1 = node_i*3 - 2;
    dof_2 = node_i*3 - 1;
    dof_3 = node_i*3;
    dof_4 = node_j*3 - 2;
    dof_5 = node_j*3 - 1;
    dof_6 = node_j*3;
    
    dof_all = [dof_1; dof_2; dof_3; dof_4; dof_5; dof_6];
    
    for k = 1:length(dof_all)
        
        for j = 1:length(dof_all)
            
            K(dof_all(k), dof_all(j)) = K(dof_all(k), dof_all(j)) + k_global{i}(k,j);
        end  
    end
end

% Define Fixed Nodes
free_dof = find(~supports);
fixed_dof = find(supports);

% Create Reduced Stiffness Matrix
%get only values corresponding to the free degrees of freedom
Kff = K(free_dof, free_dof);
Ff = external_forces(free_dof);
uf = Kff^(-1) * Ff;


%%GRAPH DISPLACED SHAPE%%
ux = uf(1:3:length(uf));
ux = [ux; zeros(15,1)];
uy = uf(2:3:length(uf));
uy = [uy; zeros(15,1)];
uz = uf(2:3:length(uf));
uz = [uz; zeros(15,1)];

x = x + ux;
y = y + uy;
z = z + uz;
node_geometry = [x,y,z];

%PLOT TOGETHER
for i = 1:length(members)
    nodei = members(i,1);
    nodej = members(i,2);
    
    xi = node_geometry(nodei, 1);
    yi = node_geometry(nodei, 2);
    zi = node_geometry(nodei, 3);
    xj = node_geometry(nodej, 1);
    yj = node_geometry(nodej, 2);
    zj = node_geometry(nodej, 3);
    
    hold on
    plot3([xi,xj],[yi,yj],[zi,zj], '-r');
    
end


%PLOT ALONE
nodei = members(1,1);
nodej = members(1,2);

xi = node_geometry(nodei, 1);
yi = node_geometry(nodei, 2);
zi = node_geometry(nodei, 3);
xj = node_geometry(nodej, 1);
yj = node_geometry(nodej, 2);
zj = node_geometry(nodej, 3);
figure()
plot3([xi,xj],[yi,yj],[zi,zj], '-r');
title('Displaced Geodesic Dome');
xlabel('x');
ylabel('y');
zlabel('z');


for i = 1:length(members)
    nodei = members(i,1);
    nodej = members(i,2);
    
    xi = node_geometry(nodei, 1);
    yi = node_geometry(nodei, 2);
    zi = node_geometry(nodei, 3);
    xj = node_geometry(nodej, 1);
    yj = node_geometry(nodej, 2);
    zj = node_geometry(nodej, 3);
    
    hold on
    plot3([xi,xj],[yi,yj],[zi,zj], '-r');
    
end






%%%%%%%  FUNCTIONS  %%%%%%%

%Get Lengths of Members
function len_members = get_member_lengths(members, node_geometry)

    len_members=zeros(length(members(:,1)),1);

    for i=1:length(members(:,1))

        iNode=members(i,1);
        jNode=members(i,2);

        iNodePosition=node_geometry(iNode,:);
        jNodePosition=node_geometry(jNode,:);

        len_members(i) = norm(iNodePosition-jNodePosition);

    end
end

%Create a Global Stiffness Matrix (K global)
function kglobal = get_global_k(E, A, L, mem_components)
    
    l11 = mem_components(1);
    l12 = mem_components(2);
    l13 = mem_components(3);
    
	kglobal = E * A/ L * [l11^2, l11*l12, l11*l13, -(l11^2), -(l11*l12), -(l11*l13)
        l11*l12, l12^2, l12*l13, -(l11*l12), -(l12^2), -(l12*l13)
        l13*l11, l13*l12, l13^2, -(l11*l13), -(l12*l13), -(l13^2)
        -(l11^2), -(l11*l12), -(l11*l13), l11^2, l11*l12, l11*l13
        -(l11*l12), -(l12^2), -(l12*l13), l11*l12, l12^2, l12*l13
        -(l13*l11), -(l13*l12), -(l13^2), l11*l13, l12*l13, l13^2];
end
