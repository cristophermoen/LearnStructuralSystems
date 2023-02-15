using InstantFrame 

#EN.560.301 Structural Systems I 
#February 2023
#HW3 solution

material = InstantFrame.Material(names=["steel"], E=[29500.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["diagonal", "chord"], A=[0.32, 0.877], Iy=[1.0, 1.0], Iz=[1.0, 1.0], J=[1.0, 1.0])

connection = InstantFrame.Connection(names=["pin"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[0.0], rz=[0.0]))

node = InstantFrame.Node(numbers=[1, 2, 3], coordinates=[(0.0, 0.0, 0.0), (24.0, -48.0, 0.0), (48.0, 0.0, 0.0)])

element = InstantFrame.Element(numbers=[1, 2, 3], nodes=[(1,2), (2,3), (1,3)], orientation=[0.0, 0.0, 0.0], connections=[("pin", "pin"), ("pin", "pin"), ("pin", "pin")], cross_section=["diagonal", "diagonal", "chord"], material=["steel", "steel", "steel"])

support = InstantFrame.Support(nodes=[1, 2, 3], stiffness=(uX=[Inf,Inf,Inf], uY=[Inf,0.0, 0.0], uZ=[Inf,Inf, Inf], rX=[Inf,Inf,Inf], rY=[Inf,Inf,Inf], rZ=[0.0,0.0,0.0]))

uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(labels = ["snow"], nodes=[3], magnitudes=(FX=[0.0], FY=[-1.875], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))

analysis_type = "first order"
model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type)

