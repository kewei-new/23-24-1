function [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_local_Guass_1D(mesh_point,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)
lower_bound=min(mesh_point);
upper_bound=max(mesh_point);
Gauss_coefficient_local_1D=(upper_bound-lower_bound)*Gauss_coefficient_reference_1D/2;
Gauss_point_local_1D=(upper_bound-lower_bound)*Gauss_point_reference_1D/2+(upper_bound+lower_bound)/2;