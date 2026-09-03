import pyvista as pv

def convert_vtk_to_custom_obj(scale, vtk_filename, obj_filename):
    # Load the binary VTK file
    mesh = pv.read(vtk_filename)

    # Extract coordinates and the number of points
    points = mesh.points
    points = scale * points  # Scale coordinates by the given scale
    num_points = mesh.n_points

    # Extract the nodal radius array
    if "radius" in mesh.point_data:
        radii = 1.0*mesh.point_data["radius"]
    else:
        radii = [0.0] * num_points

    #radii = [0.00002] * num_points

    with open(obj_filename, 'w') as f:
        f.write(f"# Extracted vertices 1-{num_points} from {vtk_filename}\n")
        f.write("# v x y z radius auxiliary0 auxiliary1\n")

        # Write vertices with the radius and auxiliary zeros
        for i in range(num_points):
            x, y, z = points[i]
            r = radii[i]
            f.write(f"v {x:.6f} {y:.6f} {z:.6f} {r:.6f} 0 0\n")

        # Write line connectivity using 1-based indexing
        # mesh.lines is a flat array: [n_points, p1, p2, ..., n_points, p1, p2, ...]
        lines = mesh.lines
        idx = 0
        while idx < len(lines):
            n_pts = lines[idx]
            # Break down polylines into individual segments
            for j in range(n_pts - 1):
                p1 = lines[idx + 1 + j] + 1
                p2 = lines[idx + 1 + j + 1] + 1
                f.write(f"l {p1} {p2}\n")
            idx += n_pts + 1

# Execute the conversion
vtkfilename = 'UIC75_01_centerline.vtk'
objfilename = 'skeleton_initial.obj'
scale = 1.0  # Set the desired scale factor
convert_vtk_to_custom_obj(scale, vtkfilename, objfilename)
