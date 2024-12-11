import numpy as np

def parse_nastran_short_field(file_path):
    """
    Parse a NASTRAN short-field .bdf file to extract nodes and tetrahedral elements.
    """
    nodes = {}
    tetrahedra = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('GRID'):
                # Parse GRID card
                node_id = int(line[8:16].strip())
                x = float(line[24:32].strip().replace('D', 'E'))
                y = float(line[32:40].strip().replace('D', 'E'))
                z = float(line[40:48].strip().replace('D', 'E'))
                nodes[node_id] = (x, y, z)
                print (nodes[node_id])
            elif line.startswith('CTETRA'):
                # Parse CTETRA card
                node_ids = [int(line[24:32].strip()),
                            int(line[32:40].strip()),
                            int(line[40:48].strip()),
                            int(line[48:56].strip())]
                tetrahedra.append(node_ids)
                #print (tetrahedra)

    # Calculate the bounding box of the mesh
    node_coords = np.array(list(nodes.values()))
    min_coords = node_coords.min(axis=0)
    max_coords = node_coords.max(axis=0)

    return nodes, tetrahedra, min_coords, max_coords

def generate_cartesian_grid(x_range, y_range, z_range, spacing):
    """
    Generate a Cartesian grid of points within the specified ranges and spacing.
    """
    x = np.arange(x_range[0], x_range[1] + spacing, spacing)
    y = np.arange(y_range[0], y_range[1] + spacing, spacing)
    z = np.arange(z_range[0], z_range[1] + spacing, spacing)
    return np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)

def signed_volume(a, b, c, d):
    """
    Compute the signed volume of a tetrahedron defined by points a, b, c, d.
    """
    mat = np.array([
        [a[0], a[1], a[2], 1],
        [b[0], b[1], b[2], 1],
        [c[0], c[1], c[2], 1],
        [d[0], d[1], d[2], 1]
    ])
    return np.linalg.det(mat) / 6.0

def is_point_inside_mesh(point, nodes, tetrahedra):
    """
    Check if a point is inside any tetrahedron in the mesh.
    """
    min = 1e4
    for tetra in tetrahedra:
        #print (tetra)
        a, b, c, d = [nodes[node_id] for node_id in tetra]
        #print ("point and abcd", point, a,b,c,d)
        v_abcd = signed_volume(a, b, c, d)
        v_pbcd = signed_volume(point, b, c, d)
        v_padc = signed_volume(point, a, d, c)
        v_pabd = signed_volume(point, a, b, d)
        v_pabc = signed_volume(point, a, b, c)
        #print("REST, ", abs(v_pbcd + v_padc + v_pabd + v_pabc - v_abcd))
        check = abs(v_pbcd + v_padc + v_pabd + v_pabc - v_abcd) 
        if check < 0.1 and \
           np.sign(v_pbcd) == np.sign(v_abcd) and \
           np.sign(v_padc) == np.sign(v_abcd) and \
           np.sign(v_pabd) == np.sign(v_abcd) and \
           np.sign(v_pabc) == np.sign(v_abcd):
            return True
        if (check < min):
          min = check
    print ("min ", min)
    return False

def write_sph_to_k_file_with_elements(points, output_file):
    """
    Write points as *NODE and *ELEMENT_SPH entries to an LS-DYNA .k file.
    """
    with open(output_file, 'w') as file:
        file.write("*KEYWORD\n")

        # Write *NODE entries
        file.write("*NODE\n")
        for i, point in enumerate(points, start=1):
            file.write(f"{i:10d}{point[0]:10.6f}{point[1]:10.6f}{point[2]:10.6f}\n")

        # Write *ELEMENT_SPH entries
        file.write("*ELEMENT_SPH\n")
        for i in range(1, len(points) + 1):
            file.write(f"{i:10d}{i:10d}\n")

        file.write("*END\n")

# Example usage
nastran_file = 'mesh.bdf'  # Replace with your NASTRAN file path
nodes, tetrahedra, min_coords, max_coords = parse_nastran_short_field(nastran_file)
print(f"Bounding box of the mesh: Min {min_coords}, Max {max_coords}")
print ("Tetra count ", len (tetrahedra))
# Define Cartesian grid range and spacing
#x_range = (0.0, 2.0)
#y_range = (0.0, 2.0)
#z_range = (0.0, 2.0)
spacing = 2.0
x_range = (min_coords[0]+spacing/2.0,max_coords[0]-spacing/2.0)
y_range = (min_coords[1]+spacing/2.0,max_coords[1]-spacing/2.0)
z_range = (min_coords[2]+spacing/2.0,max_coords[2]-spacing/2.0)



grid_points = generate_cartesian_grid(x_range, y_range, z_range, spacing)

# Check which points are inside the mesh
inside_points = []
i = 0
ip = 0
for point in grid_points:
    print ("point", point)
    print(str(i) + " of "+ str(len(grid_points)) + str(point)+ ", inside points " + str(ip))
    if is_point_inside_mesh(point, nodes, tetrahedra):
        inside_points.append(point)
        ip +=1
    i+=1
# Write the inside points to an LS-DYNA .k file with *NODE and *ELEMENT_SPH
output_file = "sph_particles.k"
write_sph_to_k_file_with_elements(inside_points, output_file)
print(f"SPH particles written to {output_file}.")
