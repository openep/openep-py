import openep

filename = "/Users/paul/github/openep-py/examples/data/new_dataset_2.mat"

case = openep.load_case(filename)
mesh = case.create_mesh()

# generate a FreeBoundary object
free_boundaries = openep.mesh.get_freeboundaries(mesh)

print(f"Perimeter lengths: {free_boundaries.calculate_lengths()}")
print(f"Cross-sectional areas: {free_boundaries.calculate_areas()}")

hsurf = openep.draw.draw_map(
    mesh,
    volt=case.fields["bip"],
    freeboundary_color="black",
    freeboundary_width=5,
    cmap="jet_r",
    minval=0,
    maxval=2,
    volt_below_color="brown",
    volt_above_color="magenta",
    nan_color="gray",
    plot=True,
)

plotter = openep.draw.draw_map(mesh)['hsurf']
plotter = openep.draw.draw_free_boundaries(free_boundaries, plotter=plotter)
plotter.background_color = "White"
plotter.show()
