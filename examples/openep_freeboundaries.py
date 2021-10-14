import openep
import openep.io
import openep.draw_routines

filename = "/Users/paul/github/openep-py/examples/data/new_dataset_2.mat"

case = openep.io.load_case(filename)
mesh = case.create_mesh(
    vertex_norms=False,
    recenter=False,
    back_faces=False
)

# generate a FreeBoundary object
free_boundaries = openep.draw_routines.get_freeboundaries(mesh)
fb = openep.draw_routines.FreeBoundary(mesh, **free_boundaries)

print(f"Perimeter lengths: {fb.calculate_lengths()}")
print(f"Cross-sectional areas: {fb.calculate_areas()}")

plotter = openep.draw_routines.draw_free_boundaries(fb)
plotter.background_color = "White"
plotter.show()
