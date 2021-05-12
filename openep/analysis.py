from matplotlib.cm import jet, rainbow, jet_r, seismic


def compute_field(
    mesh,
    fieldname,
    minval=0,
    maxval=1,
    color_map=jet_r,
    below_color=(0, 0, 0, 255),
    above_color=(255, 0, 255, 255),
    nan_color=(50, 50, 50, 255),
):
    case = mesh._kwargs["parent"]
    field = case.fields[fieldname]
    new_field = np.zeros((field.shape[0], 4), np.int32)

    for idx in range(len(field)):
        val = field[idx]

        if np.isnan(val):
            col = nan_color
        elif val < minval:
            col = below_color
        elif val > maxval:
            col = above_color
        else:
            valscaled = (val - minval) / (maxval - minval)
            col = color_map(valscaled)

            if isinstance(col[0],float):
                col=[int(c*255) for c in col]

        new_field[idx] = col

    mesh.visual.vertex_colors[:] = new_field

    return new_field
