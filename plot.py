import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

data_dir = Path('./plot_data')
out_dir = Path('./gallery')
# figure size
Nx = 2001
Ny = 2001
# colors
color_graph = np.array([0.0, 0.0, 0.0])  # graph, black
color_bg = np.array([1.0, 1.0, 1.0])  # background, white


for i in range(30):
    # read data
    array = np.fromfile(data_dir.joinpath(f'data_{i}.bin'), dtype=np.int32).astype(float)
    array = array.reshape([Nx, Ny])
    array = array[:, ::-1].T  # change axis order to make matplotlib happy
    # to RGB
    rgb_array = np.expand_dims(1.0 - array, axis=2) * color_bg + np.expand_dims(array, axis=2) * color_graph
    # save figure
    plt.imsave(out_dir.joinpath(f'{i}.png'), rgb_array)
