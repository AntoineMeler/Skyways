import numpy as np

data = np.fromfile("../data/data_11.f32", dtype=np.float32).reshape(1024, 1024, 12, 12) # (y, x, heading, altitude)

print(data)