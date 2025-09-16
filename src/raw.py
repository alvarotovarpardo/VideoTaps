from pathlib import Path
import numpy as np, cv2
import matplotlib.pyplot as plt

width, height = 1280, 1024
pxpf = width * height

raw_path = Path("input/0_DDR_1X_2Y_1280_1024_video.raw") 
mm_in = np.memmap(raw_path, dtype=np.uint16, mode='r')
nframes = mm_in.size // pxpf
video_in = mm_in.reshape((nframes, height, width))
mean_in = video_in.mean(axis=0)