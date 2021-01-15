import cv2
import numpy
import os


OUTPUT_PATH = os.path.abspath('output')
I = cv2.imread('output/bridge.jpeg')
cv2.imwrite(os.path.join(OUTPUT_PATH, 'bridge.png'), I)
