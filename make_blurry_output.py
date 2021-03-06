import cv2
import os
import numpy as np

BLURRY_OUTPUT_DIR = os.path.abspath('blurry-output')
print(BLURRY_OUTPUT_DIR)

I = cv2.imread('output/bridge.png', cv2.IMREAD_UNCHANGED)

print('Original Dimensions : ', I.shape)


def scale_image(scale_percent=100) -> np.ndarray:
    if scale_percent == 100:
        return I
    width = int(I.shape[1] * scale_percent / 100)
    height = int(I.shape[0] * scale_percent / 100)
    dim = (width, height)
    return cv2.resize(I, dim, interpolation=cv2.INTER_AREA)

# adding gaussian blur in the resized image
G = cv2.GaussianBlur(scale_image(scale_percent=100), (3, 3), 0)

cv2.imshow('original', I)
cv2.imshow("Resized image", G)
cv2.waitKey(0)
cv2.destroyAllWindows()

# save the gaussian blur image in input dir
cv2.imwrite(os.path.join(BLURRY_OUTPUT_DIR, 'bridge.png'), G)
