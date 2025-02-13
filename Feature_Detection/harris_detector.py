import numpy as np
import matplotlib.pyplot as plt

# For opening and displaying image
from skimage.data import imread
from skimage.io import imshow
# from skimage.feature import corner_harris, corner_peaks

# for convolution and gaussian function
from scipy.signal import gaussian
from scipy.signal import convolve2d

import argparse

################################################################################

def first_refine(Corner_response_Mc, threshold, local):
	"""
	Reponse above threshold and local maximum (8 neighbors) => detection
	1. Reponses below threshold * max are set to 0
	2. Keep local maxima (in (2xlocal  + 1) x (2xlocal + 1) pixels)
	:param Corner_response_Mc: matric of corners response (harris criterion)
	:param threshold: Threshold under which corner is supposed to be 0
	:param local: Width of the local window for local maxima computation
	:return: first refined corners
	"""
	# Reponses below threshold * max are set to 0
	Corner_response_Mc[Corner_response_Mc < threshold * Corner_response_Mc.max()] = 0
	# Keep local maxima
	h, w = Corner_response_Mc.shape
	detected_points = Corner_response_Mc.copy()
	for i, j in np.ndindex(detected_points.shape):
		# Don't take corner if the local window is not properly defined
		if i < local or j < local:
			continue
		local_max = Corner_response_Mc[
		max(0, i - local):min(i + local + 1, h),
		max(0, j - local):min(j + local + 1, w)
		].max()
		if detected_points[i, j] == local_max:
			detected_points[
			max(0, i - local):min(i + local + 1, h),
			max(0, j - local):min(j + local, w)].fill(0)
			detected_points[i, j] = local_max
		else:
			detected_points[i, j] = 0
	detected_points[:local, :].fill(0)
	detected_points[-max(1, local):, :].fill(0)
	detected_points[:, :local].fill(0)
	detected_points[:, -max(1, local):].fill(0)
	return detected_points

def find_harris_corners(img, sigma_d, sigma_i, kappa, threshold, local):
	"""
	Finds and returns list of corners and new image with corners drawn
	:param img: The original image
	:param window_size: The size (side length) of the sliding window
	:param k: Harris corner constant. Usually 0.04 - 0.06
	:param thresh: The threshold above which a corner is counted
	:return: refined corners as given by first_refine
	"""
	# Compute gaussian derivatives
	dx = [-0.5, 0, 0.5]
	G = np.convolve(gaussian(max(img.shape), std = sigma_d), dx)
	Gx = G.reshape(1, -1)
	Gy = Gx.T

	# Images derivatives
	Ix = convolve2d(img, Gx, mode = 'same')
	Iy = convolve2d(img, Gy, mode = 'same')

	# Add extra smoothing function
	smooth_x = gaussian(max(img.shape), std = sigma_i).reshape(1 , -1)
	smooth_y = smooth_x.T

	smooth_xy = smooth_y.dot(smooth_x)

	# Product image
	Ixy = Ix * Iy;
	Ixy_smooth = convolve2d(convolve2d(Ixy, smooth_x, mode = "same"), smooth_y, mode = "same");
	Iy_smooth = convolve2d(convolve2d(Iy, smooth_x, mode = "same"), smooth_y, mode = "same")
	Ix_smooth = convolve2d(convolve2d(Ix, smooth_x, mode = "same"), smooth_y, mode = "same")

	# Auto-corellation matrix
	A = np.array([
	[Ix_smooth**2, Ixy_smooth],
	[Ixy_smooth, Iy_smooth**2]
	])
	Corner_response_Mc = (np.linalg.det(A.T) - kappa * np.trace(A.T, axis1 = 2, axis2 = 3)**2).T

	return first_refine(Corner_response_Mc, threshold, local)

def anms(detected_points, best, c):
	"""
	Adaptative non-maximal suppresion algorithm to refine detectors.
	:param detected_points: first refined corners
	:param best: number of corners detection allowed
	:param c: Threshold to validate absolutely a detection
	:return: Second refined corners, by anms method
	"""

	# Matrix of remaining detected points
	D = detected_points.copy()
	# Processed Points
	P = np.argwhere(D == D.max())[0].reshape(1, -1)
	# Value
	F = D.max()
	# Remove P from D
	D[P[-1][0], P[-1][1]] = 0
	## Radius
	r = np.inf
	R = [np.inf]

	# Initialises m : maximal value of remaining detected points
	m = D.max()
	# Keep looking for minimal radius for every detected points
	while (m > 0):
		## Coordinates of current maximal
		p = np.argwhere(D == m)[0].reshape(1, -1)
		## Processed points large enough
		q = P[(c * F > m).reshape(-1, 1)[:, 0]]
		## Minamal radius between current maximal and Processed points
		r = np.linalg.norm(	np.repeat(p, len(q), axis = 0) - q,
							axis = 1).min() if len(q) > 0 else 0
		R.append(r)
		## Add p to processed points
		P = np.vstack((P, p))
		## Remove p from remaining detected points
		D[P[-1][0], P[-1][1]] = 0
		F = np.vstack((F, m))
		m = D.max()

	# Get the first points associated with best radius
	return P[np.argsort(R) < best]

################################################################################

if __name__ == "__main__":
	# Handling arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("image",
						help = "Indicate an image")
	parser.add_argument("-sd",
						"--sigma_d",
						type = int,
						help = "Standard deviation for gaussian kernel used in derivation",
						default = 1)
	parser.add_argument("-si",
						"--sigma_i",
						type = int,
						help = "Standard deviation for gaussian kernel used in windowing",
						default = 1)
	parser.add_argument("-k",
						"--kappa",
						type = float,
						help = "Kappa parameter for Harris Mc computation",
						default = 0.05)
	parser.add_argument("-l",
						"--local",
						type = int,
						help = "Size of local window arround an extremum (2xlocal  + 1) x (2xlocal + 1) pixels)",
						default = 3)
	parser.add_argument("-t",
						"--threshold",
						type = int,
						help = "Threshold for a valid corner detection",
						default = 0.01)
	parser.add_argument("-a",
						"--anms",
						action = "store_true")
	parser.add_argument("-c",
						"--anms_constant",
						type = int,
						help = "Constant in ANMS",
						default = 0.7)
	parser.add_argument("-b",
						"--best",
						type = int,
						help = "Best rankings in ANMS",
						default = 50)
	parser.add_argument("-s",
						"--save",
						type = str,
						help = "Output file for corners locations",
						default = None)
	parser.add_argument("-p",
						"--plot",
						type = bool,
						help = "Plot option",
						default = False)

	args = parser.parse_args()
	print ("Parameters entered :\n", vars(args))

	# Compte the Harris corners
	img = imread(args.image)
	corners = find_harris_corners(
									img,
									args.sigma_d,
									args.sigma_i,
									args.kappa,
									args.threshold,
									args.local
									)

	if args.anms == True:
		# Refine with anms
		refined_corners = anms(corners, args.best, args.anms_constant)
		print ("Number of corners detected before anms : {}".format(len(np.argwhere(corners > 0))))
		print ("Number of corners detected after anms: {}".format(len(refined_corners)))
	else:
		refined_corners = np.argwhere(corners > 0)
		print ("Number of corners detected : {}".format(len(refined_corners)))

	if args.save:
		# Save the corners locations in a txt file
		np.savetxt(args.save,
					refined_corners.astype(int),
					delimiter = '\t',
					fmt = '%d',
					header = "{} Harris corners for "
					"{}".format(len(refined_corners), args.image))

	if args.plot == True:
		# Plot the figure
		plt.figure(figsize = (20, 20))
		imshow(img);
		plt.scatter(refined_corners[:, 1],
				refined_corners[:, 0],
				color = 'red',
				s = 150)
		plt.axis('off')
		plt.show()
