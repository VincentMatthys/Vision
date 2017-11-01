// Student : Vincent Matthys
// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match
{
	float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void				algoSIFT(	Image<Color,2> I1,
								Image<Color,2> I2,
								vector<Match>& matches)
{
	// Find interest points
	SIFTDetector D;
	D.setFirstOctave(-1);
	Array<SIFTDetector::Feature> feats1 = D.run(I1);
	drawFeatures(feats1, Coords<2>(0,0));
	cout << "Im1: " << feats1.size() << flush;
	Array<SIFTDetector::Feature> feats2 = D.run(I2);
	drawFeatures(feats2, Coords<2>(I1.width(), 0));
	cout << " Im2: " << feats2.size() << flush;

	const double MAX_DISTANCE = 100.0 * 100.0;
	for(size_t i = 0; i < feats1.size(); i++)
	{
		SIFTDetector::Feature f1 = feats1[i];
		for(size_t j = 0; j < feats2.size(); j++)
		{
			double d = squaredDist(f1.desc, feats2[j].desc);
			if(d < MAX_DISTANCE)
			{
				Match m;
				m.x1 = f1.pos.x();
				m.y1 = f1.pos.y();
				m.y2 = feats2[j].pos.y();
				m.x2 = feats2[j].pos.x();
				matches.push_back(m);
			}
		}
	}
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3>	computeF(vector<Match>& matches)
{
	const float 			distMax = 1.5f; // Pixel error for inlier/outlier	discrimination
	int Niter = 100000; // Adjusted dynamically
	FMatrix<float,3,3> 		bestF;
	vector<int>				bestInliers;
	// --------------- TODO ------------
	// Vector of current inliers
	vector<Match> 			current_Inliers;
	vector<int>				current_all_inliers;
	// Counter for loops
	size_t					count;
	// Current Fundamental Matrix
	FMatrix<float, 3, 3>	current_F;
	// Matrix of points
	FMatrix<double, 9, 9> 	A;
	// Matrix for SVD decomposition
	FVector<double, 9> S;
	FMatrix<double, 9, 9> U, Vt;
	// Normalization Matrix
	FMatrix<float, 3, 3> 	N;
	// Tools for computing the distance
	FVector<float, 3>		x;
	FVector<float, 3>		x_prime;
	float 					distance;
	N(0, 0) = 0.001; N(0, 1) = 0; N(0, 2) = 0;
	N(1, 0) = 0; N(1, 1) = 0.001; N(1, 2) = 0;
	N(2, 0) = 0; N(2, 1) = 0; N(2, 2) = 1;

	// RANSAC algorithm
	// Do Niter iterations
	count = 0;
	while (count < Niter)
	{
		current_Inliers.clear();

		// Take arbitrary 8 couples (p, p') from matches
		for(size_t i = 0; i < 8; i++)
		{
			current_Inliers.push_back(matches[rand()%(matches.size())]);
			// Fill A with the equation associated to the corerspondance
			// Normalizes by 0.001 for each pixel coordinate
			A(i, 0) = 0.000001 * current_Inliers[i].x1 * current_Inliers[i].x2;
			A(i, 1) = 0.000001 * current_Inliers[i].x1 * current_Inliers[i].y2;
			A(i, 2) = 0.001 * current_Inliers[i].x1;
			A(i, 3) = 0.000001 * current_Inliers[i].y1 * current_Inliers[i].x2;
			A(i, 4) = 0.000001 * current_Inliers[i].y1 * current_Inliers[i].y2;
			A(i, 5) = 0.001 * current_Inliers[i].y1;
			A(i, 6) = 0.001 * current_Inliers[i].x2;
			A(i, 7) = 0.001 * current_Inliers[i].y2;
			A(i, 8) = 1;
			// Fill last line of A with 0 to use square SVD decomposition
			A(8, i) = 0;
		}
		A(8, 8) = 0;

		// Then compute F_tilda for these couples
		svd(A, U, S, Vt);

		current_F(0, 0) = Vt.getRow(7)[0]; current_F(0, 1) = Vt.getRow(7)[1]; current_F(0, 2) = Vt.getRow(7)[2];
		current_F(1, 0) = Vt.getRow(7)[3]; current_F(1, 1) = Vt.getRow(7)[4]; current_F(1, 2) = Vt.getRow(7)[5];
		current_F(2, 0) = Vt.getRow(7)[6]; current_F(2, 1) = Vt.getRow(7)[7]; current_F(2, 2) = Vt.getRow(7)[8];

		current_F = N*current_F*N;

		// Count the number of inliers
		current_all_inliers.clear();
		for(size_t i = 0; i < matches.size(); i++)
		{
			x[0] = matches[i].x1;
			x[1] = matches[i].y1;
			x[2] = 1;
			x_prime[0] = matches[i].x2;
			x_prime[1] = matches[i].y2;
			x_prime[2] = 1;
			x_prime = current_F*x_prime;
			distance = (x*x_prime)*(x*x_prime) / norm2(x_prime);
			// Test for beeing inlier
			if (distance <= 0.00001)
				current_all_inliers.push_back(i);
		}

		// If more outliers than before, keep it, else, repeat the loop
		if (current_all_inliers.size() > bestInliers.size())
		{
			bestInliers = current_all_inliers;
			bestF = current_F;
			// Reestimate the Niter (which is then lower than the previous one)
			Niter = min((int)(std::log(BETA) / std::log(1 - pow(((float)bestInliers.size() / matches.size()), 8))), 10000);
		}
		count++;
	}
	cout << "Number of RANSAC iterations : " << count << endl;

	// Updating matches with inliers only
	vector<Match> all=matches;
	matches.clear();
	for(size_t i = 0; i<bestInliers.size(); i++)
		matches.push_back(all[bestInliers[i]]);

	// Refine resulting F with least square minimization based on all inliers
	// Number of total inliers
	count = matches.size();
	Matrix<double> A_final(count, 8);
	Vector<double> B(count);
	while (count--)
	{
		// Fill A with the equation associated to the corerspondance
		// Normalizes by 0.001 for each pixel coordinate
		A_final(count, 0) = 0.000001 * matches[count].x1 * matches[count].x2;
		A_final(count, 1) = 0.000001 * matches[count].x1 * matches[count].y2;
		A_final(count, 2) = 0.001 * matches[count].x1;
		A_final(count, 3) = 0.000001 * matches[count].y1 * matches[count].x2;
		A_final(count, 4) = 0.000001 * matches[count].y1 * matches[count].y2;
		A_final(count, 5) = 0.001 * matches[count].y1;
		A_final(count, 6) = 0.001 * matches[count].x2;
		A_final(count, 7) = 0.001 * matches[count].y2;
		B[count] = -1;
	}
	B = linSolve(A_final, B);
	bestF(0, 0) = B[0]; bestF(0, 1) = B[1]; bestF(0, 2) = B[2];
	bestF(1, 0) = B[3]; bestF(1, 1) = B[4]; bestF(1, 2) = B[5];
	bestF(2, 0) = B[6]; bestF(2, 1) = B[7]; bestF(2, 2) = 1;

	// Put last singular value to 0 and recompose
	// Matrix for SVD decomposition
	FVector<float, 3> S_f;
	FMatrix<float, 3, 3> U_f, Vt_f;
	svd(bestF, U_f, S_f, Vt_f);
	S_f[2] = 0;
	bestF = U_f*Diagonal(S_f)*Vt_f;
	bestF = N*bestF*N;

	return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void				displayEpipolar(Image<Color> I1,
									Image<Color> I2,
									const FMatrix<float,3,3>& F)
{
	int				w = I1.width();
	int				active_image;
	FVector<float, 3>	u;
	FVector<float, 3>	x1;
	FVector<float, 3>	x2;

	while(true)
	{
		int x,y;
		if (getMouse(x, y) == 3)
			break;
		else
		{
			// --------------- TODO ------------
			active_image = (x < w) ? 0 : 1;

			// u : homogeneous coordinates of clicked point
			u[0] = x - active_image * w;
			u[1] = y;
			u[2] = 1;
			if (active_image == 1)
			// Draw the left epipolar line if clic in I2
			{
				drawCircle(x, y, 3, RED, 2);
				u = F * u;
			}
			// Or draw the right epipolar line if clic in I1
			else
			{
				drawCircle(x, y, 3, YELLOW, 2);
				u = transpose(F) * u;
			}
			// Find intersection points with images edges
			// x1 : left edge
			// x2 : right edge
			x1[0] = -u[2];
			x1[1] = 0;
			x1[2] = u[0];
			x2[0] = -u[1] * w - u[2];
			x2[1] = u[0] * w;
			x2[2] = u[0];
			x1 = x1 / x1[2];
			x2 = x2 / x2[2];

			// Draw the corersponding epipolar line
			// With the w shift for the right epipolar line
			if (active_image == 0)
				drawLine(x1[0] + w, x1[1], x2[0] + w, x2[1], YELLOW);
			else
				drawLine(x1[0], x1[1], x2[0], x2[1], RED);
		}
	}
}

int 				main(int argc, char* argv[])
{
	srand((unsigned int)time(0));

	const char* s1 = argc > 1 ? argv[1] : srcPath("im1.jpg");
	const char* s2 = argc > 2 ? argv[2] : srcPath("im2.jpg");

	// Load and display images
	Image<Color,2> I1, I2;
	if(!load(I1, s1)
	|| !load(I2, s2))
	{
		cerr << "Unable to load images" << endl;
		return 1;
	}
	int w = I1.width();
	openWindow(2 * w, I1.height());
	display(I1, 0, 0);
	display(I2, w, 0);

	vector<Match> matches;
	algoSIFT(I1, I2, matches);
	cout << " matches: " << matches.size() << endl;
	// Waits for a mouse click in active window
	click();


	FMatrix<float,3,3> F = computeF(matches);
	cout << "F = "<< endl << F;

	// Redisplay with matches
	display(I1, 0, 0);
	display(I2, w, 0);
	// Number of inliers matches
	cout << "Number of inliers matches : " << matches.size() << endl;
	for(size_t i = 0; i<matches.size(); i++)
	{
		Color c(rand()%256, rand()%256, rand()%256);
		fillCircle(matches[i].x1 + 0, matches[i].y1, 2, c);
		fillCircle(matches[i].x2 + w, matches[i].y2, 2, c);
	}
	click();

	// Redisplay without SIFT points
	display(I1, 0, 0);
	display(I2, w, 0);
	displayEpipolar(I1, I2, F);

	// Wait for user's right click
	cout << "\033[1;32m";
	endGraphics();
	cout << "\033[0m" << endl;
	return (0);
}
