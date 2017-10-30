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
	// Current Fundamental Matrix
	FMatrix<float, 3, 3>	current_F;
	// Matrix of points
	FMatrix<double, 9, 9> 	A;
	// Normalization Matrix
	FMatrix<float, 3, 3> 	N;
	// Tools for computing the distance
	FVector<float, 3>		x;
	FVector<float, 3>		x_prime;
	float 					distance;
	N(0, 0) = 0.001; N(0, 1) = 0; N(0, 2) = 0;
	N(1, 0) = 0; N(1, 1) = 0.001; N(1, 2) = 0;
	N(2, 0) = 0; N(2, 1) = 0; N(2, 2) = 1;
	// Do Niter iterations
	while (Niter--)
	{
		current_Inliers.clear();
		// Take arbitrary 8 couples (p, p') from matches
		for(size_t i = 0; i < 4; i++)
		{
			current_Inliers.push_back(matches[rand()%(matches.size())]);
			// Fill A with the equation associated to the corerspondance
			for (size_t j = 0; j < 2; j++)
			{
				// Normalizes by 0.001 for each pixel coordinate
				A(2 * i + j, 0) = 0.000001 * current_Inliers[i].x1 * current_Inliers[i].x2;
				A(2 * i + j, 1) = 0.000001 * current_Inliers[i].x1 * current_Inliers[i].y2;
				A(2 * i + j, 2) = 0.001 * current_Inliers[i].x1;
				A(2 * i + j, 3) = 0.000001 * current_Inliers[i].y1 * current_Inliers[i].x2;
				A(2 * i + j, 4) = 0.000001 * current_Inliers[i].y1 * current_Inliers[i].y2;
				A(2 * i + j, 5) = 0.001 * current_Inliers[i].y1;
				A(2 * i + j, 6) = 0.001 * current_Inliers[i].x2;
				A(2 * i + j, 7) = 0.001 * current_Inliers[i].y2;
				A(2 * i + j, 8) = 1;
				// Fill last line of A with 0 to use square SVD decomposition
				A(8, 2 * i + j) = 0;
			}
			A(8, 8) = 0;
		}
		// Then compute F_tilda for these couples
		FVector<double, 9> S;
		FMatrix<double, 9, 9> U, Vt;
		svd(A, U, S, Vt);
		// cout << "SVD check 1: "  << norm(A-U*Diagonal(S)*Vt) << endl;
		// Get the 8th right singular vectors, i.e : the 8th line of Vt.
		// cout << Vt.getRow(7) << endl;
		current_F(0, 0) = Vt.getRow(7)[0]; current_F(0, 1) = Vt.getRow(7)[1]; current_F(0, 2) = Vt.getRow(7)[2];
		current_F(1, 0) = Vt.getRow(7)[3]; current_F(1, 1) = Vt.getRow(7)[4]; current_F(1, 2) = Vt.getRow(7)[5];
		current_F(2, 0) = Vt.getRow(7)[6]; current_F(2, 1) = Vt.getRow(7)[7]; current_F(2, 2) = Vt.getRow(7)[8];

		// cout << "Before renormalization" << endl;
		// cout << current_F << endl;
		current_F = N*current_F*N;
		// cout << current_F << endl;
		cout << "F=\033[1;31m\n" << current_F << endl;
		cout << "\033[0m" << endl;

		// Count the number of inliers
		current_all_inliers.clear();
		for(size_t i = 0; i < matches.size(); i++)
		{
			// x = {matches[i].x1, matches[i].y1};
			// x_prime = {matches[i].x2, matches[i].y2};
			x[0] = matches[i].x1;
			x[1] = matches[i].y1;
			x[2] = 1;
			x_prime[0] = matches[i].x2;
			x_prime[1] = matches[i].y2;
			x_prime[2] = 1;
			x_prime = current_F*x_prime;
			distance = (x*x_prime)*(x*x_prime) / norm2(x_prime);
			// cout << "Distance" << distance << endl;
			if (distance <= 0.005)
				current_all_inliers.push_back(i);
		}
		cout <<  "Number of inliers at Niter :";
		cout << Niter << "      ->     \033[1;34m"	;
		cout << current_all_inliers.size() << endl;
		cout << "\033[0m" << endl;
		// If more outliers than before, keep it, else, repeat the loop
		if (current_all_inliers.size() > bestInliers.size())
		{
			bestInliers = current_all_inliers;
			Niter = min((int)(std::log(BETA) / std::log(1 - (float)pow(((float)bestInliers.size() / matches.size()), 8))), Niter);
			// cout << "New Niter  --- " << Niter << endl;
			// cout << std::log(1 - pow(((float)bestInliers.size() / matches.size()), 8)) << endl ;
		}
		// Reestimate the Niter (which is then lower than the previous one)
	}

	// DO NOT FORGET NORMALIZATION OF POINTS



	// Updating matches with inliers only
	vector<Match> all=matches;
	matches.clear();
	for(size_t i = 0; i<bestInliers.size(); i++)
		matches.push_back(all[bestInliers[i]]);
	return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void				displayEpipolar(Image<Color> I1,
									Image<Color> I2,
									const FMatrix<float,3,3>& F)
{
	while(true)
	{
		int x,y;
		if(getMouse(x,y) == 3)
			break;
		// --------------- TODO ------------
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
	openWindow(2*w, I1.height());
	display(I1, 0, 0);
	display(I2, w, 0);

	vector<Match> matches;
	algoSIFT(I1, I2, matches);
	cout << " matches: " << matches.size() << endl;
	// Waits for a mouse click in active window
	click();


	FMatrix<float,3,3> F = computeF(matches);
	// cout << "F = "<< endl << F;

	// // Redisplay with matches
	// display(I1,0,0);
	// display(I2,w,0);
	// for(size_t i = 0; i<matches.size(); i++)
	// {
	// 	Color c(rand()%256, rand()%256, rand()%256);
	// 	fillCircle(matches[i].x1 + 0, matches[i].y1, 2, c);
	// 	fillCircle(matches[i].x2 + w, matches[i].y2, 2, c);
	// }
	// click();
	//
	// // Redisplay without SIFT points
	// display(I1,0,0);
	// display(I2,w,0);
	// displayEpipolar(I1, I2, F);
	//

	// Wait for user's right click
	cout << "\033[1;32m";
	endGraphics();
	cout << "\033[0m" << endl;
	return (0);
}
