// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
// Vector in std : http://fr.cppreference.com/w/cpp/container/vector
void			getClicks(	Window w1,
							Window w2,
							vector<IntPoint2>& pts1,
							vector<IntPoint2>& pts2)
{
	// ------------- TODO/A completer ----------
	IntPoint2		p;
	unsigned int	button;
	int				window;

	cout << "Select similar points in both windows." << endl;
	cout << "Starting with the first window -- " << endl;
	window = 1;
	setActiveWindow(w1);
	while ((button = getMouse(p)) != 3)
	{
		drawCircle(p, 5, RED, 2);
		cout << "Pixel in window " << window << " in position " << p << endl;
		if (window == 1)
		{
			window = 2;
			pts1.push_back(p);
			setActiveWindow(w2);
		}
		else
		{
			window = 1;
			pts2.push_back(p);
			setActiveWindow(w1);
		}
	}
	cout << "Number of points selected : " << pts1.size() << endl;
	if (pts2.size() < 4)
		cerr << "Wrong number of points selected" << endl;
	else if (pts1.size() != pts2.size())
		cerr << "Missing the last similar point in window 2" << endl;
}

// Return homography compatible with point matches
Matrix<float>	getHomography(	const vector<IntPoint2>& pts1,
								const vector<IntPoint2>& pts2)
{
	size_t	n = min(pts1.size(), pts2.size());
	size_t	i;

	if (n < 4)
	{
		cout << "Not enough correspondences: " << n << endl;
		return Matrix<float>::Identity(3);
	}
	Matrix<double> A(2 * n, 8);
	Vector<double> B(2 * n);
    // ------------- TODO/A completer ----------
	i = 0;
	while (i < n)
	{
		A(2 * i, 0) = pts1.at(i)[0];
		A(2 * i, 1) = pts1.at(i)[1];
		A(2 * i, 2) = 1;
		A(2 * i, 3) = 0;
		A(2 * i, 4) = 0;
		A(2 * i, 5) = 0;
		A(2 * i, 6) = - pts1.at(i)[0] * pts2.at(i)[0];
		A(2 * i, 7) = - pts1.at(i)[0] * pts2.at(i)[1];
		A(2 * i + 1, 0) = 0;
		A(2 * i + 1, 1) = 0;
		A(2 * i + 1, 2) = 0;
		A(2 * i + 1, 3) = pts1.at(i)[0];
		A(2 * i + 1, 4) = pts1.at(i)[1];
		A(2 * i + 1, 5) = 1;
		A(2 * i + 1, 6) = - pts1.at(i)[0] * pts2.at(i)[1];
		A(2 * i + 1, 7) = - pts1.at(i)[1] * pts2.at(i)[1];
		B[2 * i] = pts2.at(i)[0];
		B[2 * i + 1] = pts2.at(i)[1];
		i++;
	}

	B = linSolve(A, B);
	Matrix<float> H(3, 3);
	H(0, 0) = B[0]; H(0, 1) = B[1]; H(0, 2) = B[2];
	H(1, 0) = B[3]; H(1, 1) = B[4]; H(1, 2) = B[5];
	H(2, 0) = B[6]; H(2, 1) = B[7]; H(2, 2) = 1;

	// Sanity check
	for(size_t i = 0; i < n; i++)
	{
		float v1[] = {(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
		float v2[] = {(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
		Vector<float> x1(v1,3);
		Vector<float> x2(v2,3);
		x1 = H * x1;
		cout	<< x1[1] * x2[2] - x1[2] * x2[1] << ' '
				<< x1[2] * x2[0] - x1[0] * x2[2] << ' '
				<< x1[0] * x2[1] - x1[1] * x2[0] << endl;
	}
	return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void			growTo(	float& x0,
						float& y0,
						float& x1,
						float& y1,
						float x,
						float y)
{
	if(x < x0)
		x0 = x;
	if(x > x1)
		x1 = x;
	if(y < y0)
		y0 = y;
	if(y > y1)
		y1 = y;
}

// Panorama construction
void			panorama(	const Image<Color,2>& I1,
							const Image<Color,2>& I2,
							Matrix<float> H)
{
	Vector<float>	v(3);
	Matrix<float>	H_inv;
	float			x0 = 0;
	float			y0 = 0;
	float			x1 = I2.width();
	float			y1 = I2.height();
	size_t			x;
	size_t			y;

	v[0] = 0; v[1] = 0; v[2] = 1;
	v = H * v; v /= v[2];
	growTo(x0, y0, x1, y1, v[0], v[1]);

	v[0] = I1.width(); v[1] = 0; v[2] = 1;
	v = H * v; v/= v[2];
	growTo(x0, y0, x1, y1, v[0], v[1]);

	v[0] = I1.width(); v[1] = I1.height(); v[2]=1;
	v = H * v; v/=v[2];
	growTo(x0, y0, x1, y1, v[0], v[1]);

	v[0] = 0; v[1] = I1.height(); v[2] = 1;
	v = H * v; v/= v[2];
	growTo(x0, y0, x1, y1, v[0], v[1]);

	cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

	Image<Color> I(int(x1 - x0), int(y1 - y0));
	setActiveWindow(openWindow(I.width(), I.height()));
	I.fill(WHITE);
	// ------------- TODO/A completer ----------
	x = 0;
	H_inv = inverse(H);
	while (x < (size_t)I.width())
	{
		y = 0;
		while (y < (size_t)I.height())
		{
			v[0] = x;
			v[1] = y;
			v[2] = 1;
			v = H_inv * v;
			// Checking if the preimage of v lies in I1
			if (0 <= v[0] && v[0] < I1.width()
			&& 0 <= v[1] && v[1] < I1.height())
			{
				// If the preimage belongs to I1, then pull the I1 value to I
				I(x, y) = I1(int(v[0]), int(v[1]));
				// cout << I(0, 0	) << endl;
				// cout << v[0] << "  " << v[1] << endl;
				// cout << I1(int(v[0]), int(v[1]))[0] << endl;
			}
			if (x < I2.width() && y < I2.height())
				I(x, y) = (I(x,y) + I2(x , y)) / 2;
			y++;
		}
		x++;
	}
	display(I,0,0);
}

// Main function
int				main(int argc, char** argv)
{
	const char* s1 = argc>1 ? argv[1] : srcPath("image0006.jpg");
	const char* s2 = argc>2 ? argv[2] : srcPath("image0007.jpg");

	// Load and display images
	Image<Color> I1, I2;
	if(!load(I1, s1) || !load(I2, s2))
	{
		cerr<< "Unable to load the images" << endl;
		return 1;
	}
	Window w1 = openWindow(I1.width(), I1.height(), s1);
	display(I1, 0, 0);
	Window w2 = openWindow(I2.width(), I2.height(), s2);
	setActiveWindow(w2);
	display(I2, 0, 0);

	// Get user's clicks in images
	vector<IntPoint2> pts1, pts2;
	// getClicks(w1, w2, pts1, pts2);
	IntPoint2 p;
	// p[0] = 533; p[1] = 162; pts1.push_back(p);
	// p[0] = 538; p[1] = 337; pts1.push_back(p);
	// p[0] = 722; p[1] = 418; pts1.push_back(p);
	// p[0] = 654; p[1] = 108; pts1.push_back(p);
	//
	// p[0] = 77; p[1] = 162; pts2.push_back(p);
	// p[0] = 78; p[1] = 332; pts2.push_back(p);
	// p[0] = 263; p[1] = 414; pts2.push_back(p);
	// p[0] = 198; p[1] = 105; pts2.push_back(p);

	// Mario points
	p[0] = 603; p[1] = 427; pts1.push_back(p);
	p[0] = 534; p[1] = 145; pts1.push_back(p);
	p[0] = 746; p[1] = 113; pts1.push_back(p);
	p[0] = 745; p[1] = 359; pts1.push_back(p);

	p[0] = 149; p[1] = 422; pts2.push_back(p);
	p[0] = 77; p[1] = 146; pts2.push_back(p);
	p[0] = 292; p[1] = 117; pts2.push_back(p);
	p[0] = 287; p[1] = 360; pts2.push_back(p);
	vector<IntPoint2>::const_iterator it;
	cout << "pts1=" << endl;
	for(it = pts1.begin(); it != pts1.end(); it++)
		cout << *it << endl;
	cout << "pts2=" << endl;
	for(it = pts2.begin(); it != pts2.end(); it++)
		cout << *it << endl;

	// Compute homography
	Matrix<float> H = getHomography(pts1, pts2);
	cout << "H=\033[1;31m" << H/H(2,2);
	cout << "\033[0m" << endl;

	// Apply homography
	panorama(I1, I2, H);

	endGraphics();
	return (0);
}
