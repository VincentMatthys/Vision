// Imagine++ project
// Project:  Seeds
// Author:   Pascal Monasse
// Student: Vincent Matthys

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <queue>
#include <iostream>
using namespace Imagine;
using namespace std;

// For debugging : use cropped image
// #define DEBUG

/// Min and max disparities
static const float dMin = -30, dMax = -7;

/// Min NCC for a seed
static const float nccSeed = 0.95;

/// Radius of patch for correlation
static const int win = (9-1) / 2;
/// To avoid division by 0 for constant patch
static const float EPS = 0.1f;

/// A seed
struct Seed
{
	Seed(int x0, int y0, int d0, float ncc0)
	: x(x0), y(y0), d(d0), ncc(ncc0) {}
	int x,y, d;
	float ncc;
};

/// Order by NCC
bool operator<(const Seed& s1, const Seed& s2)
{
	return (s1.ncc<s2.ncc);
}

/// 4-neighbors
static const int dx[] = {+1,  0, -1,  0};
static const int dy[] = { 0, -1,  0, +1};

/// Display disparity map
static void displayDisp(const Image<int> disp, Window W, int subW)
{
	Image<Color> im(disp.width(), disp.height());
	for(int j = 0; j < disp.height(); j++)
		for(int i = 0; i < disp.width(); i++)
		{
			if(disp(i,j) < dMin || disp(i,j) > dMax)
				im(i,j) = Color(0,255,255);
			else
			{
				int g = 255 * (disp(i,j) - dMin) / (dMax - dMin);
				im(i,j) = Color(g,g,g);
			}
		}
		setActiveWindow(W,subW);
		display(im);
		showWindow(W,subW);
}

/// Show 3D window
static void show3D(const Image<Color> im, const Image<int> disp)
{
#ifdef IMAGINE_OPENGL // Imagine++ must have been built with OpenGL support...
	// Intrinsic parameters given by Middlebury website
	const float f = 3740;
	const float d0 = -200; // Doll images cropped by this amount
	const float zoom = 2; // Half-size images, should double measured disparity
	const float B = 0.160; // Baseline in m
	FMatrix<float,3,3> K(0.0f);
	K(0, 0) = -f / zoom; K(0,2) = disp.width() / 2;
	K(1, 1) = f / zoom; K(1,2) = disp.height() / 2;
	K(2, 2) = 1.0f;
	K = inverse(K);
	K /= K(2,2);
	std::vector<FloatPoint3> pts;
	std::vector<Color> col;
	for(int j = 0; j < disp.height(); j++)
		for(int i = 0; i < disp.width(); i++)
			if(dMin <= disp(i,j) && disp(i,j) <= dMax)
			{
				float z = B * f / (zoom * disp(i, j) + d0);
				FloatPoint3 pt((float)i, (float)j, 1.0f);
				pts.push_back(K * pt * z);
				col.push_back(im(i,j));
			}
			Mesh mesh(&pts[0], pts.size(), 0, 0, 0, 0, VERTEX_COLOR);
			mesh.setColors(VERTEX, &col[0]);
			Window W = openWindow3D(512, 512, "3D");
			setActiveWindow(W);
			showMesh(mesh);
#else
	std::cout << "No 3D: Imagine++ not built with OpenGL support" << std::endl;
#endif
}

/// Correlation between patches centered on (i1,j1) and (i2,j2). The values
/// m1 or m2 are subtracted from each pixel value.
static float correl(const Image<byte>& im1, int i1,int j1,float m1,
					const Image<byte>& im2, int i2,int j2,float m2)
{
	float dist = 0.0f;
	// ------------- TODO -------------
	size_t			width1;
	size_t			width2;
	size_t			height1;
	size_t			height2;
	int				k;
	int				l;
	unsigned int	tmp1;
	unsigned int	tmp2;

	width1 = im1.width();
	width2 = im2.width();
	height1 = im1.height();
	height2 = im2.height();
	l = -win;
	tmp1 = 0;
	tmp2 = 0;
	while (i1 + l >= 0 && i1 + l <= height1
		&& i2 + l >= 0 && i2 + l <= height2
		&& l < win)
	{
		k = -win;
		while (j1 + k >= 0 && j1 + k <= width1
			&& j2 + k >= 0 && j2 + k <= width2
			&& k < win)
		{
			// Add the value of intensity at (i + l, j + k)
			dist += (im1(i1 + l, j1 + k) - m1) * (im2(i2 + l, j2 + k) - m2);
			tmp1 += pow(im1(i1 + l, j1 + k) - m1, 2);
			tmp2 += pow(im2(i2 + l, j2 + k) - m2, 2);
			k++;
		}
		l++;
	}
	dist /= (sqrt(tmp1) * sqrt(tmp2) + EPS);
	return dist;
}

/// Sum of pixel values in patch centered on (i,j).
static float sum(const Image<byte>& im, int i, int j)
{
	float	s = 0.0f;
	// ------------- TODO -------------
	size_t	width;
	size_t	height;
	int		k;
	int		l;

	width = im.width();
	height = im.height();
	l = -win;
	while (i + l >= 0 && i + l < height && l < win)
	{
		k = -win;
		while (j + k >= 0 && j + k < width && k < win)
		{
			// Add the value of intensity at (i + l, j + k)
			s += im(j + k, i + l);
			k++;
		}
		l++;
	}
	return s;
}

/// Centered correlation of patches of size 2*win+1.
static float ccorrel(	const Image<byte>& im1,int i1,int j1,
						const Image<byte>& im2,int i2,int j2)
{
	float m1 = sum(im1, i1, j1);
	float m2 = sum(im2, i2, j2);
	int w = 2 * win + 1;
	return correl(im1, i1, j1, m1 / (w * w), im2, i2, j2, m2 / (w * w));
}

/// Compute disparity map from im1 to im2, but only at points where NCC is
/// above nccSeed. Set to true the seeds and put them in Q.
static void find_seeds(	Image<byte> im1, Image<byte> im2,
						float nccSeed,
						Image<int>& disp, Image<bool>& seeds,
						std::priority_queue<Seed>& Q)
{
	disp.fill(dMin - 1);
	seeds.fill(false);
	while(!Q.empty())
		Q.pop();

	for(int y = win; y + win < im1.height() && y + win < im2.height(); y++)
		for(int x = win; x + win < im1.width(); x++)
		{
			// ------------- TODO -------------
			// Hint: just ignore windows that are not fully in image
			float	tmp;
			float	cur_best;
			int		d;
			int		d_best;

			d = dMin;
			d_best = dMin;
			cur_best = -1;
			while (d < dMax)
			{
				// If window is fully in image
				if (x + d - win >= 0 && x + d + win < im2.width())
				{
					tmp = ccorrel(im1, x, y, im2, x + d, y);
					if (tmp >= cur_best)
					{
						cur_best = tmp;
						d_best = d;
					}
				}
				d++;
			}
			if (cur_best > nccSeed)
			{
				disp(x, y) = d_best;
				seeds(x, y) = true;
				Seed s(x, y, d_best, cur_best);
				Q.push(s);
			}
		}
}

/// Propagate seeds
static void propagate(	Image<byte> im1, Image<byte> im2,
						Image<int>& disp, Image<bool>& seeds,
						std::priority_queue<Seed>& Q)
{
	while(!Q.empty())
	{
		Seed s = Q.top();
		Q.pop();
		// For the Neighbors
		for(int i = 0; i < 4; i++)
		{
			float x = s.x + dx[i], y = s.y + dy[i];
			if(0 <= x - win && 0 <= y - win
			&& x + win < im2.width() && y + win < im2.height()
			&& !seeds(x,y))
			{
				// ------------- TODO -------------
				float	temp;
				Seed	s_new(x, y, s.d, -1);
				temp = -1;
				for (int n = 0; n < 3; n++)
				{
					temp = ccorrel(im1, x, y, im2, x + s.d - 1 + n, y);
					if (temp >= s_new.ncc)
					{
						s_new.ncc = temp;
						s_new.d = s.d - 1 + n;
					}
				}
				disp(x, y) = s_new.d;
				seeds(x,y) = true;
				Q.push(s_new);
			}
		}
	}
}

int main(void)
{
	// Load and display images
	Image<Color> I1, I2;

#ifdef DEBUG
	if(!load(I1, srcPath("crop_1.jpg"))
	|| !load(I2, srcPath("crop_2.jpg")))
	{
		cerr<< "Unable to load images" << endl;
		return (1);
	}
#else
	if(!load(I1, srcPath("im1.jpg"))
	|| !load(I2, srcPath("im2.jpg")))
	{
		cerr<< "Unable to load images" << endl;
		return (1);
	}
#endif
	std::string names[5] = { "image 1" , "image 2" , "dense" , "seeds" , "propagation"};
	Window W = openComplexWindow(I1.width(), I1.height(),
								"Seeds propagation", 5, names);
	setActiveWindow(W, 0);
	display(I1, 0, 0);
	// display(I1,0,0);
	setActiveWindow(W, 1);
	display(I2, 0, 0);

	Image<int> disp(I1.width(), I1.height());
	Image<bool> seeds(I1.width(), I1.height());
	std::priority_queue<Seed> Q;

	// Dense disparity
	find_seeds(I1, I2, -1.0f, disp, seeds, Q);
	displayDisp(disp,W,2);

	// Only seeds
	find_seeds(I1, I2, nccSeed, disp, seeds, Q);
	displayDisp(disp,W,3);

	// Propagation of seeds
	propagate(I1, I2, disp, seeds, Q);
	displayDisp(disp,W,4);

	// Show 3D (use shift click to animate)
	show3D(I1,disp);

	// Wait for user's right click
	cout << "\033[1;32m";
	endGraphics();
	cout << "\033[0m" << endl;
	return 0;
}
