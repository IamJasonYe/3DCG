#include "stdafx.h"
#include "drawline.h"
#include <iostream>
#include <gl/glut.h>
#include <cstdlib>
#include <windows.h>
#include <iomanip>
using namespace std;
void drawLine(int x0, int x1, int y0, int y1)
{
	// cout << x0 <<" "<< y0 <<" " << x1 <<" " << y1;
	// system("pause");
	int dx, dy;
	dx = x1 - x0;
	dy = y1 - y0;
	int neg = dx * dy;
	if (dx > 0 && dy > 0)
		drawLine1(x0, x1, y0, y1, dy >= dx, neg < 0);
	else if (dx > 0 && dy < 0)
		drawLine2(x0, x1, y0, y1, dy < -dx, neg < 0);
	else if (dx < 0 && dy > 0)
		drawLine3(x0, x1, y0, y1, dy > -dx, neg < 0);
	else if (dx < 0 && dy < 0)
		drawLine4(x0, x1, y0, y1, dy <= dx, neg < 0);
	else	// horizontal or vertical line
	{
		if (dx == 0)
		{
			if (y1<y0)
				swap(y0, y1);
			for (int i = y0; i <= y1; i++)
				drawDot(x0, i, 1, 0.5, 0.5);
			//glFlush();
		}
		else if (dy == 0)
		{
			if (x1 < x0)
				swap(x0, x1);
			for (int i = x0; i <= x1; i++)
				drawDot(i, y0, 1, 0.5, 0.5);
			//glFlush();
		}
	}
}

// Draw line for dx>0 and dy>0
void drawLine1(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg)
{
	if (neg)
	{
		if (xy_interchange)
		{
			int tmpx0 = x0;
			int tmpx1 = x1;
			x0 = -y0;
			x1 = -y1;
			y0 = tmpx0;
			y1 = tmpx1;
		}
		else
		{
			y0 = -y0;
			y1 = -y1;
		}
	}
	else
	{
		if (xy_interchange)
		{
			swap(x0, y0);
			swap(x1, y1);
		}
	}
	int x = x0,
		y = y0;
	int a = y1 - y0,
		b = -(x1 - x0);
	int d = 2 * a + b;
	int IncE = 2 * a,
		IncNE = 2 * (a + b);

	while (x <= x1)
	{
		x++;
		if (d <= 0)
			d += IncE;
		else
		{
			y++;
			d += IncNE;
		}
		// line_data* tmp;
		if (neg)
		{
			if (xy_interchange)
				drawDot(y, -x, 1.0, 0.5, 0.5);
			else
				drawDot(x, -y, 1.0, 0.5, 0.5);
		}
		else
		{
			if (xy_interchange)
				drawDot(y, x, 1.0, 0.5, 0.5);
			else
				drawDot(x, y, 1.0, 0.5, 0.5);
		}
	}
	//glFlush();
}

// Draw line for dx>0 and dy<0
void drawLine2(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg)
{
	drawLine1(x0, x1, y0, y1, xy_interchange, neg);
}

// Draw line for dx<0 and dy>0
void drawLine3(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg)
{
	drawLine2(x1, x0, y1, y0, xy_interchange, neg);
}

// Draw line for dx<0 and dy>0
void drawLine4(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg)
{
	drawLine1(x1, x0, y1, y0, xy_interchange, neg);
}

void drawDot(int x, int y, float r, float g, float b)
{
	// cout << "drawDot()" << endl;
	glBegin(GL_POINTS);
	// set the color of dot
	glColor3f(r, g, b);	//color control by changeColor(R,G,B)
	// invert height because the opengl origin is at top-left instead of bottom-left
	glVertex2i(x, y);	//Here is different from Hw1 "height - y"
	glEnd();
}