#ifndef DRAWLINE_H
#define DRAWLINE_H
void drawDot(int x, int y, float r, float g, float b);
void drawLine(int x0, int x1, int y0, int y1);
void drawLine1(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg);
void drawLine2(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg);
void drawLine3(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg);
void drawLine4(int x0, int x1, int y0, int y1, bool xy_interchange, bool neg);
void changeColor(int color);
#endif