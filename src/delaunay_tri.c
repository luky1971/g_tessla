/*
 * Copyright Ahnaf Siddiqui and Sameer Varma.
 *
 * Delaunay triangulation implementation based on the algorithms presented by
 *
 * 	Der-Tsai Lee and Bruce J. Schachter, 
 	"Two Algorithms for Constructing the Delaunay Triangulation," 
	International Journal of Computer and Information Science 9(3):219-242, 1980
 *
 * and
 *
 * 	Leonidas J. Guibas and Jorge Stolfi, 
	"Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams," 
	ACM Transactions on Graphics 4(2):74-123, April 1985
 *
 */

#include "delaunay_tri.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define DTEPSILON 1e-8 // TODO: explore different epsilon values


struct vert {
	int index;
	struct vertNode *adj;
};

struct vertNode {
	struct vert *v;
	struct vertNode *prev, *next;
};


static dtreal XX(const struct vert *v);
static dtreal YY(const struct vert *v);

int compareVerts(const void *a, const void *b);

static bool ccw(const struct vert *a,
				const struct vert *b,
				const struct vert *c);
static bool rightOf(const struct vert *x, 
					const struct vert *ea, 
					const struct vert *eb);
static bool leftOf(	const struct vert *x,
					const struct vert *ea,
					const struct vert *eb);
static bool inCircle(	const struct vert *a, 
						const struct vert *b, 
						const struct vert *c, 
						const struct vert *d);

static struct vert *first(const struct vert *vi);
static struct vert *pred(	const struct vert *vi, 
							const struct vert *vj);
static struct vert *succ(	const struct vert *vi, 
							const struct vert *vj);
static void lct(const struct vert *lrightmost, 
				const struct vert *rleftmost, 
				struct vert **lctleft, 
				struct vert **lctright);

static void uct(const struct vert *lrightmost, 
				const struct vert *rleftmost, 
				struct vert **uctleft, 
				struct vert **uctright);

static void ord_dtriangulate(	struct vert *v, 
								int ia, 
								int ib, 
								struct vert **leftmost,
								struct vert **rightmost);


static struct dTriangulation *mtri;


static inline dtreal XX(const struct vert *v) {
	return mtri->points[2 * v->index];
}

static inline dtreal YY(const struct vert *v) {
	return mtri->points[2 * v->index + 1];
}


// lexicographically compares the coordinates of the given vertices
int compareVerts(const void *a, const void *b) {
	dtreal diff = XX((const struct vert*)a) - XX((const struct vert*)b);

	if(diff < DTEPSILON && diff > -DTEPSILON) {
		diff = YY((const struct vert*)a) - YY((const struct vert*)b);
	}

	return (1 - 2 * signbit(diff));
}


static bool ccw(const struct vert *a,
				const struct vert *b,
				const struct vert *c) {
	dtreal xa = XX(a);
	dtreal ya = YY(a);
	dtreal xb = XX(b);
	dtreal yb = YY(b);
	dtreal xc = XX(c);
	dtreal yc = YY(c);

	return (xa * (yb - yc) - ya * (xb - xc) + xb * yc - yb * xc) > 0;
}

static bool rightOf(const struct vert *x, 
					const struct vert *ea, 
					const struct vert *eb) {
	return ccw(x, eb, ea);
}

static bool leftOf(	const struct vert *x,
					const struct vert *ea,
					const struct vert *eb) {
	return ccw(x, ea, eb);
}

// Acknowledgments:
// Baker, M. (2015). Maths - Matrix algebra - Determinants 4D.
// http://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/fourD/index.htm
static bool inCircle(	const struct vert *a, 
						const struct vert *b, 
						const struct vert *c, 
						const struct vert *d) {
	dtreal m[4][3] = {
		{XX(a), YY(a)},
		{XX(b), YY(b)},
		{XX(c), YY(c)},
		{XX(d), YY(d)}
	};
	m[0][2] = m[0][0] * m[0][0] + m[0][1] * m[0][1];
	m[1][2] = m[1][0] * m[1][0] + m[1][1] * m[1][1];
	m[2][2] = m[2][0] * m[2][0] + m[2][1] * m[2][1];
	m[3][2] = m[3][0] * m[3][0] + m[3][1] * m[3][1];

	return (
	m[1][2] * m[2][1] * m[3][0] - m[0][2] * m[2][1] * m[3][0] -
	m[1][1] * m[2][2] * m[3][0] + m[0][1] * m[2][2] * m[3][0] +
	m[0][2] * m[1][1] * m[3][0] - m[0][1] * m[1][2] * m[3][0] -
	m[1][2] * m[2][0] * m[3][1] + m[0][2] * m[2][0] * m[3][1] +
	m[1][0] * m[2][2] * m[3][1] - m[0][0] * m[2][2] * m[3][1] -
	m[0][2] * m[1][0] * m[3][1] + m[0][0] * m[1][2] * m[3][1] +
	m[1][1] * m[2][0] * m[3][2] - m[0][1] * m[2][0] * m[3][2] -
	m[1][0] * m[2][1] * m[3][2] + m[0][0] * m[2][1] * m[3][2] +
	m[0][1] * m[1][0] * m[3][2] - m[0][0] * m[1][1] * m[3][2] -
	m[0][2] * m[1][1] * m[2][0] + m[0][1] * m[1][2] * m[2][0] +
	m[0][2] * m[1][0] * m[2][1] - m[0][0] * m[1][2] * m[2][1] -
	m[0][1] * m[1][0] * m[2][2] + m[0][0] * m[1][1] * m[2][2]) > 0;
}

static inline struct vert *first(const struct vert *vi) {
	if(vi && vi->adj)
		return vi->adj->v;
	return NULL;
}

static struct vert *pred(	const struct vert *vi, 
							const struct vert *vj) {
	if(vi && vj) {
		struct vertNode *vn = vi->adj;
		if(vn) {
			do {
				if(vn->v == vj) {
					if(vn->prev)
						return vn->prev->v;
					return NULL;
				}
				vn = vn->prev;
			} while(vn && vn != vi->adj);
		}
	}
	else {
		fprintf(stderr, 
			"Hey silly, you input null vertices in predecessor function!\n");
	}
	return NULL;
}

static struct vert *succ(	const struct vert *vi, 
							const struct vert *vj) {
	if(vi && vj) {
		struct vertNode *vn = vi->adj;
		if(vn) {
			do {
				if(vn->v == vj) {
					if(vn->next)
						return vn->next->v;
					return NULL;
				}
				vn = vn->next;
			} while(vn && vn != vi->adj);
		}
	}
	else {
		fprintf(stderr, 
			"Hey silly, you input null vertices in successor function!\n");
	}
	return NULL;
}

// Lower common tangent of two convex hulls
// TODO: what happens in edge cases? (ex. given hulls have 2 or less points)
static void lct(const struct vert *lrightmost, 
				const struct vert *rleftmost, 
				struct vert **lctleft, 
				struct vert **lctright) {
	struct vert *x = lrightmost;
	struct vert *y = rleftmost;
	struct vert *rfast = first(y);
	struct vert *lfast = pred(x, first(x));
	struct vert *temp;

	while(true) {
		if(rightOf(rfast, x, y)) {
			temp = rfast;
			rfast = succ(rfast, y);
			y = temp;
		}
		else if(rightOf(lfast, x, y)) {
			temp = lfast
			lfast = pred(lfast, x);
			x = temp;
		}
		else {
			*lctleft = x;
			*lctright = y;
			break;
		}
	}
}

// Upper common tangent of two convex hulls
static void uct(const struct vert *lrightmost, 
				const struct vert *rleftmost, 
				struct vert **uctleft, 
				struct vert **uctright) {
	struct vert *x = lrightmost;
	struct vert *y = rleftmost;
	struct vert *lfast = first(x);
	struct vert *rfast = pred(y, first(y));
	struct vert *temp;

	while(true) {
		if(leftOf(rfast, x, y)) {
			temp = rfast;
			rfast = pred(rfast, y);
			y = temp;
		}
		else if(leftOf(lfast, x, y)) {
			temp = lfast;
			lfast = succ(lfast, x);
			x = temp;
		}
		else {
			*uctleft = x;
			*uctright = y;
			break;
		}
	}
}


void dtriangulate(struct dTriangulation *tri) {
	mtri = tri;

	// construct vertex structures of points
	struct vert *v = (struct vert*)malloc(mtri->npoints * sizeof(struct vert));

	for(int i = 0; i < mtri->npoints; ++i) {
		v[i].index = i;
		v[i].adj = NULL;
	}

	// sort vertices lexicographically by point coordinates
	// TODO: remove duplicate points! (within DTEPSILON range)
	qsort(v, mtri->npoints, sizeof(struct vert), compareVerts);

	// DEBUG
	// FILE *f = fopen("points.txt", "w");

	// for(int i = 0; i < mtri->npoints; ++i) {
	// 	fprintf(f, "%f, %f\n", XX(&v[i]), YY(&v[i]));
	// }

	// fclose(f);

	// if(inCircle(&v[0], &v[1], &v[2], &v[3])) {
	// 	printf("In circle!\n");
	// }
	// else {
	// 	printf("Not in circle!\n");
	// }
	//

	// triangulate the sorted points
	struct vert *l, *r;
	ord_dtriangulate(v, 0, mtri->npoints - 1, &l, &r);

	// TODO: free vertnodes in adjacency lists!
	free(v);
}

static void ord_dtriangulate(	struct vert *v, 
								int ia, 
								int ib, 
								struct vert **leftmost,
								struct vert **rightmost) {
	if(ia >= ib)
		return;
	if(ib - ia == 1) {
		// num points = 2. Handle this base case
	}
	else if(ib - ia == 2) {
		// num points = 3. Handle this base case
	}
	else {
		struct vert *lo, *li, *ri, *ro;
		struct vert *lctl, *lctr;
		int mid = (ia + ib) / 2;

		ord_dtriangulate(v, ia, mid, &lo, &li);
		ord_dtriangulate(v, mid + 1, ib, &ri, &ro);

		lct(li, ri, &lctl, &lctr);



		*leftmost = lo;
		*rightmost = ro;
	}
}
