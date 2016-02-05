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
	struct vertNode *prev, *next; // TODO: make more memory efficient ex. xor linked list
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

static void insertNodeAfter(struct vertNode *n, struct vertNode *in);
static void insertNode(struct vert *parent, struct vert *in);
static void deleteNode(struct vert *parent, struct vert *child);

static void connectVerts(struct vert *a, struct vert *b);
static void cutVerts(struct vert *a, struct vert *b);

static struct vert *first(const struct vert *vi);
static struct vert *pred(	const struct vert *vi, 
							const struct vert *vj);
static struct vert *succ(	const struct vert *vi, 
							const struct vert *vj);

static void lct(struct vert *lrightmost, 
				struct vert *rleftmost, 
				struct vert **lctleft, 
				struct vert **lctright);

static void uct(struct vert *lrightmost, 
				struct vert *rleftmost, 
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

	return (xa * (yb - yc) - ya * (xb - xc) + xb * yc - yb * xc) > 0; // TODO: should use epsilon here?
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
	m[0][1] * m[1][0] * m[2][2] + m[0][0] * m[1][1] * m[2][2]) > 0; // TODO: should use epsilon here?
}


static inline void insertNodeAfter(struct vertNode *n, struct vertNode *in) {
	struct vertNode *temp = n->next;
	n->next = in;
	temp->prev = in;
	in->prev = n;
	in->next = temp;
}

static void insertNode(struct vert *parent, struct vert *in) {
	struct vertNode *vn = (struct vertNode*)malloc(sizeof(struct vertNode));
	vn->v = in;

	if(parent->adj) { // if parent already has neighbors, then insert in proper position
		struct vertNode *cur, *temp;
		if(rightOf(in, parent, parent->adj->v)) {
			cur = parent->adj->prev;
			while(cur != parent->adj && rightOf(in, parent, cur->v)) {
				cur = cur->prev;
			}
			if(cur == parent->adj) { // then in-vertex is convex hull successor of parent
				parent->adj = vn; // so make in-vertex "first"
				insertNodeAfter(cur->prev, vn);
			}
			else {
				insertNodeAfter(cur, vn);
			}
		}
		else {
			cur = parent->adj->next;
			while(cur != parent->adj && leftOf(in, parent, cur->v)) {
				cur = cur->next;
			}
			insertNodeAfter(cur->prev, vn);
		}
	}
	else { // if parent has no neighbors, add this node and make it a circular list
		parent->adj = vn;
		vn->prev = vn;
		vn->next = vn;
	}
}

static void deleteNode(struct vert *parent, struct vert *child) {
	struct vertNode *vn = parent->adj;
	if(vn) {
		do {
			if(vn->v == child) {
				vn->prev->next = vn->next;
				vn->next->prev = vn->prev;
				if(vn == parent->adj) {
					parent->adj = vn->next;
				}
				free(vn);
			}
			vn = vn->next;
		} while(vn && vn != parent->adj);
	}
}

static void connectVerts(struct vert *a, struct vert *b) {
	if(a && b) {
		insertNode(a, b);
		insertNode(b, a);
	}
}

static void cutVerts(struct vert *a, struct vert *b) {
	if(a && b) {
		deleteNode(a, b);
		deleteNode(b, a);
	}
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
static void lct(struct vert *lrightmost, 
				struct vert *rleftmost, 
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
			temp = lfast;
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
static void uct(struct vert *lrightmost, 
				struct vert *rleftmost, 
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

// triangulates given vertices assuming that they are lexicographically ordered
// primarily by increasing x-coordinate and secondarily by increasing y-coordinate
static void ord_dtriangulate(	struct vert *v, 
								int ia, 
								int ib, 
								struct vert **leftmost, 
								struct vert **rightmost) {
	if(ib - ia == 1) {
		// num points = 2. Handle this base case
		connectVerts(&v[ia], &v[ib]);
		*leftmost = &v[ia];
		*rightmost = &v[ib];
	}
	else if(ib - ia == 2) {
		// num points = 3. Handle this base case
		connectVerts(&v[ia], &v[ia + 1]);
		connectVerts(&v[ia + 1], &v[ib]);
		if(ccw(&v[ia], &v[ia + 1], &v[ib]) || ccw(&v[ia], &v[ib], &v[ia + 1])) {
			connectVerts(&v[ia], &v[ib]);
		} // else, the three points are collinear, so don't connect the first and third point
		*leftmost = &v[ia];
		*rightmost = &v[ib];
	}
	else if(ib - ia >= 3) { // num points >= 4
		struct vert *lo, *li, *ri, *ro;
		struct vert *lctl, *lctr, *uctl, *uctr;
		int mid = (ia + ib) / 2;

		// triangulate two halves of point set
		ord_dtriangulate(v, ia, mid, &lo, &li);
		ord_dtriangulate(v, mid + 1, ib, &ri, &ro);

		// get lower and upper common tangents between the two halves
		lct(li, ri, &lctl, &lctr);
		uct(li, ri, &uctl, &uctr);

		// merge the two halves
		li = lctl;
		ri = lctr;
		bool a, b;
		struct vert *l1, *l2, *r1, *r2;
		while(li != uctl || ri != uctr) { // connect from bottom to top
			a = false, b = false;
			connectVerts(li, ri);

			r1 = pred(ri, li);
			if(leftOf(r1, li, ri)) {
				r2 = pred(ri, r1);
				while(inCircle(r1, li, ri, r2)) {
					cutVerts(ri, r1);
					r1 = r2;
					r2 = pred(ri, r1);
				}
			}
			else {
				a = true;
			}

			l1 = succ(li, ri);
			if(rightOf(l1, ri, li)) {
				l2 = succ(li, l1);
				while(inCircle(li, ri, l1, l2)) {
					cutVerts(li, l1);
					l1 = l2;
					l2 = succ(li, l1);
				}
			}
			else {
				b = true;
			}

			if(a) {
				li = l1;
			}
			else if(b) {
				ri = r1;
			}
			else if(!inCircle(li, ri, r1, l1)) {
				ri = r1;
			}
			else {
				li = l1;
			}
		}
		connectVerts(uctl, uctr); // connect the top

		*leftmost = lo;
		*rightmost = ro;
	}
	// else, num points <=1; invalid input so do nothing
}
