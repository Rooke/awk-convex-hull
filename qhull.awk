#!/bin/gawk -f
#

# Implementation of QuickHull algorithm
#
# Input: Two points (L1, L2) defining a line, and a point P
#
# Returns: a number for distance comparison (using isLeft).
#          This is for finding the biggest distance when we don't
#          care about the magnitude (but do care about the
#          relative magnitude).
#
function cmpDistToLine(L1x,L1y, L2x,L2y, Px,Py,     tmp) {   
    tmp = isLeft(L1x, L1y, L2x, L2y, Px, Py);
    return tmp >= 0 ? tmp : -tmp;
}

# Tests if P2 is left of the base line (P0,P1)  
# Input: Three points, P0, P1, P2
#
# Returns: >0 for P2 left of the line through P0 and P1
#          =0 for P2 on the line
#          <0 for P2 right of the line
#
function isLeft(P0x,P0y, P1x,P1y, P2x,P2y){
    return (P1x - P0x)*(P2y - P0y) - (P2x - P0x)*(P1y - P0y);
}

function ajoin(array,sep) {
    r = "";
    n = 1;
    while (n in array) {
        if (n > 1) r = r sep;
        r = r array[n];
        n++;
    }
    return r;
}

# Input: Array of points P
#        Indicies of a base line (P[a], P[b])
#
# Returns: The index of the point in P which is furthest from 
#          the base line (P[a], P[b]).
#
function indexOfMax(P,a,b,    av, bv,i,v, max, maxi, dist){   

    split(P[a],av,/,/);
    split(P[b],bv,/,/);

    max = 0.0;
    maxi = a;

    for (i in P){
        split(P[i],v,/,/);
        dist = cmpDistToLine(av[1],av[2], bv[1],bv[2], v[1],v[2]);
        if (dist>max) {
            max = dist;
            maxi = i;
        }
    }
    return maxi;
}

# Input: Two arrays, P and H
#
# Returns: The union of P and H, contained in P
#
function unionOf(P,H,     i,j,S) {

    i = 0;
    while (i in P) S[P[i++]] = 1;
    for (j in H) if (!(H[j] in S)) P[i++] = H[j]; 
}
        
# Input: Array of points P
#        Indicies (l, r) of two points in the convex hull of P
#
# Returns: The convex hull of P
#
function cHull_r(P,ileft,iright,   l,r,imax,ileftH1,ileftH2,irightH1,
                                   irightH2,v,L3,H1,H2,i,j,k){


    # base case:
    if (length(P)<=2) {
        return;
    }

    # recursion
    split(P[ileft],l,/,/);
    split(P[iright],r,/,/);
    
    imax = indexOfMax(P, ileft, iright);
    split(P[imax],L3,/,/);

    # compose H1 & H2:
    j = 0;
    k = 0;

    for (i in P) {
        split(P[i],v,",");
        if (i==ileft) {
            ileftH1 = j;
            H1[j++] = P[i];
        } else if (i==imax) {
            irightH1 = j;
            H1[j++] = P[i];
            ileftH2 = k;
            H2[k++] = P[i];
        } else if (i==iright) {
            irightH2 = k;
            H2[k++] = P[i];
        } else if (isLeft(l[1],l[2], L3[1],L3[2], v[1],v[2])>0) {
            H1[j++] = P[i];
        } else if (isLeft(L3[1],L3[2], r[1],r[2], v[1],v[2])>0) {
            H2[k++] = P[i];
        }
    }

    cHull_r(H1, ileftH1, irightH1);
    cHull_r(H2, ileftH2, irightH2);
    
    unionOf(H1,H2);

    delete P; # clear out P, replace with hull (H1)

    for (i in H1) P[i] = H1[i];
}

# Input: Array of points 'P' (indexed 1..length(P))
# Returns: the convex hull of P as a polygon in standard notation
#
function cHull(P,               maxX, minX, iright, ileft,H1,H2, 
                                i,j,k,v,l,r,size,irightH1,
                                ileftH1,irightH2,ileftH2) {

    size = length(P);

    # trivial case

    if (size <= 2) return ajoin(P," ");

    # Find the extreme left and right points (minX, maxX)
    # which are definitely in the convex hull.

    split(P[1],v,/,/);

    maxX = v[1]+0;
    minX = v[1]+0;
    iright = 1;
    ileft = 1;
    
    for (i=2;i<=size;i++) {

        split(P[i],v,/,/);

        if (v[1]>maxX) {
            maxX = v[1]+0;
            iright = i;
        }

        if (v[1]<minX) {
            minX = v[1]+0;
            ileft = i;
        }
    }

    # compose H1 & H2:
    split(P[iright],r,/,/);
    split(P[ileft],l,/,/);
    j = 0;
    k = 0;

    for (i in P){
        split(P[i],v,/,/);
        if (i==ileft) {
            ileftH1 = j;
            H1[j++] = P[i];
            irightH2 = k;
            H2[k++] = P[i];
        } else if (i==iright) {
            irightH1 = j;
            H1[j++] = P[i];
            ileftH2 = k;
            H2[k++] = P[i];
        } else if (isLeft(l[1],l[2], r[1],r[2], v[1],v[2])>0) {
            H1[j++] = P[i];
        } else if (isLeft(r[1],r[2], l[1],l[2], v[1],v[2])>0) {
            H2[k++] = P[i];
        }
    }

    cHull_r(H1, ileftH1, irightH1);
    cHull_r(H2, ileftH2, irightH2);

    unionOf(H1,H2);

    delete H2;

    circular_sort(H1,H2);

    return ajoin(H2," ");
}

# Input: Array of points 'P', with index starting at 1
#
# Performs a "circular sort" on P - i.e. points follow a
# clockwise order after the sort. Puts the result in 'H'
#
# Returns: the size of 'H'
#
function circular_sort(P,H,   atans,x0,y0,Pi,ind,n,i,j){

    n = 0;
    PI = 3.1415926535897932384626433832795;

    for(i in P){
        split(P[i],Pi,/,/);
        x0 += Pi[1];
        y0 += Pi[2];
        n++;
    }

    if (n==0) return;

    x0 /= n;
    y0 /= n;

    for (i in P) {
        split(P[i],Pi,",");
        j = sprintf("%03.8f",atan2(x0-Pi[1],y0-Pi[2])+PI);
        atans[j] = P[i];
    }

    j = 1;
    for (i in atans) ind[j++] = i;

    n = asort(ind);
    
    for (i=1;i<=n;i++) H[i] = atans[ind[i]];

    H[n+1] = H[1];

    return n+1;
}

BEGIN {

    print "Some points:";

    pts[1] = "1,1";
    pts[2] = "2,2";
    pts[3] = "3,3";
    pts[4] = "2,3";
    pts[5] = "2,1";
    pts[6] = "3,1";
    pts[7] = "1,3";

    for(i in pts) print pts[i];

    print "Convex hull of those points:";

    chull = cHull(pts);
    split(chull,a,/ /);
    for(i in a) print a[i];
}
