import cv2
import numpy as np
import mep


def contours_to_cvxhull(contours):
    '''
    Params
    ------
    @contours: tuple[np.ndarray]
        tuple of contours, each contours is defined as an array of points


    Return
    ------
    cvx_hulls: list[np.ndarray]
        list of convex hull corresponding to each inputed contour
    '''

    cvx_hulls = []

    for i in range(len(contours)):
        hull = cv2.convexHull(contours[i])
        cvx_hulls.append(hull)

    return cvx_hulls


def cvxpolygon_to_mep(cvx_polygons):
    '''
    Params
    ------
    @cvx_polygons: list[np.ndarray]
        list of convex polygons, in which each polygon is defined as an array of points


    Return
    ------
    pargrams: list[np.ndarray]
        list of parallelogram, where each parallelogram is defined as an array of points 
        (each point corresponds to a vertex)
    '''
    pargrams = []

    for cvx_polygon in cvx_polygons:
        if cvx_polygon.shape[0] < 3:
            continue
        cvx_polygon = [mep.Point(pt[0, 0], pt[0, 1]) for pt in cvx_polygon]
        antipodal_evs = mep.antipodal_pairs(cvx_polygon)
        min_ep, ev1, ev2 = mep.simple_mep(antipodal_evs)
        a, b, c, d = min_ep.a.to_np(), min_ep.b.to_np(), min_ep.c.to_np(), min_ep.d().to_np()
        pts = np.array([a, b, c, d])
        pargrams.append(pts)

    return pargrams 
