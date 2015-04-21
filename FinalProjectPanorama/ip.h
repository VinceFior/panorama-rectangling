//
//  ip.h
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/20/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#ifndef __FinalProjectPanorama__ip__
#define __FinalProjectPanorama__ip__


#include "common.h"
#include "image.h"

enum SeamOrientation {
    ORIENTATION_VERTICAL,
    ORIENTATION_HORIZONTAL
};

enum NeighborDirection {
    DIRECTION_TOP_LEFT,
    DIRECTION_TOP,
    DIRECTION_TOP_RIGHT,
    DIRECTION_LEFT_TOP,
    DIRECTION_LEFT,
    DIRECTION_LEFT_BOTTOM
};

Image* ip_gray (Image* src);
Image* ip_energy(Image* src);
Image* ip_show_seam(Image* src, SeamOrientation orientation);
Image* ip_carve_seams(Image* src, SeamOrientation orientation, int numSeams);
Image* ip_rectangle(Image* src);

#endif /* defined(__FinalProjectPanorama__ip__) */
