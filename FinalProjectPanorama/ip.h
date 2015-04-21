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

Image* ip_gray (Image* src);
Image* ip_energy(Image* src);
Image* ip_show_vertical_seam(Image* src);
Image* ip_carve_vertical_seams(Image* src, int numSeams);
Image* ip_rectangle(Image* src);

#endif /* defined(__FinalProjectPanorama__ip__) */
