//
//  main.h
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/14/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//
//  This code is adapted from the Image Processing 2015 skeleton code,
//  which is copyright (c) 2015 Z Sweedyk, all rights reserved.

#ifndef FinalProjectPanorama_main_h
#define FinalProjectPanorama_main_h


#include "common.h"
#include "image.h"

extern int window_width;
extern int window_height;

extern Image* currentImage;
extern Image* originalImage;

extern bool quietMode;
extern bool textMode;


int  main (int argc, char** argv);
char* init (int argc, char** argv);
void usage ();
void display ();
void unreshape (int width, int height);
void reshape (int width, int height);


#endif
