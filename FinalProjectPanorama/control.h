//
//  control.h
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/14/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//
//  This code is adapted from the Image Processing 2015 skeleton code,
//  which is copyright (c) 2015 Z Sweedyk, all rights reserved.

#ifndef __FinalProjectPanorama__control__
#define __FinalProjectPanorama__control__

#include <stdio.h>
#include "common.h"


// init calls
int make_menu();

// glut callbacks
void menu_func (int value);
void keyboard_func (unsigned char key, int x, int y);
void mouse_click_func (int button, int state, int x, int y);
void mouse_move_func (int x, int y);
void idle_func ();

// menu calls
void menu_help ();
void image_load (const char* filename);
void image_save (const char* filename);
void image_print_info ();
void image_revert ();
void process_func (int value);

// text mode calls
void textMenuLoop ();

// helper functions
int getFilterSize();
int getInt(const char* message);
double getDouble(const char* message);
double getPositiveDouble(const char* message);




#endif /* defined(__FinalProjectPanorama__control__) */
