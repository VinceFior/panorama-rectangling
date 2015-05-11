//
//  control.cpp
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/14/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//
//  This code is adapted from the Image Processing 2015 skeleton code,
//  which is copyright (c) 2015 Z Sweedyk, all rights reserved.

#include "control.h"
#include "ip.h"
#include "main.h"
#include <stdlib.h>

enum MenuItem {
    MENU_QUIT,
    MENU_HELP,
    MENU_FILE_OPEN,
    MENU_FILE_SAVE,
    MENU_FILE_INFO,
    MENU_FILE_REVERT,
    
    MENU_PROCESS_GRAY,
    MENU_PROCESS_ENERGY,
    MENU_PROCESS_SHOW_VERTICAL_SEAM,
    MENU_PROCESS_SHOW_HORIZONTAL_SEAM,
    MENU_PROCESS_CARVE_VERTICAL_SEAM,
    MENU_PROCESS_CARVE_HORIZONTAL_SEAM,
    MENU_PROCESS_CARVE_VERTICAL_SEAMS,
    MENU_PROCESS_CARVE_HORIZONTAL_SEAMS,
    MENU_PROCESS_INSERT_VERTICAL_SEAMS,
    MENU_PROCESS_INSERT_HORIZONTAL_SEAMS,
    MENU_PROCESS_LOCAL_WARP,
    MENU_PROCESS_SHOW_INPUT_MESH,
    MENU_PROCESS_SHOW_OUTPUT_MESH,
    MENU_PROCESS_RECTANGLE,
};


int make_menu ()
{
    int file = glutCreateMenu(menu_func);
    glutAddMenuEntry( "Open..",         MENU_FILE_OPEN);
    glutAddMenuEntry( "Save..",         MENU_FILE_SAVE);
    glutAddMenuEntry( "Get Image Info",		MENU_FILE_INFO);
    glutAddMenuEntry( "Revert",		MENU_FILE_REVERT);
    
    int process = glutCreateMenu(menu_func);
    glutAddMenuEntry("Grayscale", MENU_PROCESS_GRAY);
    glutAddMenuEntry("Energy", MENU_PROCESS_ENERGY);
    glutAddMenuEntry("Show vertical seam", MENU_PROCESS_SHOW_VERTICAL_SEAM);
    glutAddMenuEntry("Show horizontal seam", MENU_PROCESS_SHOW_HORIZONTAL_SEAM);
    glutAddMenuEntry("Carve vertical seam", MENU_PROCESS_CARVE_VERTICAL_SEAM);
    glutAddMenuEntry("Carve horizontal seam", MENU_PROCESS_CARVE_HORIZONTAL_SEAM);
    glutAddMenuEntry("Carve vertical seams..", MENU_PROCESS_CARVE_VERTICAL_SEAMS);
    glutAddMenuEntry("Carve horizontal seams..", MENU_PROCESS_CARVE_HORIZONTAL_SEAMS);
    glutAddMenuEntry("Insert vertical seams..", MENU_PROCESS_INSERT_VERTICAL_SEAMS);
    glutAddMenuEntry("Insert horizontal seams..", MENU_PROCESS_INSERT_HORIZONTAL_SEAMS);
    glutAddMenuEntry("Local warp", MENU_PROCESS_LOCAL_WARP);
    glutAddMenuEntry("Show input mesh", MENU_PROCESS_SHOW_INPUT_MESH);
    glutAddMenuEntry("Show optimized mesh", MENU_PROCESS_SHOW_OUTPUT_MESH);
    glutAddMenuEntry("Rectangle", MENU_PROCESS_RECTANGLE);
    
    int main = glutCreateMenu(menu_func);
    glutAddSubMenu("File", file);
    glutAddSubMenu("Process", process);
    glutAddMenuEntry( "Help", MENU_HELP);
    glutAddMenuEntry( "Quit", MENU_QUIT);
    
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    return main;
}


static inline void checkStream (const istream& in)
{
    if (in.fail())
    {
        cerr << "Fatal error: stream failed!" << endl;
        exit(-1);
    }
}


void menu_func (int value)
{
    // variables used in the switch statement
    char filename[MAX_LINE];
    
    switch (value) {
        case MENU_QUIT:
            exit(0);
            break;
            
        case MENU_HELP:
            menu_help();
            break;
            
        case MENU_FILE_OPEN:
            if (!quietMode)
                cerr << "Open file (string - no spaces) : ";
            cin  >> filename;
            checkStream(cin);
            image_load(filename);
            break;
            
        case MENU_FILE_SAVE:
            if (!quietMode)
                cerr << "Save as (string - no spaces) : ";
            cin  >> filename;
            checkStream(cin);
            image_save(filename);
            break;
            
        case MENU_FILE_INFO:
            image_print_info();
            break;
            
        case MENU_FILE_REVERT:
            image_revert();
            break;
            
        default:
            process_func(value);
    }
    return;
}

void process_func (int value)
{
    
    Image* resultImage = NULL;
    
    //  check if we have an image to process
    if (!currentImage) {
        cerr << "Sorry, no image is loaded!" << endl;
        return;
    }
    
    switch (value) {
        case MENU_PROCESS_ENERGY:
            resultImage = ip_seam_energy(currentImage);
            break;
            
        case MENU_PROCESS_GRAY:
            resultImage = ip_gray(currentImage);
            break;
        
        case MENU_PROCESS_SHOW_VERTICAL_SEAM:
            resultImage = ip_show_seam(currentImage, ORIENTATION_VERTICAL);
            break;
            
        case MENU_PROCESS_SHOW_HORIZONTAL_SEAM:
            resultImage = ip_show_seam(currentImage, ORIENTATION_HORIZONTAL);
            break;
        
        case MENU_PROCESS_CARVE_VERTICAL_SEAM:
            resultImage = ip_carve_seams(currentImage, ORIENTATION_VERTICAL, 1);
            break;
            
        case MENU_PROCESS_CARVE_HORIZONTAL_SEAM:
            resultImage = ip_carve_seams(currentImage, ORIENTATION_HORIZONTAL, 1);
            break;
        
        case MENU_PROCESS_CARVE_VERTICAL_SEAMS:
            int numSeams;
            if (!quietMode) {
                cerr << "Number of vertical seams to carve: ";
            }
            cin >> numSeams;
            cerr << "Carving.." << endl;
            resultImage = ip_carve_seams(currentImage, ORIENTATION_VERTICAL, numSeams);
            break;
            
        case MENU_PROCESS_CARVE_HORIZONTAL_SEAMS:
            if (!quietMode) {
                cerr << "Number of horizontal seams to carve: ";
            }
            cin >> numSeams;
            cerr << "Carving.." << endl;
            resultImage = ip_carve_seams(currentImage, ORIENTATION_HORIZONTAL, numSeams);
            break;
        
        case MENU_PROCESS_INSERT_VERTICAL_SEAMS:
            if (!quietMode) {
                cerr << "Number of vertical seams to insert: ";
            }
            cin >> numSeams;
            cerr << "Inserting.." << endl;
            resultImage = ip_insert_seams(currentImage, ORIENTATION_VERTICAL, numSeams);
            break;
            
        case MENU_PROCESS_INSERT_HORIZONTAL_SEAMS:
            if (!quietMode) {
                cerr << "Number of horizontal seams to insert: ";
            }
            cin >> numSeams;
            cerr << "Inserting.." << endl;
            resultImage = ip_insert_seams(currentImage, ORIENTATION_HORIZONTAL, numSeams);
            break;

        case MENU_PROCESS_LOCAL_WARP:
            cerr << "Locally warping.." << endl;
            resultImage = ip_local_warp(currentImage);
            break;

        case MENU_PROCESS_SHOW_INPUT_MESH: {
            cerr << "Rectangling.." << endl;
            resultImage = ip_rectangle(currentImage, true, false);
            break;
        }
            
        case MENU_PROCESS_SHOW_OUTPUT_MESH: {
            cerr << "Rectangling.." << endl;
            resultImage = ip_rectangle(currentImage, false, true);
            break;
        }
            
        case MENU_PROCESS_RECTANGLE: {
            cerr << "Rectangling.." << endl;
            resultImage = ip_rectangle(currentImage, false, false);
            break;
        }
            
        default:
            break;
    }
    
    if (resultImage != NULL) {
        delete currentImage;
        currentImage = resultImage;
        
        if (currentImage->getWidth()  != window_width    ||
            currentImage->getHeight() != window_height) {
            reshape(currentImage->getWidth(), currentImage->getHeight());
        }
        
        if (!quietMode) {
            cerr << "Process done!" << endl;
        }
        
        if (!textMode) {
            glutPostRedisplay();
        }
    }
}

void keyboard_func (unsigned char key, int x, int y)
{
    switch (key) {
        case 'X':
        case 'x':
            menu_func(MENU_PROCESS_CARVE_VERTICAL_SEAM);
            break;
            
        case 'C':
        case 'c':
            menu_func(MENU_PROCESS_CARVE_HORIZONTAL_SEAM);
            break;
            
        case 'R':
        case 'r':
            menu_func(MENU_FILE_REVERT);
            break;
            
        case 'O':
        case 'o':
            menu_func(MENU_FILE_OPEN);
            break;
        
        case 'E':
        case 'e':
            menu_func(MENU_PROCESS_RECTANGLE);
            break;
            
        case 'H':
        case 'h':
            menu_help();
            break;
            
        case 'Q':
        case 'q':
            exit(0);
            break;
    }
}

void menu_help ()
{
    cerr << endl << "Open a non-rectangular panorama and apply the warp to rectangle it!" << endl;
//    cerr << endl
//    << "hmc cs155 image processor" << endl
//    << "please see the ip manual for usage and algorithm information" << endl
//    << "http://www.cs.hmc.edu/courses/2002/fall/cs155/proj1/doc/ip_manual.html"
//    << endl << endl;
}

#define MENUOP(num, tag)	cerr << " " << num << ") " << tag << endl;

void textMenuLoop ()
{
    char command[MAX_LINE];
    
    while (true) {
        
        if (!quietMode) {
            cerr << endl << "selection > " << flush;
        }
        
        cin  >> command;
        
        switch (command[0]) {
            case '\n':
            case '\0':
                //printMenu();
                break;
                
            case 'Q':
            case 'q':
                menu_func(MENU_QUIT);
                break;
                
            case 'H':
            case 'h':
                menu_func(MENU_HELP);
                break;
                
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                menu_func(atoi(command));
                break;
                
            default:
                //printMenu();
                break;
        }
    }
}

void image_load (const char* filename)
{
    if (currentImage)
        delete currentImage;
    if (originalImage)
        delete originalImage;
    currentImage  = NULL;
    originalImage = NULL;
    
    originalImage = new Image();
    originalImage->read(filename);
    
    if (originalImage->good())
    {  
        currentImage = new Image(*originalImage);
        reshape(currentImage->getWidth(), currentImage->getHeight());
    }
    else
    {
        delete originalImage;  
        originalImage = NULL;
        cerr << "Couldn't load image " << filename << "!" << endl;
        return;
    }
    
    if (!textMode)
        glutPostRedisplay();
    
    if (!quietMode)
        cerr << "Image loaded!" << endl;
}

void image_save (const char* filename)
{
    if (currentImage)
    {
        if (currentImage->write(filename) == 0)
        {
            //delete originalImage;
            //originalImage = new Image(*currentImage);
        }
    }  
    else if (originalImage)
    {
        originalImage->write(filename);
    }
    else
    {
        cerr << "No image!" << endl;
        return;
    }
    
    if (!quietMode)
        cerr << "Image saved!" << endl;
}

void image_print_info ()
{  
    if (currentImage) {
        cerr << "width:    " << currentImage->getWidth() << endl
        << "height:   " << currentImage->getHeight() << endl
        << "bits:     " << currentImage->getBits() << endl;
    }
    cerr << "Image info printed!" << endl;
}


void image_revert ()
{
    if (currentImage) {
        delete currentImage;
    }
    
    if (originalImage) {
        currentImage = new Image(*originalImage);
        
        if (window_width  != currentImage->getWidth() ||
            window_height != currentImage->getHeight())
            reshape(currentImage->getWidth(), currentImage->getHeight());
    } else {
        cerr << "No image!" << endl;
        return;
    }
    
    if (!textMode) {
        glutPostRedisplay();
    }
    
    if (!quietMode) {
        cerr << "Reverted!" << endl;
    }
}

int getFilterSize()
{
    int filtersize;
    if (!quietMode)
        cerr << "Enter filter size (positive, odd integer) : ";
    cin  >> filtersize;
    if (filtersize % 2 !=1 || filtersize<=0)
    {
        cerr << "Sorry, the filter size must be a positive, odd integer." << endl;
        filtersize=0;
    }
    checkStream(cin);
    return filtersize;
}
double getDouble(const char* message)
{
    double value;
    if (!quietMode)
        cerr << "Enter " << message << "(double): ";
    cin  >> value;
    checkStream(cin);
    return value;
}

double getPositiveDouble(const char* message)
{
    double value;
    if (!quietMode)
        cerr << "Enter positive " << message << "(double): ";
    cin  >> value;
    checkStream(cin);
    return value;
}

int getInt(const char* message)
{
    
    int value;
    if (!quietMode)
        cerr << "Enter " << message << " (integer): ";
    cin  >> value;
    checkStream(cin);
    return value;
}
