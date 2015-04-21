//
//  ip.cpp
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/20/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define INFINITY         1e10

/*
 * Applies the given effect to each pixel in the given source image.
 */
Image* ip_apply_effect_per_pixel (Image* src, Pixel (*effect)(Pixel, void*), void* args)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            Pixel outputPixel = (*effect)(srcPixel, args);
            output->setPixel(x, y, outputPixel);
        }
    }
    return output;
}


/*
 * Generates a solid-color image with the color of the given pixel.
 */
Image* ip_generate_solid_image (int width, int height, Pixel pixel)
{
    Image* output = new Image(width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            output->setPixel(x, y, pixel);
        }
    }
    return output;
}

/*
 * Calculates the psychosomatic intensity.
 */
double ip_get_intensity (Pixel srcPixel)
{
    const double redIntensityWeight = 0.2126;
    const double greenIntensityWeight = 0.7152;
    const double blueIntensityWeight = 0.0722;
    
    double red = srcPixel.getColor(RED);
    double green = srcPixel.getColor(GREEN);
    double blue = srcPixel.getColor(BLUE);
    return clamp(redIntensityWeight * red +
                 greenIntensityWeight * green +
                 blueIntensityWeight*blue,
                 0, 1);
}

/*
 * Creates a pixel with as gray of its psychosomatic intensity.
 */
Pixel ip_gray_pixel (Pixel srcPixel, void* args)
{
    double intensity = ip_get_intensity(srcPixel);
    Pixel outputPixel = Pixel(intensity, intensity, intensity);
    return outputPixel;
}


/*
 * Creates a new image with values equal to the psychosomatic intensities
 * of the source image.
 */
Image* ip_gray (Image* src)
{
    return ip_apply_effect_per_pixel(src, &ip_gray_pixel, NULL);
}

/*
 * Returns the energy function of pixel at the given coordinate in the given image,
 * defined as the sum of the change in intensity from the pixel's x-neighbor plus the change
 * in intensity from the pixel's y-neighbor (using x+1 and y+1 whenever possible, else x-1 and y-1).
 * The double returned will be in the range [0, 2].
 */
double ip_get_energy(Image* src, int x, int y)
{
    int width = src->getWidth();
    int height = src->getHeight();
    double energy = 0;
    Pixel srcPixel = src->getPixel(x, y);
    if (x + 1 < width) {
        Pixel rightNeighbor = src->getPixel(x + 1, y);
        energy += fabs(ip_get_intensity(rightNeighbor) - ip_get_intensity(srcPixel));
    } else if (x - 1 >= 0) {
        Pixel leftNeighbor = src->getPixel(x - 1, y);
        energy += fabs(ip_get_intensity(srcPixel) - ip_get_intensity(leftNeighbor));
    }
    if (y + 1 < height) {
        Pixel bottomNeighbor = src->getPixel(x, y + 1);
        energy += fabs(ip_get_intensity(bottomNeighbor) - ip_get_intensity(srcPixel));
    } else if (y - 1 >= 0) {
        Pixel topNeighbor = src->getPixel(x, y - 1);
        energy += fabs(ip_get_intensity(srcPixel) - ip_get_intensity(topNeighbor));
    }
    return energy;
}

/*
 * Shows the energy of each pixel.
 */
Image* ip_energy(Image* src)
{
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            double energy = ip_get_energy(src, x, y);
            Pixel outputPixel = Pixel(energy / 2, energy / 2, energy / 2);
            output->setPixel(x, y, outputPixel);
        }
    }
    return output;
}

/*
 * Computes the lowest-energy vertical seam in the range [topY, bottomY] and returns the
 * x-coordinates of its pixels.
 */
int* ip_get_vertical_seam_in_range(Image* src, int topY, int bottomY)
{
    int width = src->getWidth();
    //int height = src->getHeight();
    int rangeHeight = bottomY - topY + 1;
    double energyTable [width][rangeHeight]; // DP energy of pixel
    int parentTable [width][rangeHeight]; // where the previous pixel in the seam came from
    // fill out the energyTable and parentTable
    for (int y = topY; y <= bottomY; y++) {
        int tableY = y - topY;
        for (int x = 0; x < width; x++) {
            if (y == 0) {
                energyTable[x][tableY] = ip_get_energy(src, x, y);
                parentTable[x][tableY] = -1; // a pixel in the top row has no parent
            } else {
                double energyOfPathToTop = energyTable[x][tableY - 1];
                // the energies of the paths to top left or top right depend on the neighbor existing
                double energyOfPathToTopLeft = INFINITY;
                double energyOfPathToTopRight = INFINITY;
                if (x > 0) {
                    energyOfPathToTopLeft = energyTable[x - 1][tableY - 1];
                }
                if (x + 1 < width) {
                    energyOfPathToTopRight = energyTable[x + 1][tableY - 1];
                }
                double minParentEnergy = min(energyOfPathToTop, min(energyOfPathToTopLeft, energyOfPathToTopRight));
                energyTable[x][tableY] = ip_get_energy(src, x, y) + minParentEnergy;
                
                if (minParentEnergy == energyOfPathToTop) {
                    parentTable[x][tableY] = x;
                } else if (minParentEnergy == energyOfPathToTopLeft) {
                    parentTable[x][tableY] = x - 1;
                } else {
                    parentTable[x][tableY] = x + 1;
                }
            }
        }
    }
    
    int *seam = new int[rangeHeight];
    // determine the x-coordinate on the bottom row with the lowest path energy
    double minWeight = INFINITY;
    double minWeightXCoord = 0;
    for (int x = 0; x < width; x++) {
        if (energyTable[x][rangeHeight - 1] < minWeight) {
            minWeight = energyTable[x][rangeHeight - 1];
            minWeightXCoord = x;
        }
    }
    seam[rangeHeight - 1] = minWeightXCoord;
    // now that we have the end of the seam determined, backtrack up the parents to find the rest
    for (int i = rangeHeight - 1 - 1; i >= 0; i--) {
        seam[i] = parentTable[seam[i + 1]][i + 1];
    }
    return seam;
}

/*
 * Computes the lowest-energy vertical seam and returns the x-coordinates of its pixels.
 */
int* ip_get_vertical_seam(Image* src)
{
    return ip_get_vertical_seam_in_range(src, 0, src->getHeight() - 1);
}

/*
 * Creates a new image with the lowest-cost vertical seam shown in red.
 */
Image* ip_show_vertical_seam(Image* src)
{
    int *seam = ip_get_vertical_seam(src);
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    Pixel redPixel = Pixel(1, 0, 0);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            if (x == seam[y]) {
                output->setPixel(x, y, redPixel);
            } else {
                output->setPixel(x, y, srcPixel);
            }
        }
    }
    return output;
}

/*
 * Creates an image one pixel less wide with the lowest-cost vertical seam removed.
 */
Image* ip_carve_vertical_seam(Image* src)
{
    int *seam = ip_get_vertical_seam(src);
    
    int srcWidth = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image(srcWidth - 1, height);
    for (int x = 0; x < srcWidth - 1; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel;
            if (x >= seam[y]) {
                // the current pixel is to the right of the seam (or in it), so sample from the right
                srcPixel = src->getPixel(x + 1, y);
            } else {
                // this pixel is to the left of the seam and is thus unaffected by the carving
                srcPixel = src->getPixel(x, y);
            }
            output->setPixel(x, y, srcPixel);
        }
    }
    return output;
}

/*
 * Creates an image by removing numSeams vertical seams from the given source image.
 */
Image* ip_carve_vertical_seams(Image* src, int numSeams)
{
    int srcWidth = src->getWidth();
    if (numSeams > srcWidth) {
        cerr << "Cannot carve more seams than the width of the image" << endl;
        return new Image(*src);
    } else {
        Image* resultImage = new Image(*src);
        for (int i = 0; i < numSeams; i++) {
            resultImage = ip_carve_vertical_seam(resultImage);
        }
        return resultImage;
    }
}

/*
 * Fits the given src image to its rectangular boundaries, using white as transparent.
 * This is the main method of the project.
 */
Image* ip_rectangle(Image* src)
{
    1; // Currently, this method just replaces all white pixels with blue pixels to show the boundaries.
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image(width, height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            double red = srcPixel.getColor(RED);
            double green = srcPixel.getColor(GREEN);
            double blue = srcPixel.getColor(BLUE);
            if (red == 1 && green == 1 && blue == 1) {
                Pixel bluePixel = Pixel(0, 0, 1);
                output->setPixel(x, y, bluePixel);
            } else {
                output->setPixel(x, y, srcPixel);
            }
        }
    }
    return output;
}
