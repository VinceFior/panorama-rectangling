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
 * Computes the lowest-energy seam in the range [start, end] and returns the
 * coordinates of its pixels. Start is top/left, end is bottom/right.
 */
int* ip_get_seam_in_range(Image* src, SeamOrientation orientation, int start, int end)
{
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    int range = end - start + 1; // +1 because inclusive
    int outputWidth;
    int outputHeight;
    // starts and ends are all inclusive
    int xStart = 0;
    int yStart = 0;
    int xEnd = srcWidth - 1;
    int yEnd = srcHeight - 1;
    if (orientation == ORIENTATION_VERTICAL) {
        outputWidth = srcWidth;
        outputHeight = range;
        yStart = start;
        yEnd = end;
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        outputWidth = range;
        outputHeight = srcHeight;
        xStart = start;
        xEnd = end;
    }
    double energyTable [outputWidth][outputHeight]; // DP energy of pixel
    int parentTable [outputWidth][outputHeight]; // where the previous pixel in the seam came from
    // fill out the energyTable and parentTable
    if (orientation == ORIENTATION_VERTICAL) {
        for (int y = yStart; y <= yEnd; y++) {
            int tableY = y - yStart;
            for (int x = xStart; x <= xEnd; x++) {
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
                    if (x + 1 < srcWidth) {
                        energyOfPathToTopRight = energyTable[x + 1][tableY - 1];
                    }
                    double minParentEnergy = min(energyOfPathToTop, min(energyOfPathToTopLeft, energyOfPathToTopRight));
                    energyTable[x][tableY] = ip_get_energy(src, x, y) + minParentEnergy;
                    
                    if (minParentEnergy == energyOfPathToTop) {
                        parentTable[x][tableY] = x;
                    } else if (minParentEnergy == energyOfPathToTopLeft) {
                        parentTable[x][tableY] = x - 1;
                    } else { // energyOfPathToTopRight
                        parentTable[x][tableY] = x + 1;
                    }
                }
            }
        }
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        for (int x = xStart; x <= xEnd; x++) {
            int tableX = x - xStart;
            for (int y = yStart; y <= yEnd; y++) {
                if (x == 0) {
                    energyTable[tableX][y] = ip_get_energy(src, x, y);
                    parentTable[tableX][y] = -1; // a pixel in the leftmost row has no parent
                } else {
                    double energyOfPathToLeft = energyTable[tableX - 1][y];
                    // the energies of the paths to left top or left bottom depend on the neighbor existing
                    double energyOfPathToLeftTop = INFINITY;
                    double energyOfPathToLeftBottom = INFINITY;
                    if (y > 0) {
                        energyOfPathToLeftTop = energyTable[tableX - 1][y - 1];
                    }
                    if (y + 1 < srcHeight) {
                        energyOfPathToLeftBottom = energyTable[tableX - 1][y + 1];
                    }
                    double minParentEnergy = min(energyOfPathToLeft, min(energyOfPathToLeftTop, energyOfPathToLeftBottom));
                    energyTable[tableX][y] = ip_get_energy(src, x, y) + minParentEnergy;
                    
                    if (minParentEnergy == energyOfPathToLeft) {
                        parentTable[tableX][y] = y;
                    } else if (minParentEnergy == energyOfPathToLeftTop) {
                        parentTable[tableX][y] = y - 1;
                    } else { // energyOfPathToLeftBottom
                        parentTable[tableX][y] = y + 1;
                    }
                }
            }
        }
    }
    
    
    int *seam = new int[range];
    // determine the coordinate on the bottom row with the lowest path energy
    double minWeight = INFINITY;
    double minWeightCoord = 0;
    if (orientation == ORIENTATION_VERTICAL) {
        for (int x = 0; x < srcWidth; x++) {
            if (energyTable[x][range - 1] < minWeight) {
                minWeight = energyTable[x][range - 1];
                minWeightCoord = x;
            }
        }
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        for (int y = 0; y < srcHeight; y++) {
            if (energyTable[range - 1][y] < minWeight) {
                minWeight = energyTable[range - 1][y];
                minWeightCoord = y;
            }
        }
    }
    seam[range - 1] = minWeightCoord;
    // now that we have the end of the seam determined, backtrack up the parents to find the rest
    if (orientation == ORIENTATION_VERTICAL) {
        for (int i = range - 1 - 1; i >= 0; i--) {
            seam[i] = parentTable[seam[i + 1]][i + 1];
        }
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        for (int i = range - 1 - 1; i >= 0; i--) {
            seam[i] = parentTable[i + 1][seam[i + 1]];
        }
    }
    return seam;
}

/*
 * Computes the lowest-energy seam and returns the coordinates of its pixels.
 */
int* ip_get_seam(Image* src, SeamOrientation orientation)
{
    if (orientation == ORIENTATION_VERTICAL) {
        return ip_get_seam_in_range(src, orientation, 0, src->getHeight() - 1);
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        return ip_get_seam_in_range(src, orientation, 0, src->getWidth() - 1);
    } else {
        // this case should never happen
        return nullptr;
    }
}

/*
 * Creates a new image with the lowest-cost seam shown in red.
 */
Image* ip_show_seam(Image* src, SeamOrientation orientation)
{
    int *seam = ip_get_seam(src, orientation);
    
    int width = src->getWidth();
    int height = src->getHeight();
    Image* output = new Image (width, height);
    Pixel redPixel = Pixel(1, 0, 0);
    
    if (orientation == ORIENTATION_VERTICAL) {
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
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Pixel srcPixel = src->getPixel(x, y);
                if (y == seam[x]) {
                    output->setPixel(x, y, redPixel);
                } else {
                    output->setPixel(x, y, srcPixel);
                }
            }
        }
    }

    return output;
}

/*
 * Creates an image one pixel less wide or tall with the lowest-cost seam removed.
 */
Image* ip_carve_seam(Image* src, SeamOrientation orientation)
{
    int *seam = ip_get_seam(src, orientation);
    
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    int outputWidth;
    int outputHeight;
    if (orientation == ORIENTATION_VERTICAL) {
        outputWidth = srcWidth - 1;
        outputHeight = srcHeight;
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        outputWidth = srcWidth;
        outputHeight = srcHeight - 1;
    }
    Image* output = new Image(outputWidth, outputHeight);
    
    for (int x = 0; x < outputWidth; x++) {
        for (int y = 0; y < outputHeight; y++) {
            Pixel srcPixel;
            if (orientation == ORIENTATION_VERTICAL) {
                if (x >= seam[y]) {
                    // the current pixel is to the right of the seam (or in it), so sample from the right
                    srcPixel = src->getPixel(x + 1, y);
                } else {
                    // this pixel is to the left of the seam and is thus unaffected by the carving
                    srcPixel = src->getPixel(x, y);
                }
            } else if (orientation == ORIENTATION_HORIZONTAL) {
                if (y >= seam[x]) {
                    // the current pixel is to the bottom of the seam (or in it), so sample from above
                    srcPixel = src->getPixel(x, y + 1);
                } else {
                    // this pixel is above the seam and is thus unaffected by the carving
                    srcPixel = src->getPixel(x, y);
                }
            }
            output->setPixel(x, y, srcPixel);
        }
    }
    return output;
}

/*
 * Creates an image by removing numSeams vertical seams from the given source image.
 */
Image* ip_carve_seams(Image* src, SeamOrientation orientation, int numSeams)
{
    if (orientation == ORIENTATION_VERTICAL) {
        int srcWidth = src->getWidth();
        if (numSeams > srcWidth) {
            cerr << "Cannot carve more seams than the width of the image" << endl;
            return new Image(*src);
        } else {
            Image* resultImage = new Image(*src);
            for (int i = 0; i < numSeams; i++) {
                resultImage = ip_carve_seam(resultImage, orientation);
            }
            return resultImage;
        }
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        int srcHeight = src->getHeight();
        if (numSeams > srcHeight) {
            cerr << "Cannot carve more seams than the height of the image" << endl;
            return new Image(*src);
        } else {
            Image* resultImage = new Image(*src);
            for (int i = 0; i < numSeams; i++) {
                resultImage = ip_carve_seam(resultImage, orientation);
            }
            return resultImage;
        }
    } else {
        // this case should never happen
        return nullptr;
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
