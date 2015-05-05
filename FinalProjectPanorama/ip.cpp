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
#include "utils.h"
#include "lsd.h"
#include "Quadratic.h"

#define INFT         1e8

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
 * Transposes the image between x and y (replaces x with y and vice-versa).
 */
Image* ip_transpose(Image* src)
{
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    Image* output = new Image (srcHeight, srcWidth);
    for (int x = 0; x < srcHeight; x++) {
        for (int y = 0; y < srcWidth; y++) {
            Pixel srcPixel = src->getPixel(y, x);
            output->setPixel(x, y, srcPixel);
        }
    }
    return output;
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
 * Returns the forward energy given a pixel and its neighbor that are being removed.
 * This method only deals with vertical seams, since it expects horizontal images to be flipped.
 */
double ip_get_forward_energy(Image* src, int x, int y, NeighborDirection direction)
{
    int width = src->getWidth();
    double energy = 0;
    if (direction == DIRECTION_TOP_LEFT) {
        if (x <= 0 || y <= 0) {
            return energy;
        }
        Pixel bottomLeftPixel = src->getPixel(x - 1, y);
        // add top and right
        Pixel topPixel = src->getPixel(x, y - 1);
        energy += fabs(ip_get_intensity(topPixel) - ip_get_intensity(bottomLeftPixel));
        if (x + 1 < width) {
            Pixel rightPixel = src->getPixel(x + 1, y);
            energy += fabs(ip_get_intensity(rightPixel) - ip_get_intensity(bottomLeftPixel));
        }
        return energy;
    } else if (direction == DIRECTION_TOP) {
        // add left and right
        if (x - 1 >= 0 && x + 1 < width) {
            Pixel leftPixel = src->getPixel(x - 1, y);
            Pixel rightPixel = src->getPixel(x + 1, y);
            energy += fabs(ip_get_intensity(leftPixel) - ip_get_intensity(rightPixel));
        }
        return energy;
    } else if (direction == DIRECTION_TOP_RIGHT) {
        if (x >= width - 1 || y <= 0) {
            return energy;
        }
        Pixel bottomRightPixel = src->getPixel(x + 1, y);
        // add top and left
        Pixel topPixel = src->getPixel(x, y - 1);
        energy += fabs(ip_get_intensity(topPixel) - ip_get_intensity(bottomRightPixel));
        if (x - 1 >= 0) {
            Pixel leftPixel = src->getPixel(x - 1, y);
            energy += fabs(ip_get_intensity(leftPixel) - ip_get_intensity(bottomRightPixel));
        }
        return energy;
    } else {
        // this case shouldn't happen
        return 0;
    }
}

/*
 * Shows the seam carving energy of each pixel.
 */
Image* ip_seam_energy(Image* src)
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
 * Determines whether the given pixel is 'transparent'.
 */
bool ip_is_transparent(Pixel pixel) {
    return (pixel.getColor(RED) == 1 && pixel.getColor(GREEN) == 1 && pixel.getColor(BLUE) == 1);
}

/*
 * Determines whether the given pixel is part of the 'transparent' border.
 */
bool ip_is_border(Image* src, int x, int y) {
    Pixel pixel = src->getPixel(x, y);
    if (!ip_is_transparent(pixel)) {
        return false;
    }
    
    
    1;return true; //TODO: BFS
    // do a BFS to reach a border pixel
    
//    expandedPixels = set()
//    fringe = util.Queue()
//    fringe.push(startPos)
//    while not fringe.isEmpty():
//        offPixel = fringe.pop()
//        if offPixel not in expandedPixels:
//            expandedPixels.add(offPixel)
//            for neighbor in getOffNeighbors(datum, offPixel):
//                fringe.push(neighbor)
//                return expandedPixels
}

/*
 * Computes the lowest-energy seam in the range [start, end] using forward energy and returns the
 * coordinates of its pixels. Start is top/left, end is bottom/right.
 */
int* ip_get_seam_in_range(Image* src, SeamOrientation orientation, int start, int end)
{
    if (orientation == ORIENTATION_HORIZONTAL) {
        src = ip_transpose(src);
    }
    int srcWidth = src->getWidth();
    int range = end - start + 1; // +1 because inclusive
    int outputWidth = srcWidth;
    int outputHeight = range;
    // starts and ends are all inclusive
    int xStart = 0;
    int yStart = start;
    int xEnd = srcWidth - 1;
    int yEnd = end;
    
    // DP energy of pixel
    double **energyTable = new double* [outputWidth];
    for (int i = 0; i < outputWidth; i++) {
        energyTable[i] = new double[outputHeight];
    }
    // where the previous pixel in the seam came from
    double **parentTable = new double* [outputWidth];
    for (int i = 0; i < outputWidth; i++) {
        parentTable[i] = new double[outputHeight];
    }
    
    // fill out the energyTable and parentTable
    for (int y = yStart; y <= yEnd; y++) {
        int tableY = y - yStart;
        for (int x = xStart; x <= xEnd; x++) {
            // don't cut a seam through the 'transparent' border
            if (ip_is_border(src, x, y)) {
                energyTable[x][tableY] = INFT * 10;
                parentTable[x][tableY] = x;
                continue;
            }
            if (y == yStart) {
                energyTable[x][tableY] = ip_get_energy(src, x, y);
                parentTable[x][tableY] = -1; // a pixel in the top row has no parent
            } else {
                double energyOfPathFromTop = energyTable[x][tableY - 1] + ip_get_forward_energy(src, x, y, DIRECTION_TOP);
                // the energies of the paths to top left or top right depend on the neighbor existing
                double energyOfPathFromTopLeft = INFT;
                double energyOfPathFromTopRight = INFT;
                if (x > 0) {
                    energyOfPathFromTopLeft = energyTable[x - 1][tableY - 1] + ip_get_forward_energy(src, x, y, DIRECTION_TOP_LEFT);
                }
                if (x + 1 < srcWidth) {
                    energyOfPathFromTopRight = energyTable[x + 1][tableY - 1] + ip_get_forward_energy(src, x, y, DIRECTION_TOP_RIGHT);
                }
                double minParentEnergy = min(energyOfPathFromTop, min(energyOfPathFromTopLeft, energyOfPathFromTopRight));
                energyTable[x][tableY] = ip_get_energy(src, x, y) + minParentEnergy;
                
                if (minParentEnergy == energyOfPathFromTop) {
                    parentTable[x][tableY] = x;
                } else if (minParentEnergy == energyOfPathFromTopLeft) {
                    parentTable[x][tableY] = max(x - 1, 0);
                } else { // energyOfPathFromTopRight
                    parentTable[x][tableY] = x + 1;
                }
            }
        }
    }
    
    int *seam = new int[range];
    // determine the coordinate on the bottom row with the lowest path energy
    double minWeight = INFT;
    double minWeightCoord = 0;
    for (int x = 0; x < srcWidth; x++) {
        if (energyTable[x][range - 1] < minWeight) {
            minWeight = energyTable[x][range - 1];
            minWeightCoord = x;
        }
    }
    seam[range - 1] = minWeightCoord;
    // now that we have the end of the seam determined, backtrack up the parents to find the rest
    for (int i = range - 1 - 1; i >= 0; i--) {
        seam[i] = parentTable[seam[i + 1]][i + 1];
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
    if (orientation == ORIENTATION_HORIZONTAL) {
        src = ip_transpose(src);
    }
    int *seam = ip_get_seam(src, ORIENTATION_VERTICAL);
    
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
    if (orientation == ORIENTATION_HORIZONTAL) {
        return ip_transpose(output);
    } else {
        return output;
    }
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
        if (numSeams > src->getWidth()) {
            cerr << "Cannot carve more seams than the width of the image" << endl;
            return new Image(*src);
        }
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        if (numSeams > src->getHeight()) {
            cerr << "Cannot carve more seams than the height of the image" << endl;
            return new Image(*src);
        }
    } else {
        // this case should never happen
        return nullptr;
    }
    Image* resultImage = new Image(*src);
    for (int i = 0; i < numSeams; i++) {
        resultImage = ip_carve_seam(resultImage, orientation);
    }
    return resultImage;
}

/*
 * Returns a new image with the lowest-cost seam duplicated in the given range, shifting over pixels.
 */
Image* ip_insert_seam_in_range(Image* src, int* seam, SeamOrientation orientation, int start, int end, bool shiftToEnd)
{
    // this method only works with vertical seams, so we transpose horizontal images
    if (orientation == ORIENTATION_HORIZONTAL) {
        src = ip_transpose(src);
    }
    int width = src->getWidth();
    int height = src->getHeight();
    
    Image *outputImage = new Image(width, height);
    for (int x = 0; x < width; x++) {
        // above range
        for (int y = 0; y < start; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            outputImage->setPixel(x, y, srcPixel);
        }
        // below range
        for (int y = end + 1; y < height; y++) {
            Pixel srcPixel = src->getPixel(x, y);
            outputImage->setPixel(x, y, srcPixel);
        }
        // in range
        for (int y = start; y <= end; y++) {
            int seamY = y - start;
            if (x > seam[seamY] && shiftToEnd) {
                // shift
                Pixel srcPixel = src->getPixel(x - 1, y);
                outputImage->setPixel(x, y, srcPixel);
            } else if (x < seam[seamY] && !shiftToEnd) {
                // shift
                Pixel srcPixel = src->getPixel(x + 1, y);
                outputImage->setPixel(x, y, srcPixel);
            } else {
                // do nothing
                Pixel srcPixel = src->getPixel(x, y);
                outputImage->setPixel(x, y, srcPixel);
            }
        }
    }
    if (orientation == ORIENTATION_HORIZONTAL) {
        outputImage = ip_transpose(outputImage);
    }
    return outputImage;
}

/*
 * Returns a new image with the lowest-cost seam duplicated.
 */
Image* ip_insert_seam(Image* src, SeamOrientation orientation, bool shiftToEnd)
{
    if (orientation == ORIENTATION_VERTICAL) {
        int* seam = ip_get_seam_in_range(src, orientation, 0, src->getHeight() - 1);
        return ip_insert_seam_in_range(src, seam, orientation, 0, src->getHeight() - 1, shiftToEnd);
    } else if (orientation == ORIENTATION_HORIZONTAL) {
        int* seam = ip_get_seam_in_range(src, orientation, 0, src->getWidth() - 1);
        return ip_insert_seam_in_range(src, seam, orientation, 0, src->getWidth() - 1, shiftToEnd);
    } else {
        // this case should never happen
        return nullptr;
    }
}

/*
 * Inserts numSeams seams into the image by repeatedly inserting one seam. Note that this means all
 * the inserted seams will be the same.
 */
Image* ip_insert_seams(Image* src, SeamOrientation orientation, int numSeams)
{
    Image* resultImage = new Image(*src);
    for (int i = 0; i < numSeams; i++) {
        resultImage = ip_insert_seam(resultImage, orientation, false);
    }
    return resultImage;
}

/*
 * Returns the longest line of 'transparent' pixels along an edge of the image.
 */
vector<int> ip_get_longest_boundary(Image* src, BorderSide& whichSide)
{
    vector<int> boundary;
    size_t maxBoundaryLength = 0;
    int width = src->getWidth();
    int height = src->getHeight();
    
    // top
    bool isCountingBoundary = false;
    vector<int> currentBoundary;
    for (int x = 0; x < width; x++) {
        Pixel leftPixel = src->getPixel(x, 0);
        if (!ip_is_transparent(leftPixel) || (x == width - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_TOP;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(x);
            isCountingBoundary = true;
        }
    }
    
    // bottom
    isCountingBoundary = false;
    currentBoundary.clear();
    for (int x = 0; x < width; x++) {
        Pixel leftPixel = src->getPixel(x, height - 1);
        if (!ip_is_transparent(leftPixel) || (x == width - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_BOTTOM;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(x);
            isCountingBoundary = true;
        }
    }
    
    // left
    isCountingBoundary = false;
    currentBoundary.clear();
    for (int y = 0; y < height; y++) {
        Pixel topPixel = src->getPixel(0, y);
        if (!ip_is_transparent(topPixel) || (y == height - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_LEFT;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(y);
            isCountingBoundary = true;
        }
    }
    
    // right
    isCountingBoundary = false;
    currentBoundary.clear();
    for (int y = 0; y < height; y++) {
        Pixel topPixel = src->getPixel(width - 1, y);
        if (!ip_is_transparent(topPixel) || (y == height - 1)) {
            if (isCountingBoundary) {
                if (currentBoundary.size() > maxBoundaryLength) {
                    maxBoundaryLength = currentBoundary.size();
                    boundary = currentBoundary;
                    whichSide = BORDER_RIGHT;
                }
            }
            isCountingBoundary = false;
            currentBoundary.clear();
            continue;
        } else {
            currentBoundary.push_back(y);
            isCountingBoundary = true;
        }
    }
    
    return boundary;
}

/*
 * Returns a displacement map to warp the given image with an irregular 'transparent' border to fill
 * its rectangle.
 * Note that the 'transparent' border must NOT have any boundary as long as its dimension (i.e.,
 * no excessive whitespace).
 */
vector<vector<Coordinate>> ip_local_warp_displacement(Image* src)
{
    int srcWidth = src->getWidth();
    int srcHeight = src->getHeight();
    vector<vector<Coordinate>> displacementMap;
    for (int x = 0; x < srcWidth; x++) {
        vector<Coordinate> displacementRow;
        for (int y = 0; y < srcHeight; y++) {
            Coordinate coord;
            coord.x = 0;
            coord.y = 0;
            displacementRow.push_back(coord);
        }
        displacementMap.push_back(displacementRow);
    }
    vector<vector<Coordinate>> finalDisplacementMap;
    for (int x = 0; x < srcWidth; x++) {
        vector<Coordinate> displacementRow;
        for (int y = 0; y < srcHeight; y++) {
            Coordinate coord;
            coord.x = 0;
            coord.y = 0;
            displacementRow.push_back(coord);
        }
        finalDisplacementMap.push_back(displacementRow);
    }
    Image* resultImage = new Image(*src);
    while (true) {
        BorderSide whichSide;
        vector<int> longestBoundary = ip_get_longest_boundary(resultImage, whichSide);
        if (longestBoundary.size() == 0) {
            return displacementMap;
        } else {
            SeamOrientation orientation;
            bool shiftToEnd;
            if (whichSide == BORDER_TOP) {
                orientation = ORIENTATION_HORIZONTAL;
                shiftToEnd = false;
            } else if (whichSide == BORDER_BOTTOM) {
                orientation = ORIENTATION_HORIZONTAL;
                shiftToEnd = true;
            } else if (whichSide == BORDER_LEFT) {
                orientation = ORIENTATION_VERTICAL;
                shiftToEnd = false;
            } else if (whichSide == BORDER_RIGHT) {
                orientation = ORIENTATION_VERTICAL;
                shiftToEnd = true;
            } else {
                // default case, should never happen
                orientation = ORIENTATION_VERTICAL;
                shiftToEnd = true;
            }
            int start = longestBoundary[0];
            int end = longestBoundary[longestBoundary.size() - 1];
            
            int* seam = ip_get_seam_in_range(resultImage, orientation, start, end);
            resultImage = ip_insert_seam_in_range(resultImage, seam, orientation, start, end, shiftToEnd);
            
            // update displacement map based on tmp displacement map
            for (int x = 0; x < srcWidth; x++) {
                for (int y = 0; y < srcHeight; y++) {
                    // compute displacement relative to previous image
                    Coordinate mostRecentDisplacement;
                    mostRecentDisplacement.x = 0;
                    mostRecentDisplacement.y = 0;
                    if (orientation == ORIENTATION_VERTICAL && y >= start && y <= end) {
                        int seamY = y - start;
                        if (x > seam[seamY] && shiftToEnd) {
                            mostRecentDisplacement.x = -1;
                        } else if (x < seam[seamY] && !shiftToEnd) {
                            mostRecentDisplacement.x = 1;
                        }
                    } else if (orientation == ORIENTATION_HORIZONTAL && x >= start && x <= end) {
                        int seamX = x - start;
                        if (y > seam[seamX] && shiftToEnd) {
                            mostRecentDisplacement.y = -1;
                        } else if (y < seam[seamX] && !shiftToEnd) {
                            mostRecentDisplacement.y = 1;
                        }
                    }
                    // find the pixel in displacementMap that tmp points to; find the pixel that that pixel points to; point to it
                    Coordinate &finalDisplacement = finalDisplacementMap[x][y];
                    int mostRecentlyDisplacedX = x + mostRecentDisplacement.x;
                    int mostRecentlyDisplacedY = y + mostRecentDisplacement.y;
                    Coordinate displacementOfTarget = displacementMap[mostRecentlyDisplacedX][mostRecentlyDisplacedY];
                    int xInOriginalImage = mostRecentlyDisplacedX + displacementOfTarget.x;
                    int yInOriginalImage = mostRecentlyDisplacedY + displacementOfTarget.y;
                    finalDisplacement.x = xInOriginalImage - x;
                    finalDisplacement.y = yInOriginalImage - y;
                }
            }
            // update finalDisplacementMap
            for (int x = 0; x < srcWidth; x++) {
                for (int y = 0; y < srcHeight; y++) {
                    Coordinate &displacement = displacementMap[x][y];
                    Coordinate finalDisplacement = finalDisplacementMap[x][y];
                    displacement.x = finalDisplacement.x;
                    displacement.y = finalDisplacement.y;
                }
            }
            
        }
    }
    return displacementMap;
}

/*
 * Warps the given image with an irregular 'transparent' border to fill the rectangle.
 * Note that the 'transparent' border must NOT have any boundary as long as its dimension.
 */
Image* ip_local_warp(Image* src)
{
    vector<vector<Coordinate>> displacementMap = ip_local_warp_displacement(src);
    Image* result = new Image(*src);
    for (int x = 0; x < src->getWidth(); x++) {
        for (int y = 0; y < src->getHeight(); y++) {
            Coordinate displacement = displacementMap[x][y];
            Pixel srcPixel = src->getPixel(x + displacement.x, y + displacement.y);
            result->setPixel(x, y, srcPixel);
        }
    }
    return result;
}

/*
 * Performs matrix multiplication. Input vectors must be row-major order (m1[y][x]).
 */
vector<vector<double>> matmul(vector<vector<double>> m1, vector<vector<double>> m2)
{
    vector<vector<double>> result;
    for (int i = 0; i < m1.size(); i++) {
        vector<double> row;
        for (int j = 0; j < m2[0].size(); j++) {
            double val = 0;
            for (int k = 0; k < m1[0].size(); k++) {
                val += m1[i][k] * m2[k][j];
            }
            row.push_back(val);
        }
        result.push_back(row);
    }
    return result;
}

/*
 * Transposes the given matrix.
 */
vector<vector<double>> mattrans(vector<vector<double>> m)
{
    vector<vector<double>> t;
    for (int j = 0; j < m[0].size(); j++) {
        vector<double> row;
        for (int i = 0; i < m.size(); i++) {
            row.push_back(m[i][j]);
        }
        t.push_back(row);
    }
    return t;
}

/*
 * Subtracts the second matrix from the first.
 */
vector<vector<double>> matsub(vector<vector<double>> m1, vector<vector<double>> m2)
{
    vector<vector<double>> dif;
    for (int i = 0; i < m1.size(); i++) {
        vector<double> row;
        for (int j = 0; j < m1[0].size(); j++) {
            row.push_back(m1[i][j] - m2[i][j]);
        }
        dif.push_back(row);
    }
    return dif;
}

/*
 * Returns an identity matrix of the given size.
 */
vector<vector<double>> matiden(int order)
{
    vector<vector<double>> iden;
    for (int i = 0; i < order; i++) {
        vector<double> row;
        for (int j = 0; j < order; j++) {
            if (i == j) {
                row.push_back(1);
            } else {
                row.push_back(0);
            }
        }
        iden.push_back(row);
    }
    return iden;
}

/*
 * Calculates the magnitude of the given x by 1 matrix (a vector, in matrix speak).
 */
double matmag(vector<vector<double>> m)
{
    double mag = 0;
    for (int i = 0; i < m.size(); i++) {
        mag += pow(m[i][0], 2);
    }
    mag = sqrt(mag);
    return mag;
}

/*
 * Gives the "serialized" index for the given vertex coordinate dimension.
 * The ordering is x0, y0, x1, y1, ... . The vertices are ordered across and down (row-major).
 *   x0y0 x1y1 x2y2
 *   x3y3 x4y4 x5y5
 */
int variableIndexForCoordDimen(int vertexX, int vertexY, bool y, size_t numX)
{
    int vertexIndex = (int) numX * vertexY + vertexX;
    int variableIndex = vertexIndex * 2;
    if (y) {
        variableIndex += 1;
    }
    return variableIndex;
}
/*
 * Returns the quadratic terms for the translation energy for the given mesh.
 * I came up with this energy function to encourage vertices to stay anchored to their original
 * coordinates, which seems to help ensure reasonable results.
 */
Quadratic ip_energy_translation(vector<vector<CoordinateDouble>> mesh)
{
    size_t numMeshX = mesh[0].size();
    size_t numMeshY = mesh.size();
    
    int numVariables = (int) (numMeshX * numMeshY * 2);
    Quadratic energyQuadratic = Quadratic(numVariables);
    
    // Add a quadratic term about each vertex's original position in the mesh
    // (z_i - (z_i^))^2 for all vertices
    //  = (z_i)^2 - 2*(z_i)*(z_i^) + (z_i^)^2
    for (int i = 0; i < numMeshY; i++) {
        for (int j = 0; j < numMeshX; j++) {
            CoordinateDouble coord = mesh[i][j];
            //x
            int xVarIndex = variableIndexForCoordDimen(j, i, false, numMeshX);
            energyQuadratic.setCoeffForVars(xVarIndex, xVarIndex, 1);
            energyQuadratic.setCoeffForVar(xVarIndex, -2 * coord.x);
            //y
            int yVarIndex = variableIndexForCoordDimen(j, i, true, numMeshX);
            energyQuadratic.setCoeffForVars(yVarIndex, yVarIndex, 1);
            energyQuadratic.setCoeffForVar(yVarIndex, -2 * coord.y);
        }
    }
    
    return energyQuadratic;
}

/*
 * Returns the quadratic terms for the shape energy for the given mesh.
 */
Quadratic ip_energy_shape(vector<vector<CoordinateDouble>> mesh)
{
    size_t numMeshX = mesh[0].size();
    size_t numMeshY = mesh.size();
    
    int numVariables = (int) (numMeshX * numMeshY * 2);
    Quadratic energyQuadratic = Quadratic(numVariables);

    size_t numberOfQuads = (numMeshX - 1) * (numMeshY - 1);
    
    // iterate through each quad; each quad is defined by its upper-left (small y, small x) corner
    for (int row = 0; row < numMeshY - 1; row++) {
        for (int col = 0; col < numMeshX - 1; col++) {
            // 8x4 matrix aq = [[x0, -y0, 1, 0], [y0, x0, 0, 1], ..., [x3, -y3, 1, 0], [y3, x3, 0, 1]] of the input quad
            vector<vector<double>> aq;
            
            CoordinateDouble v0in = mesh[row][col]; // upper left
            
            vector<double> row0;
            row0.push_back(v0in.x);
            row0.push_back(-v0in.y);
            row0.push_back(1);
            row0.push_back(0);
            aq.push_back(row0);
            
            vector<double> row1;
            row1.push_back(v0in.y);
            row1.push_back(v0in.x);
            row1.push_back(0);
            row1.push_back(1);
            aq.push_back(row1);
            
            CoordinateDouble v1in = mesh[row][col+1]; // upper right
            
            vector<double> row2;
            row2.push_back(v1in.x);
            row2.push_back(-v1in.y);
            row2.push_back(1);
            row2.push_back(0);
            aq.push_back(row2);
            
            vector<double> row3;
            row3.push_back(v1in.y);
            row3.push_back(v1in.x);
            row3.push_back(0);
            row3.push_back(1);
            aq.push_back(row3);
            
            CoordinateDouble v2in = mesh[row+1][col]; // lower left
            
            vector<double> row4;
            row4.push_back(v2in.x);
            row4.push_back(-v2in.y);
            row4.push_back(1);
            row4.push_back(0);
            aq.push_back(row4);
            
            vector<double> row5;
            row5.push_back(v2in.y);
            row5.push_back(v2in.x);
            row5.push_back(0);
            row5.push_back(1);
            aq.push_back(row5);
            
            CoordinateDouble v3in = mesh[row+1][col+1]; // lower right
            
            vector<double> row6;
            row6.push_back(v3in.x);
            row6.push_back(-v3in.y);
            row6.push_back(1);
            row6.push_back(0);
            aq.push_back(row6);
            
            vector<double> row7;
            row7.push_back(v3in.y);
            row7.push_back(v3in.x);
            row7.push_back(0);
            row7.push_back(1);
            aq.push_back(row7);
            // end aq
            
            
            // add ((aq*(aq^T*aq)^(-1)*aq^T-I)*vq)^2
            vector<vector<double>> aqt = mattrans(aq);
            vector<vector<double>> aqtXaq = matmul(aqt, aq);
            vector<vector<double>> inv = matinvGJ(aqtXaq);
            vector<vector<double>> prdr = matmul(inv, aqt);
            vector<vector<double>> prdl = matmul(aq, prdr);
            vector<vector<double>> iden = matiden(8);
            vector<vector<double>> mat = matsub(prdl, iden);
            
            // get the indices of the four corners of this quadrilateral
            int topLeftIndexX = variableIndexForCoordDimen(col, row, false, numMeshX);
            int topLeftIndexY = variableIndexForCoordDimen(col, row, true, numMeshX);
            int topRightIndexX = variableIndexForCoordDimen(col+1, row, false, numMeshX);
            int topRightIndexY = variableIndexForCoordDimen(col+1, row, true, numMeshX);
            int bottomLeftIndexX = variableIndexForCoordDimen(col, row+1, false, numMeshX);
            int bottomLeftIndexY = variableIndexForCoordDimen(col, row+1, true, numMeshX);
            int bottomRightIndexX = variableIndexForCoordDimen(col+1, row+1, false, numMeshX);
            int bottomRightIndexY = variableIndexForCoordDimen(col+1, row+1, true, numMeshX);
            int quadCoordinateVarIndices[8] = {topLeftIndexX, topLeftIndexY,
                                               topRightIndexX, topRightIndexY,
                                               bottomLeftIndexX, bottomLeftIndexY,
                                               bottomRightIndexX, bottomRightIndexY};
            
            for (int ii = 0; ii < 8; ii++) {
                // the contributions of this row of the product (r_ii)
                
                // all terms
                for (int jj = 0; jj < 8; jj++) {
                    for (int kk = 0; kk < 8; kk++) {
                        // the squared term of quadCoordinateVarIndices[jj]
                        energyQuadratic.incrementCoeffForVars(quadCoordinateVarIndices[jj],
                                                              quadCoordinateVarIndices[kk],
                                                              mat[ii][jj] * mat[ii][kk]);
                    }
                }
            }
            
        }
    }
    
    // normalize
    energyQuadratic.scaleCoefficients(1.0 / numberOfQuads);
    
    return energyQuadratic;
}

/*
 * Returns the quadratic coefficients of the boundary energy of the given mesh on the given image. 
 * The boundary energy is zero when all the edge vertices are at the appropriate edge of the 
 * srcImage (from its dimensions).
 */
Quadratic ip_energy_boundary(Image* srcImage, vector<vector<CoordinateDouble>> mesh)
{
    int width = srcImage->getWidth();
    int height = srcImage->getHeight();
    size_t numMeshX = mesh[0].size();
    size_t numMeshY = mesh.size();
    
    int numVariables = (int) (numMeshX * numMeshY * 2);
    Quadratic energyQuadratic = Quadratic(numVariables);

    // top
    for (int i = 0; i < numMeshX; i++) {
        // (y_i)^2 for all vertices in the top row (y = 0)
        int yVarIndex = variableIndexForCoordDimen(i, 0, true, numMeshX);
        energyQuadratic.setCoeffForVars(yVarIndex, yVarIndex, 1);
    }
    // bottom
    for (int i = 0; i < numMeshX; i++) {
        // (y_i - (height - 1))^2 for all vertices in the bottom row (y = height - 1)
        //  = (y_i)^2 - 2*(y_i)*(height - 1) + (height - 1)^2
        int yVarIndex = variableIndexForCoordDimen(i, (int) numMeshY - 1, true, numMeshX);
        energyQuadratic.setCoeffForVars(yVarIndex, yVarIndex, 1);
        energyQuadratic.setCoeffForVar(yVarIndex, -2 * (height - 1));
        // neglect the constant term
    }
    // left
    for (int j = 0; j < numMeshY; j++) {
        // (x_j)^2 for all vertices in the left col (x = 0)
        int xVarIndex = variableIndexForCoordDimen(0, j, false, numMeshX);
        energyQuadratic.setCoeffForVars(xVarIndex, xVarIndex, 1);
    }
    // right
    for (int j = 0; j < numMeshY; j++) {
        // (x_j - (width - 1))^2 for all vertices in the right col (x = width - 1)
        //  = (x_j)^2 - 2*(x_j)*(width - 1) + (width - 1)^2
        int xVarIndex = variableIndexForCoordDimen((int) numMeshX - 1, j, false, numMeshX);
        energyQuadratic.setCoeffForVars(xVarIndex, xVarIndex, 1);
        energyQuadratic.setCoeffForVar(xVarIndex, -2 * (width - 1));
        // neglect the constant term
    }
    
    return energyQuadratic;
}

/*
 * Converts the given Image to a PGM image, a C-style array of doubles (to interface with LSD).
 */
double* ip_image_to_array(Image* srcImage) {
    int x, y;
    int X = srcImage->getWidth();  // x image size
    int Y = srcImage->getHeight();  // y image size
    
    double* imageArray = new double[X * Y];
    
    // the PGM image is all gray, where 0=black and 255=white
    int whiteValue = 255;
    
    for (x = 0; x < X; x++) {
        for (y = 0; y < Y; y++) {
            Pixel srcPixel = srcImage->getPixel(x, y);
            Pixel grayPixel = ip_gray_pixel(srcPixel, nullptr);
            double pixelGrayness = grayPixel.getColor(RED);
            int doubleGrayness = pixelGrayness * whiteValue;
            imageArray[x + y * X] = doubleGrayness;
        }
    }
    
    return imageArray;
}

/*
 * Determines whether the given point lies within the quad given by the vertices.
 */
bool is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight,
                CoordinateDouble bottomLeft, CoordinateDouble bottomRight)
{
    // the point must be to the right of the left line, below the top line, above the bottom line,
    // and to the left of the right line
    
    // must be right of left line
    if (topLeft.x == bottomLeft.x) {
        if (point.x < topLeft.x) {
            return false;
        }
    } else {
        double leftSlope = (topLeft.y - bottomLeft.y) / (topLeft.x - bottomLeft.x);
        double leftIntersect = topLeft.y - leftSlope * topLeft.x;
        double yOnLineX = (point.y - leftIntersect) / leftSlope;
        if (point.x < yOnLineX) {
            return false;
        }
    }
    // must be left of right line
    if (topRight.x == bottomRight.x) {
        if (point.x > topRight.x) {
            return false;
        }
    } else {
        double rightSlope = (topRight.y - bottomRight.y) / (topRight.x - bottomRight.x);
        double rightIntersect = topRight.y - rightSlope * topRight.x;
        double yOnLineX = (point.y - rightIntersect) / rightSlope;
        if (point.x > yOnLineX) {
            return false;
        }
    }
    // must be below top line
    if (topLeft.y == topRight.y) {
        if (point.y < topLeft.y) {
            return false;
        }
    } else {
        double topSlope = (topRight.y - topLeft.y) / (topRight.x - topLeft.x);
        double topIntersect = topLeft.y - topSlope * topLeft.x;
        double xOnLineY = topSlope * point.x + topIntersect;
        if (point.y < xOnLineY) {
            return false;
        }
    }
    // must be above bottom line
    if (bottomLeft.y == bottomRight.y) {
        if (point.y > bottomLeft.y) {
            return false;
        }
    } else {
        double bottomSlope = (bottomRight.y - bottomLeft.y) / (bottomRight.x - bottomLeft.x);
        double bottomIntersect = bottomLeft.y - bottomSlope * bottomLeft.x;
        double xOnLineY = bottomSlope * point.x + bottomIntersect;
        if (point.y > xOnLineY) {
            return false;
        }
    }
    // if all four constraints are satisfied, the point must be in the quad
    return true;
}

/*
 * Returns whether the given line segment intersects the line of given slope and intersect.
 * Needs to know if the line is vertical (else horizontal).
 * If there is an intersection, sets the intersectPoint to the intersection point.
 */
bool does_segment_intersect_line(LineSegment lineSegment, double slope, double intersect,
                                 bool vertical, CoordinateDouble& intersectPoint)
{
    // calculate line segment m and b
    double lineSegmentSlope = INFT;
    if (lineSegment.x1 != lineSegment.x2) {
        lineSegmentSlope = (lineSegment.y2 - lineSegment.y1) / (lineSegment.x2 - lineSegment.x1);
    }
    double lineSegmentIntersect = lineSegment.y1 - lineSegmentSlope * lineSegment.x1;
    
    // calculate intersection
    if (lineSegmentSlope == slope) {
        if (lineSegmentIntersect == intersect) {
            // same line
            intersectPoint.x = lineSegment.x1;
            intersectPoint.y = lineSegment.y1;
            return true;
        } else {
            return false;
        }
    }
    double intersectX = (intersect - lineSegmentIntersect) / (lineSegmentSlope - slope);
    double intersectY = lineSegmentSlope * intersectX + lineSegmentIntersect;
    // check if intersection is in the bounds of the line segment
    if (vertical) {
        if ((intersectY <= lineSegment.y1 && intersectY >= lineSegment.y2) ||
            (intersectY <= lineSegment.y2 && intersectY >= lineSegment.y1)) {
            intersectPoint.x = intersectX;
            intersectPoint.y = intersectY;
            return true;
        } else {
            return false;
        }
    } else {
        if ((intersectX <= lineSegment.x1 && intersectX >= lineSegment.x2) ||
            (intersectX <= lineSegment.x2 && intersectX >= lineSegment.x1)) {
            intersectPoint.x = intersectX;
            intersectPoint.y = intersectY;
            return true;
        } else {
            return false;
        }
    }
}

/*
 * Returns the point(s) at which the line segment intersects the quad.
 */
vector<CoordinateDouble> intersections_with_quad(LineSegment lineSegment, CoordinateDouble topLeft,
                                                 CoordinateDouble topRight, CoordinateDouble bottomLeft,
                                                 CoordinateDouble bottomRight)
{
    vector<CoordinateDouble> intersections;
    // find the intersections of lineSegment with each of the four quad boundaries, then add each
    // intersection that is within the appropriate dimension's endpoints of the quad boundary
    
    // left
    double leftSlope = INFT;
    if (topLeft.x != bottomLeft.x) {
        leftSlope = (topLeft.y - bottomLeft.y) / (topLeft.x - bottomLeft.x);
    }
    double leftIntersect = topLeft.y - leftSlope * topLeft.x;
    // check
    CoordinateDouble leftIntersectPoint;
    if (does_segment_intersect_line(lineSegment, leftSlope, leftIntersect, true, leftIntersectPoint)) {
        if (leftIntersectPoint.y >= topLeft.y && leftIntersectPoint.y <= bottomLeft.y) {
            intersections.push_back(leftIntersectPoint);
        }
    }
    
    // right
    double rightSlope = INFT;
    if (topRight.x != bottomRight.x) {
        rightSlope = (topRight.y - bottomRight.y) / (topRight.x - bottomRight.x);
    }
    double rightIntersect = topRight.y - rightSlope * topRight.x;
    // check
    CoordinateDouble rightIntersectPoint;
    if (does_segment_intersect_line(lineSegment, rightSlope, rightIntersect, true, rightIntersectPoint)) {
        if (rightIntersectPoint.y >= topRight.y && rightIntersectPoint.y <= bottomRight.y) {
            intersections.push_back(rightIntersectPoint);
        }
    }
    
    // top
    double topSlope = INFT;
    if (topLeft.x != topRight.x) {
        topSlope = (topRight.y - topLeft.y) / (topRight.x - topLeft.x);
    }
    double topIntersect = topLeft.y - topSlope * topLeft.x;
    // check
    CoordinateDouble topIntersectPoint;
    if (does_segment_intersect_line(lineSegment, topSlope, topIntersect, false, topIntersectPoint)) {
        if (topIntersectPoint.x >= topLeft.x && topIntersectPoint.x <= topRight.x) {
            intersections.push_back(topIntersectPoint);
        }
    }
    
    // bottom
    double bottomSlope = INFT;
    if (bottomLeft.x != bottomRight.x) {
        bottomSlope = (bottomRight.y - bottomLeft.y) / (bottomRight.x - bottomLeft.x);
    }
    double bottomIntersect = bottomLeft.y - bottomSlope * bottomLeft.x;
    // check
    CoordinateDouble bottomIntersectPoint;
    if (does_segment_intersect_line(lineSegment, bottomSlope, bottomIntersect, false, bottomIntersectPoint)) {
        if (bottomIntersectPoint.x >= bottomLeft.x && bottomIntersectPoint.x <= bottomRight.x) {
            intersections.push_back(bottomIntersectPoint);
        }
    }
    
    return intersections;
}

/*
 * Given an array of line segments and a mesh, this method returns a vector like a mesh but where the
 * element in result[i][j] is a vector of line segments where each line segment is a segment from
 * lineSegments that has been cut to fit within the quad.
 */
vector<vector<vector<LineSegment>>> ip_get_line_segments_in_mesh(double* lineSegments, int numLines,
                                                                 vector<vector<CoordinateDouble>> mesh)
{
    vector<vector<vector<LineSegment>>> lineSegmentsInMesh;
    
    // for each quad (i,j)
    // for each line segment in lineSegments
    // Find the resulting line segment cut by/into this quad - if it exists, add it to lineSegmentsInMesh[i][j]
    // If both endpoints are in the quad, just add the line segment as-is (easy)
    // if both endpoints is_in_quad
    // If exactly one endpoint is in the quad, determine where the line segment leaves the quad
    // else if one endpoint is_in_quad, call intersection_with_quad
    // If no endpoints are in the quad, determine its intersection points with the quad (if any)
    // else if no endpoints are in the quad
    // try intersections_with_quad - if two intersections, use them; if none, skip this line
    
    size_t numMeshX = mesh[0].size();
    size_t numMeshY = mesh.size();
    // TODO: fix switching x and y (it's now mesh[rol][col] = mesh[y][x])
    for (int i = 0; i < numMeshX - 1; i++) {
        vector<vector<LineSegment>> row;
        for (int j = 0; j < numMeshY - 1; j++) {
            // now we have a quad, mesh[i][j]
            CoordinateDouble topLeft = mesh[i][j];
            CoordinateDouble topRight = mesh[i+1][j];
            CoordinateDouble bottomLeft = mesh[i][j+1];
            CoordinateDouble bottomRight = mesh[i+1][j+1];
            vector<LineSegment> lineSegmentsInQuad;
            for (int lineNum = 0; lineNum < numLines; lineNum++) {
                LineSegment lineSegment;
                lineSegment.x1 = lineSegments[7*lineNum + 0];
                lineSegment.y1 = lineSegments[7*lineNum + 1];
                lineSegment.x2 = lineSegments[7*lineNum + 2];
                lineSegment.y2 = lineSegments[7*lineNum + 3];
                // now we have a lineSegment
                CoordinateDouble lineSegmentPt1;
                lineSegmentPt1.x = lineSegment.x1;
                lineSegmentPt1.y = lineSegment.y1;
                CoordinateDouble lineSegmentPt2;
                lineSegmentPt2.x = lineSegment.x2;
                lineSegmentPt2.y = lineSegment.y2;
                
                // check if the line segment cut by/into this quad exists
                bool pt1InQuad = is_in_quad(lineSegmentPt1, topLeft, topRight, bottomLeft, bottomRight);
                bool pt2InQuad = is_in_quad(lineSegmentPt2, topLeft, topRight, bottomLeft, bottomRight);
                if (pt1InQuad && pt2InQuad) {
                    lineSegmentsInQuad.push_back(lineSegment);
                } else if (pt1InQuad) {
                    vector<CoordinateDouble> intersections = intersections_with_quad(lineSegment, topLeft, topRight, bottomLeft, bottomRight);
                    if (intersections.size() != 0) {
                        LineSegment cutLineSegment;
                        cutLineSegment.x1 = lineSegmentPt1.x;
                        cutLineSegment.y1 = lineSegmentPt1.y;
                        cutLineSegment.x2 = intersections[0].x;
                        cutLineSegment.y2 = intersections[0].y;
                        lineSegmentsInQuad.push_back(cutLineSegment);
                    }
                } else if (pt2InQuad) {
                    vector<CoordinateDouble> intersections = intersections_with_quad(lineSegment, topLeft, topRight, bottomLeft, bottomRight);
                    if (intersections.size() != 0) {
                        LineSegment cutLineSegment;
                        cutLineSegment.x1 = lineSegmentPt2.x;
                        cutLineSegment.y1 = lineSegmentPt2.y;
                        cutLineSegment.x2 = intersections[0].x;
                        cutLineSegment.y2 = intersections[0].y;
                        lineSegmentsInQuad.push_back(cutLineSegment);
                    }
                } else {
                    vector<CoordinateDouble> intersections = intersections_with_quad(lineSegment, topLeft, topRight, bottomLeft, bottomRight);
                    if (intersections.size() != 0) {
                        LineSegment cutLineSegment;
                        cutLineSegment.x1 = intersections[0].x;
                        cutLineSegment.y1 = intersections[0].y;
                        cutLineSegment.x2 = intersections[1].x;
                        cutLineSegment.y2 = intersections[1].y;
                        lineSegmentsInQuad.push_back(cutLineSegment);
                    }
                }
                
            }
            row.push_back(lineSegmentsInQuad);
        }
        lineSegmentsInMesh.push_back(row);
    }

    return lineSegmentsInMesh;
}

/*
 * Returns the quadratic terms of the line energy of the mesh on the image with the given angle theta.
 */
double *ip_energy_line(Image* srcImage, vector<vector<CoordinateDouble>> mesh,
                       double *theta, int numBins)
{
    return nullptr;
//    double energy = 0;
//    int width = srcImage->getWidth();
//    int height = srcImage->getHeight();
//    
//    // get lineSegments
//    double* imageArray = ip_image_to_array(srcImage);
//    int numLines;
//    double* lineSegments = lsd(&numLines, imageArray, width, height);
//
//    vector<vector<vector<LineSegment>>> lineSegmentsInMesh = ip_get_line_segments_in_mesh(lineSegments, numLines, inputMesh);
//    
//    int numLinesOutput; // should be the same as numLines, but might not be
//    double* lineSegmentsOutput = lsd(&numLinesOutput, imageArray, width, height);
//    vector<vector<vector<LineSegment>>> lineSegmentsInMeshOutput = ip_get_line_segments_in_mesh(lineSegmentsOutput, numLinesOutput, outputMesh);
//    // get a flat copy so we can determine corresponding line segments (for the most part..)
//    vector<LineSegment> lineSegmentsInMeshOutputFlat;
//    for (int i = 0; i < lineSegmentsInMeshOutput.size(); i++) {
//        for (int j = 0; j < lineSegmentsInMeshOutput[0].size(); j++) {
//            for (int k = 0; k < lineSegmentsInMeshOutput[i][j].size(); k++) {
//                lineSegmentsInMeshOutputFlat.push_back(lineSegmentsInMeshOutput[i][j][k]);
//            }
//        }
//    }
//    
//    double thetaPerBucket = M_PI / numBins;
//    
//    size_t numMeshX = inputMesh.size();
//    size_t numMeshY = inputMesh[0].size();
//    long numLineSegments = 0;
//    for (int i = 0; i < numMeshX - 1; i++) {
//        for (int j = 0; j < numMeshY - 1; j++) {
//            
//            CoordinateDouble topLeft = inputMesh[i][j];
//            CoordinateDouble topRight = inputMesh[i+1][j];
//            CoordinateDouble bottomLeft = inputMesh[i][j+1];
//            CoordinateDouble bottomRight = inputMesh[i+1][j+1];
//            
//            vector<LineSegment> lineSegmentsInQuad = lineSegmentsInMesh[i][j];
//            for (int k = 0; k < lineSegmentsInQuad.size(); k++) {
//                
//                // e is [[a,b], [c,d]] where a is the bilinear interpolation of horiz and b is vert; c,d are for pt 2
//                LineSegment lineSegment = lineSegmentsInQuad[k];
//                vector<vector<double>> eHat;
//                vector<double> eHatRow1;
//                eHatRow1.push_back((lineSegment.x1 - topLeft.x)/(topRight.x - topLeft.x));
//                eHatRow1.push_back((lineSegment.y1 - topLeft.y)/(bottomLeft.y - topLeft.y));
//                eHat.push_back(eHatRow1);
//                vector<double> eHatRow2;
//                eHatRow2.push_back((lineSegment.x2 - topLeft.x)/(topRight.x - topLeft.x));
//                eHatRow2.push_back((lineSegment.y2 - topLeft.y)/(bottomLeft.y - topLeft.y));
//                eHat.push_back(eHatRow2);
//                
//                LineSegment lineSegmentOutput = lineSegmentsInMeshOutputFlat[numLineSegments];
//                vector<vector<double>> e;
//                vector<double> eRow1;
//                eRow1.push_back((lineSegmentOutput.x1 - topLeft.x)/(topRight.x - topLeft.x));
//                eRow1.push_back((lineSegmentOutput.y1 - topLeft.y)/(bottomLeft.y - topLeft.y));
//                e.push_back(eRow1);
//                vector<double> eRow2;
//                eRow2.push_back((lineSegmentOutput.x2 - topLeft.x)/(topRight.x - topLeft.x));
//                eRow2.push_back((lineSegmentOutput.y2 - topLeft.y)/(bottomLeft.y - topLeft.y));
//                e.push_back(eRow2);
//                
//                double lineSegmentTheta = atan((lineSegmentOutput.y2 - lineSegmentOutput.y1) /
//                                               (lineSegmentOutput.x2 - lineSegmentOutput.x1));
//                
//                int lineSegmentBucket = (lineSegmentTheta + M_PI / 2) / thetaPerBucket;
//                double targetTheta = theta[lineSegmentBucket];
//                vector<vector<double>> r;
//                vector<double> rRow1;
//                rRow1.push_back(cos(targetTheta));
//                rRow1.push_back(-sin(targetTheta));
//                r.push_back(rRow1);
//                vector<double> rRow2;
//                rRow2.push_back(sin(targetTheta));
//                rRow2.push_back(cos(targetTheta));
//                
//                // C = R*eHat*(eHat^T * eHat)^(-1)*eHat^T*R^T-I
//                vector<vector<double>> ete = matmul(mattrans(eHat), eHat);
//                vector<vector<double>> eteInv = matinv(ete); matinvMinor
//                vector<vector<double>> prdr = matmul(eteInv, mattrans(eHat));
//                vector<vector<double>> prdrr = matmul(prdr, mattrans(r));
//                vector<vector<double>> prdl = matmul(eHat, prdrr);
//                vector<vector<double>> prdll = matmul(r, prdl);
//                vector<vector<double>> c = matsub(prdll, matiden(2));
//                
//                // add (C*e)^2
//                vector<vector<double>> ce = matmul(c, e);
//                double mag = matmag(ce);
//                double magSqrd = pow(mag, 2);
//                energy += magSqrd;
//                
//                numLineSegments++;
//            }
//        }
//    }
//    
//    // normalize
//    energy = energy / numLineSegments;
//    
//    return energy;
}

/*
 * Returns the quadratic terms of the energy of the mesh on the image with the given angle theta.
 */
Quadratic ip_energy_total(Image* srcImage, vector<vector<CoordinateDouble>> mesh,
                        double *theta, int numBins)
{
    size_t numMeshX = mesh[0].size();
    size_t numMeshY = mesh.size();
    int numVariables = (int) (numMeshX * numMeshY * 2);
    Quadratic energyQuadratic = Quadratic(numVariables);
    
    // Add shape energy
    Quadratic shapeEnergyQuadratic = ip_energy_shape(mesh);
    double esWeight = 10;
    energyQuadratic.addQuadratic(shapeEnergyQuadratic, esWeight);
    
    // Add weighted line energy
//    double elWeight = 100;
//    double *lineEnergy = ip_energy_line(srcImage, mesh, theta, numBins);
//    for (int i = 0; i < numTermsInQuadratic; i++) {
//        energyQuadraticCoefficients[i] += elWeight * lineEnergy[i];
//    }
    
    // Add weighted boundary energy
    Quadratic boundaryEnergyQuadratic = ip_energy_boundary(srcImage, mesh);
    double ebWeight = INFT;
    energyQuadratic.addQuadratic(boundaryEnergyQuadratic, ebWeight);
    
    // Add my custom translation energy
    Quadratic translationEnergyQuadratic = ip_energy_translation(mesh);
    double etWeight = 1;
    energyQuadratic.addQuadratic(translationEnergyQuadratic, etWeight);
    
    return energyQuadratic;
}

/*
 * Draws the vertices of the given mesh on the given image.
 * Vertices on a boundary are drawn in blue; all other vertices are green.
 */
Image* ip_draw_vertices(Image* src, vector<vector<CoordinateDouble>> mesh)
{
    size_t numMeshX = mesh[0].size();
    size_t numMeshY = mesh.size();
    Image* result = new Image(*src);
    // set non-boundaries
    for (int i = 0; i < numMeshY; i++) {
        for (int j = 0; j < numMeshX; j++) {
            CoordinateDouble meshVertexCoord = mesh[i][j];
            Pixel greenPixel = Pixel(0, 1, 0);
            if (i == 0 || i == numMeshY - 1 || j == 0 || j == numMeshX - 1) {
                continue;
            }
            result->setPixel(clamp(meshVertexCoord.x, 0, src->getWidth()-1), clamp(meshVertexCoord.y, 0, src->getHeight() - 1), greenPixel);
        }
    }
    // set boundaries
    for (int i = 0; i < numMeshY; i++) {
        for (int j = 0; j < numMeshX; j++) {
            CoordinateDouble meshVertexCoord = mesh[i][j];
            Pixel bluePixel = Pixel(0, 0, 1);
            if (i == 0 || i == numMeshY - 1 || j == 0 || j == numMeshX - 1) {
                result->setPixel(clamp(meshVertexCoord.x, 0, src->getWidth()-1), clamp(meshVertexCoord.y, 0, src->getHeight() - 1), bluePixel);
            }
        }
    }
    return result;
}

/*
 * Fits the given src image to its rectangular boundaries, using white as transparent.
 * This is the main method of the project.
 */
Image* ip_rectangle(Image* srcImage)
{
    // === Local warp to get original mesh ===
    
    cerr << " Starting local warp" << endl;
    
    // First, use local warp to get a rectangular displacement map
    vector<vector<Coordinate>> displacementMap = ip_local_warp_displacement(srcImage);
    
    // Next, make a mesh on the rectangular image..
    int width = srcImage->getWidth();
    int height = srcImage->getHeight();
    int numMeshX = 20;
    int numMeshY = 20;
    double widthPerMesh = double(width-1) / (numMeshX - 1);
    double heightPerMesh = double(height-1) / (numMeshY - 1);
    vector<vector<CoordinateDouble>> mesh;
    for (int y = 0; y < numMeshY; y++) {
        vector<CoordinateDouble> meshRow;
        for (int x = 0; x < numMeshX; x++) {
            CoordinateDouble coord;
            coord.x = x * widthPerMesh;
            coord.y = y * heightPerMesh;
            meshRow.push_back(coord);
        }
        mesh.push_back(meshRow);
    }
    // ..and warp the mesh backward to the original image
    for (int i = 0; i < numMeshY; i++) {
        for (int j = 0; j < numMeshX; j++) {
            CoordinateDouble& meshVertexCoord = mesh[i][j];
            Coordinate vertexDisplacement = displacementMap[floor(meshVertexCoord.x)][floor(meshVertexCoord.y)];
            meshVertexCoord.x += vertexDisplacement.x;
            meshVertexCoord.y += vertexDisplacement.y;

        }
    }
    
    cerr << " Finished local warp." << endl;
    
    // === End local warp to get original mesh ===
    
    
    // === Optimize mesh ===
    
    cerr << " Starting mesh optimization." << endl;
    
    // Get equation representing the energy of an output mesh, then optimize it with math.
    
    // TO-DO: Something with this theta
    // set theta array to 50 buckets from -pi/2 to pi/2
    int numBins = 50;
    double *theta = new double[numBins];
    double thetaPerBucket = M_PI / numBins;
    for (int i = 0; i < numBins; i++) {
        theta[i] = i * thetaPerBucket;
    }
    
    // Get energy
    Quadratic energyQuadratic = ip_energy_total(srcImage, mesh, theta, numBins);
//    energyQuadratic.roundCoefficients(0.000001); // round for readability; has little impact on result
//    energyQuadratic.printEquation();
    // Get derivative
    double **derivatives = energyQuadratic.derivatives();
//    for (int i = 0; i < (numMeshX * numMeshY * 2); i++) {
//        for (int j = 0; j < (numMeshX * numMeshY * 2) + 1; j++) {
//            cerr << derivatives[i][j] << " ";
//        }
//        cerr << endl;
//    }
    // The only issue with this matrix is that it's one col too wide, due to each equation's
    // constant term. We now make a vertical matrix and a square matrix.
    int numVariables = (int) (numMeshX * numMeshY * 2);
    vector<vector<double>> colVector;
    for (int i = 0; i < numVariables; i++) {
        vector<double> row;
        row.push_back(-1 * derivatives[i][numVariables]);
        colVector.push_back(row);
    }
    vector<vector<double>>squareMatrix;
    for (int i = 0; i < numVariables; i++) {
        // i is like y (row)
        vector<double> row;
        for (int j = 0; j < numVariables; j++) {
            // j is like x (col)
            row.push_back(derivatives[i][j]);
        }
        squareMatrix.push_back(row);
    }
    // Now we can solve the system of linear equations.
    vector<vector<double>>inverseSquareMatrix = matinvGJ(squareMatrix);
    vector<vector<double>>solutionVariables = matmul(inverseSquareMatrix, colVector);
    // Finally, we can convert the list of variables to coordinates.
    vector<vector<CoordinateDouble>> outputMesh;
    for (int y = 0; y < numMeshY; y++) {
        vector<CoordinateDouble> outputMeshRow;
        for (int x = 0; x < numMeshX; x++) {
            // note that this loop repeats half the pairs, but that's okay
            CoordinateDouble coord;
            int xVariableIndex = variableIndexForCoordDimen(x, y, false, numMeshX);
            int yVariableIndex = variableIndexForCoordDimen(x, y, true, numMeshX);
            coord.x = solutionVariables[xVariableIndex][0];
            coord.y = solutionVariables[yVariableIndex][0];
            outputMeshRow.push_back(coord);
        }
        outputMesh.push_back(outputMeshRow);
    }
    
    cerr << " Finished mesh optimization." << endl;
    
//    cerr << "Output mesh:" << endl;
//    for (int y = 0; y < numMeshY; y++) {
//        for (int x = 0; x < numMeshX; x++) {
//            CoordinateDouble coord = outputMesh[y][x];
//            cerr << "(" << coord.x << ", " << coord.y << "), ";
//        }
//        cerr << endl;
//    }
    
    return ip_draw_vertices(srcImage, outputMesh);
    
    // === End optimize mesh ===
    
    // Use optimized mesh to interpolate displacement of every pixel; fill in few blank pixels
    
    // Scale the final mesh and re-interpolate to fix stretching
    
    return nullptr;
}
