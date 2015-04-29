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
 * Calculates the shape energy for the given set of meshes (input mesh and output mesh).
 */
double ip_energy_shape(vector<vector<CoordinateDouble>> inputMesh,
                       vector<vector<CoordinateDouble>> outputMesh)
{
    size_t numMeshX = inputMesh.size();
    size_t numMeshY = inputMesh[0].size();
    double energy = 0;
    size_t numberOfQuads = (numMeshX - 1) * (numMeshY - 1);
    
    // iterate through each quad; each quad is defined by its upper-left (small y, small x) corner
    for (int i = 0; i < numMeshX - 1; i++) {
        for (int j = 0; j < numMeshY - 1; j++) {
            
            // 8x4 matrix aq = [[x0, -y0, 1, 0], [y0, x0, 0, 1], ..., [x3, -y3, 1, 0], [y3, x3, 0, 1]] of the input quad
            vector<vector<double>> aq;
            
            CoordinateDouble v0in = inputMesh[i][j]; // upper left
            
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
            
            CoordinateDouble v1in = inputMesh[i+1][j]; // upper right
            
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
            
            CoordinateDouble v2in = inputMesh[i][j+1]; // lower left
            
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
            
            CoordinateDouble v3in = inputMesh[i+1][j+1]; // lower right
            
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
            
            
            // 8x1 vector vq = [x0, y0, x1, y1, x2, y2, x3, y3], the four vertices of the output quad
            vector<vector<double>> vq;
            
            CoordinateDouble v0out = outputMesh[i][j]; // upper left
            
            vector<double> vrow0;
            vrow0.push_back(v0out.x);
            vq.push_back(vrow0);
            
            vector<double> vrow1;
            vrow1.push_back(v0out.y);
            vq.push_back(vrow1);
            
            CoordinateDouble v1out = outputMesh[i+1][j]; // upper right
            
            vector<double> vrow2;
            vrow2.push_back(v1out.x);
            vq.push_back(vrow2);
            
            vector<double> vrow3;
            vrow3.push_back(v1out.y);
            vq.push_back(vrow3);
            
            CoordinateDouble v2out = outputMesh[i][j+1]; // lower left
            
            vector<double> vrow4;
            vrow4.push_back(v2out.x);
            vq.push_back(vrow4);
            
            vector<double> vrow5;
            vrow5.push_back(v2out.y);
            vq.push_back(vrow5);
            
            CoordinateDouble v3out = outputMesh[i+1][j+1]; // lower right

            vector<double> vrow6;
            vrow6.push_back(v3out.x);
            vq.push_back(vrow6);
            
            vector<double> vrow7;
            vrow7.push_back(v3out.y);
            vq.push_back(vrow7);
            // end vq
            
            
            // add ((aq*(aq^T*aq)^(-1)*aq^T-I)*vq)^2
            vector<vector<double>> aqt = mattrans(aq);
            vector<vector<double>> aqtXaq = matmul(aqt, aq);
            vector<vector<double>> inv = matinv(aqtXaq);
            vector<vector<double>> prdr = matmul(inv, aqt);
            vector<vector<double>> prdl = matmul(aq, prdr);
            vector<vector<double>> iden = matiden(8);
            vector<vector<double>> dif = matsub(prdl, iden);
            vector<vector<double>> term = matmul(dif, vq);
            double mag = matmag(term);
            double magSqrd = pow(mag, 2);
            energy += magSqrd;
        }
    }
    
    // normalize
    energy = energy / numberOfQuads;
    
    return energy;
}

/*
 * Returns the boundary energy of the given mesh on the given image. The boundary energy is zero
 * when all the edge vertices are at the appropriate edge of the srcImage (from its dimensions).
 */
double ip_energy_boundary(Image* srcImage, vector<vector<CoordinateDouble>> outputMesh)
{
    int width = srcImage->getWidth();
    int height = srcImage->getHeight();
    size_t numMeshX = outputMesh.size();
    size_t numMeshY = outputMesh[0].size();
    double energy = 0;
    // top
    for (int i = 0; i < numMeshX; i++) {
        CoordinateDouble meshVertexCoord = outputMesh[i][0];
        energy += pow(meshVertexCoord.y, 2);
    }
    // bottom
    for (int i = 0; i < numMeshX; i++) {
        CoordinateDouble meshVertexCoord = outputMesh[i][numMeshY-1];
        energy += pow((height - 1) - meshVertexCoord.y, 2);
    }
    // left
    for (int j = 0; j < numMeshY; j++) {
        CoordinateDouble meshVertexCoord = outputMesh[0][j];
        energy += pow(meshVertexCoord.x, 2);
    }
    // right
    for (int j = 0; j < numMeshY; j++) {
        CoordinateDouble meshVertexCoord = outputMesh[numMeshX-1][j];
        energy += pow((width - 1) - meshVertexCoord.x, 2);
    }
    
    return energy;
}

/*
 * Returns the line energy of the meshes on the image with the given angle theta.
 */
double ip_energy_line(Image* srcImage, vector<vector<CoordinateDouble>> inputMesh,
                      vector<vector<CoordinateDouble>> outputMesh, double theta)
{
    1; // TO-DO: implement this method
    return 0;
}

/*
 * Returns the total energy of the meshes on the image with the given angle theta.
 */
double ip_energy_total(Image* srcImage, vector<vector<CoordinateDouble>> inputMesh,
                      vector<vector<CoordinateDouble>> outputMesh, double theta)
{
    double energy = 0;
    // Add shape energy
    double shapeEnergy = ip_energy_shape(inputMesh, outputMesh);
    double esWeight = 1;
    energy += esWeight * shapeEnergy;
    // Add weighted line energy
    double elWeight = 100;
    double lineEnergy = ip_energy_line(srcImage, inputMesh, outputMesh, theta);
    energy += elWeight * lineEnergy;
    // Add weighted boundary energy
    double boundaryEnergy = ip_energy_boundary(srcImage, outputMesh);
    double ebWeight = INFT;
    energy += ebWeight * boundaryEnergy;
    return energy;
}

/*
 * Draws the vertices of the given mesh on the given image.
 */
Image* ip_draw_vertices(Image* src, vector<vector<CoordinateDouble>> mesh)
{
    size_t numMeshX = mesh.size();
    size_t numMeshY = mesh[0].size();
    Image* result = new Image(*src);
    for (int i = 0; i < numMeshX; i++) {
        for (int j = 0; j < numMeshY; j++) {
            CoordinateDouble meshVertexCoord = mesh[i][j];
            Pixel greenPixel = Pixel(0, 1, 0);
            result->setPixel(floor(meshVertexCoord.x), floor(meshVertexCoord.y), greenPixel);
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
    for (int x = 0; x < numMeshX; x++) {
        vector<CoordinateDouble> meshRow;
        for (int y = 0; y < numMeshY; y++) {
            CoordinateDouble coord;
            coord.x = x * widthPerMesh;
            coord.y = y * heightPerMesh;
            meshRow.push_back(coord);
        }
        mesh.push_back(meshRow);
    }
    // ..and warp the mesh backward to the original image
    for (int i = 0; i < numMeshX; i++) {
        for (int j = 0; j < numMeshY; j++) {
            CoordinateDouble& meshVertexCoord = mesh[i][j];
            Coordinate vertexDisplacement = displacementMap[floor(meshVertexCoord.x)][floor(meshVertexCoord.y)];
            meshVertexCoord.x += vertexDisplacement.x;
            meshVertexCoord.y += vertexDisplacement.y;

        }
    }
    
    // === End local warp to get original mesh ===
    
    
    // === Optimize mesh ===
    
    // Compute energy function for this mesh
    vector<vector<CoordinateDouble>> inputMesh = mesh;
    // TO-DO: Something with this mesh
    vector<vector<CoordinateDouble>> outputMesh = mesh;
    // TO-DO: Something with this theta
    double theta = 0;
    double energy = ip_energy_total(srcImage, inputMesh, outputMesh, theta);
    
    cerr << "Energy: " << energy << endl;
    
    // Optimize energy function to get the optimized (eventually final) mesh
    
    // === End optimize mesh ===
    
    
    // Use optimized mesh to interpolate displacement of every pixel; fill in few blank pixels
    
    // Scale the final mesh and re-interpolate to fix stretching
    
    return nullptr;
}
