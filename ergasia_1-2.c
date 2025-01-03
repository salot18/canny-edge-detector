#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// SUNFLOWER
#define N 200 // height
#define M 200 // width
#define filename "sunflower_200x200_444.yuv"

#define _PI 3.14159265359
#define KERNEL_SIZE 3
#define TILE_SIZE 100
#define CACHE_ROWS 4

#pragma arm section zidata="l1data"
int blk[3][2*KERNEL_SIZE];
int conv_sum[2][4];
#pragma arm section

#pragma arm section zidata="l2data"
int cache[CACHE_ROWS][M];
int i, j, ii, jj;
int startCacheIdx = 0;
#pragma arm section

#pragma arm section zidata="ram"
int current_y[N][M];
double gaussian_kernel[KERNEL_SIZE];
int gradX[N][M];
int gradY[N][M];
int gradMag[N][M];
int gradDir[N][M];
#pragma arm section

const int sobel_kernel_x[3][3] = {
    {-1, 0, 1},
    {-2, 0, 2},
    {-1, 0, 1},
};

const int sobel_kernel_y[3][3] = {
    {-1, -2, -1},
    {0, 0, 0},
    {1, 2, 1},
};

const float gauss_kernel_1d[KERNEL_SIZE] = {0.274069, 0.451863, 0.274069};

const enum targetArray {
    CURRENT_Y = 0,
    GRADDIR = 1,
    GRADMAG = 2,
    GRADX = 3,
    GRADY = 4
} ta;

/* FUNCTIONS */
void readImage(void);
void writeImage(char *name);

/* Canny Algorithm */
void gaussianBlur(void);
void sobel(void);
void nms(void);
void thresholding(int low, int high, int weak);
void hysteresis(int weak);

void convolution(void);
void convolutionHorizontal1D(void);
void convolutionVertical1D(void);

void thresholdCheck(int *channel, int low, int high, int weak, int strong);

void copyToL1(int start_row, int delta_rows, int start_col, int delta_cols, int skip_lines, enum targetArray target);
void copyFromL1(int start_row, int delta_rows, int start_col, int delta_cols, enum targetArray target);

void copyToL2(int start, int delta, enum targetArray target);
void copyFromL2(int start, int delta, enum targetArray target);

void copyL2ToL1(int start_row, int start_col);
void rotateL2(int row);


int main()
{
    int weak = 50;

    readImage();

    /* 1. GAUSSIAN BLUR */
    gaussianBlur();
    // writeImage("BlurredImage.yuv"); // and save it

    /* 2. SOBEL MASK */
    sobel();
    // writeImage("GradMagImage.yuv"); // and save it

    /* 3. NON-MAXIMUM SUPPRESSION */
    nms();
    // writeImage("NMSImage.yuv"); // and save it

    /* 4. HYSTERESIS THRESHOLDING */
    thresholding(5, 20, weak);
    // writeImage("ThreshImage.yuv"); // and save it

    hysteresis(weak);
    writeImage("FinalImage.yuv"); // and save it

    return 0;
}

void readImage(void)
{
    FILE *frame_c;
    if ((frame_c = fopen(filename, "rb")) == NULL)
    {
        exit(-1);
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            current_y[i][j] = fgetc(frame_c);
            current_y[i][j + 1] = fgetc(frame_c);
            current_y[i][j + 2] = fgetc(frame_c);
            current_y[i][j + 3] = fgetc(frame_c);
        }
    }

    fclose(frame_c);
}

void writeImage(char *name)
{
    FILE *frame_yuv;
    frame_yuv = fopen(name, "wb");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            fputc(current_y[i][j], frame_yuv);
            fputc(current_y[i][j + 1], frame_yuv);
            fputc(current_y[i][j + 2], frame_yuv);
            fputc(current_y[i][j + 3], frame_yuv);
        }
    }

    fclose(frame_yuv);
}

void gaussianBlur(void)
{
    convolutionHorizontal1D();
    convolutionVertical1D();
}

void sobel(void)
{
    int cache_i = 0;
    double maxG = 0.0f; // Maximum gradient magnitude value, used for normalization
    double normCoeff = 0.0f;
    int *ptrCurrentY = &current_y[0][0];
    enum targetArray ta = CURRENT_Y;
    enum targetArray ta_gradx = GRADX;
    enum targetArray ta_grady = GRADY;

    convolution();
    
    // Gradient Magnitude
    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyFromL2(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            copyToL1(i, 1, j, 4, 0, ta_gradx);
            copyToL1(i, 1, j, 4, 1, ta_grady);
        
            cache[cache_i][j + 0] = sqrt(blk[0][0] * blk[0][0] + blk[1][0] * blk[1][0]);
            cache[cache_i][j + 1] = sqrt(blk[0][1] * blk[0][1] + blk[1][1] * blk[1][1]);
            cache[cache_i][j + 2] = sqrt(blk[0][2] * blk[0][2] + blk[1][2] * blk[1][2]);
            cache[cache_i][j + 3] = sqrt(blk[0][3] * blk[0][3] + blk[1][3] * blk[1][3]);
    
        }
    }
    copyFromL2(N - CACHE_ROWS, CACHE_ROWS, ta);

    for (i = 0; i < N * M; i += 4)
    {
        if (*(ptrCurrentY + i) > maxG)
        {
            maxG = *(ptrCurrentY + i);
        }
        if (*(ptrCurrentY + i + 1) > maxG)
        {
            maxG = *(ptrCurrentY + i + 1);
        }
        if (*(ptrCurrentY + i + 2) > maxG)
        {
            maxG = *(ptrCurrentY + i + 2);
        }
        if (*(ptrCurrentY + i + 3) > maxG)
        {
            maxG = *(ptrCurrentY + i + 3);
        }
    }

    // Normalize gradient magnitude (0 - 255)
    normCoeff = 255 / maxG;
    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyFromL2(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        if (cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyToL2(i, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            cache[cache_i][j + 0] = (cache[cache_i][j + 0] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 0] * normCoeff);
            cache[cache_i][j + 1] = (cache[cache_i][j + 1] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 1] * normCoeff);
            cache[cache_i][j + 2] = (cache[cache_i][j + 2] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 2] * normCoeff);
            cache[cache_i][j + 3] = (cache[cache_i][j + 3] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 3] * normCoeff);
        }
    }
    copyFromL2(N - CACHE_ROWS, CACHE_ROWS, ta);

    ta = GRADDIR;
    // Gradient Direction
    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyFromL2(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            copyToL1(i, 1, j, 4, 0, ta_gradx);
            copyToL1(i, 1, j, 4, 1, ta_grady);

            cache[cache_i][j + 0] = atan2(blk[1][0], blk[0][0]) * 180 / _PI;
            cache[cache_i][j + 1] = atan2(blk[1][1], blk[0][1]) * 180 / _PI;
            cache[cache_i][j + 2] = atan2(blk[1][2], blk[0][2]) * 180 / _PI;
            cache[cache_i][j + 3] = atan2(blk[1][3], blk[0][3]) * 180 / _PI;
        }
    }
    copyFromL2(N - CACHE_ROWS, CACHE_ROWS, ta);
}

void nms(void)
{
    int dir, beforePixel, afterPixel;
    int cache_i = 1;
    int PI = 180;
    enum targetArray ta = CURRENT_Y;

    // Ignore the border pixels
    for (i = 1; i < N - 1; i++)
    {
        copyToL2(i - 1, CACHE_ROWS - 1, ta);
        for (j = 1; j < M - 1; j++)
        {
            dir = gradDir[i][j];

            if (((0 <= dir) || (dir < PI / 8)) || ((15 * PI / 8 <= dir) || (dir <= 2 * PI)))
            {
                beforePixel = cache[cache_i][j - 1];
                afterPixel = cache[cache_i][j + 1];
            }
            else if (((PI / 8 <= dir) || (dir < 3 * PI / 8)) || ((9 * PI / 8 <= dir) || (dir <= 11 * PI / 8)))
            {
                beforePixel = cache[cache_i + 1][j - 1];
                afterPixel = cache[cache_i - 1][j + 1];
            }
            else if (((3 * PI / 8 <= dir) || (dir < 5 * PI / 8)) || ((11 * PI / 8 <= dir) || (dir <= 13 * PI / 8)))
            {
                beforePixel = cache[cache_i - 1][j];
                afterPixel = cache[cache_i + 1][j];
            }
            else
            {
                beforePixel = cache[cache_i - 1][j - 1];
                afterPixel = cache[cache_i + 1][j + 1];
            }

            if (!(cache[cache_i][j] >= beforePixel && cache[cache_i][j] >= afterPixel))
            {
                current_y[i][j] = 0;
            }
        }
    }
}

void thresholding(int low, int high, int weak)
{
    int strong = 255;
    int cache_i = 0;

    enum targetArray ta = CURRENT_Y;

    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyFromL2(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        if (cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyToL2(i, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            thresholdCheck(&cache[cache_i][j + 0], low, high, weak, strong);
            thresholdCheck(&cache[cache_i][j + 1], low, high, weak, strong);
            thresholdCheck(&cache[cache_i][j + 2], low, high, weak, strong);
            thresholdCheck(&cache[cache_i][j + 3], low, high, weak, strong);
        }
    }
    copyFromL2(N - CACHE_ROWS, CACHE_ROWS, ta);
}

void hysteresis(int weak)
{
    int topToBottom[N][M];
    int bottomToTop[N][M];
    int rightToLeft[N][M];
    int leftToRight[N][M];
    int *ptrT2B = &topToBottom[0][0];
    int *ptrB2T = &bottomToTop[0][0];
    int *ptrL2R = &leftToRight[0][0];
    int *ptrR2L = &rightToLeft[0][0];
    int *ptrCurrentY = &current_y[0][0];
    int fp = -1;

    for (i = 0; i < N * M; i += 4)
    {
        *(ptrT2B + i) = *(ptrCurrentY + i);
        *(ptrB2T + i) = *(ptrCurrentY + i);
        *(ptrL2R + i) = *(ptrCurrentY + i);
        *(ptrR2L + i) = *(ptrCurrentY + i);

        *(ptrT2B + i + 1) = *(ptrCurrentY + i + 1);
        *(ptrB2T + i + 1) = *(ptrCurrentY + i + 1);
        *(ptrL2R + i + 1) = *(ptrCurrentY + i + 1);
        *(ptrR2L + i + 1) = *(ptrCurrentY + i + 1);

        *(ptrT2B + i + 2) = *(ptrCurrentY + i + 2);
        *(ptrB2T + i + 2) = *(ptrCurrentY + i + 2);
        *(ptrL2R + i + 2) = *(ptrCurrentY + i + 2);
        *(ptrR2L + i + 2) = *(ptrCurrentY + i + 2);

        *(ptrT2B + i + 3) = *(ptrCurrentY + i + 3);
        *(ptrB2T + i + 3) = *(ptrCurrentY + i + 3);
        *(ptrL2R + i + 3) = *(ptrCurrentY + i + 3);
        *(ptrR2L + i + 3) = *(ptrCurrentY + i + 3);
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (topToBottom[i][j] == weak)
            {
                if (topToBottom[i][j + 1] == 255 || topToBottom[i][j - 1] == 255 ||
                    topToBottom[i - 1][j] == 255 || topToBottom[i + 1][j] == 255 ||
                    topToBottom[i - 1][j - 1] == 255 || topToBottom[i + 1][j - 1] == 255 ||
                    topToBottom[i - 1][j + 1] == 255 || topToBottom[i + 1][j + 1] == 255)
                {
                    topToBottom[i][j] = 255;
                }
                else
                {
                    topToBottom[i][j] = 0;
                }
            }
        }
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (bottomToTop[i][j] == weak)
            {
                if (bottomToTop[i][j + 1] == 255 || bottomToTop[i][j - 1] == 255 ||
                    bottomToTop[i - 1][j] == 255 || bottomToTop[i + 1][j] == 255 ||
                    bottomToTop[i - 1][j - 1] == 255 || bottomToTop[i + 1][j - 1] == 255 ||
                    bottomToTop[i - 1][j + 1] == 255 || bottomToTop[i + 1][j + 1] == 255)
                {

                    bottomToTop[i][j] = 255;
                }
                else
                {
                    bottomToTop[i][j] = 0;
                }
            }
        }
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (rightToLeft[i][j] == weak)
            {
                if (rightToLeft[i][j + 1] == 255 || rightToLeft[i][j - 1] == 255 ||
                    rightToLeft[i - 1][j] == 255 || rightToLeft[i + 1][j] == 255 ||
                    rightToLeft[i - 1][j - 1] == 255 || rightToLeft[i + 1][j - 1] == 255 ||
                    rightToLeft[i - 1][j + 1] == 255 || rightToLeft[i + 1][j + 1] == 255)
                {

                    rightToLeft[i][j] = 255;
                }
                else
                {

                    rightToLeft[i][j] = 0;
                }
            }
        }
    }

    for (i = 1; i < N - 1; i++)
    {
        for (j = 1; j < M - 1; j++)
        {
            if (leftToRight[i][j] == weak)
            {
                if (leftToRight[i][j + 1] == 255 || leftToRight[i][j - 1] == 255 ||
                    leftToRight[i - 1][j] == 255 || leftToRight[i + 1][j] == 255 ||
                    leftToRight[i - 1][j - 1] == 255 || leftToRight[i + 1][j - 1] == 255 ||
                    leftToRight[i - 1][j + 1] == 255 || leftToRight[i + 1][j + 1] == 255)
                {

                    leftToRight[i][j] = 255;
                }
                else
                {

                    leftToRight[i][j] = 0;
                }
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            fp = topToBottom[i][j] + bottomToTop[i][j] + rightToLeft[i][j] + leftToRight[i][j];
            if (fp > 255)
            {
                fp = 255;
            }

            current_y[i][j] = fp;
        }
    }
}

void convolution(void)
{
    int ki, kj;
    int kernelCenter = KERNEL_SIZE / 2;
    int k, l;
    enum targetArray ta = CURRENT_Y;

    copyToL2(0, KERNEL_SIZE, ta);

    //  Convolution X
    for (i = 0; i < N; i++)
    {
        if (i != 0) {
            rotateL2(i);
        }

        for (j = 0; j < M; j += 4)
        {

            copyL2ToL1(startCacheIdx, j);
            // Fill sum array with zeros
            for (k = 0; k < 2; k++){
                for (l = 0; l < 4; l++){
                    conv_sum[k][l] = 0;
                }
            }
            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
                {
                    if (i + ki < 0 && i + ki >= N){
                        continue;
                    }
                    if (j + 0 + kj > 0 && j + 0 + kj < M)
                    {
                        conv_sum[0][0] += blk[1 + ki][1 + 0 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                        conv_sum[1][0] += blk[1 + ki][1 + 0 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 1 + kj > 0 && j + 1 + kj < M)
                    {
                        conv_sum[0][1] += blk[1 + ki][1 + 1 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                        conv_sum[1][1] += blk[1 + ki][1 + 1 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 2 + kj > 0 && j + 2 + kj < M)
                    {
                        conv_sum[0][2] += blk[1 + ki][1 + 2 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                        conv_sum[1][2] += blk[1 + ki][1 + 2 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 3 + kj > 0 && j + 3 + kj < M)
                    {
                        conv_sum[0][3] += blk[1 + ki][1 + 3 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                        conv_sum[1][3] += blk[1 + ki][1 + 3 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            gradX[i][j + 0] = conv_sum[0][0];
            gradX[i][j + 1] = conv_sum[0][1];
            gradX[i][j + 2] = conv_sum[0][2];
            gradX[i][j + 3] = conv_sum[0][3];            

            gradY[i][j + 0] = conv_sum[1][0];
            gradY[i][j + 1] = conv_sum[1][1];
            gradY[i][j + 2] = conv_sum[1][2];
            gradY[i][j + 3] = conv_sum[1][3];
        }
    }
}

void convolutionHorizontal1D(void)
{
    int ki;
    int kCenter = KERNEL_SIZE / 2;
    int output[N][M];
    int *ptrCurrentY = &current_y[0][0];
    int *ptrOutput = &output[0][0];

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            conv_sum[0][0] = 0;
            // Apply the kernel horizontally
            for (ki = -kCenter; ki <= kCenter; ki++)
            {
                if (j + ki >= 0 && j + ki < N)
                {
                    conv_sum[0][0] += current_y[i][j + ki] * gauss_kernel_1d[ki + kCenter];
                }
            }
            output[i][j] = conv_sum[0][0];
        }
    }

    for (i = 0; i < N * M; i += 4)
    {
        *(ptrCurrentY + i) = *(ptrOutput + i);
        *(ptrCurrentY + i + 1) = *(ptrOutput + i + 1);
        *(ptrCurrentY + i + 2) = *(ptrOutput + i + 2);
        *(ptrCurrentY + i + 3) = *(ptrOutput + i + 3);
    }
}

void convolutionVertical1D(void)
{
    int ki;
    int kCenter = KERNEL_SIZE / 2;
    int output[N][M];
    int *ptrCurrentY = &current_y[0][0];
    int *ptrOutput = &output[0][0];

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            conv_sum[1][0] = 0;
            // Apply the kernel horizontally
            for (ki = -kCenter; ki <= kCenter; ki++)
            {
                if (i + ki >= 0 && i + ki < M)
                {
                    conv_sum[1][0] += current_y[i + ki][j] * gauss_kernel_1d[ki + kCenter];
                }
            }
            output[i][j] = conv_sum[1][0];
        }
    }

    for (i = 0; i < N * M; i += 4)
    {
        *(ptrCurrentY + i) = *(ptrOutput + i);
        *(ptrCurrentY + i + 1) = *(ptrOutput + i + 1);
        *(ptrCurrentY + i + 2) = *(ptrOutput + i + 2);
        *(ptrCurrentY + i + 3) = *(ptrOutput + i + 3);
    }
}

void thresholdCheck(int *channel, int low, int high, int weak, int strong)
{
    if (*channel >= high)
    {
        *channel = strong;
    }
    else if (*channel >= low)
    {
        *channel = weak;
    }
    else
    {
        *channel = 0;
    }
}

void copyToL2(int start, int delta, enum targetArray target)
{
    int ci, cj;
    int(*ptrArray)[N][M] = &current_y;
    int end = start + delta;
    switch (target)
    {
    case 0:
        ptrArray = &current_y;
        break;
    case 1:
        ptrArray = &gradDir;
        break;
    case 2:
        ptrArray = &gradMag;
        break;
    default:
        ptrArray = &current_y;
    }
    for (ci = start; ci < end; ci++)
    {
        for (cj = 0; cj < M; cj++)
        {
            cache[ci - start][cj] = (*ptrArray)[ci][cj];
        }
    }
}

void copyFromL2(int start, int delta, enum targetArray target)
{
    int ci, cj;
    int(*ptrArray)[N][M] = &current_y;
    int end = start + delta;
    switch (target)
    {
    case 0:
        ptrArray = &current_y;
        break;
    case 1:
        ptrArray = &gradDir;
        break;
    case 2:
        ptrArray = &gradMag;
        break;
    default:
        ptrArray = &current_y;
    }
    for (ci = start; ci < end; ci++)
    {
        for (cj = 0; cj < M; cj++)
        {
            (*ptrArray)[ci][cj] = cache[ci - start][cj];
        }
    }
}


void copyFromL1(int start_row, int delta_rows, int start_col, int delta_cols, enum targetArray target)
{
    int ci, cj;
    int(*ptrArray)[N][M] = &current_y;
    int end_row = start_row + delta_rows;
    int end_col = start_col + delta_cols;
    switch (target)
    {
    case 0:
        ptrArray = &current_y;
        break;
    case 1:
        ptrArray = &gradDir;
        break;
    case 2:
        ptrArray = &gradMag;
        break;
    case 3:
        ptrArray = &gradX;
        break;
    case 4:
        ptrArray = &gradY;
        break;

    default:
        ptrArray = &current_y;
    }
    for (ci = start_row; ci < end_row; ci++)
    {
        if (blk[ci - start_row][0] == -1)
        {
            continue;
        }
        for (cj = start_col; cj < end_col; cj++)
        {
            (*ptrArray)[ci][cj] = blk[ci - start_row][cj - start_col];
        }
    }
}

void copyToL1(int start_row, int delta_rows, int start_col, int delta_cols, int skip_lines, enum targetArray target)
{
    int ci, cj;
    int (*ptrArray)[N][M] = &current_y;
    int end_row = start_row + delta_rows;
    int end_col = start_col + delta_cols;
    switch (target)
    {
    case 0:
        ptrArray = &current_y;
        break;
    case 1:
        ptrArray = &gradDir;
        break;
    case 2:
        ptrArray = &gradMag;
        break;
    case 3:
        ptrArray = &gradX;
        break;
    case 4:
        ptrArray = &gradY;
        break;
    default:
        ptrArray = &current_y;
    }

    for (ci = start_row; ci < end_row; ci++)
    {
        for (cj = start_col; cj < end_col; cj++)
        {
            blk[ci - start_row + skip_lines][cj - start_col] = (*ptrArray)[ci][cj];
        }
    }
}

void rotateL2(int row){
    int ri;
    
    // insert
    for (ri = 0; ri < M; ri++)
    {
        cache[startCacheIdx][ri] = current_y[row][ri];
    }
    // move pointer
    startCacheIdx = (startCacheIdx + 1) % KERNEL_SIZE;
}

void copyL2ToL1(int start_row, int start_col)
{
    int ci, cj;
    int cache_i;
    int end_col = start_col + 2*KERNEL_SIZE;

    for (ci = 0; ci < 2*KERNEL_SIZE; ci++)
    {
        cache_i = (ci + start_row) % KERNEL_SIZE;
        for (cj = start_col; cj < end_col; cj++)
        {
            blk[ci][cj - start_col] = cache[cache_i][cj];
        }
    }
}