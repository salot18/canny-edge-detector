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
#define KERNEL_SOBEL 3
#define TILE_SIZE 100
#define CACHE_ROWS 4

/* code for armulator*/
#pragma arm section zidata = "cache"
int cache[CACHE_ROWS][M];
int i, j, ii, jj;
int conv_sum0 = 0, conv_sum1 = 0, conv_sum2 = 0, conv_sum3 = 0;
#pragma arm section

#pragma arm section zidata = "ram"
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
    GRADMAG = 2
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

void copyToCache(int start, int delta, enum targetArray ta);
void copyFromCache(int start, int delta, enum targetArray ta);

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

    convolution();

    // Gradient Magnitude
    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyFromCache(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            cache[cache_i][j + 0] = sqrt(gradX[i][j + 0] * gradX[i][j + 0] + gradY[i][j + 0] * gradY[i][j + 0]);
            cache[cache_i][j + 1] = sqrt(gradX[i][j + 1] * gradX[i][j + 1] + gradY[i][j + 1] * gradY[i][j + 1]);
            cache[cache_i][j + 2] = sqrt(gradX[i][j + 2] * gradX[i][j + 2] + gradY[i][j + 2] * gradY[i][j + 2]);
            cache[cache_i][j + 3] = sqrt(gradX[i][j + 3] * gradX[i][j + 3] + gradY[i][j + 3] * gradY[i][j + 3]);
        }
    }
    copyFromCache(N - CACHE_ROWS, CACHE_ROWS, ta);

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
            copyFromCache(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        if (cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyToCache(i, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            cache[cache_i][j + 0] = (cache[cache_i][j + 0] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 0] * normCoeff);
            cache[cache_i][j + 1] = (cache[cache_i][j + 1] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 1] * normCoeff);
            cache[cache_i][j + 2] = (cache[cache_i][j + 2] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 2] * normCoeff);
            cache[cache_i][j + 3] = (cache[cache_i][j + 3] * normCoeff > 255) ? 255 : (int)(cache[cache_i][j + 3] * normCoeff);
        }
    }
    copyFromCache(N - CACHE_ROWS, CACHE_ROWS, ta);

    ta = GRADDIR;
    // Gradient Direction
    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyFromCache(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {

            cache[cache_i][j + 0] = atan2(gradY[i][j + 0], gradX[i][j + 0]) * 180 / _PI;
            cache[cache_i][j + 1] = atan2(gradY[i][j + 1], gradX[i][j + 1]) * 180 / _PI;
            cache[cache_i][j + 2] = atan2(gradY[i][j + 2], gradX[i][j + 2]) * 180 / _PI;
            cache[cache_i][j + 3] = atan2(gradY[i][j + 3], gradX[i][j + 3]) * 180 / _PI;
        }
    }
    copyFromCache(N - CACHE_ROWS, CACHE_ROWS, ta);
}

void nms(void)
{
    int dir, beforePixel, afterPixel;
    int cache_i = 1;
    int PI = 180;
    int *ptrCurrentY = &current_y[0][0];
    int *ptrGradMag = &gradMag[0][0];
    enum targetArray ta = GRADMAG;

    for (i = 0; i < N * M; i += 4)
    {
        *(ptrGradMag + i + 0) = *(ptrCurrentY + i + 0);
        *(ptrGradMag + i + 1) = *(ptrCurrentY + i + 1);
        *(ptrGradMag + i + 2) = *(ptrCurrentY + i + 2);
        *(ptrGradMag + i + 3) = *(ptrCurrentY + i + 3);
    }

    // Ignore the border pixels
    for (i = 1; i < N - 1; i++)
    {
        copyToCache(i - 1, CACHE_ROWS - 1, ta);
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
            copyFromCache(i - CACHE_ROWS, CACHE_ROWS, ta);
        }
        if (cache_i == 0 && i <= N - CACHE_ROWS)
        {
            copyToCache(i, CACHE_ROWS, ta);
        }
        for (j = 0; j < M; j += 4)
        {
            thresholdCheck(&cache[cache_i][j + 0], low, high, weak, strong);
            thresholdCheck(&cache[cache_i][j + 1], low, high, weak, strong);
            thresholdCheck(&cache[cache_i][j + 2], low, high, weak, strong);
            thresholdCheck(&cache[cache_i][j + 3], low, high, weak, strong);
        }
    }
    copyFromCache(N - CACHE_ROWS, CACHE_ROWS, ta);
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
    int kernelCenter = KERNEL_SOBEL / 2;

    int cache_i = 1;
    enum targetArray ta = CURRENT_Y;

    //  Convolution X
    for (i = 0; i < N; i++)
    {
        copyToCache(i, CACHE_ROWS - 1, ta);
        for (j = 0; j < M; j += 4)
        {
            conv_sum0 = conv_sum1 = conv_sum2 = conv_sum3 = 0;
            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
                {
                    if (j + kj >= 0 && j + kj < M)
                    {
                        conv_sum0 += cache[cache_i + ki][j + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 1 + kj >= 0 && j + 1 + kj < M)
                    {
                        conv_sum1 += cache[cache_i + ki][j + 1 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 2 + kj >= 0 && j + 2 + kj < M)
                    {
                        conv_sum2 += cache[cache_i + ki][j + 2 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 3 + kj >= 0 && j + 3 + kj < M)
                    {
                        conv_sum3 += cache[cache_i + ki][j + 3 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            gradX[i][j] = conv_sum0;
            gradX[i][j + 1] = conv_sum1;
            gradX[i][j + 2] = conv_sum2;
            gradX[i][j + 3] = conv_sum3;
        }
    }

    //  Convolution Y
    for (i = 0; i < N; i++)
    {
        copyToCache(i, CACHE_ROWS - 1, ta);

        for (j = 0; j < M; j += 4)
        {
            conv_sum0 = conv_sum1 = conv_sum2 = conv_sum3 = 0;

            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
                {
                    if (j + kj >= 0 && j + kj < M)
                    {
                        conv_sum0 += cache[cache_i + ki][j + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 1 + kj >= 0 && j + 1 + kj < M)
                    {
                        conv_sum1 += cache[cache_i + ki][j + 1 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 2 + kj >= 0 && j + 2 + kj < M)
                    {
                        conv_sum2 += cache[cache_i + ki][j + 2 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (j + 3 + kj >= 0 && j + 3 + kj < M)
                    {
                        conv_sum3 += cache[cache_i + ki][j + 3 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            gradY[i][j] = conv_sum0;
            gradY[i][j + 1] = conv_sum1;
            gradY[i][j + 2] = conv_sum2;
            gradY[i][j + 3] = conv_sum3;
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
            conv_sum0 = 0;
            // Apply the kernel horizontally
            for (ki = -kCenter; ki <= kCenter; ki++)
            {
                if (j + ki >= 0 && j + ki < N)
                {
                    conv_sum0 += current_y[i][j + ki] * gauss_kernel_1d[ki + kCenter];
                }
            }
            output[i][j] = conv_sum0;
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
            conv_sum0 = 0;
            // Apply the kernel horizontally
            for (ki = -kCenter; ki <= kCenter; ki++)
            {
                if (i + ki >= 0 && i + ki < M)
                {
                    conv_sum0 += current_y[i + ki][j] * gauss_kernel_1d[ki + kCenter];
                }
            }
            output[i][j] = conv_sum0;
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

/*
    Loads current_y rows into cache
    start: starting index row of current_y
    delta: number of rows
    target: 0->currnt_y, 1->gradDir
*/
void copyToCache(int start, int delta, enum targetArray ta)
{
    int(*ptrArray)[N][M];
    int ci, cj;
    int end = start + delta;
    switch (ta)
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
        exit(1);
    }
    for (ci = start; ci < end; ci++)
    {
        for (cj = 0; cj < M; cj++)
        {
            cache[ci - start][cj] = (*ptrArray)[ci][cj];
        }
    }
}

/*
    Loads cache rows into current_y
    start: starting index row of current_y
    delta: number of rows
*/
void copyFromCache(int start, int delta, enum targetArray ta)
{
    int(*ptrArray)[N][M];
    int ci, cj;
    int end = start + delta;
    switch (ta)
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
        exit(1);
    }
    for (ci = start; ci < end; ci++)
    {
        if (cache[ci - start][0] == -1)
        {
            continue;
        }
        for (cj = 0; cj < M; cj++)
        {
            (*ptrArray)[ci][cj] = cache[ci - start][cj];
        }
    }
}
