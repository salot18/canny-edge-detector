#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// SUNFLOWER
#define N 200 // height
#define M 200 // width
#define filename "sunflower_200x200_444.yuv"

#define _PI 3.14159265359
#define KERNEL_SIZE 7
#define KERNEL_SOBEL 3
#define TILE_SIZE 100
#define CACHE_ROWS 4

/* code for armulator*/
#pragma arm section zidata = "cache"
int cache[CACHE_ROWS][M];
#pragma arm section

#pragma arm section zidata = "ram"
int current_y[N][M];
int i, j, ii, jj;
double gaussian_kernel[KERNEL_SIZE];
int gradX[N][M];
int gradY[N][M];
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

/* FUNCTIONS */
void readImage(void);
void writeImage(char *name);

/* Canny Algorithm */
void gaussianBlur(int kernelSize, int kernelSigma);
void sobel(void);
void nms(void);
void thresholding(int low, int high, int weak);
void hysteresis(int weak);

void gaussianKernel1D(int sigma);
void convolution(void);
void convolutionHorizontal1D(void);
void convolutionVertical1D(void);

void thresholdCheck(int *channel, int low, int high, int weak, int strong);

void copyToCache(int start, int delta);
void copyFromCache(int start, int delta);
void fillCacheWithPadding(int start, int delta);

int main()
{
    int weak = 50;

    readImage();

    /* 1. GAUSSIAN BLUR */
    gaussianBlur(7, 1);
    // writeImage("BlurredImage.yuv"); // and save it

    /* 2. SOBEL MASK */
    sobel();
    writeImage("GradMagImage.yuv"); // and save it

    // /* 3. NON-MAXIMUM SUPPRESSION */
    nms();
    // writeImage("NMSImage.yuv"); // and save it

    // /* 4. HYSTERESIS THRESHOLDING */
    thresholding(5, 20, weak); // (image, low, hight, weak)
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
        printf("current frame doesn't exist\n");
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

void gaussianBlur(int kernelSize, int kernelSigma)
{
    gaussianKernel1D(kernelSigma);

    convolutionHorizontal1D();
    convolutionVertical1D();
}

void sobel(void)
{
    int cache_i = 0;
    int tempi, tempj;
    double maxG = 0.0f; // Maximum gradient magnitude value, used for normalization
    double normCoeff = 0.0f;
    int *ptrCurrentY = &current_y[0][0];

    convolution();
    fillCacheWithPadding(0, CACHE_ROWS);

    // Gradient Magnitude
    for (i = 0; i < N; i++)
    {
        cache_i = i % CACHE_ROWS;
        if (i != 0 && cache_i == 0 && i <= N - CACHE_ROWS)
        {
            printf("Wrote lines from cache %d-> %d\n", i - CACHE_ROWS, i - 1);
            copyFromCache(i - CACHE_ROWS, CACHE_ROWS);
        }
        for (j = 0; j < M; j += 4)
        {
            cache[cache_i][j + 0] = sqrt(gradX[i][j + 0] * gradX[i][j + 0] + gradY[i][j + 0] * gradY[i][j + 0]);
            cache[cache_i][j + 1] = sqrt(gradX[i][j + 1] * gradX[i][j + 1] + gradY[i][j + 1] * gradY[i][j + 1]);
            cache[cache_i][j + 2] = sqrt(gradX[i][j + 2] * gradX[i][j + 2] + gradY[i][j + 2] * gradY[i][j + 2]);
            cache[cache_i][j + 3] = sqrt(gradX[i][j + 3] * gradX[i][j + 3] + gradY[i][j + 3] * gradY[i][j + 3]);
        }
    }
    printf("Writing remaining lines %d-> %d\n", N - CACHE_ROWS, N - 1);
    copyFromCache(N - CACHE_ROWS, CACHE_ROWS);
    fillCacheWithPadding(0, CACHE_ROWS);

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
        for (j = 0; j < M; j += 4)
        {
            current_y[i][j + 0] = (current_y[i][j + 0] * normCoeff > 255) ? 255 : (int)(current_y[i][j + 0] * normCoeff);
            current_y[i][j + 1] = (current_y[i][j + 1] * normCoeff > 255) ? 255 : (int)(current_y[i][j + 1] * normCoeff);
            current_y[i][j + 2] = (current_y[i][j + 2] * normCoeff > 255) ? 255 : (int)(current_y[i][j + 2] * normCoeff);
            current_y[i][j + 3] = (current_y[i][j + 3] * normCoeff > 255) ? 255 : (int)(current_y[i][j + 3] * normCoeff);
        }
    }

    // Gradient Direction
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {

            gradDir[i][j + 0] = atan2(gradY[i][j + 0], gradX[i][j + 0]) * 180 / _PI;
            gradDir[i][j + 1] = atan2(gradY[i][j + 1], gradX[i][j + 1]) * 180 / _PI;
            gradDir[i][j + 2] = atan2(gradY[i][j + 2], gradX[i][j + 2]) * 180 / _PI;
            gradDir[i][j + 3] = atan2(gradY[i][j + 3], gradX[i][j + 3]) * 180 / _PI;
        }
    }
}

void nms(void)
{
    int dir, beforePixel, afterPixel;
    int PI = 180;
    int gradMag[N][M];
    int *ptrCurrentY = &current_y[0][0];
    int *ptrGradMag = &gradMag[0][0];

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
        for (j = 1; j < M - 1; j++)
        {
            dir = gradDir[i][j];

            // ################## REVIEW THE CONDITIONS #########################

            if (((0 <= dir) || (dir < PI / 8)) || ((15 * PI / 8 <= dir) || (dir <= 2 * PI)))
            {
                beforePixel = gradMag[i][j - 1];
                afterPixel = gradMag[i][j + 1];
            }
            else if (((PI / 8 <= dir) || (dir < 3 * PI / 8)) || ((9 * PI / 8 <= dir) || (dir <= 11 * PI / 8)))
            {
                beforePixel = gradMag[i + 1][j - 1];
                afterPixel = gradMag[i - 1][j + 1];
            }
            else if (((3 * PI / 8 <= dir) || (dir < 5 * PI / 8)) || ((11 * PI / 8 <= dir) || (dir <= 13 * PI / 8)))
            {
                beforePixel = gradMag[i - 1][j];
                afterPixel = gradMag[i + 1][j];
            }
            else
            {
                beforePixel = gradMag[i - 1][j - 1];
                afterPixel = gradMag[i + 1][j + 1];
            }

            if (!(gradMag[i][j] >= beforePixel && gradMag[i][j] >= afterPixel))
            {
                current_y[i][j] = 0;
            }
        }
    }
}

void thresholding(int low, int high, int weak)
{
    int strong = 255;

    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {
                    thresholdCheck(&current_y[ii][jj + 0], low, high, weak, strong);
                    thresholdCheck(&current_y[ii][jj + 1], low, high, weak, strong);
                    thresholdCheck(&current_y[ii][jj + 2], low, high, weak, strong);
                    thresholdCheck(&current_y[ii][jj + 3], low, high, weak, strong);
                }
            }
        }
    }
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

void gaussianKernel1D(int sigma)
{
    int x;
    int center = KERNEL_SIZE / 2;
    double res = 0.0f;
    double sum = 0.0f;

    for (i = 0; i < KERNEL_SIZE; i++) // checking purpose
    {
        x = i - center;
        res = exp(-((double)(x * x) / (2 * sigma * sigma)));
        sum += res;
        gaussian_kernel[i] = res;
    }

    // Normalize the kernel
    for (i = 0; i < KERNEL_SIZE; i += 4)
    {
        gaussian_kernel[i] = gaussian_kernel[i] / sum;
        gaussian_kernel[i + 1] = gaussian_kernel[i + 1] / sum;
        gaussian_kernel[i + 2] = gaussian_kernel[i + 2] / sum;
        gaussian_kernel[i + 3] = gaussian_kernel[i + 3] / sum;
    }
}

void convolution(void)
{
    int ki, kj;
    int kernelCenter = KERNEL_SOBEL / 2;
    int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;

    //  Convolution X
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            sum0 = sum1 = sum2 = sum3 = 0;
            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
                {
                    if (i + ki >= 0 && i + ki < N && j + kj >= 0 && j + kj < M)
                    {
                        sum0 += current_y[i + ki][j + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 1 + kj >= 0 && j + 1 + kj < M)
                    {
                        sum1 += current_y[i + ki][j + 1 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 2 + kj >= 0 && j + 2 + kj < M)
                    {
                        sum2 += current_y[i + ki][j + 2 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 3 + kj >= 0 && j + 3 + kj < M)
                    {
                        sum3 += current_y[i + ki][j + 3 + kj] * sobel_kernel_x[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            gradX[i][j] = sum0;
            gradX[i][j + 1] = sum1;
            gradX[i][j + 2] = sum2;
            gradX[i][j + 3] = sum3;
        }
    }

    //  Convolution Y
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j += 4)
        {
            sum0 = sum1 = sum2 = sum3 = 0;
            for (ki = -kernelCenter; ki <= kernelCenter; ki++)
            {
                for (kj = -kernelCenter; kj <= kernelCenter; kj++)
                {
                    if (i + ki >= 0 && i + ki < N && j + kj >= 0 && j + kj < M)
                    {
                        sum0 += current_y[i + ki][j + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 1 + kj >= 0 && j + 1 + kj < M)
                    {
                        sum1 += current_y[i + ki][j + 1 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 2 + kj >= 0 && j + 2 + kj < M)
                    {
                        sum2 += current_y[i + ki][j + 2 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                    if (i + ki >= 0 && i + ki < N && j + 3 + kj >= 0 && j + 3 + kj < M)
                    {
                        sum3 += current_y[i + ki][j + 3 + kj] * sobel_kernel_y[ki + kernelCenter][kj + kernelCenter];
                    }
                }
            }

            gradY[i][j] = sum0;
            gradY[i][j + 1] = sum1;
            gradY[i][j + 2] = sum2;
            gradY[i][j + 3] = sum3;
        }
    }
}

void convolutionHorizontal1D(void)
{
    int ki;
    int kCenter = KERNEL_SIZE / 2;
    int sum = 0;
    int output[N][M];
    int *ptrCurrentY = &current_y[0][0];
    int *ptrOutput = &output[0][0];

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            sum = 0;
            // Apply the kernel horizontally
            for (ki = -kCenter; ki <= kCenter; ki++)
            {
                if (j + ki >= 0 && j + ki < N)
                { // Boundary checki
                    sum += current_y[i][j + ki] * gaussian_kernel[ki + kCenter];
                }
            }
            output[i][j] = sum;
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
    int sum = 0;
    int output[N][M];
    int *ptrCurrentY = &current_y[0][0];
    int *ptrOutput = &output[0][0];

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            sum = 0;
            // Apply the kernel horizontally
            for (ki = -kCenter; ki <= kCenter; ki++)
            {
                if (i + ki >= 0 && i + ki < M)
                {
                    sum += current_y[i + ki][j] * gaussian_kernel[ki + kCenter];
                }
            }
            output[i][j] = sum;
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
*/
void copyToCache(int start, int delta)
{
    int ci, cj;
    int end = start + delta;
    for (ci = start; ci < end; ci++)
    {
        for (cj = 0; cj < M; cj++)
        {
            cache[ci - start][cj] = current_y[ci][cj];
        }
    }
}

/*
    Loads cache rows into current_y
    start: starting index row of current_y
    delta: number of rows
*/
void copyFromCache(int start, int delta)
{
    int ci, cj;
    int end = start + delta;
    for (ci = start; ci < end; ci++)
    {
        if (cache[ci - start][0] == -1)
        {
            printf("skipping padding line %d\n", ci);
            continue;
        }
        for (cj = 0; cj < M; cj++)
        {
            current_y[ci][cj] = cache[ci - start][cj];
        }
    }
}

/*
    Flags concecutive cache rows as invalid
    start: starting index row inside cache
    delta: number of rows
*/
void fillCacheWithPadding(int start, int delta)
{
    int ci;
    int end = start + delta;
    for (ci = start; ci < end; ci++)
    {
        cache[ci][0] = -1;
    }
}
