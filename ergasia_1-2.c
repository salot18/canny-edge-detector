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

/* code for armulator*/
#pragma arm section zidata = "ram"
int current_y[N][M];
#pragma arm section

double gaussian_kernel[KERNEL_SIZE];
int gradX[N][M];
int gradY[N][M];
int gradDir[N][M];

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

int i, j, ii, jj;

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

int main()
{
    int weak = 50;

    readImage();

    /* 1. GAUSSIAN BLUR */
    gaussianBlur(7, 1);
    // writeImage("BlurredImage.yuv"); // and save it

    /* 2. SOBEL MASK */
    sobel();
    // writeImage("GradMagImage.yuv"); // and save it

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
    double maxG = 0.0f; // Maximum gradient magnitude value, used for normalization
    double normCoeff = 0.0f;

    convolution();

    // Gradient Magnitude
    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {
                    current_y[ii][jj + 0] = sqrt(gradX[ii][jj + 0] * gradX[ii][jj + 0] + gradY[ii][jj + 0] * gradY[ii][jj + 0]);
                    current_y[ii][jj + 1] = sqrt(gradX[ii][jj + 1] * gradX[ii][jj + 1] + gradY[ii][jj + 1] * gradY[ii][jj + 1]);
                    current_y[ii][jj + 2] = sqrt(gradX[ii][jj + 2] * gradX[ii][jj + 2] + gradY[ii][jj + 2] * gradY[ii][jj + 2]);
                    current_y[ii][jj + 3] = sqrt(gradX[ii][jj + 3] * gradX[ii][jj + 3] + gradY[ii][jj + 3] * gradY[ii][jj + 3]);
                }
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            if (current_y[i][j] > maxG)
            {
                maxG = current_y[i][j];
            }
        }
    }

    // Normalize gradient magnitude (0 - 255)
    normCoeff = 255 / maxG;
    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {
                    current_y[ii][jj + 0] = (current_y[ii][jj + 0] * normCoeff > 255) ? 255 : (int)(current_y[ii][jj + 0] * normCoeff);
                    current_y[ii][jj + 1] = (current_y[ii][jj + 1] * normCoeff > 255) ? 255 : (int)(current_y[ii][jj + 1] * normCoeff);
                    current_y[ii][jj + 2] = (current_y[ii][jj + 2] * normCoeff > 255) ? 255 : (int)(current_y[ii][jj + 2] * normCoeff);
                    current_y[ii][jj + 3] = (current_y[ii][jj + 3] * normCoeff > 255) ? 255 : (int)(current_y[ii][jj + 3] * normCoeff);
                }
            }
        }
    }

    // Gradient Direction
    for (i = 0; i < N; i += TILE_SIZE)
    {
        for (j = 0; j < M; j += TILE_SIZE)
        {
            for (ii = i; ii < i + TILE_SIZE; ii++)
            {
                for (jj = j; jj < j + TILE_SIZE; jj += 4)
                {

                    gradDir[ii][jj + 0] = atan2(gradY[ii][jj + 0], gradX[ii][jj + 0]) * 180 / _PI;
                    gradDir[ii][jj + 1] = atan2(gradY[ii][jj + 1], gradX[ii][jj + 1]) * 180 / _PI;
                    gradDir[ii][jj + 2] = atan2(gradY[ii][jj + 2], gradX[ii][jj + 2]) * 180 / _PI;
                    gradDir[ii][jj + 3] = atan2(gradY[ii][jj + 3], gradX[ii][jj + 3]) * 180 / _PI;
                }
            }
        }
    }
}

void nms(void)
{
    int dir, beforePixel, afterPixel;
    int PI = 180;
    int gradMag[N][M];

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            gradMag[i][j] = current_y[i][j];
        }
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
    int fp = -1;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            topToBottom[i][j] = current_y[i][j];
            bottomToTop[i][j] = current_y[i][j];
            rightToLeft[i][j] = current_y[i][j];
            leftToRight[i][j] = current_y[i][j];
        }
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
    // double *kernel = malloc(sizeof(double *) * size);
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
    for (i = 0; i < KERNEL_SIZE; i++)
    {
        gaussian_kernel[i] = gaussian_kernel[i] / sum;
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

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            current_y[i][j] = output[i][j];
        }
    }
}

void convolutionVertical1D(void)
{
    int ki;
    int kCenter = KERNEL_SIZE / 2;
    int sum = 0;
    int output[N][M];

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

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            current_y[i][j] = output[i][j];
        }
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