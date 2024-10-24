#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

// Define constants
#define FRAME_SIZE 400
#define FRAME_SHIFT 160
#define NUM_FILTERS 26
#define NUM_MFCC 13
#define PI 3.14159265358979323846

// Pre-emphasis
void pre_emphasis(double *signal, double *pre_emphasized_signal, int length) {
    double alpha = 0.97;
    pre_emphasized_signal[0] = signal[0];
    for (int i = 1; i < length; i++) {
        pre_emphasized_signal[i] = signal[i] - alpha * signal[i - 1];
    }
}

// Hamming window
void hamming_window(double *frame, int length) {
    for (int i = 0; i < length; i++) {
        frame[i] *= 0.54 - 0.46 * cos(2 * PI * i / (length - 1));
    }
}

// FFT
void fft(double *frame, complex double *fft_result, int length) {
    for (int k = 0; k < length; k++) {
        fft_result[k] = 0;
        for (int n = 0; n < length; n++) {
            fft_result[k] += frame[n] * cexp(-I * 2 * PI * k * n / length);
        }
    }
}

// Mel filter bank
void mel_filter_bank(double *power_spectrum, double *mel_energies, int length) {
    for (int i = 0; i < NUM_FILTERS; i++) {
        mel_energies[i] = 0;
        for (int j = 0; j < length; j++) {
            // triangular filter
            if (j >= i && j <= i + 2) {
                mel_energies[i] += power_spectrum[j];
            }
        }
    }
}

// DCT
void dct(double *log_mel_energies, double *mfcc) {
    for (int k = 0; k < NUM_MFCC; k++) {
        mfcc[k] = 0;
        for (int n = 0; n < NUM_FILTERS; n++) {
            mfcc[k] += log_mel_energies[n] * cos(PI * k * (2 * n + 1) / (2.0 * NUM_FILTERS));
        }
    }
}

// Main MFCC extraction
void extract_mfcc(double *signal, double **mfcc_result, int signal_length) {
    double pre_emphasized_signal[signal_length];
    pre_emphasis(signal, pre_emphasized_signal, signal_length);

    double frame[FRAME_SIZE];
    complex double fft_result[FRAME_SIZE];
    double power_spectrum[FRAME_SIZE];
    double mel_energies[NUM_FILTERS];
    double log_mel_energies[NUM_FILTERS];

    int num_frames = (signal_length - FRAME_SIZE) / FRAME_SHIFT + 1;

    for (int i = 0; i < num_frames; i++) {
        int start = i * FRAME_SHIFT;
        for (int j = 0; j < FRAME_SIZE; j++) {
            if (start + j < signal_length) {
                frame[j] = pre_emphasized_signal[start + j];
            } else {
                frame[j] = 0.0;
            }
        }

        hamming_window(frame, FRAME_SIZE);
        fft(frame, fft_result, FRAME_SIZE);

        for (int j = 0; j < FRAME_SIZE / 2 + 1; j++) {
            power_spectrum[j] = creal(fft_result[j]) * creal(fft_result[j]) + cimag(fft_result[j]) * cimag(fft_result[j]);
        }

        mel_filter_bank(power_spectrum, mel_energies, FRAME_SIZE / 2 + 1);

        for (int j = 0; j < NUM_FILTERS; j++) {
            log_mel_energies[j] = log(mel_energies[j] + 1e-10); // Adding epsilon to avoid log(0)
        }

        dct(log_mel_energies, mfcc_result[i]);

        // Print MFCCs for the current frame
        printf("Frame %d:", i);
        for (int k = 0; k < NUM_MFCC; k++) {
            printf("%.4f ", mfcc_result[i][k]);
        }
        printf("\n");
    }
}

// Function to calculate Euclidean distance between two MFCC vectors
double euclidean_distance(double *mfcc1, double *mfcc2, int num_coefficients) {
    double distance = 0.0;
    for (int i = 0; i < num_coefficients; ++i) { // Start from 0 to include all coefficients
        distance += pow(mfcc1[i] - mfcc2[i], 2);
    }
    return sqrt(distance);
}

// Function to read samples from file
void read_samples_from_file(const char *filename, double **signal, int *signal_length) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    int count = 0;
    int capacity = 1000;
    double *temp_signal = (double *)malloc(capacity * sizeof(double));

    if (!temp_signal) {
        perror("Memory allocation error");
        exit(EXIT_FAILURE);
    }

    while (fscanf(file, "%lf", &temp_signal[count]) == 1) {
        count++;
        if (count >= capacity) {
            capacity *= 2;
            temp_signal = (double *)realloc(temp_signal, capacity * sizeof(double));
            if (!temp_signal) {
                perror("Memory reallocation error");
                exit(EXIT_FAILURE);
            }
        }
    }

    fclose(file);

    *signal_length = count;
    *signal = (double *)malloc(count * sizeof(double));
    if (!*signal) {
        perror("Memory allocation error");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < count; i++) {
        (*signal)[i] = temp_signal[i];
    }

    free(temp_signal);
}

int main() {
    const char *filenames[] = {
        "C:/Users/KIIT/Downloads/OneDrive_2024-06-18/Vowel recordings/234101011_e_8-1.txt",
        "C:/Users/KIIT/Downloads/OneDrive_2024-06-18/Vowel recordings/234101011_e_8-1.txt"
    };

    const int num_files = sizeof(filenames) / sizeof(filenames[0]);

    // Variables for MFCC calculation
    double *signal;
    int signal_length;

    // Read the first file
    read_samples_from_file(filenames[0], &signal, &signal_length);

    // Calculate MFCCs for the first file
    printf("MFCC for file 1:\n");
    int num_frames1 = (signal_length - FRAME_SIZE) / FRAME_SHIFT + 1;
    double **mfcc_result1 = (double **)malloc(num_frames1 * sizeof(double *));
    for (int i = 0; i < num_frames1; i++) {
        mfcc_result1[i] = (double *)malloc(NUM_MFCC * sizeof(double));
    }
    extract_mfcc(signal, mfcc_result1, signal_length);
    free(signal);

    // Read the second file
    read_samples_from_file(filenames[1], &signal, &signal_length);

    // Calculate MFCCs for the second file
    printf("MFCC for file 2:\n");
    int num_frames2 = (signal_length - FRAME_SIZE) / FRAME_SHIFT + 1;
    double **mfcc_result2 = (double **)malloc(num_frames2 * sizeof(double *));
    for (int i = 0; i < num_frames2; i++) {
        mfcc_result2[i] = (double *)malloc(NUM_MFCC * sizeof(double));
    }
    extract_mfcc(signal, mfcc_result2, signal_length);
    free(signal);

    // Calculate and display the distances for Frame 0 of file 1 vs all frames of file 2
    printf("\nDistances between frames:\n");
    int frame1 = 0; // Frame from file 1 to compare
    for (int j = 0; j < num_frames2; j++) {
        double distance = euclidean_distance(mfcc_result1[frame1], mfcc_result2[j], NUM_MFCC);
        printf("Frame %d of file 1 vs Frame %d of file 2: %.4f\n", frame1, j, distance);
    }

    // Free allocated memory for MFCC result of both files
    for (int i = 0; i < num_frames1; i++) {
        free(mfcc_result1[i]);
    }
    free(mfcc_result1);

    for (int i = 0; i < num_frames2; i++) {
        free(mfcc_result2[i]);
    }
    free(mfcc_result2);

    return 0;
}
