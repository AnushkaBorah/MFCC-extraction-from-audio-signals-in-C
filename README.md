# Project Title: MFCC Feature Extraction From Speech Signals in C

## Description
The primary objective of this project is to develop a robust implementation of MFCC (Mel-Frequency Cepstral Coefficients) feature extraction in C that can process speech signals and further analysis of the MFCCs obtained using the Euclidean Distance Formula.
The functions in the code implement the various steps in the MFCC Algorithm, which include:
- Pre-emphasis
- Framing
- Windowing
- FFT (Fast Fourier Transform)
- Calculation of Power Spectrum
- Mel-Energies
- Calculation of the logarithm of Mel-Energies
- DCT (Discrete Cosine Transform)
and the calculation of the Euclidean Distance.

## Table of Contents
- [Installation](#installation)
- [Execution](#execution)
- [Input Files](#input-files)
- [Output](#output)

## Installation
1. Download and install Code::Blocks 20.03
2. Download and install Cool Edit 2000

## Execution
1. Open Code::Blocks
   - Launch Code::Blocks IDE.
2. Load the Project
   - Open the project file in Code::Blocks in a workspace. 
   - Upload the input files in the 'Others' folder of workspace.
   - Make sure path for the input file has been set correctly.
3. Compile the Code
   - Compile the code to ensure there are no errors.
4. Run the Code
   - Execute the code within Code::Blocks.

## Input Files
- The input file should consist of a dataset of a vowel. 
- For the input file:
  1. Use Cool Edit 2000 to record the vowel.
  2. Save the recording in ASCII form.
  3. Slice the file to capture the steady state of the signal to obtain the correct output and avoid unnecessary frames.

## Output
- For the code:
  - 13 MFCCs per frame for both files being processed, and the Euclidean Distance between them.
