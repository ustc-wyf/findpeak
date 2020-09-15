# Find peak - Long test:
- 2 PDF files: papers about peak search algorithm. The Method used is based on the paper "Automatic program for peak detection and deconvolution of multi-overlapped chromatographic signals Part I Peak detection"
- FIRfilter (1): simple code for digital filtering, may need in the future.
- Test_noise_figure.py: Test on different smooth algorithm.
- Test_find_peak.py: Main program for peak search. The program smooth the data, calculate the 2nd deviation, find the peak based on predefined threshold and user definded threshold. 
- 4 data txt files: data used for code testing.

# Note:
- To be clear: we do NOT use peak fitting function to seach for POSITION of the peak. I though we made it clear the last time.
- With this peak seach algorithm, you can developt the windows around the peak and cut the data out and use in your peak fitting code.
- The data structure is changed
 + Before, we put the title at the beginning of the file to record important parameter of the experiment. We will NOT remove it. You can simple use numpy.genfromtxt and add skip_headers = 1 to read the file.
 + We update the data structure. If you use the txt file and read, you will output a table with 5 column. The x_axis is the 1st column and the y_axis is the 2nd column.
 + If you want to test the "vna_cont_spec_pulselike.txt" file, it will only contain 2 column.