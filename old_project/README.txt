@saamianf

PCR Plate Gene Expression Analysis

The code that I wrote takes an input file that is a PCR machine’s output file (‘plate1.xlsx’
and 'plate2.xlsx’ in the submitted folder). My current code is applicable to the standard PCR
plate setup (image below) where every 2 rows is one primer and the 2 spots in that row in the
same column are for one condition. The first 3 columns are the full controls that are used for
reference for the additional 3-column comparisons for that primer. For any blanks in the plate
reads, 40 is used to adjust the calculations. The first two rows of the first plate are the control
primer that all the other primers are compared against.

The program runs the calculations for this data using the following steps for each primer:
1. Average the 2 values for the control primer (P1)
2. Average the 2 values for the experimental primer (Pn)
3. Subtract the control primer average from the experimental primer average (Pn -
P1) to get dT
4. Average the 3 dT value from the first 3 columns to get the control dT (dT1)
5. Subtract the dT average from the other dTs and use the absolute value (|dTn - dT1|)
to get ddT
6. Calculate the Power value as 2 to the power of negative ddT (2-ddT)

The output for each primer is then a list of 12 Power values. The output file (output.csv)
lists each primer (beginning with primer 2 since primer 1 is the control) and its 12 Power values
in the corresponding order of the conditions from the input file. The Power values are then
directly used on the Prism software for the data management, processing, analysis, and
visualization.

This program is generalizable to PCR plate analysis and speeds up the calculation and
analysis steps. It can be expanded on for use in other PCR project analyses and can serve as a
starting point for projects.

* This data was tested using my lab’s PCR plate data and checked using the manually calculated
values for accuracy.


Code - project.py
Sample inputs - plate1.xlsx + plate2.xlsx
Sample output - outputs.csv
