Feb 6, 2012

copyright 2011,2012,2013, CBCSL at The Ohio State University
Use of this software is granted for basic research at non-for profit institutions only. Commercial use of this software is not permitted. Corporations or for-profit institutions interested in the use of this software should contact Prof. Aleix M. Martinez.

These codes are for doing low level feature detection using Kernel LDA.
This is the initial feature detection in my paper, so please cite for academic work:

Benitez, F., Rivera, S., Gotardo, P., & Martinez, A. (2014). Salient and Non-Salient Fiducial Detection using a Probabilistic Graphical Model. Pattern Recognition, 47 (1), 208-215.

Make sure you add the included KLDA and src folder to the MATLAB path, and cite the paper which defines the parameter tuning:

Di You, Onur C. Hamsici and Aleix M. Martinez, Kernel Optimization in Discriminant Analysis, IEEE Trans. On Pattern Analysis and Machine Intelligence, 2011.


The main function is 'detectFiducialsFull.m'

See example1.m

The function 'setCropParams' is where you can change
     the crop size (original size to crop the fiducial from the image)
     and the target size (it scales it down before doing the KLDA)
