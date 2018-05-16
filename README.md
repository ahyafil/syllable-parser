# syllable-parser
Matlab code for the PINTH neural network (theta oscillator) locking to syllables in speech signal (Hyafil et al., eLife 2015).

Main function is syllabledecoder.m. Type 'help syllabledecoder.m' to use.
Make sure to add directory to path. You will also need to download the NSL matlab toolbox (the model of pre-cortical acoustic processing in Chi et al, 2005) from: http://www.isr.umd.edu/Labs/NSL/Software.htm.

Typos from the manuscripts and missing parameters are described in "Hyafil 2015 Speech parsing typos".

The code for the entire network (with the PING module) is not provided but should easily be implemented by editing this...

Execution time can be made faster by use the mex mode. For this you will need to compile the mex code for network simulation: type "mex lif2.c" (type "mex -setup" to make sure you have a MEX compiler installed).

For any comment/bug report, please write to alexandre.hyafil (at) gmail.com. Have fun!
(Thank to Keith for his precious time debugging this)

References:
Hyafil, A., Fontolan, L., Kabdebon, C., Gutkin, B. S., & Giraud, A.-L. (2015). Speech encoding by coupled cortical theta and gamma oscillations. eLife, 4, 1–23. https://doi.org/10.7554/eLife.06213

Chi, T., Ru, P., Shamma SA. (2005) Multiresolution spectrotemporal analysis of complex sounds. The Journal of the Acoustical Society of America 118:887–906.



