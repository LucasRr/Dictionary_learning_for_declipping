# Consistent dictionary learning for signal declipping

<!--- This is a markdown README generated for the github page. For a human readable readme, see README.txt --->


This is an implementation of the algorithm proposed in:

*Consistent dictionary learning for signal declipping*, L. Rencker, F. Bach, W. Wang, M. D. Plumbley, 
Latent Variable Analysis and Signal Separation (LVA/ICA), Guildford, UK, 2018

The pre-print can be seen at `consistent_DL_LVA18_preprint.pdf`.

## Author:
Lucas Rencker,  
*Centre for Vision, Speech and Signal Processing (CVSSP)*, University of Surrey, UK  
Contact: lucas \[dot\] rencker \[at\] surrey.ac.uk  

## Quick demo:

Clipping, or saturation, is a common nonlinear distortion in signal processing. Clipping occurs when the signal reaches a maximum threshold  and the waveform is truncated:

![Clipped signal](/Figures/clipped_glockenspiel.png)

Declipping aims at recovering the clipped samples using the surrounding unclipped samples. 

This code performs declipping using 4 different approaches:
* **Iterative Hard Thresholding (IHT) for inpainting:** discards the clipped sample and performs sparse coding on the unclipped samples using IHT and a fixed DCT dictionary
* **Dictionary learning for inpainting:** discards the clipped samples and performs dictionary learning on the unclipped samples
* **Consistent IHT:** performs consistent IHT using a fixed DCT dictionary \[1\]
* **Consistent dictionary learning:** performs consistent dictionary learning using the algorithm proposed in \[2\]

You can see an example of declipping a glockspiel signal by running:
```
run Experiment/Declip_1_signal.m
```
You can see how the different algorithms perform in reconstructing a heavily clipped signal (here 73% of the samples are clipped):

![IHT](/Figures/declip_glockenspiel_IHT.png)
![DL](/Figures/declip_glockenspiel_DL.png)
![consIHT](/Figures/declip_glockenspiel_consIHT.png)
![consDL](/Figures/declip_glockenspiel_consDL.png)

## References:
\[1\]: Consistent iterative hard thresholding for signal declipping, 
   S. Kitic, L. Jacques, N. Madhu, M. P. Hopwood, A. Spriet, C. De Vleeschouwer, ICASSP, 2013

\[2\]: Consistent dictionary learning for signal declipping, 
    L. Rencker, F. Bach, W. Wang, M. D. Plumbley,
    Latent Variable Analysis and Signal Separation (LVA/ICA), Guildford, UK, 2018


