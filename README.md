# Consistent dictionary learning for signal declipping

Clipping, or saturation, is a common nonlinear distortion in signal processing. Clipping occurs when the signal reaches a maximum threshold  and the waveform is truncated:

![Clipped signal](/Figures/clipped_glockenspiel.png)

Declipping aims at recovering the clipped samples using the surrounding unclipped samples. 

This code performs declipping using 4 different approaches:
* **Iterative Hard Thresholding (IHT) for inpainting:** discards the clipped sample and performs sparse coding on the unclipped samples
* **Dictionary learning for inpainting:** discards the clipped samples and performs dictionary learning on the unclipped samples
* **Consistent IHT:** performs sparse coding while enforcing the clipped samples to be higher than the clipping threshold
* **Consistent dictionary learning:** performs dictionary learning while enforcing the clipped samples to be higher than the clipping threshold

You can run an example of declipping a glockspiel signal by running:
```
run Experiment/Declip_1_signal.m
```
You can see how the different algorithms perform in reconstructing a heavily clipped signal (here 73% of the samples are clipped):

![IHT](/Figures/declip_glockenspiel_IHT.png)
![DL](/Figures/declip_glockenspiel_DL.png)
![consIHT](/Figures/declip_glockenspiel_consIHT.png)
![consDL](/Figures/declip_glockenspiel_consDL.png)
