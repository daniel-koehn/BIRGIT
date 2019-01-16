# BIRGIT
Pre- and Postprocessing tools for full waveform inversion (FWI) field data applications written in Matlab, Python and shell scripts. Since 2013 we developed a collection of codes for field data preprocessing for SH-FWI using the [DENISE FWI code](https://github.com/daniel-koehn/DENISE-Black-Edition). These consist of 

- Seismic Unix shell scripts to distribute seismic data from a single SU file to multiple SU files for each shot

- Matlab codes covering wiggle/image plot for seismic data visualization, shot/trace-normalization, 3D-2D spreading correction, frequency filtering, trace muting, offset-windowing, time-windowing/damping, source wavelet inversion, first arrival picking, amplitude trend plots, phase-velocity frequency spectra, amplitude/phase spectra

- Python codes for visualization of FWI results from DENISE, RAJZEL and GERMAINE

BIRGIT is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.0 of the License or higher.

BIRGIT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License in LICENSE.md for more details.

If you show inversion results in a paper or presentation using BIRGIT as Pre- or Postprocessing tool please give a reference to the following papers:

- Eva Dokter, Daniel Köhn, Dennis Wilken, Denise De Nil and Wolfgang Rabbel (2017): Full-waveform inversion of SH- and Love-wave data in near-surface prospecting, Geophysical Prospecting, 65(S1), 216-236, DOI: 10.1111/1365-2478.12549 [(Link)](https://onlinelibrary.wiley.com/doi/abs/10.1111/1365-2478.12549)

- Daniel Köhn, Dennis Wilken, Denise De Nil, Tina Wunderlich, Wolfgang Rabbel, Lukas Werther, Johannes Schmidt, Christoph Zielhofer, Sven Linzen (2019): Comparison of time-domain SH waveform inversion strategies based on sequential low and bandpass filtered data for improved resolution in near-surface prospecting, Journal of Applied Geophysics, 160, 69-83, DOI: 10.1016/j.jappgeo.2018.11.001 [(Link)](https://www.sciencedirect.com/science/article/pii/S0926985118303720?via%3Dihub)

We also use Segymat by Thomas Meier Hansen to read SEGY/SU data and code snippets from SeismicLab by University of Alberta for spectral analysis. Therefore, you should also give a reference to:

[Segymat](https://github.com/cultpenguin/segymat)

[SeismicLab](http://seismic-lab.physics.ualberta.ca/index.html)

BIRGIT is developed by 

Daniel Koehn (Christian-Albrechts-University Kiel, Germany)

with contributions and bug fixes from 

Eva Dokter (University of Edinburgh, United Kingdom)

Denise De Nil (Christian-Albrechts-University Kiel, Germany)

(add future authors here)
