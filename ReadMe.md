# Cascaded noise reduction and acoustic echo cancellation based on an extended noise reduction
# License
This work is licensed under the [MIT LICENSE](LICENSE.md). By downloading and/or installing this software and associated files on your computing system you agree to use the software under the terms and conditions as specified in the license agreement.

If this code has been useful for you, please cite [[1]](#References).

# About
This repository [[2]](#References) contains the MATLAB code associated with [[1]](#References). For the problem of combined noise reduction (NR) and acoustic echo cancellation (AEC), generally, a cascaded design is adhered to where an NR precedes an AEC (NR-AEC). The NR aims at cancelling the near-end room noise (and possibly partially the echo) and operates on the microphones, and the AEC aims at cancelling the echo. On the other hand, in [[1]](#References) a cascaded design is proposed where an extended NR (NR<sub>ext</sub>) precedes an AEC (NR<sub>ext</sub>-AEC). The NRext aims at cancelling the near-end room noise and the far-end room noise component in the echo and operates on both the microphones and loudspeakers, and the AEC aims at cancelling the far-end room speech component (and residual far-end room noise component) in the echo. The NR<sub>ext</sub>-AEC design has the advantages over NR-AEC of having AEC filters that are independent of the NR<sub>ext</sub> filters, and of having NRext filters of which the degrees of freedom scale with the number of loudspeakers.

This repository provides example codes for both NR-AEC and NR<sub>ext</sub>-AEC. 

The code has been developed and tested in MATLAB R2021b.

The manuscript of [[1]](#References) has also been released as a preprint [[3]](#References).

# File structure

* [Main.m](Main.m): Main file to run the code.
* [LICENSE](LICENSE.md): License file.
* [ReadMe.md](ReadMe.md): ReadMe file.
* [Audio](Audio): Folder containing the audio files.
    - [impulse.mat](Audio/impulse.mat): File containing the impulse responses between loudspeakers and microphones.
    - [sig.mat](Audio/sig.mat): File containing the desired speech, near-end room noise, echo, and loudspeaker signals. The speech and speech component in the echo are taken from the [VCTK corpus](https://datashare.ed.ac.uk/handle/10283/3443), which is made available under the Open Data Commons Attribution License (ODC-By) v1.0.
* [Util](Util): Auxiliary code.
    - [AEC](Util/AEC): Folder the containing auxiliary code for AEC.
        + [compute_AEC.m](Util/AEC/compute_AEC.m): Calculates the AEC filters, where the filters are calculated based on the entire data.
        + [compute_AEC_adaptive.m](Util/AEC/compute_AEC_adaptive.m): Calculates the AEC filters, where the filters are computed adaptively.
    - [AEC+NR](Util/AEC+NR): Folder containing the auxiliary code for combined AEC and NR.
        + [process_NRAEC.m](Util/AEC+NR/process_NRAEC.m): NR-AEC processing, where the filters are calculated based on the entire data.
        + [process_NRAEC_adaptive.m](Util/AEC+NR/process_NRAEC_adaptive.m): NR-AEC processing, where the filters are computed adaptively.
        + [process_NRextAEC.m](Util/AEC+NR/process_NRextAEC.m): NR<sub>ext</sub>-AEC processing, where the filters are calculated based on the entire data.
        + [process_NRextAEC_adaptive.m](Util/AEC+NR/process_NRextAEC_adaptive.m): NR<sub>ext</sub>-AEC processing, where the filters are computed adaptively.
    - [Frequency_transformation](Util/Frequency_transformation): Folder containing the auxiliary code for the conversion between time- and frequency-domain.
        + [WOLA_analysis.m](Util/Frequency_transformation/WOLA_analysis.m): Weighted overlapp add (WOLA) analysis filterbank.
        + [WOLA_synthesis.m](Util/Frequency_transformation/WOLA_synthesis.m): WOLA synthesis filterbank.    
        + [WOLA2distortion.m](Util/Frequency_transformation/WOLA2distortion.m): Conversion of the short-time Fourier transform (STFT) filters to their time-domain equivalent distortion functions. 
    - [Metrics](Util/Metrics): Folder containing the auxiliary code for the evaluation of the algorithms.
        + [align_proc_unproc.m](Util/Metrics/align_proc_unproc.m): Aligns processed and unprocessed data.  
        + [compute_metrics.m](Util/Metrics/compute_metrics.m): Computes the signal-to-noise ratio (SNR), signal-to-echo ratio (SER), and speech distortion (SD).
        + [SD.m](Util/Metrics/SD.m): Computes the SD.
        + [SNR.m](Util/Metrics/SNR.m): Computes the SNR .
    - [MWF](Util/MWF): Folder containing the auxiliary code for the multichannel Wiener filters (MWFs).
        + [updateDifferenceCorrelation.m](Util/MWF/updateDifferenceCorrelation.m): Calculates the difference between two correlation matrices using a generalised eigenvalue decgit stomposition (GEVD) approximation.
        + [updateMWFGEVD.m](Util/MWF/updateMWFGEVD.m): Computes an MWF with a GEVD approximation for a particular STFT bin. 
        + [updateMWFGEVDMultichannel.m](Util/MWF/updateMWFGEVDMultichannel.m): Computes an MWF with a GEVD approximation for each STFT bin.
    - [NR](Util/NR): Folder containing the auxiliary code for NR.
        + [compute_NR.m](Util/NR/compute_NR.m): Calculates the NR<sub>ext</sub> filters, where the filters are calculated based on the entire data.
        + [compute_NR_adaptive.m](Util/NR/compute_NR_adaptive.m): Calculates NR<sub>ext</sub> filters, where the filters are computed adaptively.
        + [compute_NRext.m](Util/NR/compute_NRext.m): Calculates the NR<sub>ext</sub> filters, where the filters are calculated based on the entire data.
        + [compute_NRext_adaptive.m](Util/NR/compute_NRext_adaptive.m): Calculates the NR<sub>ext</sub> filters, where the filters are computed adaptively.
    - [VAD](Util/VAD): Folder containing the auxiliary code for the voice activity detection (VAD).
        + [VAD.m](Util/VAD/VAD.m): Computes the VAD in the STFT domain.

# Contact
Arnout Roebben, Toon van Waterschoot, and Marc Moonen\
Department of Electrical Engineering (ESAT)\
STADIUS Center for Dynamical Systems, Signal Processing and Data Analytics\
KU Leuven\
Leuven, Belgium\
E-mail: <arnout.roebben@esat.kuleuven.be>

# Acknowledgements
This research was carried out at the ESAT Laboratory of KU Leuven, in the frame of Research Council KU Leuven C14-21-0075 ”A holistic approach to the design of integrated and distributed digital signal processing algorithms for audio and speech communication devices”, and Aspirant Grant 11PDH24N (for A. Roebben) from the Research Foundation - Flanders (FWO).

# References
[1] 
```
@inproceedings{roebbenCascaded2024,
  author={Roebben, Arnout and van Waterschoot, Toon and Moonen, Marc},
  booktitle={2024 32nd European Signal Processing Conference (EUSIPCO)}, 
  title={CASCADED NOISE REDUCTION AND ACOUSTIC ECHO CANCELLATION BASED ON AN EXTENDED NOISE REDUCTION}, 
  year={2024},
  keywords={Acoustic echo cancellation (AEC), Noise reduction
  (NR), Extended noise reduction (NRext), Multichannel}
  }
```

[2]
```
@misc{roebbenGithubRepositoryCascaded2024,
  title = {Github Repository: {{Cascaded}} noise reduction and acoustic echo cancellation based on an extended noise reduction},
  author = {Roebben, Arnout},
  year = {2024},
  journal = {GitHub},
  urldate = {2024},
  howpublished = {https://github.com/Arnout-Roebben/NRAEC\_vs\_NRextAEC},
  langid = {english}
}
```

[3]
```
@misc{roebben2024cascaded,
      title={Cascaded noise reduction and acoustic echo cancellation based on an extended noise reduction}, 
      author={Arnout Roebben and Toon van Waterschoot and Marc Moonen},
      year={2024},
      eprint={2406.08974},
      archivePrefix={arXiv},
      primaryClass={id='eess.AS' full_name='Audio and Speech Processing' is_active=True alt_name=None in_archive='eess' is_general=False description='Theory and methods for processing signals representing audio, speech, and language, and their applications. This includes analysis, synthesis, enhancement, transformation, classification and interpretation of such signals as well as the design, development, and evaluation of associated signal processing systems. Machine learning and pattern analysis applied to any of the above areas is also welcome.  Specific topics of interest include: auditory modeling and hearing aids; acoustic beamforming and source localization; classification of acoustic scenes; speaker separation; active noise control and echo cancellation; enhancement; de-reverberation; bioacoustics; music signals analysis, synthesis and modification; music information retrieval;  audio for multimedia and joint audio-video processing; spoken and written language modeling, segmentation, tagging, parsing, understanding, and translation; text mining; speech production, perception, and psychoacoustics; speech analysis, synthesis, and perceptual modeling and coding; robust speech recognition; speaker recognition and characterization; deep learning, online learning, and graphical models applied to speech, audio, and language signals; and implementation aspects ranging from system architecture to fast algorithms.'}
}
```