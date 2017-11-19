# Harris Corner detector

Implementation of Harris Corner detector, with an optional adaptative non-maximal suppression method.

# Experimental review

You can find an experimental review of this implementation in the ![here](CR.pdf). You can find the images from the review in the ![images directory](images/) aswell.


# Requirements

```
scipy==0.19.1
matplotlib==2.1.0
numpy==1.13.3
scikit_image==0.13.1
skimage==0.0
```

# Usage
```
usage: harris_detector.py [-h] [-sd SIGMA_D] [-si SIGMA_I] [-k KAPPA]
                          [-l LOCAL] [-t THRESHOLD] [-a] [-c ANMS_CONSTANT]
                          [-b BEST] [-s SAVE] [-p PLOT]
                          image

positional arguments:
  image                 Indicate an image

optional arguments:
  -h, --help            show this help message and exit
  -sd SIGMA_D, --sigma_d SIGMA_D
                        Standard deviation for gaussian kernel used in
                        derivation
  -si SIGMA_I, --sigma_i SIGMA_I
                        Standard deviation for gaussian kernel used in
                        windowing
  -k KAPPA, --kappa KAPPA
                        Kappa parameter for Harris Mc computation
  -l LOCAL, --local LOCAL
                        Size of local window arround an extremum (2xlocal + 1)
                        x (2xlocal + 1) pixels)
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold for a valid corner detection
  -a, --anms
  -c ANMS_CONSTANT, --anms_constant ANMS_CONSTANT
                        Constant in ANMS
  -b BEST, --best BEST  Best rankings in ANMS
  -s SAVE, --save SAVE  Output file for corners locations
  -p PLOT, --plot PLOT  Plot option
```

# Example

![images directory](images/Figure_readme) 
