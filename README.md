# Polygonization algorithms written for the PIA2019 Munich conference
This page is devoted to the polygonization method, written for our work \[2\] presented at ISPRS PIA+MRSS 2019, (Munich, Germany), originally based on the work of Dutter \[1\]. 

It is aimed for the polygonization of building footprints. It yields very accurate polygonal outlines, without requiring any supervision or manual editting from human operators. However, they are only designed for buildings with perpendicular corners.

The contribution \[2\] (in which this polygonization participating) is awarded the Best Poster Award at the conference. You may find the poster [here](https://www.dropbox.com/s/l8xr08ksu2io0ho/poster.pdf?dl=0).


The two algorithms operate on 3 levels of shapes:
- Level 1: rectangle; 
- Level 2: L-, T- or Z-shape; 
- Level 3: U-shape. 

Author: [Thanh Huy Nguyen](mailto:nthuy190991@gmail.com).

## Citation
If you use these codes, cite the paper \[[link](https://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XLII-2-W16/181/2019/isprs-archives-XLII-2-W16-181-2019.pdf)\]:
```
@article{nguyen2019unsupervised,
	Author = {Nguyen, T. H. and Daniel, S. and Gu\'eriot, D. and Sint\`es, C. and Le Caillec, J.-M.},
	Doi = {10.5194/isprs-archives-XLII-2-W16-181-2019},
  Journal = {ISPRS - International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences},
	Pages = {181--188},
	Title = {Unsupervised Automatic Building Extraction Using Active Contour Model on Unregistered Optical Imagery and Airborne LiDAR Data},
	Volume = {XLII-2/W16},
	Year = {2019},
}
```

## Some examples
SRSM results overlapped on z-image (tile 4190) |  SRSM results overlapped on z-image (tile 4785)
:-------------------------:|:-------------------------:

## Reference
\[1\] Dutter, M. (2007). "Generalization of building footprints derived from high resolution remote sensing data", Institut für Photogrammetrie und Fernerkundung, Technische Universität Wien.

\[2\] T. H. Nguyen et al. (2019). "Unsupervised Automatic Building Extraction Using Active Contour Model on Unregistered Optical Imagery and Airborne LiDAR Data," Int. Arch. Photogramm. Remote Sens. Spatial Inf. Sci., XLII-2/W16, 181-188. DOI: 10.5194/isprs-archives-XLII-2-W16-181-2019 

## Questions/Discussions
For any other questions/issues, please open an issue on the Issues tracker.
