# Polygonization algorithm for the ISPRS PIA 2019 Munich conference
This page is devoted to the polygonization method used in our work \[2\] presented at the **ISPRS PIA+MRSS 2019 conference** (Munich, Germany), based on the work of Dutter \[1\]. 

It is aimed for the polygonization of building footprints. This algorithm works well on relatively complex boundaries,
yielding very accurate polygonal outlines, without requiring any supervision or manual editting from human operators.
However, they are only designed for buildings of perpendicular and parallel wings, as well as having perpendicular corners.

The research work \[2\] presents an unsupervised building extraction method, wrapped up by this polygonization step.
It was awarded the Best Poster Award at the PIA 2019 conference. You may find the paper [here](https://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XLII-2-W16/181/2019/isprs-archives-XLII-2-W16-181-2019.pdf).

The algorithm operates on 3 levels of shapes:
- Level 1: rectangle; 
- Level 2: L-, T- or Z-shape; 
- Level 3: U-shape. 

In this script, we used the algorithm to detect the Minimum Bounding Rectangular (MBR) written by Julien Diener contributing in Matlab Exchange \[[link](https://www.mathworks.com/matlabcentral/fileexchange/31126-2d-minimal-bounding-box)\].

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

## Code organization
The algorithm is included in the main script named **polygonization_dutter.m**. 
To perform a test run, run the script **test.m**.
Three examples are provided in folder **example**. They are real building boundaries.

### Parameters
Two threshold parameters are required:
- Vspec: standard deviation of distances between the points in a segment and its corresponding edge. (Default: 2)
- Mspec: minimum length of resulting edges. (Default: 2)

These two parameters are in **meter** units (not in pixel).

## Examples

A raw boundary is given in blue, with its MBR in black-dashed. At each level, the boundary is divided into smaller segments, color-coded for visual purposes. In the Result figures, the resulting polygons are plotted in red outlines. You may click on the figures for the full size.

### A rectangle boundary
Raw boundary | Level 1 segmentation | Result |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/rect.png" alt="Raw boundary" width="90%" height="22%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l1_rect_V1_M1.png" alt="Level 1" width="90%" height="22%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_rect_V1_M1.png" alt="Result" width="90%" height="22%">

### An L-shape boundary

Vspec=2, Mspec=2 
Raw boundary | Level 1 segmentation | Result |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/Lshape.png" alt="Raw boundary" width="90%" height="22%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l1_Lshape_V2_M2.png" alt="Level 1" width="90%" height="22%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_Lshape_V2_M2.png" alt="Result" width="90%" height="22%"/>

**Vspec=1**, Mspec=2
Level 1 segmentation | Level 2 segmentation | Result |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l1_Lshape_V2_M2.png" alt="Level 1" width="90%" height="22%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l2_Lshape_V1_M1.png" alt="Level 2" width="90%" height="22%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_Lshape_V1_M1.png" alt="Result" width="90%" height="22%"/>

### A U-shape boundary
Vspec=2, Mspec=4 
Raw boundary | Douglas-Peucker algorithm \[3\] | Level 1 segmentation |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/Ushape.png" alt="Raw boundary" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/DP_Ushape_tol2.png" alt="DP" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l1_Ushape_V2_M2.png" alt="Level 1" width="100%" height="25%"/>
Level 2 segmentation | **Level 3 segmentation** | **Result** |
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l2_Ushape_V2_M4.png" alt="Level 2" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l3_Ushape_V2_M4.png" alt="Level 3" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_Ushape_V2_M4.png" alt="Result" width="100%" height="25%"/>

Vspec=2, **Mspec=2**
Level 2 segmentation | **Level 3 segmentation** | **Result** |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l2_Ushape_V2_M2.png" alt="Level 2" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l3_Ushape_V2_M2.png" alt="Level 3" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_Ushape_V2_M2.png" alt="Result" width="100%" height="25%"/>

**Vspec=0.9**, Mspec=2
| Level 2 segmentation | **Level 3 segmentation** | **Result** |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l2_Ushape_V2_M2.png" alt="Level 2" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l3_Ushape_V09_M2.png" alt="Level 3" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_Ushape_V09_M2.png" alt="Result" width="100%" height="25%"/>

**Vspec=0.7**, Mspec=2
Level 2 segmentation | **Level 3 segmentation** | **Result** |
:-------------------------:|:-------------------------:|:-------------------------:|
<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l2_Ushape_V2_M2.png" alt="Level 2" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/l3_Ushape_V07_M2_new.png" alt="Level 3" width="100%" height="25%"/>|<img src="https://github.com/nthuy190991/polygonization_PIA2019/blob/master/figure/res_Ushape_V07_M2_new.png" alt="Result" width="100%" height="25%"/>

### Some notes
- If you put too small Vspec or Mspec value, there may be too few points (e.g. <3) in some splitted segment(s). This can make the polygon become incorrect due to the orthogonality between adjacent polygonal segments. Therefore, my advice is to begin with big Vspec and Mspec value, and then gradually reduce them.
- Despite the "surmise" that Z-shape boundaries are included in level 2, we found that they were very tricky to process, compared to the other shapes. Zigzag boundaries are even trickier. I may try to solve this later.
- One main drawback of this algorithm is that it only acts forward. Once a segment is splitted, it is splitted. There is no way to review/retroact on the previous segmentation steps. We may look into this for future improvements.


## Reference
\[1\] Dutter, M. (2007). "Generalization of building footprints derived from high resolution remote sensing data", Institut für Photogrammetrie und Fernerkundung, Technische Universität Wien.

\[2\] T. H. Nguyen et al. (2019). "Unsupervised Automatic Building Extraction Using Active Contour Model on Unregistered Optical Imagery and Airborne LiDAR Data," International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XLII-2/W16, 181-188. DOI: 10.5194/isprs-archives-XLII-2-W16-181-2019.

\[3\] Douglas, D. H., Peucker, T. K. (1973). "Algorithms for the reduction of the number of points required to represent a digitized line or its caricature". Cartographica: The International Journal for Geographic Information and Geovisualization, 10(2), 112– 122.


## Questions/Discussions
For any other questions/issues, please open an issue on the Issues tracker.
