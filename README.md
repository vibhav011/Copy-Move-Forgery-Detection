# Copy-Move-Forgery-Detection

**Topic** - ​ Copy-move Forgery Detection in Images

**Abstract** - ​ With rapid advancement in photo editing technology, image tampering has become very
common. Duplicating continuous portions of pixels in an image, after possible geometrical and
illumination adjustments is often used to forge images since it is generally hard to detect by the human
eye. In this project, we aim to detect copy-move forgery in images by accomplishing the following
sub-tasks-

  1. **Scale invariant feature extraction**​ - In order to extract certain key features from an image, we
  made use of the SIFT algorithm. 
  
  2. **Affine transformation estimation** ​ - We have used the Random Sample Consensus (RANSAC)
  algorithm which estimates the model parameters of the affine transform to account for the
  possible geometric distortions of the duplicated regions such as rotation, scaling, and shearing
  that are supported in most photo editing software. We have intoduced the best fit method for transformation.
  
  3. **Matching patch retrieval** ​ - With the estimated affine transform, we have calculated correlation
  coefficients for each matching pair and marked the regions duplicated/not-duplicated based on a
  predefined threshold.
