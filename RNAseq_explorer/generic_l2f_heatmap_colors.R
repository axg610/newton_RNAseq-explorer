generic_l2f_heatmap_colors <- function(
    mat,
    colorscale = c(
      "steelblue4", "steelblue2", "white", "firebrick2", "firebrick4"
    ),
    scaleType = c(
      "robust",      # compresses around -0.5 to 0.5 range
      "symmetric"    # symmetric range around zero
    ),
    clipFactor = 0.98
) {
  
  # Return a pretty log2 fold-change heatmap color palette given a matrix.
  
  scaleType <- match.arg(scaleType)
  
  if (scaleType == "robust") {
    
    colors <- circlize::colorRamp2(
      c(
        seq(
          quantile(mat, 1 - clipFactor),
          -0.05,
          length = 75
        ),
        seq(
          -0.49,
          0.49,
          length = 50
        ),
        seq(
          0.5,
          quantile(mat, clipFactor),
          length = 75
        )
      ),
      colorRampPalette(colorscale)(200)
    )
    
  } else if (scaleType == "symmetric") {
    
    maxVal <- quantile(abs(mat), clipFactor)
    
    colors <- circlize::colorRamp2(
      seq(
        -maxVal, 
        maxVal, 
        length = 200
      ),
      colorRampPalette(colorscale)(200)
    )
  }
  
  colors
}