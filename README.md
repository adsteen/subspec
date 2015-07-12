# subspec
## Purpose
This R package exists only to recreate the analysis from the [Steen et al paper](http://www.int-res.com/prepress/a01755.html) from Aquatic Microbial Ecology about the substrate specificity of aquatic peptidases.

## Use
If you want to see that I didn't make up the plots, use R to load the package and run the `write_paper()` function:

```library(devtools)
install_github("adsteen/subspec")
write_paper(print_plots=TRUE)
```

If you wish to use the data, load the packages: the data files listed in `/data/` will be available. (Or download the `.rda` files directly).

If you wish to use the code for similar analysis, check inside `write_paper()` - or check [http://github.com/adsteen/enzalyze], which will soon contain code to do the same thing for general purposes. Or just [contact me](mailto:asteen1@utk.edu) - I'll be happy to help you out.
