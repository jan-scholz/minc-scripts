# minc-scripts

various scripts processing minc files or minc-associated files


## Disclaimer

These scripts might be modified without notice and can change in a way that breaks backwards compatibility.
Some scripts might only be available in the developer's branch.


## Dependencies

`tag_cbind.py`, `pTFCE.py`, and `vtk_meshscale.py` depend on [numpy](http://www.numpy.org/)

`vtk_meshconvert.py` and `vtk_meshscale.py` depend on [VTK](http://www.vtk.org/download/), in particular VTK's Python bindings.
On a Mac install the [homebrew](http://brew.sh) package manager and then install VTK with `brew install vtk`.

`vtk_meshscale.py` depends on [rpy2](http://rpy.sourceforge.net). Install with `pip2 install rpy2`

`pTFCE.py` depends on [pyminc](https://github.com/Mouse-Imaging-Centre/pyminc).

`xfm_merge` depends on [minc tools](http://packages.bic.mni.mcgill.ca/minc-toolkit/)

