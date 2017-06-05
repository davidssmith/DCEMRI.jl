DCEMRI.jl NEWS
==============

Changes in v0.2.0
-----------------

- New convolution method with tremendous speed increase from [@notZaki](https://github.com/notZaki) ([#29])
- Fixed bug where not all points were fit when the number of voxels was not a multiple of the number of workers.
- Now supports Julia v0.6 and v0.7. Julia v0.4 no longer supported ([#32]).

[#29]: http://github.com/davidssmith/DCEMRI.jl/issues/29
[#32]: http://github.com/davidssmith/DCEMRI.jl/issues/32
