June 14th

software needed for consensus:
- demuxlet v2
- vireo
- scansnp
- souporcell
- dropkick 
- doubletfinder

- merge working version Davide
- remove demuxlet v1 from the pipe
- double check the value of the double prior https://github.com/GiuseppeTestaLab/demultiplexing-workflow/blob/eccc2fe17d374e73650a3374c7fde3ca8570ecff/workflow/Snakefile#L164
- add dropkick and doublefinder to demultiplexing container
- fix data for sourporcell run: remove all the working files as scratch data
- add venn.diagram to R in the container
