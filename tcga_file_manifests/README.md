# How to Download TCGA Data

These two files contain all of the files needed from the TCGA data portal.
These data were downloaded before the migration to the new data portal website.
To download the data using the current website, the easiest method is to 
use the GDC Data Transfer tool. To do this, following these steps:

1. Download the tool [here](https://gdc.cancer.gov/node/159)
2. Add the `gdc-client` path to your `$PATH` environment variable, or run it
from your downloads folder. For example `/Users/WarrenMcG/Downlaods/gdc-client -h`
3. Make sure that `GDC_CLIENT_PATH` variable is set to the correct path in 
the `step_1a.sh` script in the `step_1a` directory.
4. Run `step_1a.sh` in the `step_1a` directory

It runs the following command at the command line:
```
[gdc-client path] download -m [manifest file] -d [DIR]
``` 

If you type `gdc-client download -h` at the command line, you can see information
about other options that can facilitate download. You can modify `step_1a.sh` with these
other options.
