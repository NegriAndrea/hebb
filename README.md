# hebb
Halo Extreme Block Bootstrap package to estimate the distribution of N heaviest halo mass you can find in a survery at a given redshift.

`hebb` uses as a basedata the Uchuu simulation halo catalogue (https://skiesanduniverses.org/Simulations/Uchuu) to perform a block bootstrap by shooting boxes of a volume equal of the estimated volume of a survey, and recovering the N largest halo formed at a particular z, with an uncertainty estimate. The volume of the survey is computed given the survey's field of view and redshift depth, or manually selected. Optionally, the code can dump the list of halo found in the search, which can be used to perform a trace back in time with the Uchuu merger tree.

Please cite the following papers if you use `hebb` in your work: Negri & Belli (2026), [Ishiyama et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.506.4210I) (Uchuu Data Release 1).

## Installation
There are 2 ways to install `hebb`: via `PyPI`
```
python3 -m pip install hebb
```
Or directly from the repository `https://github.com/NegriAndrea/hebb/` and install locally
```
git clone https://github.com/NegriAndrea/hebb
cd hebb
python3 -m pip install -e .
```

#### Database Setup
In addition to installing the Python package, you must download a reduced version of the Uchuu database. Two versions are available, with a different cut in mass; a light [one](https://uses0-my.sharepoint.com/:u:/g/personal/anegri_us_es/IQCDx5T1XlTnSbM2YmYcFoYMATIk_pI1XKPzp4ycUJ6N12M?e=pDdoDg) (2.3 GB) and a more complete [one](https://uses0-my.sharepoint.com/:u:/g/personal/anegri_us_es/IQBrGP-3e0xHR6pV0YkfCb1hAYS1KLSNnjmDvVb3H6ytsUA?e=fcudjY) (~40 GB). To set the database path you can define the following environment variable in your `~/.bash_profile`
```
# bash
export HEBB_DB_PATH=/path/to/database
```
where you have to change `/path/to/database` to the path of the downloaded database. By default, the code looks first for the most complete catalugue, and if it is not found, it will look for the light one.

NOTE: In order to keep the file size manageble the database contains only log10(M200), positions and merger tree IDs; the halo positions have been binned on a 40 ckpc gridsize and stored as `uint16` integer, which is precise enough for volumes that are usually way larger than 1 cMpc.

#### Test installation
The easiest way to test if everything is correctly installed and set up is to run the command
```
hebb 100 0. -L 1000 -M 1e14
```

## Quickstart
To estimate the 5 most massive halo at z=3 (with uncertaintes) for a survey that spans from z=3 to 4, having field-of-view of 100 arcsec², using 5000 block bootstrap iterations: 
```
hebb 5000 3. --survey 3. 4. 100. -M 1e12 -n 5 --plot -t
```
where we limited the catalogue galaxies to be larger than $10^{12}~M_\odot$ for a faster database load (see later for usage of `-M`). The option `-t` creates a table with the sampled haloes, their IDs to trace them back and forth in time with the Uchuu merger tree, and a table containing a list of percentiles for every mass rank.

## Full Usage
The simplest way to use `hebb` is via command line, `hebb -h` returns the user manual:
```
usage: hebb [-h] (--survey z_min z_max fov | -L L) [-n N] [-t] [--plot] [-M M] [--lf LF] [--force-light] Nboxes z_target

Compute the N most massive dark matter haloes that you can find in a given survey with [z_min, z_max] and field-of-view by performing a non-parametric block bootstrap over the Uchuu (2/h cGpc)^3 run.

positional arguments:
  Nboxes                Number of boxes for bootstrap
  z_target              Redshift of your target

options:
  -h, --help            show this help message and exit
  -n N                  Track the N most massive haloes in each box (1 tracks only the most massive) [default: 1]
  -t                    Create a table with the sampled haloes
  --plot                Show a plot of the M200 distribution
  -M M                  OPTIMIZATION: Database mass cut, greatly speed up the database loading and search but you can incur into empty boxes [default: None]
  --lf LF               OPTIMIZATION: Leafe size for each node of the KDTree [default: 128]
  --force-light         DEBUG: force the reading of the light catalogue first

Processing mode (required, mutually exclusive):
  --survey z_min z_max fov
                        Survey min z, max z, FOV in arcmin^2 (used to compute box volume)
  -L L                  Size of the box in cMpc, alternative to the boxsize computation from the FOV and z-depth of the survey
```

## API
TBD



## Citations
```
@ARTICLE{2021MNRAS.506.4210I,
       author = {{Ishiyama}, Tomoaki and {Prada}, Francisco and {Klypin}, Anatoly A. and {Sinha}, Manodeep and {Metcalf}, R. Benton and {Jullo}, Eric and {Altieri}, Bruno and {Cora}, Sof{\'\i}a A. and {Croton}, Darren and {de la Torre}, Sylvain and {Mill{\'a}n-Calero}, David E. and {Oogi}, Taira and {Ruedas}, Jos{\'e} and {Vega-Mart{\'\i}nez}, Cristian A.},
        title = "{The Uchuu simulations: Data Release 1 and dark matter halo concentrations}",
      journal = {\mnras},
         year = 2021,
        month = sep,
       volume = {506},
       number = {3},
        pages = {4210-4231},
          doi = {10.1093/mnras/stab1755},
archivePrefix = {arXiv},
       eprint = {2007.14720},
 primaryClass = {astro-ph.CO},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021MNRAS.506.4210I},
}
```
