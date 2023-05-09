# solrindexing

Useful tools and wrappers used for indexing MMD in SolR. This software is
developed for use in the context of Arctic Data Centre, supported through
projects SIOS KC and Norwegian Scientific Data Network.


# Usage

## Usage with directory and no thumbnails

```bash
./src/indexdata.py -c etc/config.yml -d /home/johannestl/Desktop/S-ENDA/SOLR/mmd-xml-dev/arch_5/arch_9/arch_6 -n
```

## Logger object
* `SOLRINDEXER_LOGFILE` can be set to enable logging to file.
* `SOLRINDEXER_LOGLEVEL` can be set to change log level. See the Debugging section below.

# Installation

Solrindexing depends on cartopy, which cannot be installed in the regular "pip" way. On linux,
first install libproj and libgeos:

```bash
sudo apt-get install libproj-dev libgeos++-dev
```

To avoid problems with conflicting versions, we recommend using the Conda package manager. The below steps assume that Conda is installed.

Create a new environment:

```bash
conda create --name testenv python=3.11
```

Verify that the new environment is registered:

```text
conda info --envs
```

Activate (use) the new environment:

```text
conda activate testenv
```

Install the solrindexer package:

```text
pip install .
```

All dependencies should now be installed.
