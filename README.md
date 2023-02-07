# solrindexing
Useful tools and wrappers used for indexing MMD in SolR. This software is
developed for use in the context of Arctic Data Centre, supported through
projects SIOS KC and Norwegian Scientific Data Network.


# Usage

## Logger object
* `SOLRINDEXER_LOGFILE` can be set to enable logging to file.
* `SOLRINDEXER_LOGLEVEL` can be set to change log level. See the Debugging section below.

## Usage with directory and no thumbnails

```bash
./src/indexdata.py -c etc/config.yml -d /home/johannestl/Desktop/S-ENDA/SOLR/mmd-xml-dev/arch_5/arch_9/arch_6 -n
```