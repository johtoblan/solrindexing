#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
SOLR-indexer : Main script
==========================

Copyright 2021 MET Norway

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import os
import yaml
import logging
import argparse
import cartopy.crs as ccrs

from time import sleep
from requests.auth import HTTPBasicAuth
from solrindexer.indexdata import MMD4SolR, IndexMMD

logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--always_commit', action='store_true',
                        help='Specification of whether always commit or not to SolR')
    parser.add_argument('-c', '--cfg', dest='cfgfile', required=True,
                        help='Configuration file')
    parser.add_argument('-i', '--input_file',
                        help='Individual file to be ingested.')
    parser.add_argument('-l', '--list_file',
                        help='File with datasets to be ingested specified.')
    parser.add_argument('-d', '--directory',
                        help='Directory to ingest')
    parser.add_argument('-t', '--thumbnail', action='store_true',
                        help='Create and index thumbnail, do not update the main content.')
    parser.add_argument('-n', '--no_thumbnail', action='store_true',
                        help='Do not index thumbnails (done automatically if WMS available).')
    # parser.add_argument('-f','--feature_type', action='store_true',
    #                    help='Extract featureType during ingestion (to be done automatically).')
    parser.add_argument('-r', '--remove',
                        help='Remove the dataset with the specified identifier'
                        '(to be replaced by searchindex).')
    parser.add_argument('-2', '--level2', action='store_true', help='Operate on child core.')

    # Thumbnail parameters
    parser.add_argument('-m', '--map_projection', required=False,
                        help='Specify map projection for thumbnail '
                        '(e.g. Mercator, PlateCarree, PolarStereographic).')
    parser.add_argument('-t_layer', '--thumbnail_layer', required=False,
                        help='Specify wms_layer for thumbnail.')
    parser.add_argument('-t_style', '--thumbnail_style', required=False,
                        help='Specify the style (colorscheme) for the thumbnail.')
    parser.add_argument('-t_zl', '--thumbnail_zoom_level', type=float, required=False,
                        help='Specify the zoom level for the thumbnail.')
    parser.add_argument('-ac', '--add_coastlines', const=True, nargs='?', required=False,
                        help='Add coastlines too the thumbnail (True/False). Default True')
    parser.add_argument('-t_extent', '--thumbnail_extent', required=False, nargs='+', 
                        help='Spatial extent of thumbnail in lat/lon degrees like "x0 x1 y0 y1"')

    args = parser.parse_args()

    if args.cfgfile is None:
        parser.print_help()
        parser.exit()
    if not args.input_file and not args.directory and not args.list_file and not args.remove:
        parser.print_help()
        parser.exit()

    return args


def main():

    # Parse command line arguments
    try:
        args = parse_arguments()
    except Exception as e:
        logger.error("Something failed in parsing arguments: %s", str(e))
        return 1

    IDREPLS = [':', '/', '.']

    tflg = l2flg = nflg = False
    if args.level2:
        l2flg = True

    # CONFIG START
    # Read config file, can be done as a CONFIG class, such that argparser can overwrite duplicates
    with open(args.cfgfile, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

    # Specify map projection
    if args.map_projection:
        map_projection = args.map_projection
    else:
        map_projection = cfg['wms-thumbnail-projection']
    if map_projection == 'Mercator':
        mapprojection = ccrs.Mercator()
    elif map_projection == 'PlateCarree':
        mapprojection = ccrs.PlateCarree()
    elif map_projection == 'PolarStereographic':
        mapprojection = ccrs.Stereographic(central_longitude=0.0, central_latitude=90.,
                                           true_scale_latitude=60.)
    else:
        raise Exception('Map projection is not properly specified in config')

    # Enable basic authentication if configured.
    if 'auth-basic-username' in cfg and 'auth-basic-password' in cfg:
        username = cfg['auth-basic-username']
        password = cfg['auth-basic-password']
        logger.info("Setting up basic authentication")
        if username == '' or password == '':
            raise Exception('Authentication username and/or password are configured,'
                            'but have blank strings')
        else:
            authentication = HTTPBasicAuth(username, password)
    else:
        authentication = None
        logger.info("Authentication disabled")
    # Get solr server config
    SolrServer = cfg['solrserver']
    myCore = cfg['solrcore']

    # Set up connection to SolR server
    mySolRc = SolrServer+myCore
    mysolr = IndexMMD(mySolRc, args.always_commit, authentication)

    # CONFIG DONE

    # Find files to process
    # FIXME remove l2 and thumbnail cores, reconsider deletion herein
    if args.input_file:
        myfiles = [args.input_file]
    elif args.list_file:
        try:
            f2 = open(args.list_file, "r")
        except IOError as e:
            logger.error('Could not open file: %s %e', args.list_file, e)
            return
        myfiles = f2.readlines()
        f2.close()
    elif args.remove:
        mysolr.delete_level1(args.remove)
        return
    elif args.remove and args.level2:
        mysolr.delete_level2(args.remove)
        return

    # This was not implemented
    # elif args.remove and args.thumbnail:
    #    mysolr.delete_thumbnail(deleteid)
    #    return

    elif args.directory:
        try:
            myfiles = os.listdir(args.directory)
        except Exception as e:
            logger.error("Something went wrong in decoding cmd arguments: %s", e)
            return 1

    fileno = 0
    myfiles_pending = []
    for myfile in myfiles:
        myfile = myfile.strip()
        # Decide files to operate on
        if not myfile.endswith('.xml'):
            continue
        if args.list_file:
            myfile = myfile.rstrip()
        if args.directory:
            myfile = os.path.join(args.directory, myfile)

        # FIXME, need a better way of handling this, WMS layers should be interpreted
        # automatically, this way we need to know up front whether WMS makes sense or not and
        # that won't work for harvesting
        if args.thumbnail_layer:
            wms_layer = args.thumbnail_layer
        else:
            wms_layer = None
        if args.thumbnail_style:
            wms_style = args.thumbnail_style
        else:
            wms_style = None
        if args.thumbnail_zoom_level:
            wms_zoom_level = args.thumbnail_zoom_level
        else:
            wms_zoom_level = 0
        if args.add_coastlines:
            wms_coastlines = args.add_coastlines
        else:
            wms_coastlines = True
        if args.thumbnail_extent:
            thumbnail_extent = [int(i) for i in args.thumbnail_extent[0].split(' ')]
        else:
            thumbnail_extent = None

        # Index files
        logger.info('Processing file: %d - %s', fileno, myfile)

        try:
            mydoc = MMD4SolR(myfile)
        except Exception as e:
            logger.warning('Could not handle file: %s', e)
            continue
        mydoc.check_mmd()
        fileno += 1

        """ Do not search for metadata_identifier, always used id...  """
        try:
            newdoc = mydoc.tosolr()
        except Exception as e:
            logger.warning('Could not process the file: %s', e)
            continue
        if (newdoc['metadata_status'] == "Inactive"):
            continue
        if (not args.no_thumbnail) and ('data_access_url_ogc_wms' in newdoc):
            tflg = True
        # Do not directly index children unless they are requested to be children. Do always
        # assume that the parent is included in the indexing process so postpone the actual
        # indexing to allow the parent to be properly indexed in SolR.
        if 'related_dataset' in newdoc:
            # Special fix for NPI
            newdoc['related_dataset'] = newdoc['related_dataset'].replace(
                'https://data.npolar.no/dataset/', '')
            newdoc['related_dataset'] = newdoc['related_dataset'].replace(
                'http://data.npolar.no/dataset/', '')
            newdoc['related_dataset'] = newdoc['related_dataset'].replace(
                'http://api.npolar.no/dataset/', '')
            newdoc['related_dataset'] = newdoc['related_dataset'].replace('.xml', '')
            # Skip if DOI is used to refer to parent, that isn't consistent.
            if 'doi.org' in newdoc['related_dataset']:
                continue
            # Fix special characters that SolR doesn't like

            myparent = newdoc['related_dataset']
            for e in IDREPLS:
                myparent = myparent.replace(e, '-')
            myresults = mysolr.solrc.search('id:' + myparent, **{'wt': 'python', 'rows': 100})
            if len(myresults) == 0:
                logger.warning("No parent found. Staging for second run.")
                myfiles_pending.append(myfile)
                continue
            elif not l2flg:
                logger.warning('Parent found, but assumes parent will be reindexed, thus'
                               'postponing indexing of children until SolR is updated.')
                myfiles_pending.append(myfile)
                continue
        logger.info("Indexing dataset: %s", myfile)
        if l2flg:
            mysolr.add_level2(mydoc.tosolr(), addThumbnail=tflg, projection=mapprojection,
                              wmstimeout=120, wms_layer=wms_layer, wms_style=wms_style,
                              wms_zoom_level=wms_zoom_level, add_coastlines=wms_coastlines,
                              wms_timeout=cfg['wms-timeout'], thumbnail_extent=thumbnail_extent)
        else:
            if tflg:
                try:
                    mysolr.index_record(input_record=mydoc.tosolr(), addThumbnail=tflg,
                                        wms_layer=wms_layer, wms_style=wms_style,
                                        wms_zoom_level=wms_zoom_level,
                                        add_coastlines=wms_coastlines, projection=mapprojection,
                                        wms_timeout=cfg['wms-timeout'],
                                        thumbnail_extent=thumbnail_extent)
                except Exception as e:
                    logger.warning('Something failed during indexing %s', e)
            else:
                try:
                    mysolr.index_record(input_record=mydoc.tosolr(), addThumbnail=tflg)
                except Exception as e:
                    logger.warning('Something failed during indexing %s', e)
        if not args.level2:
            l2flg = False
        tflg = False

    # Now process all the level 2 files that failed in the previous
    # sequence. If the Level 1 dataset is not available, this will fail at
    # level 2. Meaning, the section below only ingests at level 2.
    fileno = 0
    if len(myfiles_pending) > 0 and not args.always_commit:
        logger.info('Processing files that were not possible to process in first take.'
                    'Waiting 20 minutes to allow SolR to update recently indexed parent datasets.')
        sleep(20*60)
    for myfile in myfiles_pending:
        logger.info('Processing L2 file: %d - %s', fileno, myfile)
        try:
            mydoc = MMD4SolR(myfile)
        except Exception as e:
            logger.warning('Could not handle file: %s', e)
            continue
        mydoc.check_mmd()
        fileno += 1
        # Do not search for metadata_identifier, always use id
        # Check if this can be used????
        newdoc = mydoc.tosolr()

        if 'data_access_resource' in newdoc.keys():
            for e in newdoc['data_access_resource']:
                if (not nflg) and "OGC WMS" in (''.join(e)):
                    tflg = True
        # Skip file if not a level 2 file
        if 'related_dataset' not in newdoc:
            continue
        logger.info("Indexing dataset: %s", myfile)
        # Ingest at level 2
        mysolr.add_level2(mydoc.tosolr(), addThumbnail=tflg, projection=mapprojection,
                          wmstimeout=120, wms_layer=wms_layer, wms_style=wms_style,
                          wms_zoom_level=wms_zoom_level, add_coastlines=wms_coastlines,
                          wms_timeout=cfg['wms-timeout'], thumbnail_extent=thumbnail_extent)
        tflg = False

    # Report status
    logger.info("Number of files processed were: %d", len(myfiles))

    # add a commit to solr at end of run
    mysolr.commit()


if __name__ == "__main__":
    main()
