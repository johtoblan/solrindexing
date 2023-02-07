"""
SOLR-indexer : Main indexer
===========================

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

PURPOSE:
    This is designed to simplify the process of indexing single or multiple datasets.

AUTHOR:
    Øystein Godøy, METNO/FOU, 2017-11-09

UPDATES:
    Øystein Godøy, METNO/FOU, 2019-05-31
        Integrated modifications from Trygve Halsne and Massimo Di Stefano
    Øystein Godøy, METNO/FOU, 2018-04-19
        Added support for level 2
    Øystein Godøy, METNO/FOU, 2021-02-19
        Added argparse, fixing robustness issues.
    Johannes Langvatn, METNO/SUV, 2023-02-07
        Started refactoring

NOTES:
    - under rewrite...
"""

import os
import math
import base64
import pysolr
import netCDF4
import logging
import warnings
import xmltodict
import dateutil.parser
import lxml.etree as ET
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import shapely.geometry as shpgeo

from collections import OrderedDict
from owslib.wms import WebMapService
from shapely.geometry import box, mapping


logger = logging.getLogger(__name__)


def getZones(lon, lat):
    "get UTM zone number from latitude and longitude"
    if lat >= 72.0 and lat < 84.0:
        if lon >= 0.0 and lon < 9.0:
            return 31
        if lon >= 9.0 and lon < 21.0:
            return 33
        if lon >= 21.0 and lon < 33.0:
            return 35
        if lon >= 33.0 and lon < 42.0:
            return 37
    if lat >= 56 and lat < 64.0 and lon >= 3 and lon <= 12:
        return 32
    return math.floor((lon + 180) / 6) + 1


class MMD4SolR:
    """ Read and check MMD files, convert to dictionary """

    def __init__(self, filename):
        logger.info('Creating an instance of IndexMMD')
        """ set variables in class """
        self.filename = filename
        try:
            with open(self.filename, encoding='utf-8') as fd:
                self.mydoc = xmltodict.parse(fd.read())
        except Exception as e:
            logger.error('Could not open file: %s; %s', self.filename, e)
            raise

    def check_mmd(self):
        """ Check and correct MMD if needed """
        """ Remember to check that multiple fields of abstract and title
        have set xml:lang= attributes... """

        """
        Check for presence of required elements
        Temporal and spatial extent are not required as of no as it will
        break functionality for some datasets and communities especially
        in the Arctic context.
        """
        # TODO add proper docstring
        mmd_requirements = {
            'mmd:metadata_version': False,
            'mmd:metadata_identifier': False,
            'mmd:title': False,
            'mmd:abstract': False,
            'mmd:metadata_status': False,
            'mmd:dataset_production_status': False,
            'mmd:collection': False,
            'mmd:last_metadata_update': False,
            'mmd:iso_topic_category': False,
            'mmd:keywords': False,
        }
        """
        Check for presence and non empty elements
        This must be further developed...
        """
        for requirement in mmd_requirements.keys():
            if requirement in self.mydoc['mmd:mmd']:
                logger.info('Checking for: %s', requirement)
                if requirement in self.mydoc['mmd:mmd']:
                    if self.mydoc['mmd:mmd'][requirement] is not None:
                        logger.info('%s is present and non empty', requirement)
                        mmd_requirements[requirement] = True
                    else:
                        logger.warning('Required element %s is missing, setting it to unknown',
                                       requirement)
                        self.mydoc['mmd:mmd'][requirement] = 'Unknown'
                else:
                    logger.warning('Required element %s is missing, setting it to unknown.',
                                   requirement)
                    self.mydoc['mmd:mmd'][requirement] = 'Unknown'

        """
        Check for correct vocabularies where necessary
        Change to external files (SKOS), using embedded files for now
        Should be collected from
            https://github.com/steingod/scivocab/tree/master/metno
        """
        mmd_controlled_elements = {
            'mmd:iso_topic_category': ['farming',
                                       'biota',
                                       'boundaries',
                                       'climatologyMeteorologyAtmosphere',
                                       'economy',
                                       'elevation',
                                       'environment',
                                       'geoscientificInformation',
                                       'health',
                                       'imageryBaseMapsEarthCover',
                                       'inlandWaters',
                                       'location',
                                       'oceans',
                                       'planningCadastre',
                                       'society',
                                       'structure',
                                       'transportation',
                                       'utilitiesCommunication'],
            'mmd:collection': ['ACCESS',
                               'ADC',
                               'AeN',
                               'APPL',
                               'CC',
                               'DAM',
                               'DOKI',
                               'GCW',
                               'NBS',
                               'NMAP',
                               'NMDC',
                               'NSDN',
                               'SIOS',
                               'SESS_2018',
                               'SESS_2019',
                               'SIOS_access_programme',
                               'YOPP'],
            'mmd:dataset_production_status': ['Planned',
                                              'In Work',
                                              'Complete',
                                              'Obsolete'],
            'mmd:quality_control': ['No quality control',
                                    'Basic quality control',
                                    'Extended quality control',
                                    'Comprehensive quality control'],
        }
        for element in mmd_controlled_elements.keys():
            logger.info('Checking %s for compliance with controlled vocabulary', element)
            if element in self.mydoc['mmd:mmd']:

                if isinstance(self.mydoc['mmd:mmd'][element], list):
                    for elem in self.mydoc['mmd:mmd'][element]:
                        if isinstance(elem, dict):
                            myvalue = elem['#text']
                        else:
                            myvalue = elem
                        if myvalue not in mmd_controlled_elements[element]:
                            if myvalue is not None:
                                logger.warning('%s contains non valid content: %s',
                                               element, myvalue)
                            else:
                                logger.warning('Discovered an empty element.')
                else:
                    if isinstance(self.mydoc['mmd:mmd'][element], dict):
                        myvalue = self.mydoc['mmd:mmd'][element]['#text']
                    else:
                        myvalue = self.mydoc['mmd:mmd'][element]
                    if myvalue not in mmd_controlled_elements[element]:
                        logger.warning('%s contains non valid content: %s', element, myvalue)

        """
        Check that keywords also contain GCMD keywords
        Need to check contents more specifically...
        """
        gcmd = False
        if isinstance(self.mydoc['mmd:mmd']['mmd:keywords'], list):
            i = 0
            # TODO: remove unused for loop
            # Switch to using e instead of self.mydoc...
            for e in self.mydoc['mmd:mmd']['mmd:keywords']:
                if str(self.mydoc['mmd:mmd']['mmd:keywords'][i]['@vocabulary']).upper() == 'GCMDSK':
                    gcmd = True
                    break
                i += 1
            if not gcmd:
                logger.warning('Keywords in GCMD are not available (a)')
        else:
            if str(self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary']).upper() == 'GCMDSK':
                gcmd = True
            else:
                # warnings.warning('Keywords in GCMD are not available')
                logger.warning('Keywords in GCMD are not available (b)')

        """
        Modify dates if necessary
        Adapted for the new MMD specification, but not all information is
        extracted as SolR is not adapted.
        FIXME and check
        """
        if 'mmd:last_metadata_update' in self.mydoc['mmd:mmd']:
            if isinstance(self.mydoc['mmd:mmd']['mmd:last_metadata_update'],
                          dict):
                for mydict in self.mydoc['mmd:mmd']['mmd:last_metadata_update'].items():
                    if 'mmd:update' in mydict:
                        for myupdate in mydict:
                            if 'mmd:update' not in myupdate:
                                mydateels = myupdate
                                # The comparison below is a hack, need to
                                # revisit later, but works for now.
                                myvalue = '0000-00-00:T00:00:00Z'
                                if isinstance(mydateels, list):
                                    for mydaterec in mydateels:
                                        if mydaterec['mmd:datetime'] > myvalue:
                                            myvalue = mydaterec['mmd:datetime']
                                else:
                                    if mydateels['mmd:datetime'].endswith('Z'):
                                        myvalue = mydateels['mmd:datetime']
                                    else:
                                        myvalue = mydateels['mmd:datetime']+'Z'

            else:
                # To be removed when all records are transformed into the
                # new format
                logger.warning('Removed D7 format in last_metadata_update')
                if self.mydoc['mmd:mmd']['mmd:last_metadata_update'].endswith('Z'):
                    myvalue = self.mydoc['mmd:mmd']['mmd:last_metadata_update']
                else:
                    myvalue = self.mydoc['mmd:mmd']['mmd:last_metadata_update']+'Z'
            mydate = dateutil.parser.parse(myvalue)
        if 'mmd:temporal_extent' in self.mydoc['mmd:mmd']:
            if isinstance(self.mydoc['mmd:mmd']['mmd:temporal_extent'], list):

                i = 0
                for item in self.mydoc['mmd:mmd']['mmd:temporal_extent']:
                    for mykey in item:
                        if (item[mykey] is None) or (item[mykey] == '--'):
                            mydate = ''
                            self.mydoc['mmd:mmd']['mmd:temporal_extent'][i][mykey] = mydate
                        else:
                            mydate = dateutil.parser.parse(str(item[mykey]))
                            self.mydoc['mmd:mmd']['mmd:temporal_extent'][i][mykey] = mydate.strftime('%Y-%m-%dT%H:%M:%SZ')
                    i += 1
            else:
                for mykey in self.mydoc['mmd:mmd']['mmd:temporal_extent']:
                    if mykey == '@xmlns:gml':
                        continue
                    if (self.mydoc['mmd:mmd']['mmd:temporal_extent'][mykey] is None) or (self.mydoc['mmd:mmd']['mmd:temporal_extent'][mykey] == '--'):
                        mydate = ''
                        self.mydoc['mmd:mmd']['mmd:temporal_extent'][mykey] = mydate
                    else:
                        try:
                            mydate = dateutil.parser.parse(str(self.mydoc['mmd:mmd']['mmd:temporal_extent'][mykey]))
                            self.mydoc['mmd:mmd']['mmd:temporal_extent'][mykey] = mydate.strftime('%Y-%m-%dT%H:%M:%SZ')
                        except Exception as e:
                            logger.error('Date format could not be parsed: %s', e)

    def tosolr(self):
        """
        Method for creating document with SolR representation of MMD according
        to the XSD.
        """

        # Defining Look Up Tables
        personnel_role_LUT = {'Investigator':'investigator',
                              'Technical contact': 'technical',
                              'Metadata author': 'metadata_author',
                              'Data center contact':'datacenter'
                             }
        related_information_LUT = {'Dataset landing page':'landing_page',
                              'Users guide': 'user_guide',
                              'Project home page': 'home_page',
                              'Observation facility': 'obs_facility',
                              'Extended metadata':'ext_metadata',
                              'Scientific publication':'scientific_publication',
                              'Data paper':'data_paper',
                              'Data management plan':'data_management_plan',
                              'Other documentation':'other_documentation',
                             }

        # Create OrderedDict which will contain all elements for SolR
        mydict = OrderedDict()

        # SolR Can't use the mmd:metadata_identifier as identifier if it contains :, replace : and other characters like / etc by _ in the id field, let metadata_identifier be the correct one.

        """ Identifier """
        idrepls = [':','/','.']
        if isinstance(self.mydoc['mmd:mmd']['mmd:metadata_identifier'],dict):
            myid = self.mydoc['mmd:mmd']['mmd:metadata_identifier']['#text']
            for e in idrepls:
                myid = myid.replace(e,'-')
            mydict['id'] = myid
            mydict['metadata_identifier'] = self.mydoc['mmd:mmd']['mmd:metadata_identifier']['#text']
        else:
            myid = self.mydoc['mmd:mmd']['mmd:metadata_identifier']
            for e in idrepls:
                myid = myid.replace(e,'-')
            mydict['id'] = myid
            mydict['metadata_identifier'] = self.mydoc['mmd:mmd']['mmd:metadata_identifier']

        """ Last metadata update """
        if 'mmd:last_metadata_update' in self.mydoc['mmd:mmd']:
            last_metadata_update = self.mydoc['mmd:mmd']['mmd:last_metadata_update']

            lmu_datetime = []
            lmu_type = []
            lmu_note = []
            # FIXME check if this works correctly
            if isinstance(last_metadata_update['mmd:update'], dict): #Only one last_metadata_update element
                    lmu_datetime.append(str(last_metadata_update['mmd:update']['mmd:datetime']))
                    lmu_type.append(last_metadata_update['mmd:update']['mmd:type'])
                    lmu_note.append(last_metadata_update['mmd:update']['mmd:note'])

            else: # multiple last_metadata_update elements
                for i,e in enumerate(last_metadata_update['mmd:update']):
                    lmu_datetime.append(str(e['mmd:datetime']))
                    lmu_type.append(e['mmd:type'])
                    if 'mmd:note' in e.keys():
                        lmu_note.append(e['mmd:note'])

            i = 0
            for myel in lmu_datetime:
                i+=1
                if myel.endswith('Z'):
                    continue
                else:
                    lmu_datetime[i-1] = myel+'Z'
            mydict['last_metadata_update_datetime'] = lmu_datetime
            mydict['last_metadata_update_type'] = lmu_type
            mydict['last_metadata_update_note'] = lmu_note

        """ Metadata status """
        if isinstance(self.mydoc['mmd:mmd']['mmd:metadata_status'],dict):
            mydict['metadata_status'] = self.mydoc['mmd:mmd']['mmd:metadata_status']['#text']
        else:
            mydict['metadata_status'] = self.mydoc['mmd:mmd']['mmd:metadata_status']
        # TODO: the string below [title, abstract, etc ...]
        #  should be comments or some sort of logging statments

        """ Collection """
        if 'mmd:collection' in self.mydoc['mmd:mmd']:
            mydict['collection'] = []
            if isinstance(self.mydoc['mmd:mmd']['mmd:collection'], list):
                i = 0
                for e in self.mydoc['mmd:mmd']['mmd:collection']:
                    if isinstance(e,dict):
                        mydict['collection'].append(e['#text'])
                    else:
                        mydict['collection'].append(e)
                    i += 1
            else:
                mydict['collection'] = self.mydoc['mmd:mmd']['mmd:collection']

        """ title """
        if isinstance(self.mydoc['mmd:mmd']['mmd:title'], list):
            i = 0
            # Switch to using e instead of self.mydoc...
            for e in self.mydoc['mmd:mmd']['mmd:title']:
                if '@xml:lang' in e:
                    if e['@xml:lang'] == 'en':
                        mydict['title'] = e['#text']
                elif '@lang' in e:
                    if e['@lang'] == 'en':
                        mydict['title'] = e['#text']
        else:
            if isinstance(self.mydoc['mmd:mmd']['mmd:title'],dict):
                if '@xml:lang' in self.mydoc['mmd:mmd']['mmd:title']:
                    if self.mydoc['mmd:mmd']['mmd:title']['@xml:lang'] == 'en':
                        mydict['title'] = self.mydoc['mmd:mmd']['mmd:title']['#text']
                if '@lang' in self.mydoc['mmd:mmd']['mmd:title']:
                    if self.mydoc['mmd:mmd']['mmd:title']['@lang'] == 'en':
                        mydict['title'] = self.mydoc['mmd:mmd']['mmd:title']['#text']
            else:
                mydict['title'] = str(self.mydoc['mmd:mmd']['mmd:title'])

        """ abstract """
        if isinstance(self.mydoc['mmd:mmd']['mmd:abstract'], list):
            for e in self.mydoc['mmd:mmd']['mmd:abstract']:
                if '@xml:lang' in e:
                    if e['@xml:lang'] == 'en':
                        mydict['abstract'] = e['#text']
                elif '@lang' in e:
                    if e['@lang'] == 'en':
                        mydict['abstract'] = e['#text']
        else:
            if isinstance(self.mydoc['mmd:mmd']['mmd:abstract'],dict):
                if '@xml:lang' in self.mydoc['mmd:mmd']['mmd:abstract']:
                    if self.mydoc['mmd:mmd']['mmd:abstract']['@xml:lang'] == 'en':
                        mydict['abstract'] = self.mydoc['mmd:mmd']['mmd:abstract']['#text']
                if '@lang' in self.mydoc['mmd:mmd']['mmd:abstract']:
                    if self.mydoc['mmd:mmd']['mmd:abstract']['@lang'] == 'en':
                        mydict['abstract'] = self.mydoc['mmd:mmd']['mmd:abstract']['#text']
            else:
                mydict['abstract'] = str(self.mydoc['mmd:mmd']['mmd:abstract'])

        """ Temporal extent """
        if 'mmd:temporal_extent' in self.mydoc['mmd:mmd']:
            if isinstance(self.mydoc['mmd:mmd']['mmd:temporal_extent'], list):
                maxtime = dateutil.parser.parse('1000-01-01T00:00:00Z')
                mintime = dateutil.parser.parse('2099-01-01T00:00:00Z')
                for item in self.mydoc['mmd:mmd']['mmd:temporal_extent']:
                    for mykey in item:
                        if item[mykey] != '':
                            mytime = dateutil.parser.parse(item[mykey])
                        if mytime < mintime:
                            mintime = mytime
                        if mytime > maxtime:
                            maxtime = mytime
                mydict['temporal_extent_start_date'] = mintime.strftime('%Y-%m-%dT%H:%M:%SZ')
                mydict['temporal_extent_end_date'] = maxtime.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                mydict["temporal_extent_start_date"] = str(
                    self.mydoc['mmd:mmd']['mmd:temporal_extent']['mmd:start_date']),
                if 'mmd:end_date' in self.mydoc['mmd:mmd']['mmd:temporal_extent']:
                    if self.mydoc['mmd:mmd']['mmd:temporal_extent']['mmd:end_date']!=None:
                        mydict["temporal_extent_end_date"] = str(
                            self.mydoc['mmd:mmd']['mmd:temporal_extent']['mmd:end_date']),

        """ Geographical extent """
        """ Assumes longitudes positive eastwards and in the are -180:180
        """
        if 'mmd:geographic_extent' in self.mydoc['mmd:mmd'] and self.mydoc['mmd:mmd']['mmd:geographic_extent'] != None:
            if isinstance(self.mydoc['mmd:mmd']['mmd:geographic_extent'],
                    list):
                logger.warning('This is a challenge as multiple bounding boxes are not supported in MMD yet, flattening information')
                latvals = []
                lonvals = []
                for e in self.mydoc['mmd:mmd']['mmd:geographic_extent']:
                    if e['mmd:rectangle']['mmd:north'] != None:
                        latvals.append(float(e['mmd:rectangle']['mmd:north']))
                    if e['mmd:rectangle']['mmd:south'] != None:
                        latvals.append(float(e['mmd:rectangle']['mmd:south']))
                    if e['mmd:rectangle']['mmd:east'] != None:
                        lonvals.append(float(e['mmd:rectangle']['mmd:east']))
                    if e['mmd:rectangle']['mmd:west'] != None:
                        lonvals.append(float(e['mmd:rectangle']['mmd:west']))

                if len(latvals) > 0 and len(lonvals) > 0:
                    mydict['geographic_extent_rectangle_north'] = max(latvals)
                    mydict['geographic_extent_rectangle_south'] = min(latvals)
                    mydict['geographic_extent_rectangle_west'] = min(lonvals)
                    mydict['geographic_extent_rectangle_east'] = max(lonvals)
                    mydict['bbox'] = "ENVELOPE("+str(min(lonvals))+","+str(max(lonvals))+","+ str(max(latvals))+","+str(min(latvals))+")"
                    #Check if we have a point or a boundingbox
                    if float(e['mmd:rectangle']['mmd:north']) == float(e['mmd:rectangle']['mmd:south']):
                        if float(e['mmd:rectangle']['mmd:east']) == float(e['mmd:rectangle']['mmd:west']):
                            point = shpgeo.Point(float(e['mmd:rectangle']['mmd:east']),float(e['mmd:rectangle']['mmd:north']))
                            #print(point.y)
                            mydict['polygon_rpt'] = point.wkt

                            print(mapping(point))
                            #mydict['geom'] = geojson.dumps(mapping(point))
                    else:
                        bbox = box(min(lonvals), min(latvals), max(lonvals), max(latvals))

                        print("First conditition")
                        print(bbox)
                        polygon = bbox.wkt
                        #p = shapely.geometry.polygon.orient(polygon, sign=1.0)
                        #print(p.exterior.is_ccw)
                        mydict['polygon_rpt'] = polygon

                        #print(mapping(polygon))
                        #pGeo = shpgeo.shape({'type': 'polygon', 'coordinates': tuple(newCoord)})
                        #mydict['geom'] = geojson.dumps(mapping(shapely.wkt.loads(polygon)))
                        #print(mydict['geom'])


                else:
                    mydict['geographic_extent_rectangle_north'] = 90.
                    mydict['geographic_extent_rectangle_south'] = -90.
                    mydict['geographic_extent_rectangle_west'] = -180.
                    mydict['geographic_extent_rectangle_east'] = 180.
            else:
                for item in self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']:
                    #print(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle'][item])
                    if self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle'][item] == None:
                        logger.warning('Missing geographical element, will not process the file.')
                        mydict['metadata_status'] = 'Inactive'
                        raise Warning('Missing spatial bounds')

                mydict['geographic_extent_rectangle_north'] = float(
                    self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:north']),
                mydict['geographic_extent_rectangle_south'] = float(
                    self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:south']),
                mydict['geographic_extent_rectangle_east'] = float(
                    self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:east']),
                mydict['geographic_extent_rectangle_west'] = float(
                    self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:west']),
                if '@srsName' in self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle'].keys():
                    mydict['geographic_extent_rectangle_srsName'] = self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['@srsName'],
                mydict['bbox'] = "ENVELOPE("+self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:west']+","+self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:east']+","+ self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:north']+","+ self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:south']+")"

                print("Second conditition")
                #Check if we have a point or a boundingbox
                if float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:south']) == float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:north']):
                    if float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:east']) == float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:west']):
                        point = shpgeo.Point(float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:east']),float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:north']))
                        #print(point.y)
                        mydict['polygon_rpt'] = point.wkt

                        print(mapping(point))

                        #mydict['geom'] = geojson.dumps(mapping(point))

                else:
                    bbox = box(float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:west']), float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:south']), float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:east']), float(self.mydoc['mmd:mmd']['mmd:geographic_extent']['mmd:rectangle']['mmd:north']), ccw=False)
                    #print(bbox)
                    polygon = bbox.wkt
                    print(polygon)
                    #p = shapely.geometry.polygon.orient(shapely.wkt.loads(polygon), sign=1.0)
                    #print(p.exterior.is_ccw)
                    mydict['polygon_rpt'] = polygon
                    #print(mapping(shapely.wkt.loads(polygon)))
                    #print(geojson.dumps(mapping(loads(polygon))))
                    #pGeo = shpgeo.shape({'type': 'polygon', 'coordinates': tuple(newCoord)})
                    #mydict['geom'] = geojson.dumps(mapping(p))
                    #print(mydict['geom'])

        logger.info('Add location element?')

        """ Dataset production status """
        if 'mmd:dataset_production_status' in self.mydoc['mmd:mmd']:
            if isinstance(self.mydoc['mmd:mmd']['mmd:dataset_production_status'],
                    dict):
                mydict['dataset_production_status'] = self.mydoc['mmd:mmd']['mmd:dataset_production_status']['#text']
            else:
                mydict['dataset_production_status'] = str(self.mydoc['mmd:mmd']['mmd:dataset_production_status'])

        """ Dataset language """
        if 'mmd:dataset_language' in self.mydoc['mmd:mmd']:
            mydict['dataset_language'] = str(self.mydoc['mmd:mmd']['mmd:dataset_language'])

        """ Operational status """
        if 'mmd:operational_status' in self.mydoc['mmd:mmd']:
            mydict['operational_status'] = str(self.mydoc['mmd:mmd']['mmd:operational_status'])

        """ Access constraints """
        if 'mmd:access_constraint' in self.mydoc['mmd:mmd']:
            mydict['access_constraint'] = str(self.mydoc['mmd:mmd']['mmd:access_constraint'])

        """ Use constraint """
        if 'mmd:use_constraint' in self.mydoc['mmd:mmd'] and self.mydoc['mmd:mmd']['mmd:use_constraint'] != None:
            # Need both identifier and resource for use constraint
            if 'mmd:identifier' in self.mydoc['mmd:mmd']['mmd:use_constraint'] and 'mmd:resource' in self.mydoc['mmd:mmd']['mmd:use_constraint']:
                mydict['use_constraint_identifier'] = str(self.mydoc['mmd:mmd']['mmd:use_constraint']['mmd:identifier'])
                mydict['use_constraint_resource'] = str(self.mydoc['mmd:mmd']['mmd:use_constraint']['mmd:resource'])
            else:
                logger.warning('Both license identifier and resource need to be present to index this properly')
                mydict['use_constraint_identifier'] =  "Not provided"
                mydict['use_constraint_resource'] =  "Not provided"
            if 'mmd:license_text' in self.mydoc['mmd:mmd']['mmd:use_constraint']:
                mydict['use_constraint_license_text'] = str(self.mydoc['mmd:mmd']['mmd:use_constraint']['mmd:license_text'])

        """ Personnel """
        if 'mmd:personnel' in self.mydoc['mmd:mmd']:
            personnel_elements = self.mydoc['mmd:mmd']['mmd:personnel']

            if isinstance(personnel_elements, dict): #Only one element
                personnel_elements = [personnel_elements] # make it an iterable list

            # Facet elements
            mydict['personnel_role'] = []
            mydict['personnel_name'] = []
            mydict['personnel_organisation'] = []
            # Fix role based lists
            for role in personnel_role_LUT:
                mydict['personnel_{}_role'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_name'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_email'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_phone'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_fax'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_organisation'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_address'.format(personnel_role_LUT[role])] = []
                # don't think this is needed Øystein Godøy, METNO/FOU, 2021-09-08 mydict['personnel_{}_address_address'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_address_city'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_address_province_or_state'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_address_postal_code'.format(personnel_role_LUT[role])] = []
                mydict['personnel_{}_address_country'.format(personnel_role_LUT[role])] = []

            # Fill lists with information
            for personnel in personnel_elements:
                role = personnel['mmd:role']
                if not role:
                    logger.warning('No role available for personnel')
                    break
                if role not in personnel_role_LUT:
                    logger.warning('Wrong role provided for personnel')
                    break
                for entry in personnel:
                    entry_type = entry.split(':')[-1]
                    if entry_type == 'role':
                        mydict['personnel_{}_role'.format(personnel_role_LUT[role])].append(personnel[entry])
                        mydict['personnel_role'].append(personnel[entry])
                    else:
                        # Treat address specifically and handle faceting elements personnel_role, personnel_name, personnel_organisation.
                        if entry_type == 'contact_address':
                            for el in personnel[entry]:
                                el_type = el.split(':')[-1]
                                if el_type == 'address':
                                    mydict['personnel_{}_{}'.format(personnel_role_LUT[role], el_type)].append(personnel[entry][el])
                                else:
                                    mydict['personnel_{}_address_{}'.format(personnel_role_LUT[role], el_type)].append(personnel[entry][el])
                        elif entry_type == 'name':
                            mydict['personnel_{}_{}'.format(personnel_role_LUT[role], entry_type)].append(personnel[entry])
                            mydict['personnel_name'].append(personnel[entry])
                        elif entry_type == 'organisation':
                            mydict['personnel_{}_{}'.format(personnel_role_LUT[role], entry_type)].append(personnel[entry])
                            mydict['personnel_organisation'].append(personnel[entry])
                        else:
                            mydict['personnel_{}_{}'.format(personnel_role_LUT[role], entry_type)].append(personnel[entry])

        """ Data center """
        if 'mmd:data_center' in self.mydoc['mmd:mmd']:

            data_center_elements = self.mydoc['mmd:mmd']['mmd:data_center']

            if isinstance(data_center_elements, dict): #Only one element
                data_center_elements = [data_center_elements] # make it an iterable list

            for data_center in data_center_elements: #elf.mydoc['mmd:mmd']['mmd:data_center']: #iterate over all data_center elements
                for key,value in data_center.items():
                    if isinstance(value,dict): # if sub element is ordered dict
                        for kkey, vvalue in value.items():
                            element_name = 'data_center_{}'.format(kkey.split(':')[-1])
                            if not element_name in mydict.keys(): # create key in mydict
                                mydict[element_name] = []
                                mydict[element_name].append(vvalue)
                            else:
                                mydict[element_name].append(vvalue)
                    else: #sub element is not ordered dicts
                        element_name = '{}'.format(key.split(':')[-1])
                        if not element_name in mydict.keys(): # create key in mydict. Repetition of above. Should be simplified.
                            mydict[element_name] = []
                            mydict[element_name].append(value)
                        else:
                            mydict[element_name].append(value)

        """ Data access """
        # NOTE: This is identical to method above. Should in future versions be implified as a method
        if 'mmd:data_access' in self.mydoc['mmd:mmd']:
            data_access_elements = self.mydoc['mmd:mmd']['mmd:data_access']

            if isinstance(data_access_elements, dict): #Only one element
                data_access_elements = [data_access_elements] # make it an iterable list

            for data_access in data_access_elements: #iterate over all data_center elements
                data_access_type = data_access['mmd:type'].replace(" ","_").lower()
                mydict['data_access_url_{}'.format(data_access_type)] = data_access['mmd:resource']

                if 'mmd:wms_layers' in data_access and data_access_type == 'ogc_wms':
                    data_access_wms_layers_string = 'data_access_wms_layers'
                    data_access_wms_layers = data_access['mmd:wms_layers']
                    mydict[data_access_wms_layers_string] = [ i for i in data_access_wms_layers.values()][0]

        """ Related dataset """
        """ TODO """
        """ Remember to add type of relation in the future ØG """
        """ Only interpreting parent for now since SolR doesn't take more
            Added handling of namespace in identifiers
        """
        self.parent = None
        if 'mmd:related_dataset' in self.mydoc['mmd:mmd']:
            idrepls = [':','/','.']
            if isinstance(self.mydoc['mmd:mmd']['mmd:related_dataset'],
                    list):
                logger.warning('Too many fields in related_dataset...')
                for e in self.mydoc['mmd:mmd']['mmd:related_dataset']:
                    if '@mmd:relation_type' in e:
                        if e['@mmd:relation_type'] == 'parent':
                            if '#text' in dict(e):
                                mydict['related_dataset'] = e['#text']
                                for e in idrepls:
                                    mydict['related_dataset'] = mydict['related_dataset'].replace(e,'-')
            else:
                """ Not sure if this is used?? """
                if '#text' in dict(self.mydoc['mmd:mmd']['mmd:related_dataset']):
                    mydict['related_dataset'] = self.mydoc['mmd:mmd']['mmd:related_dataset']['#text']
                    for e in idrepls:
                        mydict['related_dataset'] = mydict['related_dataset'].replace(e,'-')

        """ Storage information """
        if 'mmd:storage_information' in self.mydoc['mmd:mmd'] and self.mydoc['mmd:mmd']['mmd:storage_information'] != None:
            if 'mmd:file_name' in self.mydoc['mmd:mmd']['mmd:storage_information'] and self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_name'] != None:
                mydict['storage_information_file_name'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_name'])
            if 'mmd:file_location' in self.mydoc['mmd:mmd']['mmd:storage_information'] and self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_location'] != None:
                mydict['storage_information_file_location'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_location'])
            if 'mmd:file_format' in self.mydoc['mmd:mmd']['mmd:storage_information'] and self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_format'] != None:
                mydict['storage_information_file_format'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_format'])
            if 'mmd:file_size' in self.mydoc['mmd:mmd']['mmd:storage_information'] and self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_size'] != None:
                if isinstance(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_size'], dict):
                    mydict['storage_information_file_size'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_size']['#text'])
                    mydict['storage_information_file_size_unit'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:file_size']['@unit'])
                else:
                    logger.warning("Filesize unit not specified, skipping field")
            if 'mmd:checksum' in self.mydoc['mmd:mmd']['mmd:storage_information'] and self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:checksum'] != None:
                if isinstance(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:checksum'], dict):
                    mydict['storage_information_file_checksum'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:checksum']['#text'])
                    mydict['storage_information_file_checksum_type'] = str(self.mydoc['mmd:mmd']['mmd:storage_information']['mmd:checksum']['@type'])
                else:
                    logger.warning("Checksum type is not specified, skipping field")


        """ Related information """
        if 'mmd:related_information' in self.mydoc['mmd:mmd']:

            related_information_elements = self.mydoc['mmd:mmd']['mmd:related_information']

            if isinstance(related_information_elements, dict): #Only one element
                related_information_elements = [related_information_elements] # make it an iterable list

            for related_information in related_information_elements:
                for key, value in related_information.items():
                    element_name = 'related_information_{}'.format(key.split(':')[-1])

                    if value in related_information_LUT.keys():
                        mydict['related_url_{}'.format(related_information_LUT[value])] = related_information['mmd:resource']
                        if 'mmd:description' in related_information:
                            mydict['related_url_{}_desc'.format(related_information_LUT[value])] = related_information['mmd:description']

        """
        ISO TopicCategory
        """
        if 'mmd:iso_topic_category' in self.mydoc['mmd:mmd']:
            mydict['iso_topic_category'] = []
            if isinstance(self.mydoc['mmd:mmd']['mmd:iso_topic_category'], list):
                for iso_topic_category in self.mydoc['mmd:mmd']['mmd:iso_topic_category']:
                    mydict['iso_topic_category'].append(iso_topic_category)
            else:
                mydict['iso_topic_category'].append(self.mydoc['mmd:mmd']['mmd:iso_topic_category'])

        """ Keywords """
        # Added double indexing of GCMD keywords. keywords_gcmd  (and keywords_wigos) are for faceting in SolR. What is shown in data portal is keywords_keyword.
        if 'mmd:keywords' in self.mydoc['mmd:mmd']:
            mydict['keywords_keyword'] = []
            mydict['keywords_vocabulary'] = []
            mydict['keywords_gcmd'] = []
            mydict['keywords_wigos'] = [] # Not used yet
            # If there is only one keyword list
            if isinstance(self.mydoc['mmd:mmd']['mmd:keywords'], dict):
                if isinstance(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'],str):
                    if self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary'] == "GCMDSK":
                        mydict['keywords_gcmd'].append(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'])
                    mydict['keywords_keyword'].append(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'])
                    mydict['keywords_vocabulary'].append(self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary'])
                else:
                    for i in range(len(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'])):
                        if isinstance(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'][i],str):
                            if self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary'] == "GCMDSK":
                                mydict['keywords_gcmd'].append(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'][i])
                            mydict['keywords_vocabulary'].append(self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary'])
                            mydict['keywords_keyword'].append(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'][i])
            # If there are multiple keyword lists
            elif isinstance(self.mydoc['mmd:mmd']['mmd:keywords'], list):
                for i in range(len(self.mydoc['mmd:mmd']['mmd:keywords'])):
                    if isinstance(self.mydoc['mmd:mmd']['mmd:keywords'][i],dict):
                        # Check for empty lists
                        if len(self.mydoc['mmd:mmd']['mmd:keywords'][i]) < 2:
                            continue
                        if isinstance(self.mydoc['mmd:mmd']['mmd:keywords'][i]['mmd:keyword'],list):
                            for j in range(len(self.mydoc['mmd:mmd']['mmd:keywords'][i]['mmd:keyword'])):
                                if self.mydoc['mmd:mmd']['mmd:keywords'][i]['@vocabulary'] == "GCMDSK":
                                    mydict['keywords_gcmd'].append(self.mydoc['mmd:mmd']['mmd:keywords'][i]['mmd:keyword'][j])
                                mydict['keywords_vocabulary'].append(self.mydoc['mmd:mmd']['mmd:keywords'][i]['@vocabulary'])
                                mydict['keywords_keyword'].append(self.mydoc['mmd:mmd']['mmd:keywords'][i]['mmd:keyword'][j])
                        else:
                            if self.mydoc['mmd:mmd']['mmd:keywords'][i]['@vocabulary'] == "GCMDSK":
                                mydict['keywords_gcmd'].append(self.mydoc['mmd:mmd']['mmd:keywords'][i]['mmd:keyword'])
                            mydict['keywords_vocabulary'].append(self.mydoc['mmd:mmd']['mmd:keywords'][i]['@vocabulary'])
                            mydict['keywords_keyword'].append(self.mydoc['mmd:mmd']['mmd:keywords'][i]['mmd:keyword'])

            else:
                if self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary'] == "GCMDSK":
                    mydict['keywords_gcmd'].append(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'])
                mydict['keywords_vocabulary'].append(self.mydoc['mmd:mmd']['mmd:keywords']['@vocabulary'])
                mydict['keywords_keyword'].append(self.mydoc['mmd:mmd']['mmd:keywords']['mmd:keyword'])

        """ Project """
        mydict['project_short_name'] = []
        mydict['project_long_name'] = []
        if 'mmd:project' in self.mydoc['mmd:mmd']:
            if self.mydoc['mmd:mmd']['mmd:project'] == None:
                mydict['project_short_name'].append('Not provided')
                mydict['project_long_name'].append('Not provided')
            elif isinstance(self.mydoc['mmd:mmd']['mmd:project'], list):
            # Check if multiple nodes are present
                for e in self.mydoc['mmd:mmd']['mmd:project']:
                    mydict['project_short_name'].append(e['mmd:short_name'])
                    mydict['project_long_name'].append(e['mmd:long_name'])
            else:
                # Extract information as appropriate
                e = self.mydoc['mmd:mmd']['mmd:project']
                if 'mmd:short_name' in e:
                    mydict['project_short_name'].append(e['mmd:short_name'])
                else:
                    mydict['project_short_name'].append('Not provided')

                if 'mmd:long_name' in e:
                    mydict['project_long_name'].append(e['mmd:long_name'])
                else:
                    mydict['project_long_name'].append('Not provided')


        """ Platform """
        # FIXME add check for empty sub elements...
        if 'mmd:platform' in self.mydoc['mmd:mmd']:

            platform_elements = self.mydoc['mmd:mmd']['mmd:platform']
            if isinstance(platform_elements, dict): #Only one element
                platform_elements = [platform_elements] # make it an iterable list

            for platform in platform_elements:
                for platform_key, platform_value in platform.items():
                    if isinstance(platform_value,dict): # if sub element is ordered dict
                        for kkey, vvalue in platform_value.items():
                            element_name = 'platform_{}_{}'.format(platform_key.split(':')[-1],kkey.split(':')[-1])
                            if not element_name in mydict.keys(): # create key in mydict
                                mydict[element_name] = []
                                mydict[element_name].append(vvalue)
                            else:
                                mydict[element_name].append(vvalue)
                    else: #sub element is not ordered dicts
                        element_name = 'platform_{}'.format(platform_key.split(':')[-1])
                        if not element_name in mydict.keys(): # create key in mydict. Repetition of above. Should be simplified.
                            mydict[element_name] = []
                            mydict[element_name].append(platform_value)
                        else:
                            mydict[element_name].append(platform_value)

                # Add platform_sentinel for NBS
                initial_platform = mydict['platform_long_name'][0]
                if initial_platform.startswith('Sentinel'):
                    mydict['platform_sentinel'] = initial_platform[:-1]

        """ Activity type """
        if 'mmd:activity_type' in self.mydoc['mmd:mmd']:
            mydict['activity_type'] = []
            if isinstance(self.mydoc['mmd:mmd']['mmd:activity_type'], list):
                for activity_type in self.mydoc['mmd:mmd']['mmd:activity_type']:
                    mydict['activity_type'].append(activity_type)
            else:
                mydict['activity_type'].append(self.mydoc['mmd:mmd']['mmd:activity_type'])

        """ Dataset citation """
        if 'mmd:dataset_citation' in self.mydoc['mmd:mmd']:
            dataset_citation_elements = self.mydoc['mmd:mmd']['mmd:dataset_citation']

            if isinstance(dataset_citation_elements, dict): #Only one element
                dataset_citation_elements = [dataset_citation_elements] # make it an iterable list

            for dataset_citation in dataset_citation_elements:
                for k, v in dataset_citation.items():
                    element_suffix = k.split(':')[-1]
                    # Fix issue between MMD and SolR schema, SolR requires full datetime, MMD not.
                    if element_suffix == 'publication_date':
                        if v is not None:
                            v+='T12:00:00Z'
                    mydict['dataset_citation_{}'.format(element_suffix)] = v

        """ Quality control """
        if 'mmd:quality_control' in self.mydoc['mmd:mmd'] and self.mydoc['mmd:mmd']['mmd:quality_control'] != None:
            mydict['quality_control'] = str(self.mydoc['mmd:mmd']['mmd:quality_control'])

        """ Adding MMD document as base64 string"""
        # Check if this can be simplified in the workflow.
        xml_root = ET.parse(str(self.filename))
        xml_string = ET.tostring(xml_root)
        encoded_xml_string = base64.b64encode(xml_string)
        xml_b64 = (encoded_xml_string).decode('utf-8')
        mydict['mmd_xml_file'] = xml_b64

        return mydict

class IndexMMD:
    """ Class for indexing SolR representation of MMD to SolR server. Requires
    a list of dictionaries representing MMD as input.
    """

    def __init__(self, mysolrserver, always_commit=False, authentication=None):
        # Set up logging
        logger.info('Creating an instance of IndexMMD')

        # level variables

        self.level = None

        # Thumbnail variables
        self.wms_layer = None
        self.wms_style = None
        self.wms_zoom_level = 0
        self.wms_timeout = None
        self.add_coastlines = None
        self.projection = None
        self.thumbnail_type = None
        self.thumbnail_extent = None

        # Connecting to core
        try:
            self.solrc = pysolr.Solr(mysolrserver, always_commit=always_commit, timeout=1020, auth=authentication)
            logger.info("Connection established to: %s", str(mysolrserver))
        except Exception as e:
            logger.error("Something failed in SolR init: %s", str(e))
            logger.info("Add a sys.exit?")

    #Function for sending explicit commit to solr
    def commit(self):
        self.solrc.commit()
    def index_record(self, input_record, addThumbnail, level=None, wms_layer=None, wms_style=None, wms_zoom_level=0, add_coastlines=True, projection=ccrs.PlateCarree(), wms_timeout=120, thumbnail_extent=None):
        """ Add thumbnail to SolR
            Args:
                input_record() : input MMD file to be indexed in SolR
                addThumbnail (bool): If thumbnail should be added or not
                level (int): 1 or 2 depending if MMD is Level-1 or Level-2,
                            respectively. If None, assume to be Level-1
                wms_layer (str): WMS layer name
                wms_style (str): WMS style name
                wms_zoom_level (float): Negative zoom. Fixed value added in
                                        all directions (E,W,N,S)
                add_coastlines (bool): If coastlines should be added
                projection (ccrs): Cartopy projection object or name (i.e. string)
                wms_timeout (int): timeout for WMS service
                thumbnail_extent (list): Spatial extent of the thumbnail in
                                      lat/lon [x0, x1, y0, y1]
            Returns:
                bool
        """
        if level == 1 or level == None:
            input_record.update({'dataset_type':'Level-1'})
            input_record.update({'isParent':'false'})
        elif level == 2:
            input_record.update({'dataset_type':'Level-2'})
        else:
            logger.error('Invalid level given: {}. Hence terminating'.format(level))

        if input_record['metadata_status'] == 'Inactive':
            logger.warning('Skipping record')
            return False
        myfeature = None
        if 'data_access_url_opendap' in input_record:
            # Thumbnail of timeseries to be added
            # Or better do this as part of get_feature_type?
            try:
                myfeature = self.get_feature_type(input_record['data_access_url_opendap'])
            except Exception as e:
                logger.error("Something failed while retrieving feature type: %s", str(e))
                #raise RuntimeError('Something failed while retrieving feature type')
            if myfeature:
                logger.info('feature_type found: %s', myfeature)
                input_record.update({'feature_type':myfeature})

        self.id = input_record['id']
        if 'data_access_url_ogc_wms' in input_record and addThumbnail == True:
            logger.info("Checking thumbnails...")
            getCapUrl = input_record['data_access_url_ogc_wms']
            if not myfeature:
                self.thumbnail_type = 'wms'
            self.wms_layer = wms_layer
            self.wms_style = wms_style
            self.wms_zoom_level = wms_zoom_level
            self.add_coastlines = add_coastlines
            self.projection = projection
            self.wms_timeout = wms_timeout
            self.thumbnail_extent = thumbnail_extent
            thumbnail_data = self.add_thumbnail(url=getCapUrl)

            if not thumbnail_data:
                logger.warning('Could not properly parse WMS GetCapabilities document')
                # If WMS is not available, remove this data_access element from the XML that is indexed
                del input_record['data_access_url_ogc_wms']
            else:
                input_record.update({'thumbnail_data':thumbnail_data})

        logger.info("Adding records to core...")

        mmd_record = list()
        mmd_record.append(input_record)

        try:
            self.solrc.add(mmd_record)
        except Exception as e:
            logger.error("Something failed in SolR adding document: %s", str(e))
            return False
        logger.info("Record successfully added.")

        return True

    def add_level2(self, myl2record, addThumbnail=False, projection=ccrs.Mercator(), wmstimeout=120, wms_layer=None, wms_style=None, wms_zoom_level=0, add_coastlines=True, wms_timeout=120, thumbnail_extent=None):
        """ Add a level 2 dataset, i.e. update level 1 as well """
        mmd_record2 = list()

        # Fix for NPI data...
        myl2record['related_dataset'] = myl2record['related_dataset'].replace('http://data.npolar.no/dataset/','')
        myl2record['related_dataset'] = myl2record['related_dataset'].replace('https://data.npolar.no/dataset/','')
        myl2record['related_dataset'] = myl2record['related_dataset'].replace('http://api.npolar.no/dataset/','')
        myl2record['related_dataset'] = myl2record['related_dataset'].replace('.xml','')

        # Add additonal helper fields for handling in SolR and Drupal
        myl2record['isChild'] = 'true'

        myfeature = None
        if 'data_access_url_opendap' in myl2record:
            # Thumbnail of timeseries to be added
            # Or better do this as part of get_feature_type?
            try:
                myfeature = self.get_feature_type(myl2record['data_access_url_opendap'])
            except Exception as e:
                logger.error("Something failed while retrieving feature type: %s", str(e))
                #raise RuntimeError('Something failed while retrieving feature type')
            if myfeature:
                logger.info('feature_type found: %s', myfeature)
                myl2record.update({'feature_type':myfeature})

        self.id = myl2record['id']
        # Add thumbnail for WMS supported datasets
        if 'data_access_url_ogc_wms' in myl2record and addThumbnail:
            logger.info("Checking tumbnails...")
            if not myfeature:
                self.thumbnail_type = 'wms'
            self.wms_layer = wms_layer
            self.wms_style = wms_style
            self.wms_zoom_level = wms_zoom_level
            self.add_coastlines = add_coastlines
            self.projection = projection
            self.wms_timeout = wms_timeout
            self.thumbnail_extent = thumbnail_extent
            if 'data_access_url_ogc_wms' in myl2record.keys():
                getCapUrl = myl2record['data_access_url_ogc_wms']
                try:
                    thumbnail_data = self.add_thumbnail(url=getCapUrl)
                except Exception as e:
                    logger.error("Something failed in adding thumbnail: %s", str(e))
                    warnings.warning("Couldn't add thumbnail.")

        if addThumbnail and thumbnail_data:
            myl2record.update({'thumbnail_data':thumbnail_data})

        mmd_record2.append(myl2record)

        """ Retrieve level 1 record """
        myparid = myl2record['related_dataset']
        idrepls = [':','/','.']
        for e in idrepls:
            myparid = myparid.replace(e,'-')
        try:
            myresults = self.solrc.search('id:' + myparid, **{'wt':'python','rows':100})
        except Exception as e:
            logger.error("Something failed in searching for parent dataset, " + str(e))

        # Check that only one record is returned
        if len(myresults) != 1:
            logger.warning("Didn't find unique parent record, skipping record")
            return
        # Convert from pySolr results object to dict and return.
        for result in myresults:
            if 'full_text' in result:
                result.pop('full_text')
            if 'bbox__maxX' in result:
                result.pop('bbox__maxX')
            if 'bbox__maxY' in result:
                result.pop('bbox__maxY')
            if 'bbox__minX' in result:
                result.pop('bbox__minX')
            if 'bbox__minY' in result:
                result.pop('bbox__minY')
            if 'bbox_rpt' in result:
                result.pop('bbox_rpt')
            if 'ss_access' in result:
                result.pop('ss_access')
            if '_version_' in result:
                result.pop('_version_')
                myresults = result
        myresults['isParent'] = 'true'

        # Check that the parent found has related_dataset set and
        # update this, but first check that it doesn't already exists
        if 'related_dataset' in myresults:
            # Need to check that this doesn't already exist...
            if myl2record['metadata_identifier'].replace(':','_') not in myresults['related_dataset']:
                myresults['related_dataset'].append(myl2record['metadata_identifier'].replace(':','_'))
        else:
            logger.info('This dataset was not found in parent, creating it...')
            myresults['related_dataset'] = []
            logger.info('Adding dataset with identifier %s to parent %s', myl2record['metadata_identifier'].replace(':','_'),myl2record['related_dataset'])
            myresults['related_dataset'].append(myl2record['metadata_identifier'].replace(':','_'))
        mmd_record1 = list()
        mmd_record1.append(myresults)

        ##print(myresults)

        """ Index level 2 dataset """
        try:
            self.solrc.add(mmd_record2)
        except Exception as e:
            raise Exception("Something failed in SolR add level 2", str(e))
        logger.info("Level 2 record successfully added.")

        """ Update level 1 record with id of this dataset """
        try:
            self.solrc.add(mmd_record1)
        except Exception as e:
            raise Exception("Something failed in SolR update level 1 for level 2", str(e))
        logger.info("Level 1 record successfully updated.")

    def add_thumbnail(self, url, thumbnail_type='wms'):
        """ Add thumbnail to SolR
            Args:
                type: Thumbnail type. (wms, ts)
            Returns:
                thumbnail: base64 string representation of image
        """
        print(url)
        if thumbnail_type == 'wms':
            try:
                thumbnail = self.create_wms_thumbnail(url)
                return thumbnail
            except Exception as e:
                logger.error("Thumbnail creation from OGC WMS failed: %s",e)
                return None
        elif thumbnail_type == 'ts': #time_series
            thumbnail = 'TMP'  # create_ts_thumbnail(...)
            return thumbnail
        else:
            logger.error('Invalid thumbnail type: {}'.format(thumbnail_type))
            return None


    def create_wms_thumbnail(self, url):
        """ Create a base64 encoded thumbnail by means of cartopy.

            Args:
                url: wms GetCapabilities document

            Returns:
                thumbnail_b64: base64 string representation of image
        """

        wms_layer = self.wms_layer
        wms_style = self.wms_style
        wms_zoom_level = self.wms_zoom_level
        wms_timeout = self.wms_timeout
        add_coastlines = self.add_coastlines
        map_projection = self.projection
        thumbnail_extent = self.thumbnail_extent

        # map projection string to ccrs projection
        if isinstance(map_projection,str):
            map_projection = getattr(ccrs,map_projection)()

        wms = WebMapService(url,timeout=wms_timeout)
        available_layers = list(wms.contents.keys())

        if wms_layer not in available_layers:
            wms_layer = available_layers[0]
            logger.info('Creating WMS thumbnail for layer: {}'.format(wms_layer))

        # Checking styles
        available_styles = list(wms.contents[wms_layer].styles.keys())

        if available_styles:
            if wms_style not in available_styles:
                wms_style = [available_styles[0]]
            else:
                wms_style = None
        else:
            wms_style = None

        if not thumbnail_extent:
            wms_extent = wms.contents[available_layers[0]].boundingBoxWGS84
            cartopy_extent = [wms_extent[0], wms_extent[2], wms_extent[1], wms_extent[3]]

            cartopy_extent_zoomed = [wms_extent[0] - wms_zoom_level,
                    wms_extent[2] + wms_zoom_level,
                    wms_extent[1] - wms_zoom_level,
                    wms_extent[3] + wms_zoom_level]
        else:
            cartopy_extent_zoomed = thumbnail_extent

        max_extent = [-180.0, 180.0, -90.0, 90.0]

        for i, extent in enumerate(cartopy_extent_zoomed):
            if i % 2 == 0:
                if extent < max_extent[i]:
                    cartopy_extent_zoomed[i] = max_extent[i]
            else:
                if extent > max_extent[i]:
                    cartopy_extent_zoomed[i] = max_extent[i]

        subplot_kw = dict(projection=map_projection)
        logger.info(subplot_kw)

        fig, ax = plt.subplots(subplot_kw=subplot_kw)

        #land_mask = cartopy.feature.NaturalEarthFeature(category='physical',
        #                                                scale='50m',
        #                                                facecolor='#cccccc',
        #                                                name='land')
        #ax.add_feature(land_mask, zorder=0, edgecolor='#aaaaaa',
        #        linewidth=0.5)

        # transparent background
        ax.spines['geo'].set_visible(False)
        #ax.outline_patch.set_visible(False)
        ##ax.background_patch.set_visible(False)
        fig.patch.set_alpha(0)
        fig.set_alpha(0)
        fig.set_figwidth(4.5)
        fig.set_figheight(4.5)
        fig.set_dpi(100)
        ##ax.background_patch.set_alpha(1)

        ax.add_wms(wms, wms_layer,
                wms_kwargs={'transparent': False,
                    'styles':wms_style})

        if add_coastlines:
            ax.coastlines(resolution="50m",linewidth=0.5)
        if map_projection == ccrs.PlateCarree():
            ax.set_extent(cartopy_extent_zoomed)
        else:
            ax.set_extent(cartopy_extent_zoomed, ccrs.PlateCarree())

        thumbnail_fname = 'thumbnail_{}.png'.format(self.id)
        fig.savefig(thumbnail_fname, format='png', bbox_inches='tight')
        plt.close('all')

        with open(thumbnail_fname, 'rb') as infile:
            data = infile.read()
            encode_string = base64.b64encode(data)

        thumbnail_b64 = (b'data:image/png;base64,' +
                encode_string).decode('utf-8')

        # Remove thumbnail
        os.remove(thumbnail_fname)
        return thumbnail_b64

    def create_ts_thumbnail(self):
        """ Create a base64 encoded thumbnail """

    def get_feature_type(self, myopendap):
        """ Set feature type from OPeNDAP """
        logger.info("Now in get_feature_type")

        # Open as OPeNDAP
        try:
            ds = netCDF4.Dataset(myopendap)
        except Exception as e:
            logger.error("Something failed reading dataset: %s", str(e))

        # Try to get the global attribute featureType
        try:
            featureType = ds.getncattr('featureType')
        except Exception as e:
            logger.error("Something failed extracting featureType: %s", str(e))
            raise
        ds.close()

        if featureType not in ['point', 'timeSeries', 'trajectory','profile','timeSeriesProfile','trajectoryProfile']:
            logger.warning("The featureType found - %s - is not valid", featureType)
            logger.warning("Fixing this locally")
            if featureType == "TimeSeries":
                featureType = 'timeSeries'
            elif featureType == "timeseries":
                featureType = 'timeSeries'
            elif featureType == "timseries":
                featureType = 'timeSeries'
            else:
                logger.warning("The featureType found is a new typo...")

            #raise

        return(featureType)

    def delete_level1(self, datasetid):
        """ Require ID as input """
        """ Rewrite to take full metadata record as input """
        logger.info("Deleting %s from level 1.", datasetid)
        try:
            self.solrc.delete(id=datasetid)
        except Exception as e:
            logger.error("Something failed in SolR delete: %s", str(e))
            raise

        logger.info("Record successfully deleted from Level 1 core")

    def delete_level2(self, datasetid):
        """ Require ID as input """
        """ Rewrite to take full metadata record as input """
        logger.info("Deleting %s from level 2.", datasetid)
        try:
            self.solr2.delete(id=datasetid)
        except Exception as e:
            logger.error("Something failed in SolR delete: %s", str(e))
            raise

        logger.info("Records successfully deleted from Level 2 core")

    def delete_thumbnail(self, datasetid):
        """ Require ID as input """
        """ Rewrite to take full metadata record as input """
        logger.info("Deleting %s from thumbnail core.", datasetid)
        try:
            self.solrt.delete(id=datasetid)
        except Exception as e:
            logger.error("Something failed in SolR delete: %s", str(e))
            raise

        logger.info("Records successfully deleted from thumbnail core")

    def search(self):
        """ Require Id as input """
        try:
            results = pysolr.search('mmd_title:Sea Ice Extent', df='text_en', rows=100)
        except Exception as e:
            logger.error("Something failed during search: %s", str(e))

        return results

    def darextract(self, mydar):
        mylinks = {}
        for i in range(len(mydar)):
            if isinstance(mydar[i], bytes):
                mystr = str(mydar[i], 'utf-8')
            else:
                mystr = mydar[i]
            if mystr.find('description') != -1:
                t1, t2 = mystr.split(',', 1)
            else:
                t1 = mystr
            t2 = t1.replace('"', '')
            proto, myurl = t2.split(':', 1)
            mylinks[proto] = myurl

        return (mylinks)
