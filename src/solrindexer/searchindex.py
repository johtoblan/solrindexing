"""
PURPOSE:
    This searches SolR for specific records and optionally deletes them. It can also optionally
    create a list of identifiers to delete. Search is done in ID for now.

AUTHOR:
    Øystein Godøy, METNO/FOU, 2021-02-10

UPDATES:

NOTES:
    - under rewrite...

"""

import argparse
import pysolr
import yaml


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--cfg", dest="cfgfile",
                        help="Configuration file", required=True)
    parser.add_argument("-s", "--searchstringst", dest="string",
                        help="String to search for", required=True)
    parser.add_argument('-d', '--delete', action='store_true', help="Flag to delete records")
    parser.add_argument('-a', '--always_commit', action='store_true',
                        help="Flag to commit directly")

    args = parser.parse_args()

    if args.cfgfile is None or args.string is None:
        parser.print_help()
        parser.exit()

    return args


def parse_cfg(cfgfile):
    # Read config file
    print("Reading", cfgfile)
    with open(cfgfile, 'r') as ymlfile:
        cfgstr = yaml.full_load(ymlfile)

    return cfgstr


class IndexMMD:
    """ requires a list of dictionaries representing MMD as input """

    def __init__(self, mysolrserver, commit):
        """
        Connect to SolR core
        """
        try:
            self.solrc = pysolr.Solr(mysolrserver, always_commit=commit)
        except Exception as e:
            print("Something failed in SolR init", str(e))
        print("Connection established to: " + str(mysolrserver))

    def delete_item(self, datasetid, commit):
        """ Require ID as input """
        """ Rewrite to take full metadata record as input """
        print("Deleting ", datasetid, " from Level 1")
        try:
            self.solrc.delete(id=datasetid)
        except Exception as e:
            print("Something failed in SolR delete", str(e))

        print("Record successfully deleted from core")

    def search(self, myargs):
        """ Require Id as input """
        try:
            results = self.solrc.search(myargs.string, **{'wt': 'python', 'rows': 100000})
        except Exception as e:
            print("Something failed: ", str(e))

        return results


def main():

    #  Parse command line arguments
    args = parse_arguments()

    #  Parse configuration file
    cfg = parse_cfg(args.cfgfile)

    SolrServer = cfg['solrserver']
    myCore = cfg['solrcore']

    mySolRc = SolrServer+myCore

    # Search for records
    mysolr = IndexMMD(mySolRc, args.always_commit)
    myresults = mysolr.search(args)

    print('Found %d matches' % myresults.hits)
    print('Looping through matches:')
    i = 0
    for doc in myresults:
        print('\t', i, doc['id'])
        deleteid = doc['id']
        if args.delete:
            mysolr.delete_item(deleteid, commit=None)
        i += 1
    print('Found %d matches' % myresults.hits)

    return


if __name__ == "__main__":
    main()
