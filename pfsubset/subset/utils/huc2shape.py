import requests
import itertools
import urllib.parse
import xml.etree.ElementTree as ET

try:
    from osgeo import ogr
except ImportError:
    import ogr

try:
    from osgeo import osr
except ImportError:
    import osr

from shapely import geometry
'''
date: 04.17.2020
author: Tony Castronova
affiliation: Consortium of Universties for the Advancement of Hydrologic
             Sciences Inc.
contact: acastronova@cuahsi.org
description: This script demonstrates how to query a WFS server using
             hydrologic unit codes (HUC) and build a GIS shapefile from
             the response. This uses a CUAHSI ArcServer that is only intended
             for testing purposes and shouldn't be used by any
             non-prototype application.
installation: conda install requests gdal shapely
usage: python3 build-shapefile-from-wfs.py 050600030504 050902010302
'''



# disable warnings
requests.packages.urllib3.disable_warnings()


def build_shape_object(hucs):

    watershed = WatershedBoundary()

    # loop through each HUC and query the geometry
    # from arcgis.cuahsi.org
    for huc in hucs:
        print(f'Querying geometry of HUC:{huc}')
        try:
            host_url = 'https://arcgis.cuahsi.org/arcgis/services/US_WBD/HUC_WBD/MapServer/WFSServer?'

            huc_filter = "<ogc:Filter>" \
                         "<ogc:PropertyIsEqualTo>" \
                         "<ogc:PropertyName>HUC12</ogc:PropertyName>" \
                         f"<ogc:Literal>{huc}</ogc:Literal>" \
                         "</ogc:PropertyIsEqualTo>" \
                         "</ogc:Filter>"

            defaultParameters = {'service': 'WFS',
                                 'request': 'GetFeature',
                                 'typeName': 'HUC_WBD:HUC12_US',
                                 'SrsName': 'EPSG:4326',
                                 'Filter': huc_filter}

            params = urllib.parse.urlencode(defaultParameters)
            request_url = host_url + params

            # make the http request
            response = requests.get(request_url, verify=False)

            # parse the response
            xmlbody = response.text
            tree = ET.ElementTree(ET.fromstring(xmlbody))
            root = tree.getroot()
            namespaces = {'gml': 'http://www.opengis.net/gml/3.2'}

            # add this geometry to the watershed object
            if watershed.srs is None:
                watershed.set_srs_from_gml_polygon(root.findall(
                    './/gml:Polygon', namespaces)[0])
            polygons_gml = root.findall('.//gml:Polygon//gml:posList',
                                        namespaces)
            watershed.add_boundary_from_gml(polygons_gml, huc)

        except Exception as e:
            print(f'Error: {e}')

    return watershed


class WatershedBoundary(object):
    """
    Helper class for organizing WFS response
    """

    def __init__(self, srs=None):
        self.srs = srs
        self.polygons = {}
        self.polygons_object = None

    def set_srs_from_gml_polygon(self, gml):

        longsrs = gml.attrib['srsName']
        longsrs = longsrs.split('::')
        self.srs = (longsrs[0].split(':')[-1], longsrs[-1])

    def get_polygons(self):
        return self.polygons

    def get_polygon_object(self):
        return self.polygons_object

    def get_polygon_wkt(self):
        wkt = []
        for hucid, poly_shape in self.polygons.items():
            wkt.append(poly_shape.wkt)
        return wkt

    def add_boundary_from_gml(self, gml, hucname):
        """
        Creates Shapely polygon objects from GML returned via ArcServer
        """

        polys = []
        for polygon in gml:
            vertex_list = polygon.text.split(' ')
            vertex_list = [float(v) for v in vertex_list]
            coords = list(zip(vertex_list[0::2], vertex_list[1::2]))
            polys.append(geometry.Polygon(coords))

        if len(polys) > 1:
            # create multipolygon
            plist = list(itertools.chain.from_iterable(self.polygons.values()))
            shape = geometry.MultiPolygon(plist)
        else:
            shape = polys[0]

        if hucname not in self.polygons.keys():
            self.polygons[hucname] = shape
        else:
            raise Exception(f'Huc already exists in WatershedBoundary: {hucname}')

    def write_shapefile(self, out_shp):
        """
        Creates a Shapefile using the watershed geometry
        """
        print(f'Creating Shapefile: {out_shp}')

        driver = ogr.GetDriverByName('Esri Shapefile')
        ds = driver.CreateDataSource(out_shp)
        ref = None
        if self.srs is not None and self.srs[0] == 'EPSG':
            ref = osr.SpatialReference()
            ref.ImportFromEPSG(int(self.srs[1]))

        layer = ds.CreateLayer('', ref, ogr.wkbPolygon)

        # Add one attribute
        layer.CreateField(ogr.FieldDefn('ID', ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn('HUC_ID', ogr.OFTString))
        defn = layer.GetLayerDefn()
        shapeid = 0
        for hucid, poly in self.polygons.items():
            shapeid += 1

            # Create a new feature (attribute and geometry)
            feat = ogr.Feature(defn)
            feat.SetField('ID', shapeid)
            feat.SetField('HUC_ID', hucid)

            geom = ogr.CreateGeometryFromWkb(poly.wkb)
            feat.SetGeometry(geom)

            layer.CreateFeature(feat)
            feat = geom = None  # destroy these

        # Save and close everything
        ds = layer = feat = geom = None
