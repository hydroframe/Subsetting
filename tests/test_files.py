from pathlib import Path

test_inputs = Path(__file__).parent / 'test_inputs'
huc10190004_inputs = test_inputs / 'HUC10190004'

conus1_dem = test_inputs / 'CONUS1_Inputs' / 'CONUS2.0_RawDEM_CONUS1clip.tif'
conus1_dem_pfb = test_inputs / 'CONUS1_Inputs'/'CONUS2.0_RawDEM_CONUS1clip.pfb'
conus2_dem = test_inputs / 'CONUS2_Inputs'/'CONUS2.0_RawDEM.tif'
conus1_mask = test_inputs / 'CONUS1_Inputs'/'Domain_Blank_Mask.tif'
conus2_mask = test_inputs / 'CONUS2_Inputs'/'conus_1km_PFmask2.tif'
forcings_pfb = test_inputs / 'NLDAS.Temp.000001_to_000024.pfb'
forcings_sa = test_inputs / 'NLDAS.Temp.000001_to_000024.sa'
forcings_tif = test_inputs / 'NLDAS.Temp.000001_to_000024.tif'
conus1_latlon = test_inputs / 'CONUS1_Inputs'/'conus1_Grid_Centers_Short_Deg.format.sa'
conus1_landcover = test_inputs / 'CONUS1_Inputs'/'conus1_landcover.sa'
conus2_landcover = test_inputs / 'CONUS2_Inputs'/'1km_CONUS2_landcover_IGBP.tif'
conus2_subsurface = test_inputs / 'CONUS2_Inputs'/'3d-grid.v3.tif'
test_bbox_input = test_inputs / 'bbox_test_file.txt'
partial_shape_file = test_inputs / 'partial_shape.shp'

regression_truth_tif = test_inputs / 'test_truth.tif'

test_domain_manifest = test_inputs / 'test_domain_manifest.yaml'

test_domain_inputs = test_inputs / 'testdom_inputs'

huc10190004 = {'conus1_mask': huc10190004_inputs/'conus1'/'WBDHU8_conus1_mask.tif',
               'conus1_sol': huc10190004_inputs/'conus1'/'WBDHU8_conus1.pfsol',
               'conus1_vtk': huc10190004_inputs/'conus1'/'WBDHU8_conus1_ref.vtk',
               'conus2_mask': huc10190004_inputs/'conus2'/'WBDHU8_conus2_mask.tif',
               'conus2_sol': huc10190004_inputs/'conus2'/'WBDHU8_conus2.pfsol',
               'conus2_vtk': huc10190004_inputs/'conus2'/'WBDHU8_conus2_ref.vtk',
               'conus1_dem': huc10190004_inputs/'conus1'/'WBDHU8_conus1_dem.tif',
               'conus1_dem_box': huc10190004_inputs/'conus1'/'WBDHU8_conus1_dem_box.pfb',
               'conus1_dem_padded_box': huc10190004_inputs/'conus1'/'WBDHU8_conus1_dem_padded.pfb',
               'conus1_inset': huc10190004_inputs/'conus1/WBDHU8_conus1_mask_crop.tif',
               'conus1_bbox': [1040, 717, 85, 30],
               'conus2_dem': huc10190004_inputs/'conus2/WBDHU8_conus2_dem.tif',
               'conus2_inset': huc10190004_inputs/'conus2/WBDHU8_conus2_mask_crop.tif',
               'conus2_bbox': [1469, 1666, 82, 29],
               'shapefile': huc10190004_inputs/'WBDHU8.shp',
               'conus1_mask_-9999999': huc10190004_inputs/'conus1/WBDHU8_conus1_mask_-9999999.tif',
               'conus1_vegm': huc10190004_inputs/'conus1/WBDHU8_conus1_vegm.dat',
               'conus1_latlon': huc10190004_inputs/'conus1/WBDHU8_conus1_latlon.sa',
               'conus1_subsurface': huc10190004_inputs/'conus1/WBDHU8_conus1_subsurface.pfb',
               'conus1_slopex': huc10190004_inputs/'conus1/WBDHU8_conus1_slopex.pfb',
               'conus1_slopey': huc10190004_inputs/'conus1/WBDHU8_conus1_slopey.pfb'}
