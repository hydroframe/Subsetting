"""Represent a ParFlow domain
"""
import os
import logging
import yaml
import errno
import pfsubset.subset.data as data
import pfsubset.subset.utils.io as file_io_tools


class ParflowDomain:
    """A ParFlow domain described by a YAML file"""
    def __init__(self, name, local_path, manifest_file=None, version=1):
        """Information about the ParFlow domain and dataset we are working with

        Parameters
        ----------
        name : str
            the name of the Domain to create
        local_path : str
             path on system where ALL domain inputs live
        manifest_file : str, optional
            path on system to domain_manifest.yml file
        version : int, optional
            the version of the domain to load, as defined in the yaml file (default=1)
        Returns
        -------
        ParflowDomain
        """
        self.name = name
        self.version = version
        self.local_path = local_path
        self.manifest_file = manifest_file
        self.required_files = {}
        self.optional_files = {}
        if self.manifest_file is not None:
            self._read_manifest(self.required_files, self.optional_files)
        self.check_inputs_exist()
        self.check_destination()
        self.mask_tif = None
        self.mask_array = None

    def _load_domain_tif(self, domain_mask_key='DOMAIN_MASK'):
        """Load the domain raster
        Parameters
        ----------
        domain_mask_key : str, optional
             Key value to load from domain definition (Default value = 'DOMAIN_MASK')

        Returns
        -------
        None
        """
        tif_filename = os.path.join(self.local_path, self.required_files.get(domain_mask_key))
        self.mask_tif = file_io_tools.read_geotiff(tif_filename)
        self.mask_array = file_io_tools.read_file(tif_filename)

    def get_domain_mask(self, domain_mask_key='DOMAIN_MASK'):
        """get the domain full_dim_mask array

        Parameters
        ----------
        domain_mask_key : str, optional
            key in yaml file which defines the domain full_dim_mask (Default value = 'DOMAIN_MASK')

        Returns
        -------
        mask_array : ndarray
            numpy array for the domain full_dim_mask

        """
        if self.mask_array is None:
            self._load_domain_tif(domain_mask_key)
        return self.mask_array

    def get_domain_tif(self, domain_mask_key='DOMAIN_MASK'):
        """get the domain full_dim_mask tif as a gdal geotif object

        Parameters
        ----------
        domain_mask_key : str, optional
            key in yaml file which defines the domain full_dim_mask (Default value = 'DOMAIN_MASK')

        Returns
        -------
        mask_tif
            gdal object for the domain full_dim_mask

        """
        if self.mask_tif is None:
            self._load_domain_tif(domain_mask_key)
        return self.mask_tif

    def check_inputs_exist(self):
        """Look for each input file to see if it exists

        Returns
        -------
        None

        """
        # check for required files
        required_missing = self._identify_missing_inputs(self.required_files)
        if len(required_missing) > 0:
            logging.error(f'could not locate required model input file(s) {required_missing}')
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), required_missing)
        # check for optional files
        optional_missing = self._identify_missing_inputs(self.optional_files)
        if len(optional_missing) > 0:
            logging.warning(f'could not locate optional model input file(s) {optional_missing}')
        return None

    def _identify_missing_inputs(self, file_dict):
        """Identify any missing files from the file dictionary

        Parameters
        ----------
        file_dict : dict
            dictionary mapping the input file and filenames for the model

        Returns
        -------
        missing : list
            list of missing files

        """
        missing = []
        for name, file_name in file_dict.items():
            file_path = os.path.join(self.local_path, file_name)
            if not os.path.isfile(file_path):
                missing.append(file_path)
        return missing

    def check_destination(self):
        """make sure the local folder to store inputs exists

        Returns
        -------
        None
            None if folder exists
        Raises
        ------
        FileNotFoundError
            if a destination folder is not found
        """
        if not os.path.isdir(self.local_path):
            msg = self.local_path
            logging.exception(msg)
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), msg)
        return None

    def _read_manifest(self, required_file_dict, optional_file_dict):
        """read a manifest file in yaml format like so:
        <model>:
            <ver>:
                required_files:
                    <file mappings>
                optional_files:
                    <file mappings>

        Parameters
        ----------
        required_file_dict : dict
            domain file dict to update
        optional_file_dict : dict
            optional_files data dict to update

        Returns
        -------
        None

        Raises
        ------
        AttributeError
            if required_files key is missing

        """
        with open(self.manifest_file, 'r') as manifest_file:
            yml = yaml.safe_load_all(manifest_file)
            for doc in yml:
                model_dict = doc.get(self.name).get(self.version)
                try:
                    for input_name, file_name in model_dict.get('required_files').items():
                        required_file_dict.update({input_name: file_name})
                except AttributeError as err:
                    logging.error(f'required file definitions not found in manifest {self.manifest_file} for model '
                                  f'{self.name}, version {self.version}')
                    raise err
                try:
                    for input_name, file_name in model_dict.get('optional_files').items():
                        optional_file_dict.update({input_name: file_name})
                except AttributeError:
                    logging.warning(f'optional file definitions not found in manifest {self.manifest_file} for model '
                                    f'{self.name}, version {self.version}')

    def __repr__(self):
        return f"{self.__class__.__name__}(name:{self.name!r}, version:{self.version!r}, " \
               f"manifest:{self.manifest_file!r}, local_path:{self.local_path!r}, " \
               f"required_files:{self.required_files!r}, optional_files:{self.optional_files!r})"


class Conus(ParflowDomain):

    def __init__(self, local_path, manifest_file=None, version=1):
        """Information about the CONUS domain we are working with

        Parameters
        ----------
        local_path : str
            path on system where conus inputs live
        manifest_file : str, optional
            a file containing the keys and values for CONUS input required_files
        version : int, optional
            the conus version to create

        Returns
        -------
        Conus
        """
        if manifest_file is None:
            manifest_file = data.conus_manifest
        super().__init__('conus', local_path, manifest_file, version)
        # had to do this because conus1 full_dim_mask is all 0's
        if self.version == 1:
            self.mask_array = self.get_domain_mask() + 1
            self.mask_array = self.mask_array * 3
            self._patch_map = {1: 'land', 3: 'top', 5: 'sink', 6: 'bottom'}
        elif self.version == 2:
            self._patch_map = {1: 'ocean', 2: 'land', 3: 'top', 4: 'lakes', 5: 'sink', 6: 'bottom'}
        self._border_mask = None

    def get_border_mask(self):
        if self._border_mask is None:
            tif_filename = os.path.join(self.local_path, self.required_files.get('CELL_TYPES'))
            self._border_mask = file_io_tools.read_file(tif_filename)
            if self.version == 1:
                self._border_mask = self._border_mask + 1
        return self._border_mask

    def get_patch_name(self, patch_id):
        return self._patch_map.get(patch_id)

    def __repr__(self):
        return f"{self.__class__.__name__}(name:{self.name!r}, version:{self.version!r}, " \
               f"manifest:{self.manifest_file!r}, local_path:{self.local_path!r}, mask_array:{self.mask_array!r}), " \
               f"required_files:{self.required_files!r}, optional_files:{self.optional_files!r})"
