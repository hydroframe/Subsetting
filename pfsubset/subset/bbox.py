"""Bounding Box to clip data with

"""


class BBox:

    """Bounding box to identify data mask and optional padding"""

    def __init__(self, x_1=1, y_1=1, nx=0, ny=0, pad=(0, 0, 0, 0)):
        """A point (x,y) plus extents (nx, ny) and optional padding values
        _1 indexes are 1's based
        _0 are zero based
        human bbox creates bbox using (1,1) as lower left
        bbox translates to 0's based needed for data processing
        padding order specified as (top, right, bottom, left)

        Parameters
        ----------
        x_1: int, optional
           starting X location (default=1)
        y_1: int, optional
            starting Y location (default=1)
        nx : int, optional
            number of X cells to include (default=0)
        ny : int, optional
            number of Y cells to include (default=0)
        pad : tuple, optional
            number of no_data cells to add around data region.
            Specified clockwise, starting from top (top,right,bot,left) (CSS Style)
            default=(0,0,0,0)
        Returns
        -------
        BBox
        """
        self.pad_top = pad[0]
        self.pad_right = pad[1]
        self.pad_bottom = pad[2]
        self.pad_left = pad[3]
        self.x_1 = x_1
        self.y_1 = y_1
        self.nx = nx
        self.ny = ny
        self.x_0 = self.x_1 - 1
        self.y_0 = self.y_1 - 1

    def __repr__(self):
        return f"{self.__class__.__name__}(x_1:{self.x_1}, y_1:{self.y_1}, nx:{self.nx}, ny:{self.ny}, " \
               f"x_0:{self.x_0}, y_0:{self.y_0}, pad_top:{self.pad_top}, pad_right:{self.pad_right}, " \
                f"pad_bottom:{self.pad_bottom}, pad_left:{self.pad_left})"

    def get_inner_extents(self):
        """get the extent range of data without padding

        Returns
        -------
        int
            minimum Y value containing data
        int
            maximum Y value containing data
        int
            minimum X value containing data
        int
            maximum X value containing data
        """
        # return min_y, max_y, min_x, max_x without padding
        return self.y_0, self.y_0 + self.ny, self.x_0, self.x_0 + self.nx

    def get_padded_extents(self):
        """get the extent range of the data with padding

        Returns
        -------
        int
            minimum Y value including padding
        int
            maximum Y value including padding
        int
            minimum X value including padding
        int
            maximum X value including padding
        """

        x_0 = self.x_0 - self.pad_left
        x_end = x_0 + self.pad_left + self.pad_right + self.nx
        y_0 = self.y_0 - self.pad_bottom
        y_end = y_0 + self.pad_top + self.pad_bottom + self.ny
        return y_0, y_end, x_0, x_end

    def get_padding(self):
        """ get the padding values for the bbox

        Returns
        -------
        tuple of ints
            padding values in (top, right, bot, left) format
        """
        return tuple([self.pad_top, self.pad_right, self.pad_bottom, self.pad_left])

    def get_system_bbox(self):
        """get 0's based bbox values
        """
        return self.x_0, self.y_0, self.nx, self.ny

    def get_human_bbox(self):
        """get 1's based bbox values
        """
        return self.x_1, self.y_1, self.nx, self.ny

