# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2012-2019 GEM Foundation
#
# OpenQuake is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenQuake is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with OpenQuake. If not, see <http://www.gnu.org/licenses/>.

"""
Module :mod:`openquake.hazardlib.scalerel.allen_hayes2017` implements
:class:AllenHayes2017Interface`
:class:`AllenHayes2017Intraslab`
"""

from numpy import log10
from openquake.hazardlib.scalerel.base import BaseMSRSigma, BaseASRSigma


class AllenHayes2017Interface(BaseMSRSigma, BaseASRSigma):
    """
    Allen & Hayes magnitude -- rupture area relationships for
    interface events.
    
    Allen, T. I., and G. P. Hayes (2017). Alternative rupture-scaling relationships 
    for subduction interface and other offshore environments, 
    Bull. Seism. Soc. Am. 107, doi: 10.1785/0120160255.

    Implements both magnitude-area and area-magnitude bi-linear scaling relationships.
    """

    def get_median_area(self, mag, rake):
        """
        Calculates median fault area from magnitude.
        """
        if mag <= 8.63:
            return 10 ** (-5.62 + 1.22 * mag)
        else:
            return 10 ** (2.23 + 0.31 * mag)

    def get_std_dev_area(self, mag, rake):
        """
        Returns std
        """
        return 0.256

    def get_median_mag(self, area, rake):
        """
        Returns magnitude for a given fault area
        
        :param area:
            Area in square km.
        """
        if area < 74000:
            return (log10(area)+5.62)/1.22
        else:
            return (log10(area)-2.23)/0.31

    def get_std_dev_mag(self, area, rake):
        """
        Returns std
        """
        return 0.267
        
class AllenHayes2017Intraslab(BaseMSRSigma, BaseASRSigma):
    """
    Allen & Hayes magnitude -- rupture area relationships for
    intraslab events.

    Allen, T. I., and G. P. Hayes (2017). Alternative rupture-scaling relationships 
    for subduction interface and other offshore environments, 
    Bull. Seism. Soc. Am. 107, doi: 10.1785/0120160255.

    Implements both magnitude-area and area-magnitude scaling relationships.
    """
    def get_median_area(self, mag, rake):
        """
        Calculates median areas as `10** (a + b*mag)`.
        The values are a function of magnitude. Rake is ignored.

        """
        return 10 ** (-3.89 + 0.96 * mag)
        
    def get_std_dev_area(self, mag, rake):
        """
        Standard deviation for Strasser et al 2010. Magnitude is ignored.
        """
        return 0.19

    def get_median_mag(self, area, rake):
        """
        Return magnitude (Mw) given the area Rake is ignored.

        :param area:
            Area in square km.
        """
        return (log10(area)+3.89)/0.96

    def get_std_dev_mag(self, rake):
        """
        Standard deviation on the magnitude for the Allen & Hayes (2017)
        area relation.
        """
        return 0.19

class AllenHayes2017StrikeSlip(BaseMSRSigma, BaseASRSigma):
    """
    Allen & Hayes magnitude -- rupture area relationships for
    offshore strike-slip events.

    Allen, T. I., and G. P. Hayes (2017). Alternative rupture-scaling relationships 
    for subduction interface and other offshore environments, 
    Bull. Seism. Soc. Am. 107, doi: 10.1785/0120160255.

    Implements both magnitude-area and area-magnitude scaling relationships.
    """
    def get_median_area(self, mag, rake):
        """
        Calculates median areas as `10** (a + b*mag)`.
        The values are a function of magnitude. Rake is ignored.

        """
        return 10 ** (-4.04 + 0.96 * mag)
        
    def get_std_dev_area(self, mag, rake):
        """
        Standard deviation for Strasser et al 2010. Magnitude is ignored.
        """
        return 0.20

    def get_median_mag(self, area, rake):
        """
        Return magnitude (Mw) given the area Rake is ignored.

        :param area:
            Area in square km.
        """
        return (log10(area)+4.04)/0.96

    def get_std_dev_mag(self, rake):
        """
        Standard deviation on the magnitude for the Allen & Hayes (2017)
        area relation.
        """
        return 0.20

class AllenHayes2017OuterRise(BaseMSRSigma, BaseASRSigma):
    """
    Allen & Hayes magnitude -- rupture area relationships for
    outer-rise events.

    Allen, T. I., and G. P. Hayes (2017). Alternative rupture-scaling relationships 
    for subduction interface and other offshore environments, 
    Bull. Seism. Soc. Am. 107, doi: 10.1785/0120160255.

    Implements both magnitude-area and area-magnitude scaling relationships.
    """
    def get_median_area(self, mag, rake):
        """
        Calculates median areas as `10** (a + b*mag)`.
        The values are a function of magnitude. Rake is ignored.

        """
        return 10 ** (-3.89 + 0.96 * mag)
        
    def get_std_dev_area(self, mag, rake):
        """
        Standard deviation for Strasser et al 2010. Magnitude is ignored.
        """
        return 0.11

    def get_median_mag(self, area, rake):
        """
        Return magnitude (Mw) given the area Rake is ignored.

        :param area:
            Area in square km.
        """
        return (log10(area)+3.89)/0.96

    def get_std_dev_mag(self, rake):
        """
        Standard deviation on the magnitude for the Allen & Hayes (2017)
        area relation.
        """
        return 0.11