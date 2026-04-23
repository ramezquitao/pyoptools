# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     SimpleDMD surface definition module
# Symbols Defined: SimpleDMD
# ------------------------------------------------------------------------------

'''Module that defines a simplified DMD (Digital Micromirror Device) surface class

SimpleDMD models a DMD as a flat reflective surface with a state-dependent tilted
normal vector. This is suitable for system-level optical design where pixel-level
effects are negligible.
'''

from libc.math cimport sin, cos

from pyoptools.raytrace.surface.plane cimport Plane
from pyoptools.misc.cmisc.eigen cimport Vector3d


cdef class SimpleDMD(Plane):
    '''Simplified DMD surface with state-dependent tilted normal.

    SimpleDMD represents a Digital Micromirror Device as a geometrically flat
    surface (z=0 plane) with a tilted normal vector that depends on the DMD state.
    This provides a computationally efficient model for system-level optical design
    without modeling individual micromirrors.

    The surface supports three states:
    - "on": Normal tilted by tilt_angle in on_direction_angle
    - "off": Normal tilted by tilt_angle in off_direction_angle
    - "flat": Normal perpendicular to surface (0, 0, 1)

    Parameters
    ----------
    tilt_angle : float
        Tilt angle in radians. This is the angle between the surface normal and
        the perpendicular direction (Z-axis). When tilt_angle=0, the normal is
        perpendicular to the surface. Use math.radians() to convert from degrees.
    on_direction_angle : float
        Direction angle in radians defining where the surface normal is rotated
        to for the ON state. Measured counter-clockwise from the +X axis in the
        XY plane when viewed from +Z. For example, 0° points toward +X, 90° points
        toward +Y, 180° points toward -X, and 270° points toward -Y.
    off_direction_angle : float
        Direction angle in radians defining where the surface normal is rotated
        to for the OFF state. Measured counter-clockwise from the +X axis in the
        XY plane when viewed from +Z. For example, 0° points toward +X, 90° points
        toward +Y, 180° points toward -X, and 270° points toward -Y.
    state : str, optional
        Initial state: "on", "off", or "flat". Default is "flat".
    shape : Shape
        Aperture shape (e.g., Rectangular, Circular) defining the active area.
    **kwargs
        Additional arguments passed to Surface base class (except reflectivity,
        which is hardcoded to 1.0).

    Attributes
    ----------
    state : str
        Current DMD state ("on", "off", or "flat"). Can be changed at runtime.

    Notes
    -----
    - All angles are in radians (consistent with pyoptools conventions)
    - The surface is always fully reflective (reflectivity = 1.0)
    - Ray-surface intersection is geometrically flat (z=0 plane)
    - Reflection uses the state-dependent tilted normal vector
    - Normal vectors are pre-computed for performance and automatically normalized

    The normal vector calculation uses spherical coordinates:
        normal_x = sin(tilt_angle) * cos(direction_angle)
        normal_y = sin(tilt_angle) * sin(direction_angle)
        normal_z = cos(tilt_angle)

    This is equivalent to rotating (0,0,1) using:
        compute_rotation_matrix(Rx=0, Ry=tilt_angle, Rz=direction_angle)

    Limitations
    -----------
    - This is a simplified model suitable for system-level design
    - Does not model individual pixels, diffraction, or fill factor effects
    - Surface is geometrically flat; only the optical normal is tilted
    - For pixel-level accuracy, use PixelatedDMD (future implementation)

    Examples
    --------
    >>> from math import radians
    >>> from pyoptools.raytrace.surface import SimpleDMD
    >>> from pyoptools.raytrace.shape import Rectangular
    >>>
    >>> # Create a DMD with 12° tilt magnitude
    >>> # ON state: normal rotates 12° toward 45° direction (northeast)
    >>> # OFF state: normal rotates 12° toward 225° direction (southwest)
    >>> dmd = SimpleDMD(
    ...     tilt_angle=radians(12),           # 12° tilt magnitude
    ...     on_direction_angle=radians(45),   # Normal tilts toward northeast
    ...     off_direction_angle=radians(225), # Normal tilts toward southwest
    ...     state="on",
    ...     shape=Rectangular(size=(10.8, 10.8))  # Mirror size in mm
    ... )
    >>>
    >>> # Switch state at runtime
    >>> dmd.state = "off"
    >>> dmd.state = "flat"
    '''

    # C-level attributes for pre-computed normal vectors
    cdef Vector3d _normal_on
    cdef Vector3d _normal_off
    cdef Vector3d _normal_flat
    cdef str _state
    cdef double _tilt_angle
    cdef double _on_direction_angle
    cdef double _off_direction_angle

    def __init__(self,
                 double tilt_angle,
                 double on_direction_angle,
                 double off_direction_angle,
                 str state="flat",
                 shape=None,
                 **kwargs):
        '''Initialize SimpleDMD surface.

        Parameters defined in class docstring.
        '''
        # Store parameters for __repr__
        self._tilt_angle = tilt_angle
        self._on_direction_angle = on_direction_angle
        self._off_direction_angle = off_direction_angle

        # Validate state parameter
        if state not in ("on", "off", "flat"):
            raise ValueError(
                f"Invalid state '{state}'. Must be 'on', 'off', or 'flat'."
            )

        # Remove reflectivity from kwargs if present (we hardcode it to 1.0)
        kwargs.pop('reflectivity', None)

        # Initialize parent Plane with reflectivity=1.0 (always reflective)
        Plane.__init__(self, shape=shape, reflectivity=1.0, **kwargs)

        # Pre-compute normal vectors for each state
        # Using direct trigonometric formula (spherical coordinates)
        # This is equivalent to:
        # compute_rotation_matrix(Rx=0, Ry=tilt_angle, Rz=direction_angle) @ (0,0,1)

        # ON state normal - use pointer syntax like in Plane
        (<double*>(&self._normal_on(0)))[0] = sin(tilt_angle) * cos(on_direction_angle)
        (<double*>(&self._normal_on(1)))[0] = sin(tilt_angle) * sin(on_direction_angle)
        (<double*>(&self._normal_on(2)))[0] = cos(tilt_angle)
        # Normalize to ensure unit length
        self._normal_on.normalize()

        # OFF state normal
        (<double*>(&self._normal_off(0)))[0] = \
            sin(tilt_angle) * cos(off_direction_angle)
        (<double*>(&self._normal_off(1)))[0] = \
            sin(tilt_angle) * sin(off_direction_angle)
        (<double*>(&self._normal_off(2)))[0] = cos(tilt_angle)
        # Normalize to ensure unit length
        self._normal_off.normalize()

        # FLAT state normal (perpendicular to surface)
        (<double*>(&self._normal_flat(0)))[0] = 0.0
        (<double*>(&self._normal_flat(1)))[0] = 0.0
        (<double*>(&self._normal_flat(2)))[0] = 1.0
        # Already unit length, but normalize for consistency
        self._normal_flat.normalize()

        # Set initial state
        self._state = state

    @property
    def state(self):
        '''Get current DMD state.

        Returns
        -------
        str
            Current state: "on", "off", or "flat"
        '''
        return self._state

    @state.setter
    def state(self, str value):
        '''Set DMD state.

        Parameters
        ----------
        value : str
            New state: "on", "off", or "flat"

        Raises
        ------
        ValueError
            If value is not a valid state.
        '''
        if value not in ("on", "off", "flat"):
            raise ValueError(
                f"Invalid state '{value}'. Must be 'on', 'off', or 'flat'."
            )
        self._state = value

    cdef void _calculate_normal(self,
                                Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        '''Calculate state-dependent surface normal.

        Returns the pre-computed normal vector based on current state.

        Parameters
        ----------
        intersection_point : Vector3d
            Point on surface (not used, included for interface compatibility)
        normal : Vector3d
            Output parameter for the normal vector
        '''
        # Return appropriate pre-computed normal based on state
        # Note: String comparison in nogil context - use simple if-else chain
        if self._state == "on":
            normal = self._normal_on
        elif self._state == "off":
            normal = self._normal_off
        else:  # "flat"
            normal = self._normal_flat

    def __repr__(self):
        '''Return string representation of SimpleDMD surface.

        Returns
        -------
        str
            String representation showing shape and current state
        '''
        return (f"SimpleDMD(shape={self.shape}, state='{self._state}', "
                f"reflectivity={self.reflectivity}, "
                f"tilt_angle={self._tilt_angle:.5f}, "
                f"on_direction_angle={self._on_direction_angle:.5f}, "
                f"off_direction_angle={self._off_direction_angle:.5f})")
