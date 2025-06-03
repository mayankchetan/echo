# This script will calculate structural properties of a simple turbine
import numpy as np
import pandas as pd

STEEL_DENSITY = 7850 # kg/m^3
STEEL_YOUNGS_MODULUS = 210e9 # Pa
GFRP_DENSITY = 2000 # kg/m^3
GFRP_YOUNGS_MODULUS = 45e9 # Pa

class SimpleTurbine:
    """
    A simple turbine model for calculating structural properties. The tower and blade are modeled as cylinders
    The dimentions best try to match a real turbine, but are not exact.

    """

    def __init__(self, hub_height = 120, # m 
                        tower_height = 117, # m
                        tower_mass  = 256000, # kg
                        tower_OD = 4, # m
                        blade_length = 63, # m
                        blade_mass = 13700, # kg
                        balde_OD = 2, # m
                        rotor_radius = 64, # m,
                        overhang = 4, # m
                        towerNodes = 14, # #
                        bladeNodes = 8, # #
                        ):
        
        self.hub_height = hub_height
        self.tower_height = tower_height
        self.tower_mass = tower_mass
        self.tower_OD = tower_OD
        self.blade_length = blade_length
        self.blade_mass = blade_mass
        self.blade_OD = balde_OD
        self.rotor_radius = rotor_radius
        self.overhang = overhang
        self.towerNodes = towerNodes
        self.bladeNodes = bladeNodes

        # calc thicknesses
        self.tower_thickness = self._calc_thickness_annulus(self.tower_OD/2, self.tower_height, self.tower_mass, STEEL_DENSITY)
        self.blade_thickness = self._calc_thickness_annulus(self.blade_OD/2, self.blade_length, self.blade_mass, GFRP_DENSITY)


    # def buildTurbine(self):

    #     self.tower = pd.DataFrame(columns=[


    def _calc_thickness_annulus(self, rad, len, mass, density):
        """
        Calculate the thickness of an annulus given the radius, length, mass and density.
        
        Parameters:
        rad (float): Outer radius of the annulus (m)
        len (float): Length of the annulus (m)
        mass (float): Mass of the annulus (kg)
        density (float): Density of the material (kg/m^3)
        
        Returns:
        float: Thickness of the annulus (m)
        """        
        # Volume calculation from mass and density
        volume = mass / density
        # Solve the quadratic equation:
        # or: t² - 2*rad*t + volume/(π*length) = 0
        
        # Solve using quadratic formula with NumPy
        a = 1.0, b = -2.0 * rad, c = volume / (np.pi * len)
        
        # Calculate discriminant
        discriminant = b**2 - 4*a*c
        
        if discriminant < 0:
            raise ValueError("No real solution exists. Check your input values.")
        
        # Find the two possible solutions
        t = np.roots([a, b, c])
        
        # Select the physically meaningful solution (positive and less than radius)
        valid_solutions = t[(t > 0) & (t < rad)]
        
        if len(valid_solutions) == 0:
            raise ValueError("No valid thickness found. Check your input values.")
        
        # Return the smallest valid solution
        return np.min(valid_solutions)

    def _sectional_modulus(self, rad, thickness): # of an annulus
        # returns area, polar moment of inertia, and Ixx_Iyy
        area = np.pi * (rad**2 - (rad - thickness)**2)
        polar_moment_of_inertia = (np.pi / 2) * (rad**4 - (rad - thickness)**4)
        Ixx_Iyy = (np.pi / 4) * (rad**4 - (rad - thickness)**4)
        return area, polar_moment_of_inertia, Ixx_Iyy
