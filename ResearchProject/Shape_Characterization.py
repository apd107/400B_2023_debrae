

class Shape_Char:
    def __init__(self, x_radii, y_radii, z_radii):
        '''
        This function uses the major, minor, and intermediate
        axes to calculate the shape of the halo. 

        INPUT
        -----
        x_radii: 'float'
            The x radius of the ellipse.

        y_radii: 'float'
            The y radius of the ellipse.
        
        z_radii: 'float'
            The z radius of the ellipse.
        
        OUTPUT
        ------
        T: 'float'
            The float value which defines the shape of the halo 
            according to the derivations from Garavito-Camargo 2021
        '''
        # Calculating Shape! 
        # c = minor 
        # a = major 
        # b = intermediate
        self.a = x_radii
        self.b = y_radii
        self.c = z_radii


    def T_calc(self):
        # prolate - T >= 0.67
        # triaxial - 0.67 >= T >= 0.33
        # oblate - T <= 0.33
        s = self.c/self.a
        q = self.b/self.a
        T = (1-q**2)/(1-s**2)

        if T >= 0.67:
            print('Prolate')
        elif T <= 0.67 and T >= 0.33:
            print('Triaxial')
        elif T <= 0.33:
            print('Oblate')
        return T


