import math

class Sphere3d:
    def __init__(self, x, y, z, radius):
        self.x = x
        self.y = y
        self.z = z
        self.radius = radius

    def volume(self):
        return (4/3) * math.pi * self.radius**3
    
    def surface_area(self):
        return 4 * math.pi * self.radius**2
    
    def contains_points(self, x, y, z, atom_radius):
        distance = math.sqrt((x - self.x)**2 + (y - self.y)**2 + (z - self.z)**2)
        return distance <= (self.radius + atom_radius)