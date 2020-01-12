class Particle(object):
    __instance__ = None

    def __init__(self):
        Particle.__instance__ = self

        self.x = None
        self.y = None
        self.z = None


        self.vx = None
        self.vy = None
        self.vz = None

        self.m = None

    def isNone(self):
        if(self.x == None or self.y == None or self.z == None
                or self.vx == None or self.vy == None or self.vz == None
                or self.m == None):
            return True
        else:
            return False

        pass


