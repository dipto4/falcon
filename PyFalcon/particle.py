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



