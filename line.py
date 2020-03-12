

class Line:
    def __init__(self,x1,y1,x2,y2):
        super().__init__()
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.k = (x1-x2)/(y1-y2)


    def gety(self,x):
        return 1/self.k*(x-self.x2)+self.y2

    def getx(self,y):
        return self.k*(y-self.y2)+self.x2 



linea = Line(-282.463, -65.1037, 87.0396, 18.5869)

print(linea.gety(-1))
print(linea.gety(0))
print(linea.gety(1))