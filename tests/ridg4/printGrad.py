import sympy as sp
from sympy.printing.str import StrPrinter

class CustomStrPrinter(StrPrinter):
    """
    Examples of how to customize the StrPrinter for both a SymPy class and a
    user defined class subclassed from the SymPy Basic class.
    """

    def _print_Pow(self, expr):
        return "(" + str(expr.args[0]) + ")^(" + str(expr.args[1]) + ")"

sp.Basic.__str__ = lambda self: CustomStrPrinter().doprint(self)
x,y,z = sp.symbols("x y z ")

x0 = x*5.0
y0 = y*5.0
z0 = z*5.0
#f = sp.exp(-(x*x + 3*y*y))
#f = (y-(x*x*x - 2*x))*(y-(x*x*x - 2*x)) #-- gordon
#f = z0*z0 - (y0-(x0*x0*x0 - 2*x0))**2 #2.5 - 0.1 -- gordon--scale it up by 5 to shrink them togeather.
f = y*y*x + z*z #-- eberly
#f = sp.exp(y*y*x)*z--really strange
#f = z*z * sp.sin(sp.sin(x*x+y*y+z*z)) # -- hight strenght - concentric spheres -- fstrength=20 - circle from 4 to 5
#f = sp.sin(4*3.14*x)*sp.cos(4*3.14*y) disaster
#f = sp.sin(4*3.14*x)*sp.cos(4*3.14*y)*sp.cos(4*3.14*z)
#f = y*y*x*z*z--pretty cool but not a surface really # same if you ahd +z^2
#f = y*y*x*z #could be good for intersection with sphere
#f = y*y*x
grad = [sp.diff(f, x), sp.diff(f, y), sp.diff(f,z)]

print("vec3 g = [{0}, {1}, {2}];".format(*[str(x) for x in grad]))

hess = [[sp.diff(g, x), sp.diff(g, y), sp.diff(g, z)] for g in grad]
hessStr = [[str(x) for x in y] for y in hess]
hessStrFlat = [y for x in hessStr for y in x]
print("tensor[3,3] hess = [[{0}, {1}, {2}], [{3}, {4}, {5}], [{6}, {7}, {8}]];".format(*hessStrFlat))
