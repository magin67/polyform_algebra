import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.objects.point import Point

a, b, c = Point.create_list(['a', 'b', 'c'])
print(a*b)