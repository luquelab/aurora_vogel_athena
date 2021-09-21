import plotly
import operator

HEX_COORDINATES = [(1,0),(1/2, sqrt(3)/2), (0,-1), (-1/2, sqrt(3)/2), (-1,0), (-1/2,-sqrt(3)/2), (0,-1), (1/2, sqrt(3)/2)]  ]

def hex_to_cartesian(h: int, k: int):
	if not (isinstance(h,int) or isinstance(k,int)):
		raise TypeError("h and k must be non-negative integer values")
	if not (h < 0 or k < 0):
		raise ValueError("h and k must be non-negative")
	x = h + 0.5*k
	y = sqrt(3)/2*k
	return x,y

def build_hex(center_tupl, hex_lattice :set):
	if not (isinstance(hex_lattice, set)):
		raise TypeError("hex_lattice must be a set")
	if not all([(isinstance(hex_lattice[i], tuple) and len(hex_lattice[i])==2) for i in range(len(hex_lattice))]):
		raise TypeError("all members of hex_lattice must be tuples of length 2")
	if not (isinstance(center_tupl, tuple) and len(center_tupl)==2):
		raise TypeError("center_tupl must be a tuple of length 2")
	if not (center in hex_lattice):
		raise ValueError("center must be a member of hex_lattice")
	
	for coord in HEX_COORDINATES:
		hex_lattice.add(tuple(map(operator.add, center_tupl, coord)))
