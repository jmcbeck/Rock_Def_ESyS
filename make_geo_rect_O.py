# --- geometry setup script for block with clustered particles ---
from gengeo import *
from esys.lsm.geometry import *
import math

# make layers using boxes with lines

# - input parameters --
# block dimensions
# mm
xdim = 20.0
ydim = 40.0
zdim = 20.0

xmin = 0.0
ymin = 0.0
zmin = 0.0

# thickness of layers
thickness = 5.0

# orientation of layers w.r.t ydim (max compressive stress)
#theta = 90 #0, 15, 30, 45, 60, 75, 90
arg_list = sys.argv
theta = int(arg_list[1])
print theta

# amount of overlap between layers
# perpendicular overlap distance between layers
overlap = 0.5 #0.2, 0.5, 1.0 # do these produce different roughnesses?  

# for different particle packing
geoIter = 1 # if >1, then seed bool = True if want different packing

geo_file = "rectO_G"+str(geoIter)+"_L"+str(thickness) +"_T"+str(theta)+"_O"+str(overlap)
print geo_file
#exit()

# horizontal layers
if theta == 0:
    thick_y = thickness
# vertical layers
elif theta == 90:
    thick_y = 0
# layers oriented at some intermediate angle 0-90 degrees
else:
    # with horizontal layers, number of layers = height (y) / thickness
    thick_y = abs(thickness/math.sin(math.radians(90-theta)))
    nlayers = ydim/thick_y

# particle size range
# wider range of particle radius leads to low porosity (~25%)
minRadius = 0.1 
maxRadius = 1.0

# ---------------------
# corner points of total domain
minPoint = Vector3(xmin,ymin,zmin)
maxPoint = Vector3(xdim,ydim,zdim)

# set larger area of mntable to ensure some settling of particles does
# not move outside
minPointTbl = Vector3(xmin-5.*maxRadius,ymin-5.*maxRadius,zmin-5.*maxRadius) 
maxPointTbl = Vector3(xdim+5.*maxRadius,ydim+5.*maxRadius,zdim+5.*maxRadius)

# neighbour table
mntable = MNTable3D(
    minPoint=minPointTbl, 
    maxPoint=maxPointTbl,
    gridSize=2.5*maxRadius,
    numGroups=2)

# -- setup packer --
# iteration parameters
insertFails = 1000 # 1e4
maxIter = 1000
tol = 1.0e-6

# packer
# with random seeding that changes
if geoIter > 1:
    packer = InsertGenerator3D(minRadius, maxRadius, insertFails, maxIter, tol, True)
else:
    # set as False for calibration, achieve same packing seed
    packer = InsertGenerator3D(minRadius, maxRadius, insertFails, maxIter, tol, False)

# boundary planes
bottomPlane=Plane(minPoint,Vector3(0.0,1.0,0.0))
leftPlane=Plane(minPoint,Vector3(1.0,0.0,0.0))
frontPlane=Plane(minPoint,Vector3(0.0,0.0,1.0))
topPlane=Plane(maxPoint,Vector3(0.0,-1.0,0.0))
rightPlane=Plane(maxPoint,Vector3(-1.0,0.0,0.0))
backPlane=Plane(maxPoint,Vector3(0.0,0.0,-1.0))

# get normals to layers
norm_x = math.cos(math.radians(90-theta))
norm_y = math.sin(math.radians(90-theta))

# y coordinate of lower point of layers
bottom_y = ymin

# thickness in y dimension
#overlap_y = overlap/math.cos(math.radians(90-theta)) #######################
overlap_y = overlap/math.cos(math.radians(theta))

#overlap_y = overlap
#print overlap_y

#exit()

# normals to planes when
# when layers fully inside the volume
top_normal_inside = Vector3(-1*norm_x, -1*norm_y, 0.0)
bottom_normal_inside = Vector3(norm_x, norm_y, 0.0)

# normals to top and bottoms of total model
top_normal_outside = Vector3(0.0,-1.0,0.0)
bottom_normal_outside = Vector3(0.0,1.0,0.0)

bottomPt_poly = Vector3(xmin, ymin, zmin)
topPt_poly = Vector3(xdim, ydim, zdim)

# if horizontal layers
if theta==0:
    overlap_y = overlap
    top_normal_inside = top_normal_outside
    bottom_normal_inside = bottom_normal_outside
    ydist = xdim/math.tan(math.radians(90-theta))
elif theta==90:
    ydist = 0
else:
    ydist = xdim/math.tan(math.radians(90-theta))

print "overlap perpen separation"
print overlap_y

# make polyhedrons in model
tag_i =1
top_y = 0
xL = xmin
bottom_y_right = bottom_y - ydist
while bottom_y_right < ydim and xL < xdim:

    # vertical layers
    if theta==90:
        bottomPlane_poly = bottomPlane
        topPlane_poly = topPlane

        leftPt = Vector3(xL, ymin, zmin)
        leftPlane_poly = Plane(leftPt, Vector3(1.0, 0.0, 0.0)) 

        xR = xL + thickness + overlap
        if xR > xdim:
            xR = xdim

        rightPt = Vector3(xR, ymin, zmin)
        rightPlane_poly = Plane(rightPt, Vector3(-1.0, 0.0, 0.0)) 

        xL = xL + thickness
    else:
        # normal to layer planes
        bottom_normal = bottom_normal_inside
        top_normal = top_normal_inside

        # if bottom layer
        if tag_i == 1:
            bottom_normal = bottom_normal_outside
            #top_normal = top_normal_outside
        
        # max y coordinate of upper wall of layer
        top_y = thick_y+bottom_y+overlap_y
            
        # make points to define top and bottom layer planes
        bottomPt_layer = Vector3(xmin, bottom_y, zmin)
        topPt_layer = Vector3(xmin, top_y, zmin)
            
        # define points for each polyhedra
        bottomPlane_poly = Plane(bottomPt_layer, bottom_normal)
        topPlane_poly = Plane(topPt_layer, top_normal)

        leftPlane_poly = leftPlane
        rightPlane_poly = rightPlane

        y_coors = str(bottom_y)+' '+str(top_y)
        print y_coors

        bottom_y = bottom_y+thick_y # next bottom layer point
        bottom_y_right = bottom_y - ydist
                
    # over number of non-side intersecting polyhedra
    poly = ConvexPolyhedron(bottomPt_poly, topPt_poly) 
                
    # add planes to the polyhedron
    # planes are agnostic to each other: okay if overlapping
    poly.addPlane(bottomPlane_poly)
    poly.addPlane(leftPlane)
    poly.addPlane(frontPlane)
    poly.addPlane(topPlane_poly)
    poly.addPlane(rightPlane)
    poly.addPlane(backPlane)
            
    # add plane for top and bottom of total boundary
    poly.addPlane(topPlane) 
    poly.addPlane(bottomPlane) 

    # add left and right plane, that differ when theta == 90
    poly.addPlane(leftPlane_poly)
    poly.addPlane(rightPlane_poly)

    # pack volume
    # pack particles into volume
    packer.generatePacking(volume=poly,
 		               ntable=mntable,
 		               groupID=0,
 		               tag=tag_i)
	
	
    tag_i=tag_i+1
	#exit()

max_tag = tag_i-1                   

# tag particles based on particle tags, which depend on
# what polyhedoron it is in
#ONLY WORKS WHEN PARTICLES TAGGED SEQUENTIALLY AND ONLY NEIGHBORED BY
#PARTICLE TAGS AT +/- SELF.PARTICLE TAG
layer_i =1
while layer_i <= max_tag:
    # set bond tags for same values within each layer
    mntable.generateBondsTagged(
        groupID=0,
        tolerance = 1.0e-4, # 1.0e-5
        bondID = layer_i, 
        particleTag1 = layer_i, # for particles with same tag, set bond tag equal to particle tag+1
        particleTag2 = layer_i        
    )

    # set bond tag=0 for values between layers
    mntable.generateBondsTagged(
        groupID=0,
        tolerance=1.0e-4,
        bondID = 0, 
        particleTag1 = layer_i, 
        particleTag2 = layer_i+1        
    )

    layer_i = layer_i+1
    

# calculate and print the porosity:
volume = xdim*ydim*zdim
porosity = (volume - mntable.getSumVolume(groupID=0))/volume
print "Porosity: ", porosity

# write a geometry file
mntable.write("geo_"+geo_file +".geo", 1)
mntable.write("vtk_"+geo_file+".vtu", 2)

print "geo_"+geo_file

