# Distance: Find distance between sites in cell, in quoted Cartesian coordinates

# Andy Paul Chen, Monday, 9 August 2021, Little Italy, Cleveland, Ohio (National Day of Singapore)
# 2 Nov 2023 I can't believe I am back again

from Cartesian import * # Cartesian coordinates package
from CQuery import * # coordinates query package

def distance(c1, c2):
    # Read Cartesian coordinates 1 and 2, derive distance (in Angstroms)
    # c1, c2 are arrays of three float values
    dist=math.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)
    return dist

def isNearest(data, c):
    # read c (array of three floats)
    # return index of closest site to coordinates
    # display distance
    dcindex = 8 if isSeldyn(data) else 7 # Index of Direct/Cartesian line
    
    # Switch to Cartesian coordinates
    if isCart(data):
        print("Coordinates form: Cartesian")
    else:
        print("Coordinates form: Fractional; converting to Cartesian")
        data=switchCart(data)
        
    # Line-by-line analysis
    for line in range(len(data)):
        if line > dcindex:
            dirs = re.findall(r"-?\d+\.\d+", data[line].strip()) #coordinates extracted from line in data
            dist = distance([float(dirs[0]),float(dirs[1]),float(dirs[2])], c)
            # print("Index: ", line-dcindex, "Distance: ", "{:11f}".format(dist)) # test line
            if line == dcindex+1:
                mindex = 1
                mindist = dist
            else:
                if dist < mindist:
                    mindex = line-dcindex
                    mindist = dist
                    # print("NEW LOW") # test line
                    
    print("Minimum distance agreement: ", "{:11f}".format(mindist), ", Index #", mindex)
    return mindex