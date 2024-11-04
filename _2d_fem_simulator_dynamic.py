from random import triangular
from re import X
from ssl import ALERT_DESCRIPTION_DECOMPRESSION_FAILURE
import pygame
import numpy as np
import math
import sys
import pygame.gfxdraw
import random

pygame.init()

win = pygame.display.set_mode((800,800))
pygame.display.set_caption("PYTHON TEST APPLICATION")
font = pygame.font.SysFont("courier", 12, True)
fontsmall = pygame.font.SysFont("courier", 10, True)

gamestate = 0
clock = pygame.time.Clock() 
user_text = '' 

play = 0

#buttons = [[0, 20, "draw polygon"], [0, 60, "draw void"], [0, 100, "clear screen"], [0, 140, "construct mesh:"], [0, 180, ""]]


#add a type of material (youngs modulus, poisson ratio, ...)
#draw polygon (type of material (0 = void))
#enter polygon coordanites manually
#clear screen
#construct mesh
#add point load
#add node displacement restriction
#run simulation
#turn on/off scale grid

#pan (rclick), zoom (mouse wheel)

buttons = np.array([[0, 20, "reset everything", 0],
                    [0, 60, "create material polygon", 1],
                    [0, 100, "create void polygon", 2],
                    [0, 140, "construct mesh", 3],
                    [0, 180, "set bulk properties (displacement x, displacement y, load x, load y, youngs modulus, poisson ratio, density)", 4],
                    [0, 220, "set point properties (displacement x, displacement y, load x, load y, point 1, <point 2, point 3...>)", 5],
                    [0, 260, "run simulation", 6],
                    [0, 300, "set # of iterations per program tick (0 = pause, 1 = 'normal speed')", 7]])

# textbox_coordsx = [0, 0, 0, 0, 0, 0, 0]
# textbox_coordsy = [20, 60, 100, 140, 180, 220, 260]

simulationspeed = 1
timestep = 0.05

renderingstartx = 0
renderingstarty = 0

renderingoriginx = 0
renderingoriginy = 0

renderingboxsize = 240
renderingboxscale = 1

color_active = pygame.Color('lightskyblue3') 
color_passive = pygame.Color('chartreuse4') 
color = color_passive 
  
active = -1
  
def checkCollision(colx, coly, colrectx1, colrecty1, colrectx2, colrecty2):
    if colx < colrectx2 and colx > colrectx1 and coly < colrecty2 and coly > colrecty1:
        return True
    else:  
        return False

# def checkTriangleCollision(x1, y1, x2, y2, x3, y3, xp, yp):
#     #x1, y1, x2, y2, x3, y3, xp, yp = map(float, input().split())

#     # Calculate the cross products (c1, c2, c3) for the point relative to each edge of the triangle
#     c1 = (x2 - x1) * (yp - y1) - (y2 - y1) * (xp - x1)
#     c2 = (x3 - x2) * (yp - y2) - (y3 - y2) * (xp - x2)
#     c3 = (x1 - x3) * (yp - y3) - (y1 - y3) * (xp - x3)

#     # Check if all cross products have the same sign (inside the triangle) or different signs (outside the triangle)
#     if (c1 < 0 and c2 < 0 and c3 < 0) or (c1 > 0 and c2 > 0 and c3 > 0):
#         return True
#     else:
#         return False




points = []
points2 = []
elements = []
#elements2 = []


polygon_points = []
polygon_voids = []
width = 3
height = 3



def drawPolygon(dppoints, dpcolor):
    if len(dppoints) > 2:
        pygame.draw.polygon(win, dpcolor, dppoints, 0)

# def circumcircle(ax, ay, bx, by, cx, cy):
#     d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
#     ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d
#     uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d

#     a = math.sqrt((ax - bx)*(ax-bx) + (ay-by)*(ay-by))
#     b = math.sqrt((cx - bx)*(cx-bx) + (cy-by)*(cy-by))
#     c = math.sqrt((cx - ax)*(cx-ax) + (cy-ay)*(cy-ay))

#     R = (a*b*c) / math.sqrt((a + b + c)*(b + c - a)*(c + a - b)*(a + b - c))

#     return (ux, uy, R)

# def distance(point1,point2):
#     return math.sqrt((point1[0]-point2[0])**2 + (point1[1]-point2[1])**2)

# def compareEdges(ce1, ce2, cevals):

#     ce_i = 0
#     while ce_i < 3:
#         ce_j = 0 
#         while ce_j < 3:
#             if ce1[ce_i][0][0] == ce2[ce_j][0][0] and ce1[ce_i][0][1] == ce2[ce_j][0][1] and ce1[ce_i][1][0] == ce2[ce_j][1][0] and ce1[ce_i][1][1] == ce2[ce_j][1][1]: #0 = 0 and 1 = 1
#                 cevals[ce_i] = 1
#                 #this edge is the same
#             elif ce1[ce_i][1][0] == ce2[ce_j][0][0] and ce1[ce_i][1][1] == ce2[ce_j][0][1] and ce1[ce_i][0][0] == ce2[ce_j][1][0] and ce1[ce_i][0][1] == ce2[ce_j][1][1]: #1 = 0 and 0 = 1
#                 #this edge is the same
#                 cevals[ce_i] = 1


#             ce_j = ce_j + 1
#         ce_i = ce_i + 1

#     return cevals




        
# def BowyerWatson(pointList):
#     #// pointList is a set of coordinates defining the points to be triangulated
#     #print(pointList)

#     supertriangle = [[10000, 10000], [-10000, 10000], [0, -10000]] #:= empty triangle mesh data structure
#     triangulation = [supertriangle]
#     #add super-triangle to triangulation // must be large enough to completely contain all the points in pointList
    
#     i = 0
#     while i < len(pointList):
#         badTriangles = []
#         badTrianglesindex = []
#         j = 0
#         while j < len(triangulation):
#                 circleresult = circumcircle(triangulation[j][0][0], triangulation[j][0][1], triangulation[j][1][0], triangulation[j][1][1], triangulation[j][2][0], triangulation[j][2][1])
#                 if math.sqrt((pointList[i][0] - circleresult[0])*(pointList[i][0]-circleresult[0]) + (pointList[i][1] - circleresult[1])*(pointList[i][1]-circleresult[1])) < circleresult[2]:
#                     badTriangles.append(triangulation[j])
#                     badTrianglesindex.append(j)
#                 j = j + 1
#         validedges = []
#         j = 0
#         while j < len(badTriangles):
#             edges = [[badTriangles[j][0], badTriangles[j][1]], [badTriangles[j][1], badTriangles[j][2]], [badTriangles[j][2], badTriangles[j][0]]]
#             k = j+1
#             edgevals = [0, 0, 0]
#             while k < len(badTriangles):
#                 #if k != j:
#                     #compare each edge
#                 compareresult = compareEdges(edges, [[badTriangles[k][0], badTriangles[k][1]], [badTriangles[k][1], badTriangles[k][2]], [badTriangles[k][2], badTriangles[k][0]]], edgevals)
#                 #print(compareresult)
#                 edgevals = compareresult
#                 k = k + 1
#             if edgevals[0] == 0:
#                 validedges.append(edges[0])
#                 #add 1st edge to polygon
#             if edgevals[1] == 0:
#                 validedges.append(edges[1])
#                 #add 2nd edge to polygon
#             if edgevals[2] == 0:
#                 validedges.append(edges[2])
#                 #add 3rd edge to polygon
               
#             j = j + 1
#         newtriangulation = []
#         j = 0
#         while j < len(triangulation):
#             if j in badTriangles == False:
#                 newtriangulation.append(triangulation[j])
#             j = j + 1
#         triangulation = newtriangulation

#         j = 0
#         while j < len(validedges):
#             triangulation.append([pointList[i], validedges[j][0], validedges[j][1]])
#             j = j + 1

#         i = i + 1

#     finalnodes = []
#     finaltriangulation = []
#     j = 0
#     while j < len(triangulation):
#         nogood = 0
#         k = 0
#         while k < 3:
#             l = 0
#             while l < 3:
#                 if triangulation[j][k][0] == supertriangle[l][0] and triangulation[j][k][1] == supertriangle[l][1]:
#                     l = 2
#                     k = 2
#                     nogood = 1
#                 l = l + 1
#             k = k + 1
#         if nogood == 0:

#             thistriangle = []
#             l = 0
#             while l < 3:
#                 m = 0
#                 nogood = -1
#                 while m < len(finalnodes):
#                     if triangulation[j][l][0] == finalnodes[m][0] and triangulation[j][l][1] == finalnodes[m][1]:
#                         nogood = m
#                         m = len(finalnodes)-1
#                     m = m + 1

#                 if nogood != -1: #triangulation[j][l] in finalnodes == True:
#                     #k = 0
#                     #while k < len(finalnodes):
#                     #    if finalnodes[k] == triangulation[j][l]:
#                     thistriangle.append(nogood)
#                     #        k = len(finalnodes)-1
#                     #    k = k + 1
#                 else:

#                     finalnodes.append([triangulation[j][l][0], triangulation[j][l][1], 0, 0, 0, 0, 0, 0, 0])
#                     thistriangle.append(len(finalnodes)-1)
#                 l = l + 1
#             finaltriangulation.append(thistriangle)
#         j = j + 1

#     print([finalnodes, finaltriangulation])
#     return [finalnodes, finaltriangulation]

                   #if checkTriangleCollision(triangulation[0][0], triangulation[0][1], triangulation[1][0], triangulation[1][1], triangulation[2][0], triangulation[2][1], pointList[i][0], pointList[i][1]) == True:


    #for each point in pointList do // add all the points one at a time to the triangulation
    #    badTriangles := empty set
    #    for each triangle in triangulation do // first find all the triangles that are no longer valid due to the insertion
            #if point is inside circumcircle of triangle
             #   add triangle to badTriangles
        #polygon := empty
        #for each triangle in badTriangles do // find the boundary of the polygonal hole
            #for each edge in triangle do
            #    if edge is not shared by any other triangles in badTriangles
        #            add edge to polygon
        #for each triangle in badTriangles do // remove them from the data structure
        #    remove triangle from triangulation
        #for each edge in polygon do // re-triangulate the polygonal hole
        #    newTri := form a triangle from edge to point
        #    add newTri to triangulation
    # for each triangle in triangulation // done inserting points, now clean up
    #     if triangle contains a vertex from original super-triangle
    #         remove triangle from triangulation
    # return triangulation




# def constructMeshFromPolygonsquare(fineness):
#     #for all points from 0 - 800 (size of screen)
#     meshpoints = []
#     meshcolor = (0, 0, 255)

#     meshsize = math.ceil(800/fineness)
#     meshgrid = np.zeros((800, 800))
#     meshpoints = []
#     meshelements = []

#     k = 0
#     i = 0
#     while i < 800:
        
#         j = 0
#         while j < 800:
#             if win.get_at((i, j)) == meshcolor:
#                 k = k + 1
#                 meshgrid[int(i/fineness)][int(j/fineness)] = k
#                 meshpoints.append([i, j])
#                 #add to points

#             j = j + fineness
#         i = i + fineness

#     if len(meshpoints) > 0:

#         i = 0
#         while i < len(meshpoints):
            
#             meshpx = int(meshpoints[i][0]/fineness)
#             meshpy = int(meshpoints[i][1]/fineness)
#             if meshgrid[meshpx+1][meshpy] != 0 and meshgrid[meshpx+1][meshpy+1] != 0 and meshgrid[meshpx][meshpy+1] != 0:
#                 meshelements.append([i, int(meshgrid[meshpx+1][meshpy]) - 1, int(meshgrid[meshpx+1][meshpy+1]) - 1, int(meshgrid[meshpx][meshpy+1]) - 1])



#             i = i + 1

#     return (meshpoints, meshelements)







def constructMeshFromPolygon(fineness):
    #for all points from 0 - 800 (size of screen)
    meshpoints = []
    meshcolor = (0, 0, 255)

    meshsize = math.ceil(800/fineness)
    meshsize2 = math.ceil(800/(fineness*(math.sqrt(3)/2)))
    meshgrid = np.zeros((800, 800))
    meshpoints = []
    meshelements = []

    k = 0
    i = 0
    type = 0
    
    while i < meshsize2:
        if type == 0:
            type = 0.5*fineness
        else:
            type = 0

        j = 0
        while j < meshsize:
            if win.get_at((int(j*fineness+type), int(i*fineness*(math.sqrt(3)/2)))) == meshcolor:
                k = k + 1
                meshgrid[j][i] = k
                meshpoints.append([(j*fineness+type), i*fineness*(math.sqrt(3)/2), j, i, type, 0, 0, 0, 0])
                #add to points

            j = j + 1#fineness
        i = i + 1#fineness*(math.sqrt(3)/2)

    if len(meshpoints) > 0:
        i = 0

        while i < len(meshpoints):


            meshpx = meshpoints[i][2]
            meshpy = meshpoints[i][3]
            type = meshpoints[i][4]

            if type != 0:
                #print ([i, int(meshgrid[meshpx][meshpy+1]) - 1, int(meshgrid[meshpx+1][meshpy+1]) - 1, int(meshgrid[meshpx+1][meshpy]) - 1])

                if meshgrid[meshpx+1][meshpy+1] != 0 and meshgrid[meshpx+1][meshpy] != 0:
                    meshelements.append([i, int(meshgrid[meshpx+1][meshpy+1]) - 1, int(meshgrid[meshpx+1][meshpy]) - 1])
                    

                if meshgrid[meshpx+1][meshpy+1] != 0 and meshgrid[meshpx][meshpy+1] != 0:
                    meshelements.append([i, int(meshgrid[meshpx+1][meshpy+1]) - 1, int(meshgrid[meshpx][meshpy+1]) - 1])


                k = 0
            else:
                if meshgrid[meshpx][meshpy+1] != 0 and meshgrid[meshpx+1][meshpy] != 0:
                    meshelements.append([i, int(meshgrid[meshpx][meshpy+1]) - 1, int(meshgrid[meshpx+1][meshpy]) - 1])

                if meshgrid[meshpx][meshpy+1] != 0 and meshgrid[meshpx-1][meshpy+1] != 0:
                    meshelements.append([i, int(meshgrid[meshpx][meshpy+1]) - 1, int(meshgrid[meshpx-1][meshpy+1]) - 1])

                k = 0
            i = i + 1

    return (meshpoints, meshelements)


def triangularMeshFEM(tm_nodes, tm_elements):

    num_nodes = 0
    num_elements = 0#len(elements)

    coords = []#np.zeros((num_nodes, 2))
    supports = []#np.zeros((num_nodes, 2))
    loads = []#np.zeros((num_nodes, 2))
    node_directions = []
    node_index = []

    #mass_local_template = area*density*[[2, 1], [1, 2]]/6

    rawmassmatrix = np.array([[2, 0, 1, 0, 1, 0], 
                              [0, 2, 0, 1, 0, 1],
                              [1, 0, 2, 0, 1, 0],
                              [0, 1, 0, 2, 0, 1],
                              [1, 0, 1, 0, 2, 0],
                              [0, 1, 0, 1, 0, 2]])


    
  
    i = 0
    j = 0
    while i < len(tm_nodes):
        if len(tm_nodes[i]) > 0:
            num_nodes = num_nodes + 1
            coords.append([tm_nodes[i][0], tm_nodes[i][1]])
            loads.append([tm_nodes[i][5], tm_nodes[i][6]])
            supports.append([tm_nodes[i][7], tm_nodes[i][8]])
            node_directions.append([j*2, j*2+1])
            node_index.append(j)
            j = j + 1
        i = i + 1

    num_dofs = 2*num_nodes
    dofs = []
    stiffness_global = np.zeros((num_dofs, num_dofs))
    mass_global = np.zeros((num_dofs, num_dofs))
    connectivity = []
    forces_global = np.zeros(num_dofs)

    #k = tA(BT)DB
    t = 1 #assume thickness of 1

    i = 0
    while i < len(tm_elements):
        if len(tm_nodes[tm_elements[i][0]]) > 0 and len(tm_nodes[tm_elements[i][1]]) > 0 and len(tm_nodes[tm_elements[i][2]]) > 0:

            connectivity.append([node_index[tm_elements[i][0]], node_index[tm_elements[i][1]], node_index[tm_elements[i][2]]])
        
            nodes_for_this_element = (connectivity[len(connectivity)-1][0], connectivity[len(connectivity)-1][1], connectivity[len(connectivity)-1][2])

            dofs.append([node_directions[nodes_for_this_element[0]][0], node_directions[nodes_for_this_element[0]][1], node_directions[nodes_for_this_element[1]][0], node_directions[nodes_for_this_element[1]][1], node_directions[nodes_for_this_element[2]][0], node_directions[nodes_for_this_element[2]][1]])

            
            


            x13 = coords[nodes_for_this_element[0]][0] - coords[nodes_for_this_element[2]][0]
            x23 = coords[nodes_for_this_element[1]][0] - coords[nodes_for_this_element[2]][0]
            y13 = coords[nodes_for_this_element[0]][1] - coords[nodes_for_this_element[2]][1]
            y23 = coords[nodes_for_this_element[1]][1] - coords[nodes_for_this_element[2]][1]

            y31 = coords[nodes_for_this_element[2]][1] - coords[nodes_for_this_element[0]][1]
            y12 = coords[nodes_for_this_element[0]][1] - coords[nodes_for_this_element[1]][1]
            x32 = coords[nodes_for_this_element[2]][0] - coords[nodes_for_this_element[1]][0]
            x21 = coords[nodes_for_this_element[1]][0] - coords[nodes_for_this_element[0]][0]

            detJ = x13*y23 - y13*x23
            A = abs(detJ)/2
            B = np.array([[y23, 0, y31, 0, y12, 0], 
                          [0, x32, 0, x13, 0, x21],
                          [x32, y23, x13, y31, x21, y12]])/detJ

            BT = B.transpose() 
            v = poissonratio
            E = youngsmod
            

            D = np.array([[1, v, 0], 
                          [v, 1, 0],
                          [0, 0, (1-v)/2]])*(E/(1-v*v))
    
            #print(B)
            #print(D)
            #print(BT)
            #print(np.dot(B, D))

            k_local = t*A*np.dot(np.dot(BT, D), B)

            m_local = A*density*1*rawmassmatrix/12

            #mass_local = distance_for_this_element*mass_local_template
            #mass_local_inv = numpy.linalg.inv(mass_local)


            #print(k_local)

            j = 0
            k = len(dofs)-1 #for this node...
            while j < 6:
                stiffness_global[int(dofs[k][j])][int(dofs[k][0])] += k_local[j][0]
                stiffness_global[int(dofs[k][j])][int(dofs[k][1])] += k_local[j][1]
                stiffness_global[int(dofs[k][j])][int(dofs[k][2])] += k_local[j][2]
                stiffness_global[int(dofs[k][j])][int(dofs[k][3])] += k_local[j][3]
                stiffness_global[int(dofs[k][j])][int(dofs[k][4])] += k_local[j][4]
                stiffness_global[int(dofs[k][j])][int(dofs[k][5])] += k_local[j][5]

                mass_global[int(dofs[k][j])][int(dofs[k][0])] += m_local[j][0]
                mass_global[int(dofs[k][j])][int(dofs[k][1])] += m_local[j][1]
                mass_global[int(dofs[k][j])][int(dofs[k][2])] += m_local[j][2]
                mass_global[int(dofs[k][j])][int(dofs[k][3])] += m_local[j][3]
                mass_global[int(dofs[k][j])][int(dofs[k][4])] += m_local[j][4]
                mass_global[int(dofs[k][j])][int(dofs[k][5])] += m_local[j][5]


                j = j + 1
        i = i + 1

    kt_global = stiffness_global*1 #why *1??? is it to establish kt_global as its own thing instead of copying the id to the array?? weird.
    m_global = mass_global*1
    i = 0

    while i < num_nodes:

        j = 0
        while j < 2:
            #print(int(supports[i][j]))
            if int(supports[i][j]) == 1:
                #print("success")
                kt_global[int(node_directions[i][j]), :] = 0
                kt_global[:, int(node_directions[i][j])] = 0
                kt_global[int(node_directions[i][j]), int(node_directions[i][j])] = 1

                m_global[int(node_directions[i][j]), :] = 0
                m_global[:, int(node_directions[i][j])] = 0
                m_global[int(node_directions[i][j]), int(node_directions[i][j])] = 1
                forces_global[int(node_directions[i][j])] = 0
            else:
                forces_global[int(node_directions[i][j])] = loads[i][j]
            j = j + 1
        i = i + 1

    #print("connectivity:")
    #print(connectivity)
    #print("forces_global:")
    #print(forces_global)


    displacement = np.zeros(num_dofs)
    velocities = np.zeros(num_dofs)
    accelerations = np.zeros(num_dofs)
    forces = forces_global
    beta = 1/4
    gamma = 1/2


    
    invmass = np.linalg.inv(m_global)
    # print('a')
    # print(forces) #18x1
    # print(invmass) #18x18
    # print(kt_global) #18x18
    # print(displacement) #18x1

    # print(kt_global*displacement)
    # print(invmass*(forces - kt_global*displacement))
    # print('b')
    #return values are:
    #displacement, velocity, acceleration, force, global K, global M, beta, gamma, const_a

    return [displacement, velocities, np.dot(invmass, (forces - np.dot(kt_global, displacement))), forces, kt_global, m_global, invmass, beta, gamma, m_global/(beta*timestep*timestep)]




    #need: mass_local & mass_global


    #displacement = np.linalg.solve(kt_global, forces_global)
    #reaction_forces = np.dot(stiffness_global, displacement)

    #num_of_steps_to_do = 100
    #timestep = 0.05

    #step = 1
    #while step < num_of_steps_to_do:
        #accelerations = 
        #const_a = mass_global/(beta*timestep*timestep)




        #step = step + 1





    # nodes_next = []
    # i = 0
    # while i < len(node_index):
    #     index = node_index[i]
    #     nodes_next.append([tm_nodes[index][0] + displacement[2*i]/1000, tm_nodes[index][1] + displacement[2*i + 1]/1000])
    #     i = i + 1

    # print("Displacements (mm):")
    # print(displacement * 1000)
    # print("Reaction Forces (N):")
    # print(reaction_forces)


    # #print (nodes_next)
    # return nodes_next


    





    







# def doDynamicalFEMStep():


#                         #     FEMdisplacement = bFresult[6][0]
#                         # FEMvelocity = bFresult[6][1]
#                         # FEMacceleration = bFresult[6][2]
#                         # FEMforces = bFresult[6][3]
#                         # FEMk = bFresult[6][4]
#                         # FEMm = bFresult[6][5]
#                         # FEMinvm = bFresult[6][6]
#                         # FEMbeta = bFresult[6][7]
#                         # FEMgamma = bFresult[6][8]
#                         # FEMalpha = bFresult[6][9]


#     kprime = kt_global + const_a
#     fprime = forces + const_a*(displacement + timestep*velocities + (1/2 - beta)*(timestep*timestep)*accelerations)

#     newdisplacements = np.linalg.solve(kprime, fprime) #update displacements

#     newaccelerations = (newdisplacements - displacements - timestep*velocities - timestep*timestep*(1/2 - beta)*accelerations)/(beta*timestep*timestep)

#     velocities = velocities + timestep*((1 - gamma)*accelerations + gamma*newaccelerations)

#     displacements = newdisplacements
#     accelerations = newaccelerations

#     print ("a")















def drawPoints(dpoints, delements, dpoints_color):
    if len(dpoints) > 0:
        if len(delements) > 0:
            i = 0
            while i < len(delements):
                #print ("b" + str(i))
                if len(dpoints[delements[i][0]]) == 0 or len(dpoints[delements[i][1]]) == 0:
                    #delements[i] = []
                    j = 0 #unneeded
                else:
                    if len(delements[i]) == 3:
                        drawx = (dpoints[delements[i][0]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy = (dpoints[delements[i][0]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        drawx2 = (dpoints[delements[i][1]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy2 = (dpoints[delements[i][1]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        drawx3 = (dpoints[delements[i][2]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy3 = (dpoints[delements[i][2]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        #print(([drawx, drawy], [drawx2, drawy2], [drawx3, drawy3]))




                        x13 = drawx - drawx3
                        x23 = drawx2 - drawx3
                        y13 = drawy - drawy3
                        y23 = drawy2 - drawy3

                        y31 = drawy3 - drawy
                        y12 = drawy - drawy2
                        x32 = drawx3 - drawx2
                        x21 = drawx2 - drawx

                        detJ = x13*y23 - y13*x23
                        A = abs(detJ)/2

                        area_ratio = 255*(A/(30*30*math.sqrt(3)/4))*(A/(30*30*math.sqrt(3)/4))*(A/(30*30*math.sqrt(3)/4))

                        if area_ratio > 255:
                            area_ratio = 255
                        elif area_ratio < 0:
                            area_ratio = 0

                        #colornew = pygame.color.hsva()

                        #newcolor = pygame.color(255, 255, 0)

                        # hprime = area_ratio/60
                        # Xitalic = (1 - abs(float(hprime % 2) - 1))

                        # if hprime < 1:
                        #     ctriangle = (255, 0, Xitalic)
                        # elif hprime < 2:
                        #     ctriangle = (Xitalic, 0, 255)
                        # elif hprime < 3:
                        #     ctriangle = (0, 255, Xitalic)
                        # elif hprime < 4:
                        #     ctriangle = (0, Xitalic, 255)
                        # elif hprime < 5:
                        #     ctriangle = (Xitalic, 255, 0)
                        # else:
                        #     ctriangle = (255, Xitalic, 0)

                        #print(hsv_to_rgb(area_ratio, 1, 1, 1))

                        pygame.draw.polygon(win, (255-area_ratio, 0, area_ratio), ([drawx, drawy], [drawx2, drawy2], [drawx3, drawy3]), 0)
                        #pygame.draw.line(win, (255, 255, 255), [drawx, drawy], [drawx2, drawy2], 3)
                        #pygame.draw.line(win, (255, 255, 255), [drawx3, drawy3], [drawx2, drawy2], 3)
                        #pygame.draw.line(win, (255, 255, 255), [drawx3, drawy3], [drawx, drawy], 3)
                        
                        

                        #color = dpoints[]
                        
                        
                    else:
                        drawx = (dpoints[delements[i][0]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy = (dpoints[delements[i][0]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        drawx2 = (dpoints[delements[i][1]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy2 = (dpoints[delements[i][1]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        drawx3 = (dpoints[delements[i][2]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy3 = (dpoints[delements[i][2]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        drawx4 = (dpoints[delements[i][3]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                        drawy4 = (dpoints[delements[i][3]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                        #print(([drawx, drawy], [drawx2, drawy2], [drawx3, drawy3]))

                        pygame.draw.polygon(win, dpoints_color, ([drawx, drawy], [drawx2, drawy2], [drawx3, drawy3], [drawx4, drawy4]), 0)

                    #pygame.draw.line(win, dpoints_color, [drawx, drawy], [drawx2, drawy2], 3)
                i = i + 1
        i = 0
        while i < len(dpoints):
            #print ("A" + str(i))
            if len(dpoints[i]) > 0:

                drawx = (dpoints[i][0] + renderingoriginx)*renderingboxscale + renderingstartx
                drawy = (dpoints[i][1] + renderingoriginy)*renderingboxscale + renderingstarty

                pygame.draw.rect(win, dpoints_color, (drawx, drawy, width, height))
                text = fontsmall.render(str(i), 1, (255,255,255))
                #text = fontsmall.render(str(i) + ": " + str(points[i][0]) + ", " + str(points[i][1]) + ": " + str(points[i][2]) + ", " + str(points[i][3]), 1, (255,255,255))
                win.blit(text, (drawx+2, drawy-2))
            i = i + 1


scalar = float #??????
def hsv_to_rgb( h:scalar, s:scalar, v:scalar, a:scalar ) -> tuple:
    a = int(255*a)
    if s:
        if h > 0.99: h = 0.0
        i = int(h*6.0); f = h*6.0 - i
        
        w = v * (255.0 - s)
        q = v * (255.0 - s * f)
        t = v * (255.0 - s * (255.0 - f))
        v = int(255*v)
        if i==0: return (v, t, w, a)
        if i==1: return (q, v, w, a)
        if i==2: return (w, v, t, a)
        if i==3: return (w, q, v, a)
        if i==4: return (t, w, v, a)
        if i==5: return (v, w, q, a)
    else:
        v = int(255*v)
        return (v, v, v, a)

def redrawGameWindow():
    win.fill((0,0,0))



    drawPoints(points, elements, (255, 0, 0))

    drawPoints(points2, elements, (0, 255, 0))

    if len(polygon_points) > 0:
        i = 0
        while i < len(polygon_points):
            drawPolygon(polygon_points[i], (0, 0, 255))
            i = i + 1

    if len(polygon_voids) > 0:
        i = 0
        while i < len(polygon_voids):
            drawPolygon(polygon_voids[i], (0, 0, 0))
            i = i + 1

    i = 0
    while i < len(buttons):
        text = font.render(str(buttons[i][2]), 1, (255,255,255))

        win.blit(text, (int(buttons[i][0]), int(buttons[i][1])-20))

        if active == int(buttons[i][3]):
            
            color = color_active 
            input_rect = pygame.Rect(int(buttons[i][0]), int(buttons[i][1]), 200, 20)
            pygame.draw.rect(win, color, input_rect) 
            text_surface = font.render(user_text, True, (255, 255, 255))
            win.blit(text_surface, (int(buttons[i][0])+3, int(buttons[i][1])+3)) 
        else: 
            color = color_passive 
            input_rect = pygame.Rect(int(buttons[i][0]), int(buttons[i][1]), 200, 20)
            pygame.draw.rect(win, color, input_rect) 
        #clock.tick(60) 
        i = i + 1

    pygame.display.update() 


def buttonFunctions(bFnum, bFpoints, bFpoints2, bFelements, bFpolygon_points, bFpolygon_voids, bFactive):
    bFextra = 0

    user_text = ''
    if bFnum == 0:
         print ('a')
         bFpoints = []
         bFpoints2 = []
         bFelements = []
         bFpolygon_points = []
         bFpolygon_voids = []
    if bFnum == 1:
        bFpolygon_points.append([])
        bFactive = bFnum
    if bFnum == 2:
        bFpolygon_voids.append([])
        bFactive = bFnum
    if bFnum == 4 or bFnum == 5 or bFnum == 7:
        bFactive = bFnum
    if bFnum == 3:
        #combine all points from polygons and voids into shared list
        combinedpointlist = []
        i =0
        while i < len(bFpolygon_points):
            combinedpointlist = combinedpointlist + bFpolygon_points[i]
            i = i + 1
        i =0
        while i < len(bFpolygon_voids):
            combinedpointlist = combinedpointlist + bFpolygon_voids[i]
            i = i + 1


        #result = BowyerWatson(combinedpointlist)##constructMeshFromPolygon(30)

        result = constructMeshFromPolygon(30)
                
        bFpoints = result[0]
        bFelements = result[1]

        bFpolygon_points = []
        bFpolygon_voids = []
        bFpolygon_points.append([])
        bFpolygon_voids.append([])
    if bFnum == 6:
        bFextra = triangularMeshFEM(bFpoints, bFelements)


        #[displacement, velocities, invmass*(forces - kt_global*displacement), forces, kt_global, m_global, invmass, beta, gamma, m_global/(beta*timestep*timestep)]



        

    return (bFpoints, bFpoints2, bFelements, bFpolygon_points, bFpolygon_voids, bFactive, bFextra)











# buttons = np.array([[0, 20, "reset everything", 0],
#                     [0, 60, "create material polygon", 1],
#                     [0, 100, "create void polygon", 2],
#                     [0, 140, "construct mesh", 3]
#                     [0, 180, "set bulk properties (displacement x, displacement y, load x, load y, youngs modulus, poisson ratio)", 4]
#                     [0, 220, "set point properties (displacement x, displacement y, load x, load y, point 1, <point 2, point 3...>)", 5],
#                     [0, 260, "run simulation", 6]])


run = True

while run:
    #pygame.time.delay(100)
    if play == 1:
        #doDynamicalFEMStep()

        i = 0
        while i < simulationspeed:

            #print('f')
            kprime = FEMk + FEMalpha
            fprime = FEMforces + np.dot(FEMalpha, (FEMdisplacement + timestep*FEMvelocity + (1/2 - FEMbeta)*(timestep*timestep)*FEMacceleration))

            newdisplacement = np.linalg.solve(kprime, fprime) #update displacements

            newacceleration = (newdisplacement - FEMdisplacement - timestep*FEMvelocity - timestep*timestep*(1/2 - FEMbeta)*FEMacceleration)/(FEMbeta*timestep*timestep)

            FEMvelocity = FEMvelocity + timestep*((1 - FEMgamma)*FEMacceleration + FEMgamma*newacceleration)
            FEMdisplacement = newdisplacement
            FEMacceleration = newacceleration 
            #print(FEMdisplacement[0])
            #print(i)

            j = 0
            while j < len(points):
                #print([FEMdisplacement[2*i], FEMdisplacement[2*i+1]])
                points[j][0] = points_bookmark[j][0] + FEMdisplacement[2*j]
                points[j][1] = points_bookmark[j][1] + FEMdisplacement[2*j+1]
                #index = node_index[i]
                #nodes_next.append([tm_nodes[index][0] + displacement[2*i]/1000, tm_nodes[index][1] + displacement[2*i + 1]/1000])
                j = j + 1

            i = i + 1





        

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run = False
        #if event.type == pygame.MOUSEBUTTONUP:
        #    pos = pygame.mouse.get_pos()
        #    points.append(pos)
        if event.type == pygame.MOUSEBUTTONDOWN: 
            pos = pygame.mouse.get_pos()
            #points.append(pos)
            
            buttonclicked = 0
            i = 0
            while i < len(buttons): #and i != -7:
                
                if checkCollision(pos[0], pos[1], int(buttons[i][0]), int(buttons[i][1]), int(buttons[i][0])+200, int(buttons[i][1])+20):
                    #
                    
                    bFresult = buttonFunctions(i, points, points2, elements, polygon_points, polygon_voids, active)


                    if bFresult[6] != 0:

                        print('g')
                        FEMdisplacement = bFresult[6][0]
                        FEMvelocity = bFresult[6][1]
                        FEMacceleration = bFresult[6][2]
                        FEMforces = bFresult[6][3]
                        FEMk = bFresult[6][4]
                        FEMm = bFresult[6][5]
                        FEMinvm = bFresult[6][6]
                        FEMbeta = bFresult[6][7]
                        FEMgamma = bFresult[6][8]
                        FEMalpha = bFresult[6][9]
                        points_bookmark = points
                        play = 1


                    points = bFresult[0]
                    points2 = bFresult[1]
                    elements = bFresult[2]
                    polygon_points = bFresult[3]
                    polygon_voids = bFresult[4]
                    active = int(bFresult[5])
                    #print(active)
                    buttonclicked = 1
                i = i + 1

            if buttonclicked == 0:
                if active == 1:
                    #print(polygon_points)
                    pos = pygame.mouse.get_pos()
                    polygon_points[len(polygon_points)-1].append(pos)
                elif active == 2:
                    pos = pygame.mouse.get_pos()
                    polygon_voids[len(polygon_voids)-1].append(pos)


        if event.type == pygame.KEYDOWN and active != -1: 
            if event.key == pygame.K_BACKSPACE: 
                user_text = user_text[:-1] 
            elif event.key == pygame.K_RETURN:
                if active == 7:
                    simulationspeed = int(user_text)
                else:
                    newtuple = tuple(user_text.split(", "))
                    if active == 4: #(displacement x, displacement y, load x, load y, youngs modulus, poisson ratio, density)
                        i = 0
                        while i < len(points):
                            points[i][5] = float(newtuple[2])
                            points[i][6] = float(newtuple[3])
                            points[i][7] = float(newtuple[0])
                            points[i][8] = float(newtuple[1])
                            i = i + 1
                        youngsmod = float(newtuple[4])
                        poissonratio = float(newtuple[5])
                        density = float(newtuple[6])
                        #points.append([float(newtuple[0]), float(newtuple[1]), float(newtuple[2]), float(newtuple[3]), float(newtuple[4]), float(newtuple[5])])
                    elif active == 5: #(displacement x, displacement y, load x, load y, point 1, <point 2, point 3...>)
                        i = 4
                        while i < len(newtuple):
                            points[int(newtuple[i])][5] = float(newtuple[2])
                            points[int(newtuple[i])][6] = float(newtuple[3])
                            points[int(newtuple[i])][7] = float(newtuple[0])
                            points[int(newtuple[i])][8] = float(newtuple[1])
                            i = i + 1




                    #elements.append([int(newtuple[0]), int(newtuple[1]), float(newtuple[2]), float(newtuple[3])])
                user_text = ''
            else:
                user_text += event.unicode



    redrawGameWindow()


    
pygame.quit()


