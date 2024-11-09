from random import triangular
from re import X
from ssl import ALERT_DESCRIPTION_DECOMPRESSION_FAILURE #dont even remember writing these 'from' statements
import pygame
import numpy as np
import math
import time
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

buttons = np.array([[0, 20, "reset everything", 0],
                    [0, 60, "create material polygon", 1],
                    [0, 100, "create void polygon", 2],
                    [0, 140, "construct mesh", 3],
                    [0, 180, "set bulk properties (displacement x, displacement y, load x, load y, youngs modulus, poisson ratio, density)", 4],
                    [0, 220, "set point properties (displacement x, displacement y, load x, load y, point 1, <point 2, point 3...>)", 5],
                    [0, 260, "run simulation", 6],
                    [0, 300, "set # of iterations per program tick (0 = pause, 1 = 'normal speed')", 7]])


simulationspeed = 1
timestep = 0.0005

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

points = []
points2 = []
elements = []
#elements2 = []


polygon_points = []
polygon_voids = []
width = 3
height = 3



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

                        pygame.draw.polygon(win, (255-area_ratio, 0, area_ratio), ([drawx, drawy], [drawx2, drawy2], [drawx3, drawy3]), 0)
                        #pygame.draw.line(win, (255, 255, 255), [drawx, drawy], [drawx2, drawy2], 3)
                        #pygame.draw.line(win, (255, 255, 255), [drawx3, drawy3], [drawx2, drawy2], 3)
                        #pygame.draw.line(win, (255, 255, 255), [drawx3, drawy3], [drawx, drawy], 3)

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
        i = i + 1

    pygame.display.update() 


def buttonFunctions(bFnum, bFpoints, bFpoints2, bFelements, bFpolygon_points, bFpolygon_voids, bFactive):
    bFextra = 0

    user_text = ''
    if bFnum == 0:
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


    return (bFpoints, bFpoints2, bFelements, bFpolygon_points, bFpolygon_voids, bFactive, bFextra)











# buttons = np.array([[0, 20, "reset everything", 0],
#                     [0, 60, "create material polygon", 1],
#                     [0, 100, "create void polygon", 2],
#                     [0, 140, "construct mesh", 3]
#                     [0, 180, "set bulk properties (displacement x, displacement y, load x, load y, youngs modulus, poisson ratio)", 4]
#                     [0, 220, "set point properties (displacement x, displacement y, load x, load y, point 1, <point 2, point 3...>)", 5],
#                     [0, 260, "run simulation", 6]])



def drawPolygon(dppoints, dpcolor):
    if len(dppoints) > 2:
        pygame.draw.polygon(win, dpcolor, dppoints, 0)

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

    rawmassmatrix = np.array([[2, 0, 1, 0, 1, 0], 
                              [0, 2, 0, 1, 0, 1],
                              [1, 0, 2, 0, 1, 0],
                              [0, 1, 0, 2, 0, 1],
                              [1, 0, 1, 0, 2, 0],
                              [0, 1, 0, 1, 0, 2]]) #this is a template of the local mass matrix without being multiplied by anything

    i = 0
    j = 0
    while i < len(tm_nodes):
        if len(tm_nodes[i]) > 0:
            num_nodes = num_nodes + 1
            coords.append([tm_nodes[i][0], tm_nodes[i][1]]) #node coordanites
            loads.append([tm_nodes[i][5], tm_nodes[i][6]]) #node loads
            supports.append([tm_nodes[i][7], tm_nodes[i][8]]) #if the node is supported or not
            node_directions.append([j*2, j*2+1]) #reference to where the x and y direction of each node is (e.g. in a global matrix)
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

            #add the nodes for this triangle to a connectivity matrix
            connectivity.append([node_index[tm_elements[i][0]], node_index[tm_elements[i][1]], node_index[tm_elements[i][2]]])
            #get the index of each node in this triangle
            nodes_for_this_element = (connectivity[len(connectivity)-1][0], connectivity[len(connectivity)-1][1], connectivity[len(connectivity)-1][2]) 
            #add degrees of freedom for each node in this triangle to dofs
            dofs.append([node_directions[nodes_for_this_element[0]][0], node_directions[nodes_for_this_element[0]][1], node_directions[nodes_for_this_element[1]][0], node_directions[nodes_for_this_element[1]][1], node_directions[nodes_for_this_element[2]][0], node_directions[nodes_for_this_element[2]][1]])

            #get distances between nodes for calculation purposes (e.g. determinant)
            #for coords[nodes_for_this_element[A]][B], A is the node (0, 1, 2) and B is whether or not you want the x (0) or y (1) value
            x13 = coords[nodes_for_this_element[0]][0] - coords[nodes_for_this_element[2]][0]
            x23 = coords[nodes_for_this_element[1]][0] - coords[nodes_for_this_element[2]][0]
            y13 = coords[nodes_for_this_element[0]][1] - coords[nodes_for_this_element[2]][1]
            y23 = coords[nodes_for_this_element[1]][1] - coords[nodes_for_this_element[2]][1]
            y31 = coords[nodes_for_this_element[2]][1] - coords[nodes_for_this_element[0]][1]
            y12 = coords[nodes_for_this_element[0]][1] - coords[nodes_for_this_element[1]][1]
            x32 = coords[nodes_for_this_element[2]][0] - coords[nodes_for_this_element[1]][0]
            x21 = coords[nodes_for_this_element[1]][0] - coords[nodes_for_this_element[0]][0]

            detJ = x13*y23 - y13*x23 #determinant
            A = abs(detJ)/2 #area
            B = np.array([[y23, 0, y31, 0, y12, 0], 
                          [0, x32, 0, x13, 0, x21],
                          [x32, y23, x13, y31, x21, y12]])/detJ
            
            BT = B.transpose() 
            v = poissonratio
            E = youngsmod

            D = np.array([[1, v, 0], 
                          [v, 1, 0],
                          [0, 0, (1-v)/2]])*(E/(1-v*v))
    
            k_local = t*A*np.dot(np.dot(BT, D), B) #local stiffness matrix
            m_local = t*A*density*rawmassmatrix/12 #local mass matrix

            j = 0
            k = len(dofs)-1 #this adds values from local stiffness matrix and local mass matrix values to their global counterparts using dofs
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

    kt_global = stiffness_global*1
    m_global = mass_global*1
    i = 0

    while i < num_nodes: #this removes stiffness and mass matrix rows/columns from supported nodes
        j = 0
        while j < 2:
            if int(supports[i][j]) == 1:
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

    #set up values for dynamical FEM; assume displacement and velocity starts at 0, and forces are from input forces. (basically 'step 1' of newmark's method)
    #beta and gamma are set to ensure simulation is 'stable' (although the simulation isn't stable at the moment, for some reason)

    displacement = np.zeros(num_dofs)
    velocities = np.zeros(num_dofs)
    gamma = 1/2
    beta = 1/4 #cannot equal 0
    invmass = np.linalg.inv(m_global)

    #acceleration is 2nd entry in this list (np.dot(invmass, (forces - np.dot(kt_global, displacement)))) ('step 2' of newmark's method)
    

    return [displacement, velocities, np.dot(invmass, (forces_global - np.dot(kt_global, displacement))), forces_global, kt_global, m_global, invmass, beta, gamma, m_global/(beta*timestep*timestep), num_dofs]


run = True

while run:
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

                        #take result from FEM and put it into variables that are not confined to one function
                        FEMdisplacement = np.copy(bFresult[6][0])
                        FEMvelocity = np.copy(bFresult[6][1])
                        FEMacceleration = np.copy(bFresult[6][2])
                        FEMforces = np.copy(bFresult[6][3])
                        FEMk = np.copy(bFresult[6][4])
                        FEMm = np.copy(bFresult[6][5])
                        FEMinvm = np.copy(bFresult[6][6])
                        FEMbeta = bFresult[6][7]
                        FEMgamma = bFresult[6][8]
                        FEMalpha = np.copy(bFresult[6][9])
                        FEMdofs = bFresult[6][10]
                        
                        
                        points_bookmark = np.copy(points)

                        #points[0][0] = 772

                        #print(points[0][0])
                        #print(points_bookmark[0][0])





                        #print(str(points_bookmark))
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



    if play == 1:

        i = 0
        while i < simulationspeed: #this simulationspeed loop allows the simulation to iterate multiple times over 1 program tick (essentially a 'fast forward' function)

            starttime = time.process_time()

            #step 3 of newmark's method
            kprime = FEMk + FEMalpha
            fprime = FEMforces + np.dot(FEMalpha, (FEMdisplacement + timestep*FEMvelocity + ((1/2) - FEMbeta)*(timestep*timestep)*FEMacceleration))
            newdisplacement = np.linalg.solve(kprime, fprime) #update displacements
            
            time1 = time.process_time()

            #step 4 of newmark's method
            newacceleration = (newdisplacement - FEMdisplacement - timestep*FEMvelocity - timestep*timestep*((1/2) - FEMbeta)*FEMacceleration)/(FEMbeta*timestep*timestep)
            time2 = time.process_time()
            #step 5 of newmark's method
            FEMvelocity = FEMvelocity + timestep*((1 - FEMgamma)*FEMacceleration + FEMgamma*newacceleration)
            time3 = time.process_time()
            #update displacement and acceleration (previous values of acceleration/displacement/velocity are not stored anywhere)
            FEMdisplacement = np.copy(newdisplacement)
            FEMacceleration = np.copy(newacceleration)
            time4 = time.process_time()


            print(str([time1-starttime, time2-time1, time3-time2, time4-time3]))
            i = i + 1

            

        j = 0
        #points = 0#np.zeros((8, len(points)))
        while j < len(points): #update position of nodes; use 'points_bookmark' as the initial position before simulation and 'FEMdisplacement' as displacement from initial
            #points[j][0] = 0
            #points[j][1] = 0
            points[j][0] = points_bookmark[j][0] + FEMdisplacement[2*j]
            points[j][1] = points_bookmark[j][1] + FEMdisplacement[2*j+1]
            #print(str([points_bookmark[j], points[j], FEMdisplacement[2*j], FEMdisplacement[2*j+1]]))
            j = j + 1

        












    redrawGameWindow()


    
pygame.quit()



