import pygame
import numpy as np
import math
import sys

pygame.init()

win = pygame.display.set_mode((800,800))
pygame.display.set_caption("PYTHON TEST APPLICATION")
font = pygame.font.SysFont("courier", 12, True)
fontsmall = pygame.font.SysFont("courier", 10, True)

gamestate = 0
clock = pygame.time.Clock() 
user_text = '' 

textbox_coordsx = [0, 0, 0, 0, 0, 0, 0]
textbox_coordsy = [20, 60, 100, 140, 180, 220, 260]

renderingstartx = 200
renderingstarty = 200

renderingoriginx = 0
renderingoriginy = 0

renderingboxsize = 240
renderingboxscale = 20

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
width = 3
height = 3

def drawPolygon(dppoints, dpcolor):
    if len(dppoints) > 2:
        pygame.draw.polygon(win, dpcolor, dppoints, 0)








def constructMeshFromPolygon(fineness):
    #for all points from 0 - 800 (size of screen)
    meshpoints = []
    meshcolor = (0, 0, 255)

    meshsize = math.ceil(800/fineness)
    meshgrid = np.zeros((800, 800))
    meshpoints = []
    meshelements = []

    k = 0
    i = 0
    while i < 800:
        
        j = 0
        while j < 800:
            if win.get_at((i, j)) == meshcolor:
                k = k + 1
                meshgrid[int(i/fineness)][int(j/fineness)] = k
                meshpoints.append([i, j])
                #add to points

            j = j + fineness
        i = i + fineness

    if len(meshpoints) > 0:

        i = 0
        while i < len(meshpoints):
            
            meshpx = int(meshpoints[i][0]/fineness)
            meshpy = int(meshpoints[i][1]/fineness)
            if meshgrid[meshpx+1][meshpy] != 0:
                meshelements.append([i, int(meshgrid[meshpx+1][meshpy]) - 1])

            if meshgrid[meshpx][meshpy+1] != 0:
                meshelements.append([i, int(meshgrid[meshpx][meshpy+1]) - 1])

            if meshgrid[meshpx+1][meshpy+1] != 0:
                meshelements.append([i, int(meshgrid[meshpx+1][meshpy+1]) - 1])

            if meshgrid[meshpx+1][meshpy-1] != 0:
                meshelements.append([i, int(meshgrid[meshpx+1][meshpy-1]) - 1])

            i = i + 1

    return (meshpoints, meshelements)

def drawPoints(dpoints, delements, dpoints_color):
    if len(dpoints) > 0:
        i = 0
        while i < len(dpoints):
            if len(dpoints[i]) > 0:

                drawx = (dpoints[i][0] + renderingoriginx)*renderingboxscale + renderingstartx
                drawy = (dpoints[i][1] + renderingoriginy)*renderingboxscale + renderingstarty

                pygame.draw.rect(win, dpoints_color, (drawx, drawy, width, height))
                text = fontsmall.render(str(i), 1, (255,255,255))
                #text = fontsmall.render(str(i) + ": " + str(points[i][0]) + ", " + str(points[i][1]) + ": " + str(points[i][2]) + ", " + str(points[i][3]), 1, (255,255,255))
                win.blit(text, (drawx+2, drawy-2))
            i = i + 1
        if len(delements) > 0:
            i = 0
            while i < len(delements):
                if len(dpoints[delements[i][0]]) == 0 or len(dpoints[delements[i][1]]) == 0:
                    #delements[i] = []
                    j = 0 #unneeded
                else:
                    drawx = (dpoints[delements[i][0]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                    drawy = (dpoints[delements[i][0]][1] + renderingoriginy)*renderingboxscale + renderingstarty
                    drawx2 = (dpoints[delements[i][1]][0] + renderingoriginx)*renderingboxscale + renderingstartx
                    drawy2 = (dpoints[delements[i][1]][1] + renderingoriginy)*renderingboxscale + renderingstarty

                    pygame.draw.line(win, dpoints_color, [drawx, drawy], [drawx2, drawy2], 3)
                i = i + 1



def redrawGameWindow():
    win.fill((0,0,0))
    #win.blit(bg, (0,0))
    text = font.render('Create node at coords "x, y, forcex, forcey, displacementx, displacementy" (0 for unknown/unconfined):', 1, (255,255,255))
    win.blit(text, (0, 0))

    text = font.render('Create element at nodes "node1, node2, E, A":', 1, (255,255,255))
    win.blit(text, (0, 40))

    text = font.render('Delete node "node":', 1, (255,255,255))
    win.blit(text, (0, 80))

    text = font.render('Delete element "element":', 1, (255,255,255))
    win.blit(text, (0, 120))

    text = font.render('Run Simulation:', 1, (255,255,255))
    win.blit(text, (0, 160))

    text = font.render('Switch to polygon mode:', 1, (255,255,255))
    win.blit(text, (0, 200))

    text = font.render('Convert polygon to mesh:', 1, (255,255,255))
    win.blit(text, (0, 240))

    i = 0
    while i < len(textbox_coordsx):
        if active == i:
            color = color_active 
            input_rect = pygame.Rect(textbox_coordsx[i], textbox_coordsy[i], 200, 20)
            pygame.draw.rect(win, color, input_rect) 
            text_surface = font.render(user_text, True, (255, 255, 255))
            win.blit(text_surface, (textbox_coordsx[i]+3, textbox_coordsy[i]+3)) 
        else: 
            color = color_passive 
            input_rect = pygame.Rect(textbox_coordsx[i], textbox_coordsy[i], 200, 20)
            pygame.draw.rect(win, color, input_rect) 
        clock.tick(60) 
        i = i + 1



    drawPoints(points, elements, (255, 0, 0))
    drawPoints(points2, elements, (0, 255, 0))
    drawPolygon(polygon_points, (0, 0, 255))


    pygame.display.update() 






def simulateFEM():
    
    num_nodes = 0
    num_elements = 0#len(elements)

    coords = []#np.zeros((num_nodes, 2))
    supports = []#np.zeros((num_nodes, 2))
    loads = []#np.zeros((num_nodes, 2))
    node_directions = []
    node_index = []


  
    i = 0
    j = 0
    while i < len(points):
        if len(points[i]) > 0:
            num_nodes = num_nodes + 1
            coords.append([points[i][0], points[i][1]])
            loads.append([points[i][2], points[i][3]])
            supports.append([points[i][4], points[i][5]])
            node_directions.append([j*2, j*2+1])
            node_index.append(j)
            j = j + 1
        i = i + 1

    connectivity = []#np.zeros((num_elements, 2))
    stiffness = []
    dofs = []
    num_dofs = 2*num_nodes
    stiffness_global = np.zeros((num_dofs, num_dofs))

    forces_global = np.zeros(num_dofs)
    #stiffness_inputs = []#np.zeros(num_elements)

    i = 0
    while i < len(elements):
        if len(points[elements[i][0]]) == 0 or len(points[elements[i][1]]) == 0:
            j = 0 #unnecessary
        else:
            num_elements = num_elements + 1
            #find leftmost point of element
            if points[elements[i][0]][0] > points[elements[i][1]][0]:
                connectivity.append([node_index[elements[i][1]], node_index[elements[i][0]]])
            else:
                connectivity.append([node_index[elements[i][0]], node_index[elements[i][1]]])
        
            nodes_for_this_element = (connectivity[len(connectivity)-1][0], connectivity[len(connectivity)-1][1])

            
            xdistance_for_this_element = coords[nodes_for_this_element[0]][0] - coords[nodes_for_this_element[1]][0]
            ydistance_for_this_element = coords[nodes_for_this_element[0]][1] - coords[nodes_for_this_element[1]][1]
            distance_for_this_element = np.sqrt(pow(xdistance_for_this_element, 2) + pow(ydistance_for_this_element, 2))
            stiffness.append(elements[i][2]*elements[i][3]/distance_for_this_element)

            #print([distance_for_this_element, xdistance_for_this_element, ydistance_for_this_element, coords[int(connectivity[i][0])][0], coords[int(connectivity[i][0])][1], coords[int(connectivity[i][1])][0], coords[int(connectivity[i][1])][1]])
            #if loading_state == 1:
                #if simple_stiffness == 0:
                #    cross_sectional_area_for_this_element = float(input("What is the cross sectional area for element " + str(i) + "?"))
                #    youngs_modulus_for_this_element = float(input("What is the young's Modulus for element " + str(i) + "?"))
                #    stiffness_inputs[i] = cross_sectional_area_for_this_element*youngs_modulus_for_this_element/distance_for_this_element
                #else:

                    #stiffness_inputs[i] = simple_stiffness/distance_for_this_element
            
            dofs.append([node_directions[nodes_for_this_element[0]][0], node_directions[nodes_for_this_element[0]][1], node_directions[nodes_for_this_element[1]][0], node_directions[nodes_for_this_element[1]][1]])

            c = xdistance_for_this_element/distance_for_this_element
            s = ydistance_for_this_element/distance_for_this_element

            k_local = stiffness[len(stiffness)-1]*np.array([
              [c ** 2, c * s, -c ** 2, -c * s],
              [c * s, s ** 2, -c * s, -s ** 2],
              [-c ** 2, -c * s, c ** 2, c * s],
              [-c * s, -s ** 2, c * s, s ** 2]
             ])
            #print(k_local)

    
            #add to global matrix
            #coded from first principles instead of using .ix_
            j = 0
            k = len(stiffness)-1
            while j < 4:
                stiffness_global[int(dofs[k][j])][int(dofs[k][0])] += k_local[j][0]
                stiffness_global[int(dofs[k][j])][int(dofs[k][1])] += k_local[j][1]
                stiffness_global[int(dofs[k][j])][int(dofs[k][2])] += k_local[j][2]
                stiffness_global[int(dofs[k][j])][int(dofs[k][3])] += k_local[j][3]
                j = j + 1
    
                
                
        i = i + 1

    kt_global = stiffness_global*1 #why *1??? is it to establish kt_global as its own thing instead of copying the id to the array?? weird.
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
                forces_global[int(node_directions[i][j])] = 0
            else:
                forces_global[int(node_directions[i][j])] = loads[i][j]
            j = j + 1
        i = i + 1

    print("kt_global:")
    print(kt_global)
    print("forces_global:")
    print(forces_global)
    displacement = np.linalg.solve(kt_global, forces_global)
    reaction_forces = np.dot(stiffness_global, displacement)

    i = 0
    while i < len(node_index):
        index = node_index[i]
        points2.append([points[index][0] + displacement[2*i]/1000, points[index][1] + displacement[2*i + 1]/1000])
        i = i + 1

    print("Displacements (mm):")
    print(displacement * 1000)
    print("Reaction Forces (N):")
    print(reaction_forces)


run = True

while run:
    pygame.time.delay(100)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run = False
        #if event.type == pygame.MOUSEBUTTONUP:
        #    pos = pygame.mouse.get_pos()
        #    points.append(pos)
        if event.type == pygame.MOUSEBUTTONDOWN: 
            pos = pygame.mouse.get_pos()
            #points.append(pos)

            i = 0
            while i < len(textbox_coordsx) and i != -7:
                #input_rect = pygame.Rect(textbox_coordsx,textbox_coordsy, textbox_coordsx+60, textbox_coordsy+14)

                if checkCollision(pos[0], pos[1], textbox_coordsx[i],textbox_coordsy[i], textbox_coordsx[i]+200, textbox_coordsy[i]+20):
                    if i == 5:
                        points = []
                        points2 = []
                        elements = []
                    active = i
                    user_text = ''
                    i = -8

                i = i + 1

            if i != -7:
                if active == 5:

                    pos = pygame.mouse.get_pos()
                    polygon_points.append(pos)
                else:
                    active = -1
            if active == 4:
                simulateFEM()
            if active == 6:
                result = constructMeshFromPolygon(10)
                
                points = result[0]
                elements = result[1]

                polygon_points = []


        if event.type == pygame.KEYDOWN and active != -1: 
            if event.key == pygame.K_BACKSPACE: 
                user_text = user_text[:-1] 
            elif event.key == pygame.K_RETURN:
                newtuple = tuple(user_text.split(", "))
                #print (newtuple)
                #user_text = str(newtuple)
                if active == 0:
                    points.append([float(newtuple[0]), float(newtuple[1]), float(newtuple[2]), float(newtuple[3]), float(newtuple[4]), float(newtuple[5])])
                    


                    
                    
                    
                    
                elif active == 1:
                    elements.append([int(newtuple[0]), int(newtuple[1]), float(newtuple[2]), float(newtuple[3])])
                elif active == 2:
                    points[int(newtuple[0])] = []
                elif active == 3:
                    elements[int(newtuple[0])] = []


                # elif active == 2:

                # elif active == 3:

                # elif active == 4:
                user_text = ''
            



            else:
                user_text += event.unicode






    #keys = pygame.key.get_pressed()
    
    # if keys[pygame.K_LEFT]:
    #     x -= vel

    # if keys[pygame.K_RIGHT]:
    #     x += vel

    # if keys[pygame.K_UP]:
    #     y -= vel

    # if keys[pygame.K_DOWN]:
    #     y += vel

    #print (points)
    #print (elements)
    redrawGameWindow()


    
pygame.quit()
