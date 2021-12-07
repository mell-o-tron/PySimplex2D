import pygame
from time import sleep
import numpy as np
import math
from pygame.locals import *

import itertools



def simpView(A, B, b, c):
    mb = max(b)
    maxRoot = max(mb / (1 if i == 0 else i) for i in A[b.index(mb)])
    print("doot: "+str(maxRoot))
    winsize = 500
    scale = 200 / maxRoot
    point_size = 5
    constraint_thickness = 2
    axis_thickness = 3
    limit = (3/2) * maxRoot


    constr_colours = [
        (255, 226, 0),
        (254, 48, 111),
        (185, 30, 225),
        (167, 215, 7),
        (221, 39, 129),
        (255, 146, 2)
    ]

    started = False
    done = 0


    solution = [0,0]

    # a1x1 + a2x2 = bi
    # x2 = (bi - a1x1)/a2


    ########################################### PLOTTING FUNCTIONS ###########################################

    def constraint(a1, a2, bi, scale):
        if a2 == 0:
            start = [bi * scale + winsize/2, -limit * -scale + winsize/2]
            end = [bi * scale + winsize/2, limit * -scale + winsize/2]
            return [start, end]
        
        start = [-limit * scale + winsize/2, -scale * (bi-a1 * (-limit))/(a2)+winsize/2]
        end = [limit * scale + winsize/2, -scale * (bi-a1 * (limit))/(a2)+winsize/2]
        return [start, end]

    pygame.init()


    def drawConstraints(A, b, scale):
        for i in range(len(b)):
            v = constraint(A[i][0], A[i][1], b[i], scale)
            
            pygame.draw.line(screen, 
                            constr_colours[i % (len(constr_colours ))],
                            v[0],v[1],
                            constraint_thickness)


    def drawFOb(c, scale):
        pygame.draw.line(screen, 
                        (100, 100, 255), 
                        (winsize/2,winsize/2), 
                        (winsize/2 + c[0] * scale/2, winsize/2 - c[1] * scale/2),
                        constraint_thickness)
        


    def drawAxis():
        p = [[0, winsize/2], 
            [winsize, winsize/2], 
            [winsize/2, 0], 
            [winsize/2, winsize]]
        
        pygame.draw.line(screen, (255, 255, 255),p[0],p[1], axis_thickness)
        pygame.draw.line(screen, (255, 255, 255),p[2],p[3], axis_thickness)

    def drawPolygon(A, b, scale):
        intersections = []
        filteredsections = []
        for v1 in range(len(A)):
            for v2 in range(len(A)):
                if v1 != v2:
                    
                    a1 = A[v1][0]
                    a2 = A[v2][0]
                    b1 = A[v1][1]
                    b2 = A[v2][1]
                    c1 = -b[v1]
                    c2 = -b[v2]
                    
                    if a1 * b2 - a2 * b1 != 0 and a1*b2 - a2*b1 != 0:
                        inters = ((b1 * c2 - b2 * c1)/(a1 * b2 - a2 * b1), (c1*a2 - c2*a1)/(a1*b2 - a2*b1))
                        
                        if inters not in intersections:
                            intersections.append(inters)
        
        for i in intersections:
            if all(A[v][0] * i[0] + A[v][1] *i[1] <= b[v] for v in range(len(A))):
                filteredsections.append((i[0] * scale + winsize/2, i[1] * -scale + winsize/2))


        # So this is a terrible approach, a decent one could be visiting the vertices by moving along the edges.
        for i in itertools.permutations(filteredsections):
            pygame.draw.polygon(screen, (0, 100, 100), i, width = 0)
        
    
    
    ########################################### SIMPLEX ITERATION ###########################################

    def simp_primale(A, B, b, c, first):
        
        print("\n\n")
        print("PRIMAL SIMPLEX")
        if first:
            for i in range(len(B)):
                B[i] -= 1
        AB = []
        AN = []
        bB = []
        N = []
        
        # Riempi AB, bB, AN, N
        for i in range(len(b)):
            if i in B:
                AB.append(A[i])
                bB.append(b[i])
            else:
                AN.append(A[i])
                N.append(i)

        
        print("B: " + str([i+1 for i in B]) + ("\t\tN: " + str([i+1 for i in N])))
        print()
        print("AB:")
        print(np.matrix(AB))
        
        
        ABinv = np.linalg.inv(np.matrix(AB))
        
        x = np.asarray(np.matmul(ABinv, bB))[0]
        yB = np.asarray(np.matmul(c, ABinv))[0]
        
        
        y = []
        j = 0
        for i in range(0, len(b)):
            if i in B:
                y.append(yB[j])
                j+=1
            else:
                y.append(0)
        print()
        print("x: " + str(x) + "\t\ty: " + str(y))

        if all(i >= 0 for i in y):
            print("optimal solutions.")
            pygame.draw.circle(screen, 
                            (0,255,0), 
                            (scale * x[0]+ winsize/2, 
                                -scale * x[1]+ winsize/2), 
                            point_size)
            global solution
            solution = [scale * x[0]+ winsize/2, -scale * x[1]+ winsize/2]
            return -2
        
        print("y is not admissible")
        pygame.draw.circle(screen, 
                        (255,0,0), 
                        (scale * x[0] + winsize/2, 
                            -scale * x[1]+ winsize/2), 
                        point_size)
        h = 100
        for i in range(len(y)):
            if i in B and y[i] < 0:
                h = i
                break
        
        print("Exiting index h: " + str(h + 1))
        
        xi = np.matrix(np.subtract(np.zeros(shape = (2,2)), np.matrix(ABinv)))[:, B.index(h)]

        print("xi: " + str(xi))
        if all(np.matmul(AN, xi) <= 0):
            print("Problem is unlimited")
            return -1
        
        lam = []
        for i in range(len(y)):
            if i in N and all(np.matmul(A[i], xi) > 0):
                lam.append((b[i] - np.matmul(A[i], x))/np.matmul(A[i], xi)[0,0])
            else:
                lam.append(math.inf)
        k = lam.index(min(lam))
        print("Entering index k: "+ str(k + 1))
        B.remove(h)
        B.append(k)
        B.sort()
        print("new base: \n" + str([i+1 for i in B]))
        return B


    ########################################### PYGAME WINDOW ###########################################

    # Set up the drawing window
    screen = pygame.display.set_mode([winsize, winsize])

    running = True
    # Run until the user asks to quit
    if len(c) != 2 or (any( len(i) != 2 for i in A)):
        print("Wrong problem format. Problem must be in R^2")
        running = False

    while running:
        # Did the user click the window close button?
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        # Fill the background with white\
        screen.fill((10, 10, 10))

        # Draw a solid blue circle in the center
        
        drawPolygon(A, b, scale)
        drawAxis()
        drawConstraints(A, b, scale)
        drawFOb(c, scale)
        res = 0
        if done == 0:
            res = simp_primale(A, B, b, c, not started)
            started = True
            
            if res == -1 or res == -2:
                done = 1
            else:
                B = res;
        else:
            pygame.draw.circle(screen, 
                            (0,255,0), 
                            (solution[0],
                                solution[1]), 
                            point_size)
            
        

        # Flip the display
        pygame.display.flip()
        
        while True:
            event = pygame.event.wait()
            if event.type == KEYDOWN and event.key == K_SPACE:
                if done == 1:
                    running = False
                break
            if event.type == pygame.QUIT:
                running = False
                break
        
        
        
        
        
    # Done! Time to quit.
    pygame.quit()
