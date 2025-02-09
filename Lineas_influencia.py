import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from math import *
from itertools import accumulate
import plotly.graph_objects as go
from plotly.subplots import make_subplots


class Influence_Lines():
    def __init__(self, Type, nbeams, L, EI, loc_sec, dist_sec):
        self.Type = Type
        self.nbeams = nbeams
        self.L = L
        self.EI = EI
        self.loc_sec = loc_sec
        self.dist_sec = dist_sec
    
    def Local_Matrix(self, EI, L):
        K = np.array([[4, 2],
                      [2, 4]])*(EI/L)
        return K
    
    def General_Matrix(self, nbeams, EI, L):
        ndof = nbeams + 1
        Kg = np.zeros((ndof,ndof))
        for i in range(nbeams):
            K = self.Local_Matrix(EI[i], L[i])
            dof = [i, i+2]
            Kg[dof[0]:dof[1], dof[0]:dof[1]] += K
        return Kg
    
    def General_Forces(self, nbeams, L, xloc, xdist):
        ndof = nbeams + 1
        F = np.zeros((ndof, 1))
        M1 = xdist*((L[xloc-1]-xdist)/L[xloc-1])**2
        M2 = -(L[xloc-1]-xdist)*(xdist/L[xloc-1])**2
        Fu = np.array([[M1],
                       [M2]])
        dof = [xloc-1, xloc+1]
        F[dof[0]:dof[1]] -= Fu
        return F, Fu
    
    def Displacements(self, nbeams, Kg, F):
        u = np.linalg.solve(Kg, F)
        ue = []
        for i in range(nbeams):
            ue.append(np.array(u[i:i+2]))
        return ue
    
    def Local_Forces(self, nbeams, EI, L, ue, xloc, Fu):
        Fe = []
        for i in range(nbeams):
            K = self.Local_Matrix(EI[i], L[i])
            if i == xloc-1:
                Fe.append(Fu + K@ue[i])
            else:
                Fe.append(K@ue[i])
        return Fe
    
    def Shear_Forces(self, nbeams, L, xloc, xdist, Fe):
        ndof = nbeams + 1
        Re = np.zeros((ndof, 1))
        Ve = []
        for i in range(nbeams):
            dof = [i, i+2]
            if i == xloc-1:
                V1 = (np.sum(Fe[i]) - xdist)/L[i] + 1
                V2 = 1 - V1
            else:
                V1 = (np.sum(Fe[i]))/L[i]
                V2 = -V1
            V = np.array([[V1],
                          [V2]])
            Ve.append(V)
            Re[dof[0]:dof[1]] += V
        return Ve, Re
    
    def main(self, Type, nbeams, EI, L, loc_sec, dist_sec):
        Lt = list(accumulate(L))
        M = []; V = []; R = []; loc = []
        Kg = self.General_Matrix(nbeams, EI, L)
        npoints = 40

        # Generador de lista de ubicación
        for i in range(nbeams):
            for j in np.linspace(0, 1, npoints):
                if i == 0:
                    loc.append(L[i]*j)
                else:
                    loc.append(Lt[i-1] + L[i]*j)
        loc = [float(x) for x in loc]
        Vloc = [x for x in loc]
        

        # Generador de linea de influencia
        for i in range(nbeams):
            for j in np.linspace(0, 1, npoints):
                F, Fu= self.General_Forces(nbeams, L, i+1, j*L[i])
                ue = self.Displacements(nbeams, Kg, F)
                Fe = self.Local_Forces(nbeams, EI, L, ue, i+1, Fu)
                Ve, Re = self.Shear_Forces(nbeams, L, i+1, j*L[i], Fe)

                # Lista para momentos flectores y fuerzas cortantes
                Ve_loc = float(Ve[loc_sec-1][0, 0])
                Fe_loc = float(Fe[loc_sec-1][0, 0])
                if i == loc_sec-1 and j < dist_sec/L[i]:
                    Mx = (Ve_loc - 1)*dist_sec - Fe_loc + j*L[i]
                    M.append(Mx)
                    if dist_sec != 0:
                        V.append((Fe_loc + Mx - j*L[i])/dist_sec)
                    else:
                        V.append(Ve_loc)
                    
                else:
                    Mx = Ve_loc*dist_sec - Fe_loc
                    M.append(Mx)
                    if dist_sec != 0:
                        V.append((Fe_loc + Mx)/dist_sec)
                    else:
                        V.append(Ve_loc)

                M = [float(x) for x in M]
                V = [float(x) for x in V]
              
                # Lista para reacciones
                R.append(Re[loc_sec-1][0])
    
        # Fuerzas para la sección solicitada
        if dist_sec != 0 and Type == 'Internal Forces':

            id_find = min(range(len(V) - 1), key=lambda i: abs(V[i] - V[i+1] + 1))
            
            F, Fu= self.General_Forces(nbeams, L, loc_sec, dist_sec)
            ue = self.Displacements(nbeams, Kg, F)
            Fe = self.Local_Forces(nbeams, EI, L, ue, loc_sec, Fu)
            Ve, Re = self.Shear_Forces(nbeams, L, loc_sec, dist_sec, Fe)
            Ve_loc = float(Ve[loc_sec-1][0, 0])
            Fe_loc = float(Fe[loc_sec-1][0, 0])

            Msec = Ve_loc*dist_sec - Fe_loc
            Vsec = (Fe_loc + Msec) / dist_sec
            Lt_loc = [0] + Lt
            Lsec = Lt_loc[loc_sec-1] + dist_sec
            
            M.insert(id_find + 1, Msec)
            loc.insert(id_find + 1, Lsec)

            V[id_find+1:id_find+1] = [Vsec - 1, Vsec]
            Vloc[id_find+1:id_find+1] = [Lsec, Lsec]


        # Ploteo de gráfica de línea de influencia
        fig = plt.figure(figsize=(20,15), constrained_layout=True)
        spec = gridspec.GridSpec(ncols=3, nrows=3, figure=fig)

        ax1 = fig.add_subplot(spec[0, :])
        ax1.plot([0] + Lt, (nbeams + 1)*[0], linewidth=4, color='k')
        ax1.plot([0] + Lt, (nbeams + 1)*[0], marker=6, markersize=20, color='k')
        ax1.set_xlim(-1, Lt[-1]+1)
        # ax1.set_title('Influences Lines of Continuous beam of ' + str(nbeams) + ' spans', size=16)
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        if Type == 'Internal Forces':
            ax1.plot(sum(L[:loc_sec-1]) + dist_sec, [0], marker='|', markersize=40, color='r')
            ax2 = fig.add_subplot(spec[1, :])
            ax3 = fig.add_subplot(spec[2, :])

            ax2.plot(loc, M, 'b')
            ax2.plot([0] + Lt, (nbeams + 1)*[0], 'o',  markersize=8, color='k')
            ax2.grid(color = 'silver')
            ax2.set_xlim(-1, Lt[-1]+1)
            ax2.set_ylabel('Bending Moment (t.m)', size = 20)

            ax3.plot(Vloc, V, 'g')
            ax3.plot([0] + Lt, (nbeams + 1)*[0], 'o',  markersize=8, color='k')
            ax3.grid(color = 'silver')
            ax3.set_xlim(-1, Lt[-1]+1)
            ax3.set_ylabel('Shear Force (t)', size = 20)
            ax3.set_xlabel('Location (m)', size = 20)
        
        elif Type == 'Reaction':
            ax1.plot(sum(L[:loc_sec-1]), [0], marker='|', markersize=40, color='r')
            ax2 = fig.add_subplot(spec[1, :])

            ax2.plot(loc, R, 'b')
            ax2.plot([0] + Lt, (nbeams + 1)*[0], 'o',  markersize=8, color='k')
            ax2.grid(color = 'silver')
            ax2.set_xlim(-1, Lt[-1]+1)
            ax2.set_ylabel('Reaction (t)')
            ax2.set_xlabel('Location (m)')

        return fig