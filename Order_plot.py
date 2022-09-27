from plotly.figure_factory import create_trisurf
from colorsys import hsv_to_rgb
import pandas as pd
import plotly.graph_objects as go
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


class Order_plot():
    
    def __init__(self, hkl,I):
        h=np.array([x[0] for x in hkl])
        k=np.array([x[1] for x in hkl])
        l=np.array([x[2] for x in hkl])
        I=np.array(I)
        
        ###Scale?
        I=np.array([x**1.5 for x in I])*50
        
        
        sigma, n = 0.2, 20000
        xyzs = self.fibonacci_sphere(n)

        grids = np.zeros([n, 3])
        grids[:, :2] = self.xyz2sph(xyzs)
        
        #####EDITED
        # data = np.loadtxt('stereograph.txt', skiprows=1)
        # h, k, l, I = data[:, 0], data[:, 1], data[:, 2], data[:, 8]
        # df=pd.read_csv('data.csv')
        # h=np.array(df['h'])
        # k=np.array(df['k'])
        # l=np.array(df['l'])
        # I=np.array(df['I'])
        # I=np.array([x**1.5 for x in I])*50


        pts = []
        for i in range(len(h)):
            p, r = self.hkl2tp(h[i], k[i], l[i])
            pts.append([p, r, I[i]])
        pts = np.array(pts)

        vals = self.calculate_density(pts, xyzs, sigma=sigma)


        valss=(vals+100)

        valss/=valss.max()
        for i,x in enumerate(valss):
            xyzs[i]*=np.abs(x)

        phi=[]
        rho=[]
        for row in xyzs:
            r,p=self.hkl2tp(row[0],row[1],row[2])
            phi.append(p)
            rho.append(r)
        phi=np.array(phi)
        rho=np.array(rho)
        # phi=pts[:,1]
        # rho=pts[:,0]
        # I=pts[:,2]
        # print(len(simplices))
        x=xyzs[:,0]
        y=xyzs[:,1]
        z=xyzs[:,2]
        points2D=np.vstack([phi,rho]).T
        tri=Delaunay(points2D)
        simplices=tri.simplices
        trisurf=create_trisurf(x=x,y=y,z=z,colormap='Jet', simplices=simplices,plot_edges=False,color_func=self.color_func)

        self.colorscale = [
            [0, "rgb(84,48,5)"],
            [1, "rgb(84,48,5)"],
        ]

        xy_a=[x_ for i,x_ in enumerate(x) if np.abs(z[i])<2e-2]
        xy_b=[z_ for i,z_ in enumerate(y) if np.abs(z[i])<2e-2]
        xz_a=[x_ for i,x_ in enumerate(x) if np.abs(y[i])<2e-2]
        xz_b=[z_ for i,z_ in enumerate(z) if np.abs(y[i])<2e-2]
        yz_a=[x_ for i,x_ in enumerate(y) if np.abs(x[i])<2e-2]
        yz_b=[z_ for i,z_ in enumerate(z) if np.abs(x[i])<2e-2]
        
        fig,ax=plt.subplots(nrows=1,ncols=3,figsize=(15,5))
        # plt.gca().set_aspect('equal', adjustable='box')
        ax[0].scatter(xy_a,xy_b,color='r')
        ax[0].set_xlim(-1,1)
        ax[0].set_ylim(-1,1)
        ax[0].set_title('X-Y')
        ax[0].set_xlabel('X')
        ax[0].set_ylabel('Y')
        ax[1].scatter(xz_a,xz_b,color='r')
        ax[1].set_title('X-Z')
        ax[1].set_xlabel('X')
        ax[1].set_ylabel('Z')
        ax[1].set_xlim(-1,1)
        ax[1].set_ylim(-1,1)
        ax[2].scatter(yz_a,yz_b,color='r')
        ax[2].set_title('Y-Z')
        ax[2].set_xlabel('Y')
        ax[2].set_ylabel('Z')
        ax[2].set_xlim(-1,1)
        ax[2].set_ylim(-1,1)
        plt.show()




        layout = go.Layout(scene=dict(aspectmode='data',annotations=self.get_axis_names()))
        # fig=go.Figure(data=marker_data1, layout=layout)
        fig=go.Figure(data=trisurf, layout=layout)
        fig.add_trace(go.Scatter3d(x = [1.1,0,0], y = [0,1.1,0], z=[0,0,1.1], mode="text", text = ['a','b','c']))

        self.add_axis_arrows(fig)

        fig.update_scenes(camera_projection_type='orthographic')
        fig.show()
    
    def color_func(self,x,y,z):
        mag=np.sqrt(x**2 + y**2 + z**2)
        return np.floor(mag*255.9999)

    def calculate_density(self,pts, xyzs, sigma=0.01):
        """
        calculate the projected density on the unit sphere
        """
        vals = np.zeros(len(xyzs))
        pi = np.pi
        for pt in pts:
            t0, p0, h = pt
            x0, y0, z0 = np.sin(t0)*np.cos(p0), np.sin(t0)*np.sin(p0), np.cos(t0)
            dst = np.linalg.norm(xyzs - np.array([x0, y0, z0]), axis=1)
            vals += h*np.exp(-(dst**2/(2.0*sigma**2)))
        return vals

    def hkl2tp(self,h, k, l):
        #convert hkl to theta and phi
        mp = [h,k,l]
        r = np.linalg.norm(mp)

        theta = np.arctan2(mp[1],mp[0])
        phi = np.arccos(mp[2]/r)

        #return theta, phi
        return phi, theta

    def fibonacci_sphere(self,samples=1000):
        """
        Sampling the sphere grids

        Args:
            samples: number of pts to generate

        Returns:
            3D points array in Cartesian coordinates
        """
        points = []
        phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
            radius = np.sqrt(1 - y * y)  # radius at y
            theta = phi * i  # golden angle increment
            x = np.cos(theta) * radius
            z = np.sin(theta) * radius
            points.append((x, y, z))

        return np.array(points)


    def xyz2sph(self,xyzs, radian=True):
        """
        convert the vectors (x, y, z) to the sphere representation (theta, phi)

        Args:
            xyzs: 3D xyz coordinates
            radian: return in radian (otherwise degree)
        """
        pts = np.zeros([len(xyzs), 2])   
        for i, r_vec in enumerate(xyzs):
            r_mag = np.linalg.norm(r_vec)
            theta0 = np.arccos(r_vec[2]/r_mag)
            if abs((r_vec[2] / r_mag) - 1.0) < 10.**(-8.):
                theta0 = 0.0
            elif abs((r_vec[2] / r_mag) + 1.0) < 10.**(-8.):
                theta0 = np.pi

            if r_vec[0] < 0.:
                phi0 = np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif 0. < r_vec[0] and r_vec[1] < 0.:
                phi0 = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif 0. < r_vec[0] and 0. <= r_vec[1]:
                phi0 = np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] == 0. and 0. < r_vec[1]:
                phi0 = 0.5 * np.pi
            elif r_vec[0] == 0. and r_vec[1] < 0.:
                phi0 = 1.5 * np.pi
            else:
                phi0 = 0.
            pts[i, :] = [theta0, phi0]
        if not radian:
            pts = np.degree(pts)

        return pts



    def get_arrow(self,axisname="x"):

        # Create arrow body
        body = go.Scatter3d(
            marker=dict(size=1, color=self.colorscale[0][1]),
            line=dict(color=self.colorscale[0][1], width=3),
            showlegend=False,  # hide the legend
        )

        head = go.Cone(
            sizeref=0.1,
            autocolorscale=None,
            colorscale=self.colorscale,
            showscale=False,  # disable additional colorscale for arrowheads
            hovertext=axisname,
        )
        for ax, direction in zip(("x", "y", "z"), ("u", "v", "w")):
            if ax == axisname:
                body[ax] = -1,1
                head[ax] = [1]
                head[direction] = [1]
            else:
                body[ax] = 0,0
                head[ax] = [0]
                head[direction] = [0]

        return [body, head]


    def add_axis_arrows(self,fig):
        for ax in ("x", "y", "z"):
            for item in self.get_arrow(ax):
                fig.add_trace(item)

    def get_annotation_for_ax(self,ax):
        d = dict(showarrow=False, text=ax, xanchor="left", font=dict(color="#1f1f1f"))

        if ax == "a":
            d["x"] = 1.1
            d["y"] = 0
            d["z"] = 0
        elif ax == "b":
            d["x"] = 0
            d["y"] = 1.1 
            d["z"] = 0
        else:
            d["x"] = 0
            d["y"] = 0
            d["z"] = 1.1 

        if ax in {"a", "b"}:
            d["xshift"] = 15

        return d


    def get_axis_names(self):
        return [self.get_annotation_for_ax(ax) for ax in ("a", "b", "c")]
        
if __name__ == "__main__":
    df=pd.read_csv('data.csv')
    hkl=[]
    for i in range(len(df['h'])):
        hkl.append([df['h'][i],df['k'][i],df['l'][i]])
    I=df['I']
    p=Order_plot(hkl,I)