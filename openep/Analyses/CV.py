import sys
import os
import math
import pyvista as pv
import numpy as np
import scipy.io
import scipy
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
import pdb
from scipy.interpolate import Rbf
from vedo import *
import pyvistaqt 
from sklearn import metrics

def func(coef, x1, x2, x3, m):
        ''' This the function we use to fit a surface to the activations in plane 
        fitting method'''
        a = coef[0]
        b = coef[1]
        c = coef[2]
        d = coef[3]
        e = coef[4]
        f = coef[5]
        g = coef[6]
        h = coef[7]
        q = coef[8]
        l = coef[9]
        #print(np.linalg.norm(m-(a*x1**2 + b*x2**2 + c*x3**2 + d*x1*x2 + e*x1*x3 + f*x2*x3 + g*x1 +h*x2 + q*x3 + l)))
        return np.linalg.norm(m-(a*x1**2 + b*x2**2 + c*x3**2 + d*x1*x2 + e*x1*x3 + f*x2*x3 + g*x1 +h*x2 + q*x3 + l))     
class CV:
    
    def Triangulation(self,egmX,LAT,data_type=''):
        ''' This function calculates local atrial conduction velocity from recorded 
                electro-anatomical data, using triangulation method. (C.D. Cantwell et.al. 2015
                ''Techniques for automated local activation time annotation and conduction velocity estimation in cardiac mapping")
            egmX: Recorded Electrode locations (x,y,z) 
            LAT : Recorded local activation local activation time'''
        ##### initiation ######
        cv = []
        min_theta = 30
        max_elec_distance = 10
        min_elec_distance = 1.5
        min_lat_difference = 2
        cv_centers_x = []
        cv_centers_y = []
        cv_centers_z = []
        cv_centers_id = []
        alpha = 5
        
        ######## Create a triangulation mesh from the recording electrodes
        egmX_point_cloud = pv.PolyData(egmX)
        print('---> Triangulation process')
        surf = egmX_point_cloud.delaunay_3d(alpha=alpha)
        print('---> Triangulation process')
        print('---> CV calculation (Triangulation method started ...)')
        for num_point in range(len(surf.cells_dict[5])):
            
            vtx_id = []
            lat = []
            
            for i in range (3):
                vtx_id.append(surf.cells_dict[5][num_point][i])
                lat.append(LAT[surf.cells_dict[5][num_point][i]])
            
            id_lat_sorted  = np.argsort(lat)

            O = [egmX[int(vtx_id[id_lat_sorted[0]])],lat[id_lat_sorted[0]]]
            A = [egmX[int(vtx_id[id_lat_sorted[1]])],lat[id_lat_sorted[1]]]
            B = [egmX[int(vtx_id[id_lat_sorted[2]])],lat[id_lat_sorted[2]]]
            
            OA = np.sqrt(sum(np.power(np.subtract(O[0],A[0]),2)))
            OB = np.sqrt(sum(np.power(np.subtract(O[0],B[0]),2)))
            AB = np.sqrt(sum(np.power(np.subtract(A[0],B[0]),2)))
            
            tOA = A[1] - O[1]
            tOB = B[1] - O[1]
            
            theta = np.arccos((np.power(OA,2)+ np.power(OB,2) - np.power(AB,2))/(2 * OA * OB))
            # check if the conditions are meet to accept the triangle set as a viable acceptable one
            if (math.degrees(theta) >= min_theta and OA >= min_elec_distance and OA <= max_elec_distance and
                OB >= min_elec_distance and OB <= max_elec_distance and tOA >= min_lat_difference and tOB >= min_lat_difference):

                alpha = np.arctan((tOB * OA - tOA * OB * np.cos(theta)) / (tOA * OB * np.sin(theta)))
                cv_temp = (OA/tOA) * np.cos(alpha)
                cv.append(cv_temp)
                cv_centers_x.append(O[0][0])
                cv_centers_y.append(O[0][1])
                cv_centers_z.append(O[0][2])
                cv_centers_id.append(vtx_id[id_lat_sorted[0]])
                #if cv_temp >= 0.2 and cv_temp <=2:
                #    cv.append(cv_temp)
                #    cv_centers_x.append(O[0][0])
                #    cv_centers_y.append(O[0][1])
                #    cv_centers_z.append(O[0][2])
                #    cv_centers_id.append(vtx_id[id_lat_sorted[0]])
            #else:
            #    cv_temp = (OA/tOA) * np.cos(alpha)
            #    cv_centers_x.append(O[0][0])
            #    cv_centers_y.append(O[0][1])
            #    cv_centers_z.append(O[0][2])
            #    cv_centers_id.append(vtx_id[id_lat_sorted[0]])
            #    cv.append(cv_temp)
            ####### Create a triangulation mesh from egmX pointcloud #####
        cv_centers = [cv_centers_x, cv_centers_y, cv_centers_z]
        print('---> CV calculation ended')
        return cv, cv_centers, cv_centers_id
    
    def plane_fitting(self,egmX, LAT):
            
        cv = []
        cv_centers_x = []
        cv_centers_y = []
        cv_centers_z = []
        cv_centers_id = []
        cent = np.random.random((len(egmX), 3))
        direction = np.random.random((len(egmX), 3))
        
        tree = KDTree(egmX, leaf_size=5) 
        all_nn_indices = tree.query_radius(egmX, r=10)
        all_nns = [[egmX[idx] for idx in nn_indices] for nn_indices in all_nn_indices]
        
        print(len(egmX))
        for k in range(len(egmX)):
            print(k,end='\r')
            if len(all_nns[k]) >= 5:
                cv_data = np.zeros(shape= (len(all_nns[k]) + 1,4))
                cv_data[0][:3] = egmX[k]
                cv_data[0][3] = LAT[k]
                for i in range(len(all_nns[k])):
                    cv_data[i+1][:3] = all_nns[k][i]
                    cv_data[i+1][3] = LAT[all_nn_indices[k][i]]
                    
                x1 = cv_data[:,0]
                x2 = cv_data[:,1]
                x3 = cv_data[:,2]
                m = cv_data[:,3]
                #mesh_original = pv.PolyData(Points,Face)
                #p2.add_points(cv_data[:,:3], render_points_as_spheres=True, point_size=5.0, color = 'yellow')
                #p2.add_mesh(mesh_original,opacity=0.4)  
                #p2.show()
                
                res = minimize(func, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], args=(x1, x2, x3, m))
                a, b, c, d,e, f, g, h, q, l= res.x
                #print(a*x1**2+b*x2**2+c*x3**2+d*x1*x2+e*x1*x3+f*x2*x3+g*x1+h*x2+q*x3+l -m)
                
                
                x = egmX[k,0]
                y = egmX[k,1]
                z = egmX[k,2]
                
                tx = 2*a*x + d*y + e*z + g
                ty = 2*b*y + d*x + f*z + h
                tz = 2*c*z + e*x + f*y + q
                
                vx = tx/(tx**2 + ty**2 + tz**2)
                vy = ty/(tx**2 + ty**2 + tz**2)
                vz = tz/(tx**2 + ty**2 + tz**2) 
                
                    
                cent[k,0] = x
                cent[k,1] = y
                cent[k,2] = z
                
                direction[k,0] = vx
                direction[k,1] = vy
                direction[k,2] = vz
                
                cv_temp = np.sqrt(np.sum(np.power(vx,2)+np.power(vy,2)+np.power(vz,2)))
                direction[k,0] = vx 
                direction[k,1] = vy
                direction[k,2] = vz 
                
                if cv_temp >=0.2 and cv_temp <=2:
                    cv.append(cv_temp)
                    cv_centers_x.append(x) 
                    cv_centers_y.append(y)
                    cv_centers_z.append(z)
                    cv_centers_id.append(k)
                else:
                    
                    cv.append(cv_temp)
                    cv_centers_x.append(x) 
                    cv_centers_y.append(y)
                    cv_centers_z.append(z)
            else:
                x = egmX[k,0]
                y = egmX[k,1]
                z = egmX[k,2]
                cv.append(np.nan)
                cv_centers_x.append(x) 
                cv_centers_y.append(y)
                cv_centers_z.append(z)
                cv_centers_id.append(k)
            cv_centers = [cv_centers_x, cv_centers_y, cv_centers_z]    

        return cv, cv_centers, cv_centers_id

    def RBF_method(self, egmX, LAT, points, face,Fiber_data='',):
        
        #rbfi = Rbf(egmX[:,0],egmX[:,1],egmX[:,2],LAT,kernel='gaussian',epsilon=0.5)
        rbfi = Rbf(egmX[:,0],egmX[:,1],egmX[:,2],LAT)
        lat_i = []
        for i in range(len(points)):
            lat_i.append(rbfi(points[i][0], points[i][1], points[i][2]))
        
        mesh_original_RBF = pv.PolyData(points,face)
        mesh_original_RBF['values'] = lat_i 
        deriv = mesh_original_RBF.compute_derivative('values')
        df = deriv['gradient']
        
        cvx =  deriv['gradient'][:,0]/np.sum(np.power( deriv['gradient'],2),axis=1)
        cvy =  deriv['gradient'][:,1]/np.sum(np.power( deriv['gradient'],2),axis=1)
        cvz =  deriv['gradient'][:,2]/np.sum(np.power( deriv['gradient'],2),axis=1)
        cv = np.sqrt(np.power(cvx,2) + np.power(cvy,2) + np.power(cvz,2))
        
        cv_data = []
        cv_x = []
        cv_y = []
        cv_z = []
        cv_center_id =[]
        cv_centers_x = []
        cv_centers_y = []
        cv_centers_z = []
        tree = KDTree(points, leaf_size=5) 
        dist, ind = tree.query(egmX, k=1)
        for k in range(len(ind)):
            
            if cv[ind[k][0]] > 0 and cv[ind[k][0]] <= 10:
                cv_data.append(cv[ind[k][0]])
                cv_x.append(cvx[ind[k][0]])
                cv_y.append(cvy[ind[k][0]])
                cv_z.append(cvz[ind[k][0]])
                cv_center_id.append(ind[k][0])
                #cv_centers_x.append(points[ind[k][0]][0])
                #cv_centers_y.append(points[ind[k][0]][1])
                #cv_centers_z.append(points[ind[k][0]][2])
                cv_centers_x.append(egmX[k][0])
                cv_centers_y.append(egmX[k][1])
                cv_centers_z.append(egmX[k][2])
                
                
            else:
                cv_data.append(np.nan) 
                cv_centers_x.append(points[ind[k][0]][0])
                cv_centers_y.append(points[ind[k][0]][1])
                cv_centers_z.append(points[ind[k][0]][2])
                if 0:
                    cv_data.append(cv[ind[k][0]])
                    cv_center_id.append(ind[k][0])
                    cv_centers_x.append(points[ind[k][0]][0])
                    cv_centers_y.append(points[ind[k][0]][1])
                    cv_centers_z.append(points[ind[k][0]][2])
                    cv_x.append(cvx[ind[k][0]])
                    cv_y.append(cvy[ind[k][0]])
                    cv_z.append(cvz[ind[k][0]])
                
        direction = [cv_x,cv_y,cv_z]
        cv_centers = [cv_centers_x,cv_centers_y,cv_centers_z]
        
            
        
        
        return cv_data, cv_centers, direction, cv_center_id

    def visualization(self,lat_data = '',
                    cv_data = '',
                    mesh = '', 
                    egmX = '', 
                    method = '',
                    file_name = '',
                    direction='', 
                    direction_centers = '', 
                    show_histogram= False,
                    show_cv_fibrosis = True,
                    interpolation_radius = 5):
        
        boring_cmap = plt.cm.get_cmap("jet", 256)
        boring_cmap_r = plt.cm.get_cmap("jet_r", 256)
        boring_cmap_2 = plt.cm.get_cmap("Paired", 4)
        
        boring_cmap_2.colors[0][0]= 34/255 #,237,222#1
        boring_cmap_2.colors[0][1]=94/255
        boring_cmap_2.colors[0][2]=168/255
        boring_cmap_2.colors[1][0]= 62/255 #,190,133#0
        boring_cmap_2.colors[1][1]=182/255 #1
        boring_cmap_2.colors[1][2]=196/255 #1
        boring_cmap_2.colors[2][0]=253/255
        boring_cmap_2.colors[2][1]=141/255
        boring_cmap_2.colors[2][2]= 60/255
        boring_cmap_2.colors[3][0]=217/255
        boring_cmap_2.colors[3][1]=71/255
        boring_cmap_2.colors[3][2]=1/255
       
        



        pv.rcParams['transparent_background'] = True
        
        
        
        if method == 'LAT':
            #pl = pv.Plotter(off_screen=True)
            pl = pv.Plotter()
            egmX_point_cloud = pv.PolyData(egmX)
            egmX_point_cloud['LAT[ms]'] = lat_data
            interpolated_mesh_lat = mesh.interpolate(egmX_point_cloud, radius=interpolation_radius)
            pl.add_mesh(egmX_point_cloud,scalars="LAT[ms]",render_points_as_spheres=True, 
                        point_size=20, cmap = boring_cmap_r)
            pl.add_mesh(interpolated_mesh_lat,scalars="LAT[ms]", clim = [np.min(lat_data), np.max(lat_data)], 
                        below_color='brown', above_color = 'magenta', cmap = boring_cmap_r)
            if direction != '' and direction_centers != '':
                pl.add_arrows(direction_centers, direction, mag=2, color=[90, 90, 90])
            #here
            pl.show(screenshot=f'/Users/aligharaviri/Desktop/CV_Data/Figures_picsLAT_{file_name}.png')
        if method == 'CV':
            pl = pv.Plotter(shape=(1,2),off_screen=False)
            sargs = dict(
                    title_font_size=32,
                    label_font_size=32,
                    shadow=True,
                    n_labels=11,
                    italic=True,
                    fmt="%.1f",
                    font_family="arial",
                    color='black'
                    )      
            egmX_point_cloud = pv.PolyData(egmX)
            egmX_point_cloud['CV[m/s]'] = cv_data
            ######## Ver 1
            #interpolated_mesh_cv = mesh.interpolate(egmX_point_cloud,radius=5 ,n_points=2 ,sharpness=2)
            ######## Ver 1
            interpolated_mesh_cv = mesh.interpolate(egmX_point_cloud,radius=10 ,n_points=None,null_value=np.nan)
            
            pl.subplot(0,0)
            
            pl.add_mesh(egmX_point_cloud,scalars='CV[m/s]',render_points_as_spheres=True, 
                        point_size=20, cmap = boring_cmap,clim=[0,1.5],below_color='navy', above_color = 'darkred')
            pl.add_mesh(interpolated_mesh_cv,scalars='CV[m/s]', clim=[0,1.5],
                        below_color='navy', above_color = 'darkred', cmap = boring_cmap,scalar_bar_args=sargs)
            
            if direction != '':
                pl.add_arrows(egmX, np.transpose(direction), mag=4, color='black')#color=[90, 90, 90])
            pl.camera.roll += 20
            pl.camera.elevation -=10
            pl.camera.azimuth = 10
            pl.camera.zoom (1.2)
            pl.subplot(0,1)
            
            pl.add_mesh(egmX_point_cloud,scalars='CV[m/s]',render_points_as_spheres=True, 
                        point_size=20, cmap = boring_cmap,clim=[0,1.2],below_color='navy', above_color = 'darkred')
            pl.add_mesh(interpolated_mesh_cv,scalars='CV[m/s]', clim=[0,1.5],
                        below_color='navy', above_color = 'darkred', cmap = boring_cmap,scalar_bar_args=sargs)
            if direction != '':
                pl.add_arrows(egmX, np.transpose(direction), mag=4, color='black')#color=[90, 90, 90])
            #here
            pl.camera.azimuth += 190
            pl.camera.roll -=30
            pl.camera.elevation -=40
            pl.camera.zoom (1.2)
            #pl.show(screenshot=f'/Users/aligharaviri/Desktop/CV_Data/Figures_pics/{file_name}.png')
            pl.show()
            #pl.close()

        if method== 'CV' and show_histogram:
            plt.rcParams['font.size'] = '14'
            plt.rcParams['font.family'] = 'arial'
            bin_width = 0.01
            plt.hist(cv_data, bins=np.arange(np.nanmin(cv_data), np.nanmax(cv_data) + bin_width, bin_width))
            
            plt.xlabel('CV (m/s)')
            plt.ylabel('Number of Events')
            
            #plt.savefig(f'/Volumes/Simulations/CV_Data/Figures_pics/hist_{file_name}.pdf',bbox_inches='tight',transparent=True, dpi=400)
            #here
            plt.show()
            #plt.close()
        if method == 'CV' and show_cv_fibrosis:
            
            interpolated_mesh_cv_face = interpolated_mesh_cv.point_data_to_cell_data()
            correlation_data = []
            
            low_cv_mask = interpolated_mesh_cv_face['CV[m/s]'] <= 0.3
            high_cv_mask = interpolated_mesh_cv_face['CV[m/s]'] > 0.3
            fibrosis_mask =  mesh['regions'] == 34
            normal_mask =  mesh['regions'] != 34
            low_cv_fibrosis = np.multiply(low_cv_mask,fibrosis_mask)
            high_cv_no_fibrosis = np.multiply(high_cv_mask,normal_mask)
            low_cv_no_fibrosis = np.multiply(low_cv_mask,normal_mask)
            high_cv_fibrosis = np.multiply(high_cv_mask,fibrosis_mask)
            correlation_data = (np.multiply(low_cv_fibrosis,3) + np.multiply(high_cv_fibrosis,2) +  np.multiply(low_cv_no_fibrosis,1) + np.multiply(high_cv_no_fibrosis,0))
            accuracy =  np.sum(np.multiply(low_cv_fibrosis,1) +  np.multiply(high_cv_no_fibrosis,1)) / len(interpolated_mesh_cv_face['CV[m/s]'])
            print('Accuracy= ',accuracy)
            mesh['correlation'] = correlation_data
            pl = pv.Plotter(off_screen=True)
            pl.add_mesh(mesh,scalars='correlation', cmap = boring_cmap_2)
            pl.show(screenshot=f'/Users/aligharaviri/Desktop/CV_Data/Figures_pics/correlation_{file_name}.png')
            pl.close()
            cv_not_nan = interpolated_mesh_cv_face['CV[m/s]'][~np.isnan(interpolated_mesh_cv_face['CV[m/s]'])]
            region_not_nan = fibrosis_mask[~np.isnan(interpolated_mesh_cv_face['CV[m/s]'])]
            reg_id = np.multiply(region_not_nan,1)
        elif method == 'CV' and not show_cv_fibrosis:
            reg_id = []
            cv_not_nan = []
            #cv_corr = []
            #for i in range(len(mesh['regions'])):
            #    print(i,'end:\r')
            #    if mesh['regions'][i] == 34: 
            #        if ~np.isnan(interpolated_mesh_cv_face['CV[m/s]'][i]):
            #            reg_id.append(1)
            #            cv_corr.append(interpolated_mesh_cv_face['CV[m/s]'][i])
            #    else:
            #        if ~np.isnan(interpolated_mesh_cv_face['CV[m/s]'][i]):
            #            reg_id.append(0)
            #            cv_corr.append(interpolated_mesh_cv_face['CV[m/s]'][i])
                        
           
        if method == 'CV':
            return  interpolated_mesh_cv['CV[m/s]'], reg_id, cv_not_nan
    
    def gradient_method(self, egmX, LAT, points, face, Fiber_data = ''):
        '''#sampling_points = Points[sampling_points_ID]'''
        
        mesh_original = pv.PolyData(points,face)
        mesh_original['LAT[ms]'] = LAT
        deriv = mesh_original.compute_derivative('LAT[ms]')
        df = deriv['gradient']
        d_lat_x =  deriv['gradient'][:,0]#/np.sum(np.power( deriv['gradient'],2),axis=1)
        d_lat_y =  deriv['gradient'][:,1]#/np.sum(np.power( deriv['gradient'],2),axis=1)
        d_lat_z =  deriv['gradient'][:,2]#/np.sum(np.power( deriv['gradient'],2),axis=1)
        #cv_data = 1 / (np.power(d_lat_x,2) + np.power(d_lat_y,2) +  np.power(d_lat_z,2))
        cvx = d_lat_x/np.sum(np.power( deriv['gradient'],2),axis=1)
        cvy = d_lat_y/np.sum(np.power( deriv['gradient'],2),axis=1)
        cvz = d_lat_z/np.sum(np.power( deriv['gradient'],2),axis=1)
        cv_data = np.sqrt(np.power(cvx,2) + np.power(cvy,2) + np.power(cvz,2))#/np.sum(np.power( deriv['gradient'],2),axis=1)
        cv_data_sampled = cv_data[egmX]

        return cv_data, cv_data_sampled
        
    def calculate_divergence (self, egmX, LAT, points, face, sampling_data,collision_threshold=-1,focal_threshold=1,interpolation_radius=5):
        boring_cmap = plt.cm.get_cmap("jet", 256)
        boring_cmap_r = plt.cm.get_cmap("jet_r", 256)
        egmX_point_cloud = pv.PolyData(egmX)
        egmX_point_cloud['LAT[ms]'] = LAT 
        mesh_original = pv.PolyData(points,face)
        interpolated = mesh_original.interpolate(egmX_point_cloud, radius=interpolation_radius,n_points=None)
        pl = pv.Plotter()
        pl.add_mesh(interpolated, cmap=boring_cmap_r)
        #pl.show()
        pl.close()
        deriv = interpolated.compute_derivative('LAT[ms]')
        
        cv_x =  deriv['gradient'][:,0]/np.sum(np.power( deriv['gradient'],2),axis=1)
        cv_y =  deriv['gradient'][:,1]/np.sum(np.power( deriv['gradient'],2),axis=1)
        cv_z =  deriv['gradient'][:,2]/np.sum(np.power( deriv['gradient'],2),axis=1)
        
        direction = np.ndarray(shape=(len(cv_x),3))
        cv_direction = np.ndarray(shape=(len(cv_x),3))
        for i in range(len(cv_x)):
            direction[i][0]=cv_x[i] /np.sqrt(np.sum(np.power(cv_x[i],2) + np.power(cv_y[i],2) + np.power(cv_z[i],2)))
            direction[i][1]=cv_y[i] /np.sqrt(np.sum(np.power(cv_x[i],2) + np.power(cv_y[i],2) + np.power(cv_z[i],2)))
            direction[i][2]=cv_z[i] /np.sqrt(np.sum(np.power(cv_x[i],2) + np.power(cv_y[i],2) + np.power(cv_z[i],2)))
            cv_direction[i][0]=cv_x[i]
            cv_direction[i][1]=cv_y[i]
            cv_direction[i][2]=cv_z[i]
            
        mesh_original['activation_direction'] = cv_direction
        div = mesh_original.compute_derivative(scalars='activation_direction',divergence=True)
        sargs = dict(
                    title_font_size=32,
                    label_font_size=32,
                    shadow=True,
                    n_labels=11,
                    italic=True,
                    fmt="%.1f",
                    font_family="arial",
                    color='black'
                    )      
        pv.rcParams['transparent_background'] = True         
        pl = pv.Plotter(shape=(1,2))
        pl.subplot(0,0)
        pl.add_mesh(div.cell_data_to_point_data(), scalars='divergence',cmap=plt.cm.get_cmap("PRGn", 9),clim=[-1.5,1.5],below_color=[63, 1, 44],above_color=[1, 50,32],
                    scalar_bar_args=sargs)
        pl.add_arrows(points[sampling_data],direction[sampling_data],mag=2,color=[90, 90, 90])
        pl.set_background('white')
        pl.camera.roll += 20
        pl.camera.elevation -=10
        pl.camera.azimuth = 10
        pl.camera.zoom (1.5)
        
        pl.subplot(0,1)
        pl.add_mesh(div.cell_data_to_point_data(), scalars='divergence',cmap=plt.cm.get_cmap("PRGn", 9),clim=[-1.5,1.5],below_color=[63, 1, 44],above_color=[1, 50,32],
                    scalar_bar_args=sargs)
        pl.add_arrows(points[sampling_data],direction[sampling_data],mag=2,color=[90, 90, 90])
        pl.set_background('white')
        pl.camera.azimuth += 190
        pl.camera.roll -=30
        pl.camera.elevation -=40
        pl.camera.zoom (1.5)
        #pl.show(screenshot='/Volumes/Simulations/CV_Data/Figures_pics/foo.png')
        #pl.show()
        pl.close()
    
        div_value = []
        for i in range(len(div['divergence'])):
            if div['divergence'][i] < collision_threshold  or div['divergence'][i] > focal_threshold :
                div_value.append(1)
                
            else:
                div_value.append(0)
        mesh_divergence = pv.PolyData(points,face)
        mesh_divergence ['div'] = div_value
        divergence_per_point = mesh_divergence.cell_data_to_point_data()
        pl = pv.Plotter()
        pl.add_mesh(divergence_per_point,cmap=plt.cm.get_cmap("jet", 2))
        pl.add_arrows(points[sampling_data],direction[sampling_data],mag=2,color=[90, 90, 90])
        #here
        #pl.show()
        pl.close()
        
        
        return direction, mesh_original, div_value

    def exclude_collision_points(self, conduction_velocity_centers = '', conduction_velocity= '', 
                                divergence_value = '', upper_threshold = '',mesh_points = ''):
        tree = KDTree(mesh_points, leaf_size=2) 
        ind = tree.query_radius(conduction_velocity_centers, r=0.5)    
        non_colision_points = []
        collision_detection = []   
        cv_data = []
        for i in range(len(ind)):
            print(i,end='\r')
            for k in range(len(ind[i])):
                if divergence_value[ind[i][k]] == 1:
                    collision_detection.append(1)
                else:
                    collision_detection.append(0)
            
            if np.sum(collision_detection) or (conduction_velocity[i] > upper_threshold):
                non_colision_points.append(False)
                cv_data.append(np.nan)
            else:
                non_colision_points.append(True)
                cv_data.append(conduction_velocity[i])
            collision_detection = [] 
        
        #bin_width = 0.01
        #plt.hist(cv_data, bins=np.arange(np.nanmin(cv_data), np.nanmax(cv_data) + bin_width, bin_width))
            
        #plt.xlabel('CV (m/s)')
        #plt.ylabel('Number of Events')
        #here
        #plt.show()
        return non_colision_points, cv_data
        #pl = pv.Plotter()
        #pl.add_mesh(interpolated,cmap = boring_cmap_r)
        #pl.add_arrows(points,direction,mag=2,color=[90, 90, 90])
        #pl.show()
        
       
        
        
        #mesh_original['divergence'] = div_value
        #pl.add_mesh(mesh_original, scalars='divergence',cmap=plt.cm.get_cmap("jet", 2))
        

