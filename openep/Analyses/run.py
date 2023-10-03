import openep 
from CV import CV
import pyvista as pv
import numpy as np
def run_cv (openep_file_name,cv_method, visualsiation):
    case =  openep.load_openep_mat(openep_file_name)
    egmX = case.electric.bipolar_egm.points
    voltage = case.electric.bipolar_egm.voltage
    print('---> read electrode data')
    LAT = case.electric.annotations.local_activation_time - case.electric.annotations.reference_activation_time
    print('---> read bi-polar LAT data')
    surface = case.create_mesh()

    valid_id = np.where(LAT!=-10000)
    invalid_id = np.where(LAT==-10000) 
    egmX_clean = np.delete(egmX, np.where(LAT == -10000),0)
    LAT_clean = np.delete(LAT, np.where(LAT == -10000))
    Points = surface.points
    faces = surface.faces
    Face = faces.reshape(-1,4)
    temp_face = Face[:,1:4]
    mesh_original = pv.PolyData(Points,Face)
    cal_cv = CV()
    if cv_method == 'Plane':
        cv_plane, cv_centers_plane, cv_centers_id_plane = cal_cv.plane_fitting(egmX_clean, LAT_clean)
        cv_total = np.zeros((len(egmX),1))
        for i in range(len(cv_plane)):
            cv_total[valid_id[0][i]] = cv_plane[i]
        
        for i in range(len(invalid_id)):
            cv_total[invalid_id[i]] = np.nan
        
        if visulaisation:
            cal_cv.visualization(lat_data = '', cv_data = cv_plane, mesh = mesh_original, egmX = egmX_clean, method = 'CV',file_name = '', 
                                show_histogram= True, direction_centers = Points,show_cv_fibrosis=False)
        cv = cv_plane
            
    if cv_method == 'RBF':
        cv_rbf, cv_centers_rbf, direction_rbf, cv_centers_id_rbf = cal_cv.RBF_method(egmX_clean, LAT_clean,Points, Face)
        if visulaisation:
            cal_cv.visualization(lat_data = '', cv_data = cv_rbf, mesh = mesh_original, egmX = np.transpose(cv_centers_rbf), method = 'CV',file_name = '', 
                                show_histogram= True, direction_centers = Points,show_cv_fibrosis=False)
        cv = cv_rbf

    if cv_method == 'Tri':
        cv_tri, cv_centers_tri,  cv_centers_id_tri = cal_cv.Triangulation(egmX_clean,LAT_clean)
        if visulaisation:
            cal_cv.visualization(lat_data = '', cv_data = cv_tri, mesh = mesh_original, egmX = np.transpose(cv_centers_tri), method = 'CV',file_name = '', 
                                show_histogram= True, direction_centers = Points,show_cv_fibrosis=False)
        cv = cv_tri
    ######## Save the CV data as openep #######
    case.fields['conduction_velocity'] = cv
    return case
    

    
    
    


if __name__=="__main__":
    file_path = '/Users/aligharaviri/Downloads/LAWT_CartoAF_ EP_Studies/'
    file_name = '521.mat'
    openep_file_name = f"{file_path}{file_name}"
    out_put_file_name = f"{file_name[:-4]}_cv.mat"
    method = "Plane" ######## methods are Plane, RBF, and Tri
    visulaisation = True
    case = run_cv (openep_file_name,method, visulaisation)
    openep.export_openep_mat(case,out_put_file_name)
    
