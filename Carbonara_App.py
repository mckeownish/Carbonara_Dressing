# Library Imports

import streamlit as st
import CarbonaraDataTools as CDT
import biobox
from glob import glob
import re
import os
import shutil
from stmol import showmol,obj_upload
import py3Dmol
import numpy as np
import subprocess

import writhe as wr
import plotly.graph_objects as go


# Config streamlit + plotly figures
st.set_page_config(layout="wide")
config = {'displayModeBar': False}


# Tabs for 
tab_n1 = 'Carbonara Setup'
tab_n2 = 'Carbonara Build Run-File'
tab_n3 = 'Carbonara Analysis'
tab_n4 = 'Predictions'
tab1, tab2, tab3, tab4 = st.tabs([tab_n1, tab_n2, tab_n3, tab_n4])


with tab1:


    st.title('Carbonara Setup')

    setup_c1, setup_c2 = st.columns(2)

    PDB_file = setup_c1.file_uploader("Upload PDB file")

    SAXS_file = setup_c2.file_uploader("Upload SAXS profile")

    with setup_c1:
        if PDB_file is not None:
            obj = obj_upload(PDB_file)
            obj.setStyle({'cartoon':{'color':'spectrum'}})

            showmol(obj,width=800, height=600)

    if SAXS_file is not None:
        scattering_preview = CDT.SAXS_selection_plotter(SAXS_file, 0, 10)
        scattering_preview.update_layout(height=500)
        setup_c2.plotly_chart(scattering_preview, use_container_width=True, config=config)
   
    molecule_name = st.text_input('Name your molecule!', value="Test_Molecule")

    setup_dir = st.button('Setup Run')

    if setup_dir:

        molecule_path = CDT.setup_molecule_directory(molecule_name=molecule_name)
        
        if PDB_file == None:
            st.write('Upload a PDB file!')

        else:

            pdb_file_new = os.path.join(molecule_path, PDB_file.name)
            with open(pdb_file_new,"wb") as f: 
                f.write(PDB_file.getbuffer())         
            st.success("Saved PDB File")

            M = CDT.pdb_2_biobox(pdb_file_new)

            coords = CDT.extract_CA_coordinates(M)
            sequence = CDT.extract_sequence(M)
            secondary = CDT.DSSP_structure_extractor(pdb_file_new)

            CDT.write_fingerprint_file(1, sequence, secondary, molecule_path)
            CDT.write_coordinates_file(coords, molecule_path)

        if SAXS_file == None:
            st.write('Upload a SAXS file!')

        else:
            with open(os.path.join(molecule_path, SAXS_file.name[:-4]+'_SAXS.dat'),"wb") as f: 
                f.write(SAXS_file.getbuffer())         
                st.success("Saved SAXS File")
        

with tab2:

    st.title(tab_n2)

    molecule_folders = glob('newFitData/*')
    molecule_folder_lst = []
    molecule_name_lst = [] 
    for d in molecule_folders:
        if os.path.isdir(d):
            molecule_folder_lst.append(d)
            molecule_name_lst.append(d.split('/')[-1])
        
    molecule_dict = dict(zip(molecule_name_lst, molecule_folder_lst))

    user_select_molecule_folder = st.selectbox('Select Molecule Folder', molecule_name_lst)
    molecule_path = molecule_dict[user_select_molecule_folder]+'/'
    # st.write(molecule_path)

    user_fit_name = st.text_input('Name your run!', value="Test_run")

    user_coord_file = molecule_path+'coordinates1.dat'
    user_fingerprint_file = molecule_path+'fingerPrint1.dat'

    user_coords = CDT.read_coords(user_coord_file)
    user_secondary = CDT.get_secondary(user_fingerprint_file)

    build_c1, build_c2 = st.columns(2)


    run_auto = True
    @st.cache_data
    def auto_varying_function(run_auto, user_coord_file, user_fingerprint_file):
        if run_auto:
            return CDT.find_non_varying_linkers(user_coord_file, user_fingerprint_file)


    allowed_linker, linker_indices = auto_varying_function(run_auto, user_coord_file, user_fingerprint_file)

    selected_linkers_lst = build_c1.multiselect('Select allowed varying linkers',
                                    linker_indices,
                                    allowed_linker,
                                    label_visibility="collapsed")

    x_lst, y_lst, z_lst, color_lst = CDT.smooth_me_varying(user_coords, user_secondary, selected_linkers_lst, oversample=5)
    varying_fig = CDT.line_plotly(x_lst, y_lst, z_lst, color_lst, outline=True)
    varying_fig.update_scenes(aspectratio=dict(x=1.5,y=1.5,z=1.5))
    varying_fig.update_layout(height=550)

    build_c1.write('')
    build_c1.write('')
    build_c1.plotly_chart(varying_fig, use_container_width=True, config=config)


    user_SAXS_file = glob(molecule_path+'*_SAXS.dat')[0] # might want a check here

    q_exp_min, q_exp_max, q_Q1, q_Q3 = CDT.get_minmax_q(user_SAXS_file)

    # st.write(user_SAXS_file)

    build_c2.write(user_SAXS_file)
    min_q, max_q = list( build_c2.slider( 'Select a range of q values', 0.0, q_exp_max+0.01, (q_exp_min, q_exp_max*.9), step=0.001 ))
    scattering_fig = CDT.SAXS_selection_plotter(user_SAXS_file, min_q, max_q)
    scattering_fig.update_layout(height=600)
    build_c2.plotly_chart(scattering_fig, use_container_width=True, config=config)


    build2_c1, build2_c2, build2_c3, build2_c4 = st.columns(4)

    fit_n_times = build2_c1.number_input('Choose number of Carbonara repeats', min_value=1, value=5)
    max_fit_steps = build2_c2.number_input('Choose naximum number of fit steps in each repeat', min_value=10, value=5000, step=100)


    @st.cache_data
    def finalise_run_file(molecule_path, user_SAXS_file, selected_linkers_lst, user_select_molecule_folder, user_fit_name, fit_n_times, min_q, max_q, max_fit_steps):
        # if execute:
            
            

        CDT.write_mixture_file(molecule_path)

        CDT.write_saxs(user_SAXS_file, molecule_path)
        CDT.write_varysections_file(selected_linkers_lst, molecule_path)
        CDT.write_sh_qvary_file(working_path=molecule_path, mol_name=user_select_molecule_folder, fit_name=user_fit_name, fit_n_times=fit_n_times, min_q=min_q, max_q=max_q, max_fit_steps=max_fit_steps)
        st.write('sh file written to: '+ '/RunMe_'+ user_select_molecule_folder + '_' + user_fit_name + '.sh')

            # result = subprocess.run(['bash', 'RunMe.sh'], capture_output=True, text=True)

    
    build2_c3.write('')
    build2_c3.write('')


    overwrite = build2_c3.checkbox('Overwrite existing folder')

    if build2_c3.button('Build Carbonara Run-File'):

        if os.path.exists(molecule_path+user_fit_name):
                
            # st.error("Do you really want to overwrite" + user_fit_name + " ?")
            if overwrite:
                st.write("Overwriting...")
                shutil.rmtree(molecule_path+user_fit_name)
                os.mkdir(molecule_path+user_fit_name)
                finalise_run_file(molecule_path, user_SAXS_file, selected_linkers_lst, user_select_molecule_folder, user_fit_name, fit_n_times, min_q, max_q, max_fit_steps)
           
            else:
                st.error("Run name exists - You need to check overwrite")
                
                

        else:
            st.write("Writing...")
            os.mkdir(molecule_path+user_fit_name)

            finalise_run_file(molecule_path, user_SAXS_file, selected_linkers_lst, user_select_molecule_folder, user_fit_name, fit_n_times, min_q, max_q, max_fit_steps)

        

with tab3:

    st.title(tab_n3)

    a_select_c1, a_select_c2 = st.columns(2)

    # User select molecule file
    molecule_folders_A = glob('newFitData/*')
    molecule_folder_lst_A = []
    molecule_name_lst_A = [] 
    for d in molecule_folders_A:
        if os.path.isdir(d):
            molecule_folder_lst_A.append(d)
            molecule_name_lst_A.append(d.split('/')[-1])
        
    molecule_dict_A = dict(zip(molecule_name_lst_A, molecule_folder_lst_A))
    user_molecule_folder_analysis = a_select_c1.selectbox('Select Molecule Folder for Analysis', molecule_name_lst_A)
    molecule_path_A = molecule_dict_A[user_molecule_folder_analysis]+'/'

    # Load the relavent fingerprint
    fingerprint_file = molecule_path_A+'fingerPrint1.dat'
    secondary = np.asarray( list(open(fingerprint_file, 'r').readlines()[-1][:-1]) )

    # User select run file
    test_folders_A = glob(molecule_path_A+'*')
    test_folder_lst_A = []
    test_name_lst_A = [] 
    for d in test_folders_A:
        if os.path.isdir(d):
            test_folder_lst_A.append(d)
            test_name_lst_A.append(d.split('/')[-1])
        

    test_dict_A = dict(zip(test_name_lst_A, test_folder_lst_A))
    user_test_folder_analysis = a_select_c2.selectbox('Select Test Folder for Analysis', test_name_lst_A)
    spec_folder = test_dict_A[user_test_folder_analysis]+'/'

    

    # Just a sanity check for any file existing in chosen directory!
    sub_file_size = len( glob(spec_folder+'*') )
    if sub_file_size < 1:

        st.error('Uh oh.. No files in selected directory')

    else:
        regex = re.compile(r'\d+')
        tot_run_num = regex.findall(CDT.sort_by_creation(glob(spec_folder+'fitLog*'))[-1].split('/')[-1])[0]
        if tot_run_num == '1':
            tot_run_num = regex.findall(CDT.sort_by_creation(glob(spec_folder+'fitLog*'))[-2].split('/')[-1])[0]


        

        sel_c1, sel_c2, sel_c3 = st.columns(3)


        run_num = sel_c1.selectbox('Select run number', np.arange(2, int(tot_run_num)+1))

        log_file = glob(spec_folder+'fitLog' + str(run_num) + '.dat')

        if len(log_file)!=1:

            if len(log_file)==0:
                st.error("Log file doesn't exist! :o  **"+'fitLog' + str(run_num) + '.dat**')

            else:
                st.error("Multiple log files with the same name! :o")

        df_log = CDT.log2df(log_file[0])

        # st.table(df_log)

        tot_structs = df_log.shape[0]-1
        # struc_num = sel_c2.selectbox('Select run number', np.arange(1, int(tot_structs)), key='unique_structure')
        struc_num = sel_c2.select_slider('Select run number', np.arange(0, int(tot_structs)), key='unique_structure')

        sel_c3.write('')
        full_q = sel_c3.checkbox('Show all q')
        outline = sel_c3.checkbox('Contrast outline', True)


        SAXS_file = glob(molecule_path_A+'*_SAXS.dat')[0]
        # SAXS_file = working_path+'human_SMARCAL1Saxs.dat'
        # fit_file = 'smarcalData/withTagTest1/mol2SubstepScatter_27.dat'
        
        if struc_num > 0:

            # auto select Scatter file with mol-# and step-#
            fit_file = glob(spec_folder+'mol' + str(run_num) + 'SubstepScatter_' + str(struc_num) + '.dat')[0]

            # auto select coordinates with mol-# and step-#
            coord_file = glob(spec_folder+'mol'+ str(run_num) + 'Substep_' + str(struc_num) + '_*')[0]

        elif struc_num == 0:

            fit_file = glob(spec_folder+'mol' + str(run_num) + 'origScatter.dat')[0]
            coord_file = molecule_path_A+'coordinates1.dat'


        @st.cache_data
        def plot_log_st(df_log, struc_num):

            # Plot log file
            log_fig = CDT.df2plot(df_log, highlight=struc_num)
            
            return log_fig

        @st.cache_data
        def plot_saxs_st(SAXS_file, fit_file, full_q):
            
            # Plot scatter fit
            st.write(SAXS_file)
            st.write(fit_file)
            fit_fig = CDT.SAXS_fit_plotter(SAXS_file, fit_file, full_q)
            fit_fig.update_layout(height=600)
            
            return fit_fig


        @st.cache_data
        def plot_coords_st(coord_file, secondary, outline):
        
            coords = CDT.extract_coords(coord_file)
            x_lst, y_lst, z_lst, color_lst = CDT.smooth_me(coords, secondary, oversample=5)
            structure_fig = CDT.line_plotly(x_lst, y_lst, z_lst, color_lst, outline=outline)
            structure_fig.update_layout(height=600)

            structure_fig.update_scenes(aspectratio=dict(x=1.5,y=1.5,z=1.5))
            
            return structure_fig


        log_fig = plot_log_st(df_log, struc_num)
        st.plotly_chart(log_fig, use_container_width=True, config=config)

        # st.table(df_log)
        fig_c1, fig_c2 = st.columns(2)

        fit_fig = plot_saxs_st(SAXS_file, fit_file, full_q)
        fig_c1.plotly_chart(fit_fig, use_container_width=True, config=config)

        structure_fig =plot_coords_st(coord_file, secondary, outline)
        fig_c2.plotly_chart(structure_fig, use_container_width=True, config=config)



with tab4:

    st.title(tab_n3)

    # p_select_c1, p_select_c2 = st.columns(2)

    # User select molecule file
    molecule_folders_P = glob('newFitData/*')
    molecule_folder_lst_P = []
    molecule_name_lst_P = [] 
    for d in molecule_folders_P:
        if os.path.isdir(d):
            molecule_folder_lst_P.append(d)
            molecule_name_lst_P.append(d.split('/')[-1])
        
    molecule_dict_P = dict(zip(molecule_name_lst_P, molecule_folder_lst_P))
    user_molecule_folder_analysis = st.selectbox('Select Molecule Folder for Analysis', molecule_name_lst_P, key='Molecule select global')
    molecule_path_P = molecule_dict_P[user_molecule_folder_analysis]+'/'

    # Load the relavent fingerprint
    fingerprint_file = molecule_path_P+'fingerPrint1.dat'
    secondary = np.asarray( list(open(fingerprint_file, 'r').readlines()[-1][:-1]) )

    # # User select run file

    test_folders_P = glob(molecule_path_P+'*')
    test_folder_lst_P = []
    test_name_lst_P = [] 
    for d in test_folders_P:
        if os.path.isdir(d):
            test_folder_lst_P.append(d)
            test_name_lst_P.append(d.split('/')[-1])
        


    log_files_all = []
    for folder in test_name_lst_P:
        log_files_all = log_files_all + glob(molecule_path_P+folder+ '/fitLog*')


    best_fit_files = []
    max_q_lst = []
    best_fit_scatter_files = []

    for f in log_files_all:
    #     print(f)

        try:
            log_df = CDT.log2df(f)[1:-1]
            
            tmp = log_df.loc[log_df.groupby('max_fitting_k')['total_penalty'].idxmin()]
            
            best_fits = list( tmp['file_name'].values )
            max_qs = list( tmp['max_fitting_k'].values )
            
            for idx, bf in enumerate(best_fits):
        
                name = log_df['file_name'].values[-1].split('/')[-1]
                nums = re.findall(r'\d+', name)

                path = '/'.join(log_df['file_name'].values[-1].split('/')[:-1]) 
                file = glob(path+'/mol'+nums[0]+'Substep_'+nums[1]+'*')
                
                if len(file) == 1:
                    
                    best_fit_scatter_files.append(best_fits[idx])
                    best_fit_files.append(file[0])
                    max_q_lst.append(max_qs[idx])
        except:
            pass


    if len(best_fit_files) < 1:

        st.error('Yikes.. No best fitting files in selected molecule')


    else:
        coord_tensor = CDT.extract_coords(best_fit_files[0])


        for f in best_fit_files[1:]:
            
            tmp = CDT.extract_coords(f)
            coord_tensor = np.dstack((coord_tensor, tmp))  


        @st.cache_data
        def cluster_run(coord_tensor, min_cluster_size=8):

            rmsd_arr = CDT.coord_tensor_pairwise_rmsd(coord_tensor)   


            labels, probabilities = CDT.cluster(rmsd_arr, min_cluster_size)
            
            return labels, probabilities


        @st.cache_data
        def cluster_plot(labels):
            cluster_fig, bf_names_sort = CDT.visualise_clusters(coord_tensor, labels, best_fit_files)
            return cluster_fig, bf_names_sort
        

        labels, probabilities = cluster_run(coord_tensor, min_cluster_size=8)

        if len( np.unique(labels) ) > 4:
            st.error('Too many clusters (try increasing the min cluster size!)')

        else:
            cluster_fig, bf_names_sort = cluster_plot(labels)

            sub_sel1, sub_sel2 = st.columns(2)
            user_A_idx = sub_sel1.select_slider('Select structure A', np.arange(len(bf_names_sort[1:])))
            user_B_idx = sub_sel2.select_slider('Select structure B', np.arange(len(bf_names_sort[1:]))) # explore later - double counting of first file?

            user_struct_A = bf_names_sort[1:][user_A_idx]
            user_struct_B = bf_names_sort[1:][user_B_idx]

            
            # st.write(bf_names_sort)

            # cluster_fig.update_layout(yaxis=dict(scaleanchor='x'))
            @st.cache_data
            def update_cluster_plot(cluster_fig, bf_names_sort, user_A_idx, user_B_idx):

                # idx_i = np.where(bf_names_sort[1:]==user_struct_A)[0].item()
                # idx_j = np.where(bf_names_sort[1:]==user_struct_B)[0].item()

                cluster_fig.add_trace(go.Scatter(x=[user_A_idx], y=[user_B_idx],  mode='markers', line=dict(color="#FF0000")))
                cluster_fig.add_trace(go.Scatter(x=[user_B_idx], y=[user_A_idx],  mode='markers', line=dict(color="#FF0000")))
                cluster_fig.update_layout(xaxis_range=[0,len(bf_names_sort)], yaxis_range=[0,len(bf_names_sort)], showlegend=False)
                return cluster_fig


            
            up_cluster_fig = update_cluster_plot(cluster_fig, bf_names_sort, user_A_idx, user_B_idx)
            st.plotly_chart(up_cluster_fig, config=config)

            if st.button('Calculate Writhe Fingerprints'):


                user_coords_A = CDT.extract_coords(molecule_path_P+user_struct_A+'.dat')
                user_coords_B = CDT.extract_coords(molecule_path_P+user_struct_B+'.dat')

                writhe_fp_A = wr.writhe_fingerprint(segments=user_coords_A[::4])
                writhe_fp_B = wr.writhe_fingerprint(segments=user_coords_B[::4])

                writhe_fig_A = wr.writhe_heatmap_plot(writhe_fp_A)
                writhe_fig_A.update_layout(yaxis=dict(scaleanchor='x'))
                writhe_fig_B = wr.writhe_heatmap_plot(writhe_fp_B)
                writhe_fig_B.update_layout(yaxis=dict(scaleanchor='x'))

                writhe_fig_diff = wr.writhe_heatmap_plot(writhe_fp_B-writhe_fp_A)
                writhe_fig_diff.update_layout(yaxis=dict(scaleanchor='x'))

                wr_c1, wr_c2, wr_c3 = st.columns(3)
                wr_c1.plotly_chart(writhe_fig_A, config=config)
                wr_c2.plotly_chart(writhe_fig_B, config=config)
                wr_c3.plotly_chart(writhe_fig_diff, config=config)

            if st.button('Show overlay'):
                

                unique_labels = np.sort( np.unique(labels) )[1:]
                cols_overlay = st.columns(len(unique_labels))

                for il, l in enumerate(unique_labels):
                    conds = (probabilities > 0.8) & (labels==l)

                    tensor = coord_tensor[:,:,conds]
                    new_tensor = CDT.align_coords(tensor)
                    overlay_fig = CDT.overlay_coords(new_tensor)
                    cols_overlay[il].plotly_chart(overlay_fig, config=config)
