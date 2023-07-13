def check_gene_in_patient(gene_name,celltype5_id,ylabel_cell_name):
    ylabel='Normalized '+str(gene_name)+' expression in '+str(ylabel_cell_name)
    dict_use={
        'AHLD054':'M','AFBY380':'M','AHBL363':'F','AIEB037':'M','AID4428':'M','AICO097':'M','AIEL007':'F','AIEU059':'M','AHBQ224':'F','AGFY330':'M','AFCD398':'F','AGAK419':'M','BMC21673':'M','BMC34582':'F','BMC48768':'F','BMC54678':'F','BxT0DCD9':'F','BxT0DCD14':'M','BxT0LD2':'M','BxT0LD3':'F','BxT0LD4':'M','BxT0LD5':'F','BxT1LD2':"NA",'BxT1LD5':"NA"
    }
    dict_use2={
        'AHLD054':'CKD','AFBY380':'AKI','AHBL363':'Health','AIEB037':'CKD','AID4428':'AKI','AICO097':'AKI','AIEL007':'AKI','AIEU059':'CKD','AHBQ224':'AKI','AGFY330':'AKI','AFCD398':'AKI','AGAK419':'AKI','BMC21673':'Health','BMC34582':'Health','BMC48768':'Health','BMC54678':'Health','BxT0DCD9':'Health','BxT0DCD14':'Health','BxT0LD2':'Health','BxT0LD3':'Health','BxT0LD4':'CKD','BxT0LD5':'Health','BxT1LD2':"NA",'BxT1LD5':"NA"
    }
    cell_use=[]
    custom_lines = [
    Line2D([0], [0], marker='o', color='limegreen', markersize=8),
    Line2D([0], [0], marker='o', color='red', markersize=8),
    Line2D([0], [0], marker='o', color='dodgerblue', markersize=8)]

    custom_lines2 = [
    Line2D([0], [0], color='black', lw=1.5, linestyle='-'),
    Line2D([0], [0], color='black', lw=1.5, linestyle='--')]
    
    for i in adata.obs['celltype5']:
        if i in celltype5_id:
            cell_use.append(1)
        else:
            cell_use.append(0)
    adata.obs['cell_use']=cell_use
    adata_use = adata[adata.obs['cell_use'] == 1, :]

    gene_exp={}
    for i in set(adata_use.obs['patient_id']):
        test=adata_use[adata_use.obs['patient_id'].isin([i]),:].X
        test_use=pd.DataFrame(test.T)
        idx=int(list(adata_use.var_names).index(gene_name))
        gene_exp[i]=float(np.mean(test_use.iloc[idx]))
        
    xlabel=['Patient age','Patient BMI','Patient BUN (mg/dL)','Patient creatinine (mg/dL)','Patient GS%']
    count=0
    for all_dict in [[age_dict,age_dict2],[bmi_dict,bmi_dict2],[bun_dict,bun_dict2],[crea_dict,crea_dict2],[gs_dict,gs_dict2]]:
        for xxx in [all_dict[0],all_dict[1]]:
            xaxis=[]
            yaxis=[]
            scatter_color=[]
            edge_style=[]
            for i,j in xxx.items():
                if isinstance(j,int)==True or isinstance(j,float)==True:
                    if i in gene_exp.keys():
                        xaxis.append(xxx[i])
                        yaxis.append(gene_exp[i])
                        if dict_use2[i]=='AKI':
                            scatter_color.append('red')
                        elif dict_use2[i]=='CKD':
                            scatter_color.append('dodgerblue')
                        elif dict_use2[i]=='Health':
                            scatter_color.append('limegreen')
                        else:
                            scatter_color.append('dimgray')
                        
                        if dict_use[i]=='M':
                            edge_style.append('-')
                        elif dict_use[i]=='F':
                            edge_style.append('--')
            r, p = stats.pearsonr(xaxis, yaxis)
            if p < 0.1:
                fig, ax = plt.subplots(figsize=(3,3), dpi=300)
                sns.regplot(xaxis,yaxis,
                           line_kws={"color": "orange"},
                            scatter_kws={'color':scatter_color,"s": 30,'linestyle':edge_style,'edgecolors':'black'})
                ax.text(0.65, .15, 'r={:.2f}'.format(r),transform=ax.transAxes)
                ax.text(0.65, .05, 'p={:.2f}'.format(p),transform=ax.transAxes)
                ax.grid(True)
                plt.ylabel(ylabel,fontsize=9.5)
                plt.xlabel(xlabel[count])
                plt.tight_layout()
                #plt.title('title')
                first_legend = plt.legend(custom_lines,[' Health',' AKI',' CKD'], bbox_to_anchor=(1.3, 0.7),loc='center')
                plt.gca().add_artist(first_legend)
                plt.legend(custom_lines2,[' Male',' Female'],handlelength=2, bbox_to_anchor=(1.37, 0.35),loc='center')
                plt.show()
                plt.close()
                print(r,p)
            else:
                fig, ax = plt.subplots(figsize=(3,3), dpi=100)
                sns.regplot(xaxis,yaxis,
                           line_kws={"color": "orange"},
                            scatter_kws={'color':scatter_color,"s": 30,'linestyle':edge_style,'edgecolors':'black'})
                ax.text(0.65, .15, 'r={:.2f}'.format(r),transform=ax.transAxes)
                ax.text(0.65, .05, 'p={:.2f}'.format(p),transform=ax.transAxes)
                ax.grid(True)
                plt.ylabel(ylabel,fontsize=9.5)
                plt.xlabel(xlabel[count])
                plt.tight_layout()
                #plt.title('title')
                first_legend = plt.legend(custom_lines,[' Health',' AKI',' CKD'], bbox_to_anchor=(1.3, 0.7),loc='center')
                plt.gca().add_artist(first_legend)
                plt.legend(custom_lines2,[' Male',' Female'],handlelength=2, bbox_to_anchor=(1.37, 0.35),loc='center')
                plt.show()
                plt.close()
                print(r,p)
        count+=1


#run like this
check_gene_in_patient('VCAM1',['PT','PT_VCAM1','PT_dediff'],'PT')

def check_gene_in_patient_region(gene_name,region_name):
    region_name_use={'C':'cortex','M':'medulla','P':'papilla','RA':'renal artery','U':'ureter'}
    ylabel='Normalized '+str(gene_name)+' expression in '+str(region_name_use[region_name])
    adata_use = adata[adata.obs['renal_region_new'] == region_name, :]
    dict_use={
        'AHLD054':'M','AFBY380':'M','AHBL363':'F','AIEB037':'M','AID4428':'M','AICO097':'M','AIEL007':'F','AIEU059':'M','AHBQ224':'F','AGFY330':'M','AFCD398':'F','AGAK419':'M','BMC21673':'M','BMC34582':'F','BMC48768':'F','BMC54678':'F','BxT0DCD9':'F','BxT0DCD14':'M','BxT0LD2':'M','BxT0LD3':'F','BxT0LD4':'M','BxT0LD5':'F','BxT1LD2':"NA",'BxT1LD5':"NA"
    }
    dict_use2={
        'AHLD054':'CKD','AFBY380':'AKI','AHBL363':'Health','AIEB037':'CKD','AID4428':'AKI','AICO097':'AKI','AIEL007':'AKI','AIEU059':'CKD','AHBQ224':'AKI','AGFY330':'AKI','AFCD398':'AKI','AGAK419':'AKI','BMC21673':'Health','BMC34582':'Health','BMC48768':'Health','BMC54678':'Health','BxT0DCD9':'Health','BxT0DCD14':'Health','BxT0LD2':'Health','BxT0LD3':'Health','BxT0LD4':'CKD','BxT0LD5':'Health','BxT1LD2':"NA",'BxT1LD5':"NA"
    }
    custom_lines = [
    Line2D([0], [0], marker='o', color='limegreen', markersize=8),
    Line2D([0], [0], marker='o', color='red', markersize=8),
    Line2D([0], [0], marker='o', color='dodgerblue', markersize=8)]

    custom_lines2 = [
    Line2D([0], [0], color='black', lw=1.5, linestyle='-'),
    Line2D([0], [0], color='black', lw=1.5, linestyle='--')]

    gene_exp={}
    for i in set(adata_use.obs['patient_id']):
        test=adata_use[adata_use.obs['patient_id'].isin([i]),:].X
        test_use=pd.DataFrame(test.T)
        idx=int(list(adata_use.var_names).index(gene_name))
        gene_exp[i]=float(np.mean(test_use.iloc[idx]))
        
    xlabel=['Patient age','Patient BMI','Patient BUN (mg/dL)','Patient creatinine (mg/dL)','Patient GS%']
    count=0
    for all_dict in [[age_dict,age_dict2],[bmi_dict,bmi_dict2],[bun_dict,bun_dict2],[crea_dict,crea_dict2],[gs_dict,gs_dict2]]:
        for xxx in [all_dict[0],all_dict[1]]:
            xaxis=[]
            yaxis=[]
            scatter_color=[]
            edge_style=[]
            for i,j in xxx.items():
                if isinstance(j,int)==True or isinstance(j,float)==True:
                    if i in gene_exp.keys():
                        xaxis.append(xxx[i])
                        yaxis.append(gene_exp[i])
                        if dict_use2[i]=='AKI':
                            scatter_color.append('red')
                        elif dict_use2[i]=='CKD':
                            scatter_color.append('dodgerblue')
                        elif dict_use2[i]=='Health':
                            scatter_color.append('limegreen')
                        else:
                            scatter_color.append('dimgray')
                        
                        if dict_use[i]=='M':
                            edge_style.append('-')
                        elif dict_use[i]=='F':
                            edge_style.append('--')
            r, p = stats.pearsonr(xaxis, yaxis)
            if p < 0.1:
                fig, ax = plt.subplots(figsize=(3,3), dpi=300)
                sns.regplot(xaxis,yaxis,
                           line_kws={"color": "orange"},
                            scatter_kws={'color':scatter_color,"s": 30,'linestyle':edge_style,'edgecolors':'black'})
                ax.text(0.65, .15, 'r={:.2f}'.format(r),transform=ax.transAxes)
                ax.text(0.65, .05, 'p={:.2f}'.format(p),transform=ax.transAxes)
                ax.grid(True)
                plt.ylabel(ylabel,fontsize=9.5)
                plt.xlabel(xlabel[count])
                plt.tight_layout()
                #plt.title('title')
                first_legend = plt.legend(custom_lines,[' Health',' AKI',' CKD'], bbox_to_anchor=(1.3, 0.7),loc='center')
                plt.gca().add_artist(first_legend)
                plt.legend(custom_lines2,[' Male',' Female'],handlelength=2, bbox_to_anchor=(1.37, 0.35),loc='center')
                plt.show()
                plt.close()
                print(r,p)
            else:
                fig, ax = plt.subplots(figsize=(3,3), dpi=100)
                sns.regplot(xaxis,yaxis,
                           line_kws={"color": "orange"},
                            scatter_kws={'color':scatter_color,"s": 30,'linestyle':edge_style,'edgecolors':'black'})
                ax.text(0.65, .15, 'r={:.2f}'.format(r),transform=ax.transAxes)
                ax.text(0.65, .05, 'p={:.2f}'.format(p),transform=ax.transAxes)
                ax.grid(True)
                plt.ylabel(ylabel,fontsize=9.5)
                plt.xlabel(xlabel[count])
                plt.tight_layout()
                #plt.title('title')
                first_legend = plt.legend(custom_lines,[' Health',' AKI',' CKD'], bbox_to_anchor=(1.3, 0.7),loc='center')
                plt.gca().add_artist(first_legend)
                plt.legend(custom_lines2,[' Male',' Female'],handlelength=2, bbox_to_anchor=(1.37, 0.35),loc='center')
                plt.show()
                plt.close()
                print(r,p)
        count+=1

#run like this
check_gene_in_patient_region('ACTA2','C')