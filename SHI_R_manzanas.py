# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 13:49:34 2020
MDF extension on Jan 4 12:00 2021

@author: andr3 + agodoy
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import os.path
import geopandas as gpd

np.set_printoptions(suppress=True, #prevent numpy exponential 
    formatter={'float_kind':'{:f}'.format})                        #notation on print, default False

MXC_population = 19438576
pumped_water = 430430000 #annual water pumped in from surrounding catchment.

debug=1
in_dir="Input_files/"#C:\\Modelling\\Mexico\\Python_SHI\\"
con_dir="Control_files/"

#Zero
read_fail=0

#====================================================
#=====Get data, weightings and derive parameters=====
#====================================================
file_dir=(con_dir+"Input_files.txt")
indicator_data = pd.read_csv(file_dir, header=0, delim_whitespace=True)
file_dir=(con_dir+"ahp_weighting_new.txt")
ahp_weight = pd.read_csv(file_dir, header=0, delim_whitespace=True)


#======Universal data ====== 
#===========================
#shapefile with area of polygons:
fp='shapefiles/manzanas/manzanas_cdmx.shp'
map_df = gpd.read_file(fp)
# check data type so we can see that this is not a normal dataframe, but a GEOdataframe
print map_df.crs
map_df["area"] = map_df['geometry'].area*1000000 #m*m

map_df['ID2']=map_df['ID']
map_df2=map_df[['ID','ID2','area']]
print map_df2.head() 
#create dataframe wiht all data:
inputfile = 'fakedataset_manzanas.csv'
alldata = pd.read_csv(inputfile, header=0)

alldata = alldata.set_index('ID').join(map_df2.set_index('ID2'))
del map_df['CVEGEO']
print alldata.head()
print(alldata.columns.tolist())


#-----Raw Data-----
#------------------
data_name="Population"
print("Getting "+data_name+" data.....")
find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
file_name=find_data.iloc[0][2] #Find file name
file_dir=(in_dir+file_name) #Create file directory string
alldata[file_name] = alldata[file_name].replace(-4444,np.nan)
name_colum=alldata.loc[:,file_name]
population_array = name_colum.values
print 'Zeros in Population', np.count_nonzero(population_array==0)
population_array[population_array==-3333]=np.nan 
WSI_array = np.zeros(len(population_array))
ACI_array = np.zeros(len(population_array)) 
 
  
data_name="Recharge"
print("Getting "+data_name+" data.....")
find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
file_name=find_data.iloc[0][2] #Find file name
file_dir=(in_dir+file_name) #Create file directory string
name_colum=alldata.loc[:,file_name]
Rech_array = name_colum.values
name_colum=alldata.loc[:,'area']
area_array = name_colum.values
Rech_array = (Rech_array/1000)*365*area_array  #annual volume m3 
Rech_array[Rech_array==-4444]=np.nan
print 'Nans in Recharge', np.isnan(Rech_array).sum()      
   
data_name="Water_useage_pc"
print("Getting "+data_name+" data.....")
find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
file_name=find_data.iloc[0][2] #Find file name
file_dir=(in_dir+file_name)
name_colum=alldata.loc[:,file_name]
WUpc_array = name_colum.values
WUpc_array = WUpc_array.astype('float64')
WUpc_array[ WUpc_array==-4444]=np.nan
print 'Nans in Water_useage_pc', np.isnan(WUpc_array).sum() 
     
         
#-----Derived data-----
#----------------------
if read_fail!=1: #Check to see if all files read correctly 
    WRpc_array = (Rech_array/population_array)+(pumped_water/MXC_population)   #Water resources per person
    #======Water Variation (WV)====== 
    #================================
    print("Processing WV =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Water_variation')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
        
    data_name="Rainfall_variation"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    indicator_data[indicator_data['Data'].str.match(data_name)]
    file_dir=(in_dir+file_name) #Create file directory string
    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-4444]=np.nan  
    print 'Nans in Rainfall_variation', np.isnan(indicator_array_1).sum()  
  
    #-----undertake data calculations-----
    #-------------------------------------
    indicator_array_1[indicator_array_1>=0.4]=1  
            
    #-----Combine data into parameter-----
    #-------------------------------------
    #Jaramillo Eq.2 and Eq.17
    indicator_array_X=indicator_array_1
    WSI_array=WSI_array+(indicator_array_X*ahp)
                 

    if debug==1:   
	WV_data = alldata[['ID','CVEGEO']]
	WV_data['WV_out'] = indicator_array_X
	WV_data.to_csv('WV_out.csv',index=False) 
        
	map_df2 = map_df.set_index('ID').join(WV_data.set_index('ID'))


	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='WV_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
        #size of colorbars:
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
        # add the colorbar to the figure
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'WV Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('WV_out.png', dpi=fig.dpi)
	
	del WV_data
      
    #======Water Scarcity (WS)====== 
    print("Processing WS =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Water_scarcity')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
    
    data_name="Pop_poor_access"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-4444]=np.nan 
    indicator_array_1[indicator_array_1==-3333]=np.nan 
    print 'Nans in Pop_poor_access', np.isnan(indicator_array_1).sum()      
    indicator_array_1[indicator_array_1>=1]=1 
  
    #-----undertake data calculations-----
    #-------------------------------------
    indicator_array_1 = indicator_array_1 #percentage of the population with poor access to water
    indicator_array_2 = 1700/WRpc_array #Falkenmark water stress indicator Jaramillo Eq.5 
    indicator_array_2[indicator_array_2>=1]=1  
    indicator_array_2[indicator_array_2==-4444]=np.nan 

    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=np.sqrt((indicator_array_1**2)+(indicator_array_2**2))/np.sqrt(2)
    WSI_array=WSI_array+(indicator_array_X*ahp) #Jaramillo Eq.4 without adaquacy (a)
  
        
    if debug==1:   
	WS_data = alldata[['ID','CVEGEO']]
	WS_data['WS_out'] = indicator_array_X
	WS_data.to_csv('WS_out.csv',index=False)
 	
	#add the corresponding column to the geopandas dataframe:
	map_df2 = map_df.set_index('ID').join(WS_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='WS_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'WS Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('WS_out.png', dpi=fig.dpi)

        del WS_data

    #======Water Exploitation (WE)====== 
    print("Processing WE =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Water_resource_exploitation')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
            
    #-----undertake data calculations-----
    #-------------------------------------
    indicator_array_1 = (WUpc_array/WRpc_array) #resource development rate
    indicator_array_1[indicator_array_1>=0.4]=1  
    indicator_array_1[indicator_array_1==-4444]=np.nan  
           
    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=indicator_array_1
    WSI_array=WSI_array+(indicator_array_X*ahp)
    if debug==1:    
	WE_data = alldata[['ID','CVEGEO']]
	WE_data['WE_out'] = indicator_array_X
	WE_data.to_csv('WE_out.csv',index=False)
        
	map_df2 = map_df.set_index('ID').join(WE_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='WE_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'WE Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('WE_out.png', dpi=fig.dpi)

	del WE_data

     
    #======Water Polution (WP)====== 
    print("Processing WP =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Water_pollution')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
    
    #Find required 
    data_name="Wastewater_pc"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-3333]=np.nan 
    indicator_array_1[indicator_array_1==-4444]=np.nan 
    print 'Nans in Wastewater_pc', np.isnan(indicator_array_1).sum()      
                
    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=indicator_array_1
    WSI_array=WSI_array+(indicator_array_X*ahp)
       
    if debug==1:   
	WP_data = alldata[['ID','CVEGEO']]
	WP_data['WP_out'] = indicator_array_X
	WP_data.to_csv('WP_out.csv',index=False)

	map_df2 = map_df.set_index('ID').join(WP_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='WP_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'WP Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('WP_out.png', dpi=fig.dpi)
	
        del WP_data

    #=====Natural Capacity (NC)=====
    #Capacity to provide a hydrological service
    print("Processing NC =======") 
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Natural_capacity')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter


    data_name="Land_class"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-4444]=np.nan
        
    #-----undertake data calculations-----
    #-------------------------------------
    indicator_array_1[indicator_array_1>=1]=1 
    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=indicator_array_1
    ACI_array=ACI_array+(indicator_array_X*ahp) 
        
    if debug==1:  
	NC_data = alldata[['ID','CVEGEO']]
	NC_data['NC_out'] = indicator_array_X
	NC_data.to_csv('NC_out.csv',index=False)
  
	map_df2 = map_df.set_index('ID').join(NC_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='NC_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'NC Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('NC_out.png', dpi=fig.dpi)

	del NC_data
   

    #=====Physical Capacity (PC)=====
    #integrated capacity of farmers, local community members, and government entities to divert water to meet water demand
    print("Processing PC =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Physical_capacity')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
    
    data_name="Land_irrigation"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-4444]=np.nan  #any, mostly zeros
    indicator_array_1[indicator_array_1>=1]=1 
      
    data_name="Drinking_water"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
    name_colum=alldata.loc[:,file_name]
    indicator_array_2 = name_colum.values
    indicator_array_2[indicator_array_2==-4444]=np.nan	
    indicator_array_2[indicator_array_2==-3333]=np.nan
    print 'Nans in Drinking_water', np.isnan(indicator_array_2).sum()      
    indicator_array_2[indicator_array_2>=1]=1 
  

    data_name="Surface_water_storage"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name) #Create file directory string
    name_colum=alldata.loc[:,file_name]
    indicator_array_3 = name_colum.values
    indicator_array_3 = indicator_array_3.astype('float64')
    indicator_array_3[indicator_array_3==-3333]=np.nan  
    indicator_array_3[indicator_array_3==-4444]=np.nan  
    #-----undertake data calculations-----
    #-------------------------------------
    indicator_array_1=indicator_array_1*0.53 #irrigation multiplied by fraction of water used for irrigation/industry
    indicator_array_2=indicator_array_2*0.47 #drinkable water multiplied by fraction of water used for drinking
    indicator_array_3=indicator_array_3/WUpc_array #surface water 
                        
    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=indicator_array_1+indicator_array_2+indicator_array_3 #Jaramillo Eq.18
    ACI_array=ACI_array+(indicator_array_X*ahp)
        
    if debug==1:   
	PC_data = alldata[['ID','CVEGEO']]
	PC_data['PC_out'] = indicator_array_X
	PC_data.to_csv('PC_out.csv',index=False) 

	map_df2 = map_df.set_index('ID').join(PC_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='PC_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'PC Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('PC_out.png', dpi=fig.dpi)
	
	del PC_data
       
    #=====Human Resource Capacity (HC)=====
    print("Processing HC =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Human_resource_capacity')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
    
    data_name="Literacy_rate"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)

    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-4444]=np.nan 
    indicator_array_1[indicator_array_1==-3333]=np.nan
    print 'Nans in Literacy_rate', np.isnan(indicator_array_1).sum()      
    indicator_array_1[indicator_array_1>=1]=1     
        
    data_name="Econ_act_pop"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
    name_colum=alldata.loc[:,file_name]
    indicator_array_2 = name_colum.values
    indicator_array_2[indicator_array_2==-4444]=np.nan
    indicator_array_2[indicator_array_2==-3333]=np.nan  
    print 'Nans in Econ_act_pop', np.isnan(indicator_array_2).sum()      
    indicator_array_2[indicator_array_2>=1]=1 
    #-----undertake data calculations-----
    #-------------------------------------
    #Not needed
    
    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=(indicator_array_1+indicator_array_2)/2 #Jaramillo Eq.18
    ACI_array=ACI_array+(indicator_array_X*ahp)
        
    if debug==1:   
	HC_data = alldata[['ID','CVEGEO']]
	HC_data['HC_out'] = indicator_array_X
	HC_data.to_csv('HC_out.csv',index=False) 

	map_df2 = map_df.set_index('ID').join(HC_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='HC_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'HC Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('HC_out.png', dpi=fig.dpi)

	del HC_data
        
    #=====Economic Capacity (EC)=====
    print("Processing EC =======")
    find_weight=ahp_weight[ahp_weight['Sub_index'].str.match('Economic_capacity')] #Find info from table
    ahp=find_weight.iloc[0][1] #Get weighting once for parameter
    
    data_name="Average_income"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)

    name_colum=alldata.loc[:,file_name]
    indicator_array_1 = name_colum.values
    indicator_array_1[indicator_array_1==-4444]=np.nan

    data_name="Min_income"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)

    name_colum=alldata.loc[:,file_name]
    indicator_array_2 = name_colum.values
    #indicator_array_2[indicator_array_2==-4444]=np.nan 
            
    data_name="Max_income"
    print("Getting "+data_name+" data.....")
    find_data=indicator_data[indicator_data['Data'].str.match(data_name)] #Find info from table
    file_name=find_data.iloc[0][2] #Find file name
    file_dir=(in_dir+file_name)
 
    name_colum=alldata.loc[:,file_name]
    indicator_array_3 = name_colum.values    
    #indicator_array_3[indicator_array_3==-4444]=np.nan 
    
    #-----undertake data calculations-----
    #-------------------------------------
    indicator_array_1=np.log(indicator_array_1)
    indicator_array_2=np.log(indicator_array_2)
    indicator_array_3=np.log(indicator_array_3)
                        
    #-----Combine data into parameter-----
    #-------------------------------------
    indicator_array_X=(indicator_array_1-indicator_array_2)/(indicator_array_3-indicator_array_2) #Jaramillo Eq.18
    ACI_array=ACI_array+(indicator_array_X*ahp)
                
        
    if debug==1:  
        EC_data = alldata[['ID','CVEGEO']]
	EC_data['EC_out'] = indicator_array_X
	EC_data.to_csv('EC_out.csv',index=False)

	map_df2 = map_df.set_index('ID').join(EC_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='EC_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]),norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]), norm=plt.Normalize(vmin=np.nanmin(indicator_array_X), vmax=np.nanmax(indicator_array_X)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'EC Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('EC_out.png', dpi=fig.dpi)

	del EC_data

#===========================  
#=====Index and Outputs=====
#=========================== 
if read_fail!=1: #Check to see if all files read correctly 
 
    if debug==1:    
	WSI_data = alldata[['ID','CVEGEO']]
	WSI_data['WSI_out'] = WSI_array
	WSI_data.to_csv('WSI_out.csv',index=False)

	map_df2 = map_df.set_index('ID').join(WSI_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='WSI_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]),norm=plt.Normalize(vmin=np.nanmin(WSI_array), vmax=np.nanmax(WSI_array)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","yellow","red"]), norm=plt.Normalize(vmin=np.nanmin(WSI_array), vmax=np.nanmax(WSI_array)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'WSI Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('WSI_out.png', dpi=fig.dpi)

	del WSI_data
        
        ACI_data = alldata[['ID','CVEGEO']]
	ACI_data['ACI_out'] = ACI_array
	ACI_data.to_csv('ACI_out.csv',index=False)

	map_df2 = map_df.set_index('ID').join(ACI_data.set_index('ID'))
	fig=plt.figure()
	ax = plt.subplot(111, aspect='equal')
	variable='ACI_out'

	map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]),norm=plt.Normalize(vmin=np.nanmin(ACI_array), vmax=np.nanmax(ACI_array)),ax=ax)
	sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["red","yellow","green"]), norm=plt.Normalize(vmin=np.nanmin(ACI_array), vmax=np.nanmax(ACI_array)))
	sm._A = []
	cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
	ax.patch.set(hatch='x', edgecolor='black')

	cbar = fig.colorbar(sm,cax=cax)
	cbar.set_label(r'ACI Index', rotation=270,labelpad=18, fontsize=11)
	ax.axis('off')
	fig.savefig('ACI_out.png', dpi=fig.dpi)

	del ACI_data
    
    #=====Calculate Index=====
    #========================= 
    print("Calculating SHI.....")
    ACI_array[ACI_array==0]=np.nan
    SHI_index_array=(WSI_array/ACI_array) #Jaramillo Eq.1
    print 'Nans in WSI_array', np.isnan(WSI_array).sum()      
    print 'Nans in ACI_array', np.isnan(ACI_array).sum()      
    print 'Nans in SHI_index_array', np.isnan(SHI_index_array).sum()   
    print 'WSI_array min max', np.nanmin(WSI_array),np.nanmax(WSI_array)
    print 'ACI_array', np.nanmin(ACI_array),np.nanmax(ACI_array),np.count_nonzero(ACI_array==0)
    print 'SHI', np.nanmin(SHI_index_array),np.nanmax(SHI_index_array)
    sortSHI=np.sort(-SHI_index_array)

    #=====Outputs======
    #==================
    
    #-----Output SHI as ASCII grid-----
    print("Writing SHI to gridded ASCII.....")
    SHI_data = alldata[['ID','CVEGEO']]
    SHI_data['SHI_out'] = SHI_index_array
    SHI_data.to_csv('SHI_out.csv',index=False)

    map_df2 = map_df.set_index('ID').join(SHI_data.set_index('ID'))
    fig=plt.figure()
    ax = plt.subplot(111, aspect='equal')
    variable='SHI_out'

    map_df2.plot(column=variable,cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","orange","orangered","red"]),norm=plt.Normalize(vmin=np.nanmin(SHI_index_array), vmax=np.nanmax(SHI_index_array)),ax=ax)
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("", ["green","orange","orangered","red"]), norm=plt.Normalize(vmin=np.nanmin(SHI_index_array), vmax=np.nanmax(SHI_index_array)))
    sm._A = []
    cax = fig.add_axes([ax.get_position().x1+0.0,ax.get_position().y0,0.02,ax.get_position().height])
    ax.patch.set(hatch='x', edgecolor='black')

    cbar = fig.colorbar(sm,cax=cax)
    cbar.set_label(r'SHI Index', rotation=270,labelpad=18, fontsize=11)
    ax.axis('off')
    fig.savefig('SHI_out.png', dpi=fig.dpi)

    del SHI_data
    
    print("Finished :)")
    
else:
    print ("FILE READ FAILED!!!!!")
 

