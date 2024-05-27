import sys
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import geopandas as gpd
from datetime import datetime
from goes_tools import bin_ndarray,get_like_index
from dotenv import dotenv_values,load_dotenv
import  os

load_dotenv()
inputPath=os.getenv("DOWNLOADPATH")
outPath=os.getenv("PROD_PATH")

def prod_lst(date):
    print("lst :",date) 
    prodPath="lst/"
    dominio_1="rep"
    dominio_2="came"
    title_dom_1="LST - México "
    title_dom_2="LST - Centro México" 
    out_file_dom_1=outPath+dominio_1+"/"+prodPath+"lst_"+date+"_"+dominio_1+".png"
    out_file_dom_2=outPath+dominio_2+"/"+prodPath+"lst_"+date+"_"+dominio_2+".png"

    red_file=inputPath+"goes16.abi-"+date+"-C02_0.5km.tif"
    green_file=inputPath+"goes16.abi-"+date+"-C03_1km.tif"
    blue_file=inputPath+"goes16.abi-"+date+"-C01_1km.tif"
    lst_file=inputPath+"goes16.abi-"+date+"-C01_1km.tif"

    file = inputPath+'goes16.abi-'+date+'-LST_LST_2km.tif'
    df = xarray.open_dataset(file)

    date_time_obj = datetime.strptime(date, '%Y.%m%d.%H%M')
    B_dataset=xarray.open_dataset(blue_file, engine="rasterio")
    matriz=B_dataset.Rad.sel(band=1).values
    nx=matriz.shape[1]
    ny=matriz.shape[0]
    B=matriz[:ny,:nx-1]
    low_result=B.shape
    print(B.shape)

    R_dataset=xarray.open_dataset(red_file, engine="rasterio")
    matriz=R_dataset.Rad.sel(band=1).values
    nx=matriz.shape[1]
    ny=matriz.shape[0]
    R_aux=matriz[:ny,:nx-1]

    R= bin_ndarray(R_aux, low_result, operation='avg') 
    print(R.shape)

    G_dataset=xarray.open_dataset(green_file, engine="rasterio")
    matriz=G_dataset.Rad.sel(band=1).values
    nx=matriz.shape[1]
    ny=matriz.shape[0]
    G=matriz[:ny,:nx-1]
    print(G.shape)

    Esun_Ch_01 = 441.868715 #726.721072
    Esun_Ch_02 = 663.274497 # 663.274497
    Esun_Ch_03 = 726.721072
    d2 = .3

    R_reflec = (R * np.pi * d2) / Esun_Ch_02
    G_reflec = (G * np.pi * d2) / Esun_Ch_03
    B_reflec = (B * np.pi * d2) / Esun_Ch_01

    R_color = np.clip(R_reflec, 0, 1)
    G_color = np.clip(G_reflec, 0, 1)
    B_color= np.clip(B_reflec, 0, 1)

    gamma = 2.2
    R_bright = np.power(R_color, 1/gamma)
    G_bright = np.power(G_color, 1/gamma)
    B_bright = np.power(B_color, 1/gamma)

    G_true = 0.45706946 * R_bright +0.06038137* G_bright + 0.48358168   * B_bright
    G_true = np.clip(G_true, 0, 1)
    RGB = np.dstack([R_bright, G_true, B_bright])
    bounds = (B_dataset.x.data[0],B_dataset.x.data[-1], B_dataset.y.data[-1], B_dataset.y.data[0])


    matriz=df.LST.sel(band=1)-273
    gdf = gpd.read_file('archivos_extra/destdv1gw/destdv1gw.shp')


    bounds = (df.x.data[0],df.x.data[-1], df.y.data[-1], df.y.data[0])

    fig = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.1, color='white',linestyle='--')
    gl.top_labels=False   # suppress top labels
    gl.right_labels=False
    gl.ylabel_style = {'size': 6, 'color': 'white'}
    gl.xlabel_style = {'size': 6, 'color': 'white'}
    gl.ypadding = -10 
    gl.xpadding = -10 
    ax.coastlines(color='white', linewidth=0.2)
    ax.imshow(RGB, extent = bounds, origin = 'upper')
    imgd=ax.imshow(matriz, extent = bounds, origin = 'upper', cmap="coolwarm", vmin = -15,
                                vmax = 55)
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.1)
    ##cbar=fig.colorbar(imgd, extend='both', orientation='horizontal', pad=0.25,fraction=0.05)
    ##cbar.ax.set_title('(°C) Celsius',fontsize = 6)
    plt.title('LST '+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9, x=0.43, y=0.95)
    plt.savefig(out_file_dom_1,bbox_inches = 'tight', pad_inches = 0 )
    print("salida lst_d1: ", out_file_dom_1)
    
    latitudes=B_dataset.y.values
    longitudes=B_dataset.x.values
    max_lat= get_like_index(latitudes,21.95 )
    min_lat=get_like_index(latitudes,15.3468056)
    min_lon=get_like_index(longitudes,-102.96770000)
    max_lon=get_like_index(longitudes,-93.19)

    bounds2 = (longitudes[min_lon],longitudes[max_lon],latitudes[min_lat],latitudes[max_lat])

    latitudes1=df.y.values
    longitudes1=df.x.values

    max_lat1= get_like_index(latitudes1,21.95 )
    min_lat1=get_like_index(latitudes1,15.3468056)
    min_lon1=get_like_index(longitudes1,-102.96770000)
    max_lon1=get_like_index(longitudes1,-93.19)
    bounds3 = (longitudes1[min_lon1],longitudes1[max_lon1],latitudes1[min_lat1],latitudes1[max_lat1])

    fig = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(),extent = bounds2)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.1, color='white',linestyle='--')
    gl.top_labels=False   # suppress top labels
    gl.right_labels=False
    gl.ylabel_style = {'size': 6, 'color': 'white'}
    gl.xlabel_style = {'size': 6, 'color': 'white'}
    gl.ypadding = -10 
    gl.xpadding = -10 
    ax.coastlines(color='white', linewidth=0.2)

    ax.imshow(RGB[max_lat:min_lat,min_lon:max_lon], extent = bounds2, origin = 'upper')
    imgd=ax.imshow(matriz[max_lat1:min_lat1,min_lon1:max_lon1], extent = bounds3, origin = 'upper', cmap="coolwarm", vmin = -15, vmax = 55)

    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.1)
    #cbar=fig.colorbar(imgd,extend='both', orientation='horizontal', pad=0.1,fraction=0.05)
    #cbar.ax.set_title('(°C) Celsius',fontsize = 6)
    plt.title('LST '+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9, x=0.43, y=0.95)
    plt.savefig(out_file_dom_2,bbox_inches = 'tight', pad_inches = 0 )
    print("salida lst_d2: ", out_file_dom_2)
