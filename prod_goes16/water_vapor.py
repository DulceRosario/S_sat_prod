import sys
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import geopandas as gpd
from datetime import datetime
from goes_tools import radiance_to_bt,get_like_index
from dotenv import dotenv_values,load_dotenv
import  os

load_dotenv()
inputPath=os.getenv("DOWNLOADPATH")
outPath=os.getenv("PROD_PATH")

def prod_water_vapor(date):
    prodPath="water_vapor/"

    dominio_1="rep"
    dominio_2="came"
    title_dom_1="Water Vapor RGB - México "
    title_dom_2="Water Vapor RGB - Centro México "
    out_file_dom_1=outPath+dominio_1+"/"+prodPath+"water-vapor_"+date+"_"+dominio_1+".png"
    out_file_dom_2=outPath+dominio_2+"/"+prodPath+"water-vapor_"+date+"_"+dominio_2+".png"
    blue_file=inputPath+"goes16.abi-"+date+"-C10_2km.tif"
    red_file=inputPath+"goes16.abi-"+date+"-C13_2km.tif"
    green_file=inputPath+"goes16.abi-"+date+"-C08_2km.tif"


    date_time_obj = datetime.strptime(date, '%Y.%m%d.%H%M')

    R_dataset= xarray.open_dataset(red_file, engine="rasterio")
    R=(R_dataset.Rad.sel(band=1).values)
    fk1=1.07364e+04
    fk2=1.38986e+03 
    bc1=0.13445
    bc2=0.99955
    RT=radiance_to_bt(R,fk1,fk2,bc1,bc2)-273.15

    G_dataset= xarray.open_dataset(green_file, engine="rasterio")
    G=G_dataset.Rad.sel(band=1).values
    fk1=5.03614e+04
    fk2=2.32657e+03 
    bc1=2.12504 
    bc2=0.99541
    GT=radiance_to_bt(G,fk1,fk2,bc1,bc2)-273.15 

    B_dataset= xarray.open_dataset(blue_file, engine="rasterio")
    B=B_dataset.Rad.sel(band=1).values
    fk1=3.00925e+04
    fk2=1.95961e+03
    bc1=0.06984
    bc2=0.99983
    BT=radiance_to_bt(B,fk1,fk2,bc1,bc2)-273.15

    # Minimuns and Maximuns
    Rmin =  5.81
    Rmax = -70.86
    Gmin = -30.48
    Gmax = -58.49
    Bmin = -12.12
    Bmax = -28.03

    RT[RT<Rmax] = Rmax
    RT[RT>Rmin] = Rmin
    GT[GT<Gmax] = Gmax
    GT[GT>Gmin] = Gmin
    BT[BT<Bmax] = Bmax
    BT[BT>Bmin] = Bmin

    # Choose the gamma
    gamma_R = 1.5
    gamma_G = 1.5
    gamma_B = 1.5

    # Normalize the data
    R_bright = ((RT - Rmin) / (Rmax - Rmin)) ** (1/gamma_R)
    G_bright = ((GT - Gmin) / (Gmax - Gmin)) ** (1/gamma_G)
    B_bright = ((BT - Bmin) / (Bmax - Bmin)) ** (1/gamma_B) 

    RGB = np.dstack([R_bright, G_bright, B_bright])

    bounds = (B_dataset.x.data[0],B_dataset.x.data[-1],B_dataset.y.data[-1],B_dataset.y.data[0])
    gdf = gpd.read_file('archivos_extra/destdv1gw/destdv1gw.shp')
    mpl.rc('text', color='white')
    fig = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.3,color='white',linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'size': 6, 'color': 'white'}
    gl.xlabel_style = {'size': 6, 'color': 'white'}
    gl.ypadding = -10 
    gl.xpadding = -10 
    ax.coastlines(color='white', linewidth=0.2)
    ax.imshow(RGB, extent = bounds, origin = 'upper')
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.08)
    plt.title(title_dom_1+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9,x=0.43, y=0.95)
    plt.savefig(out_file_dom_1,bbox_inches = 'tight', pad_inches = 0)

    latitudes=B_dataset.y.values
    longitudes=B_dataset.x.values

    max_lat= get_like_index(latitudes,21.95 )
    min_lat=get_like_index(latitudes,15.3468056)
    min_lon=get_like_index(longitudes,-102.96770000)
    max_lon=get_like_index(longitudes,-93.19)

    bounds2 = (longitudes[min_lon],longitudes[max_lon],latitudes[min_lat],latitudes[max_lat])


    fig = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.imshow(RGB[max_lat:min_lat,min_lon:max_lon], extent = bounds2, origin = 'upper')
    im2=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.15)
    ax.set_xlim(longitudes[min_lon], longitudes[max_lon])
    ax.set_ylim(latitudes[min_lat], latitudes[max_lat])
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.2,color='white',linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabel_style = {'size': 6, 'color': 'white'}
    gl.xlabel_style = {'size': 6, 'color': 'white'}    
    gl.ypadding = -10 
    gl.xpadding = -10 
    plt.title(title_dom_2+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9,x=0.43, y=0.95)

    plt.savefig(out_file_dom_2,bbox_inches = 'tight', pad_inches = 0 )