import sys
import xarray
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import geopandas as gpd
from datetime import datetime
from goes_tools import creapaleta, radiance_to_bt, temptoindex,get_like_index, get_wind_file
from dotenv import dotenv_values,load_dotenv
import  os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

load_dotenv()
inputPath=os.getenv("DOWNLOADPATH")
inputNoaPath=os.getenv("DOWINDPATH")
outPath=os.getenv("PROD_PATH")

def prod_wind_vec(date):
    prodPath="wind/"
    dominio_1="rep"
    dominio_2="came"

    title_dom_1="Dirección del viento - México"
    title_dom_2="Dirección del viento - Centro México "
    out_file_dom_1=outPath+dominio_1+"/"+prodPath+"wind_"+date+"_"+dominio_1+".png"
    out_file_dom_2=outPath+dominio_2+"/"+prodPath+"wind_"+date+"_"+dominio_2+".png"

    date_time_obj = datetime.strptime(date, '%Y.%m%d.%H%M')
    d4_file=inputPath+"goes16.abi-"+date+"-C13_2km.tif"
    D4_dataset= xarray.open_dataset(d4_file, engine="rasterio")
    D4=D4_dataset.Rad.sel(band=1).values
    min_lon=D4_dataset.y.data[-1]
    max_lon=D4_dataset.y.data[0]
    bounds = (D4_dataset.x.data[0],D4_dataset.x.data[-1], D4_dataset.y.data[-1], D4_dataset.y.data[0])
    
    fk1=1.07364e+04
    fk2=1.38986e+03
    bc1=0.13445
    bc2=0.99955

    D4T=radiance_to_bt(D4,fk1,fk2,bc1,bc2)-273

    print("DATE $$$$$$$$$   ", date)
    dmwf= get_wind_file(date,"Noaa/")
    print("///////////////PATH",dmwf)
    
    dmw =  xarray.open_dataset(dmwf)
    croped=dmw.where((dmw.lon >= bounds[0]) & (dmw.lon <= bounds[1]) & (dmw.lat >= bounds[2]) & (dmw.lat <= bounds[3]))
    croped['wind_direction']-=180
    u = croped.wind_speed * np.sin(croped.wind_direction * np.pi/180)
    v = croped.wind_speed * np.cos(croped.wind_direction * np.pi/180)

    gdf = gpd.read_file('archivos_extra/destdv1gw/destdv1gw.shp')

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
    ax.imshow(D4T, extent=bounds,cmap="Greys")
    barbs = ax.barbs(croped['lon'], croped['lat'], u, v, croped['pressure'], linewidth=0.2, length=3.5,cmap='rainbow')
    axins = inset_axes(ax, 
                    width = '100%', 
                    height = '3%', 
                    loc = 'lower center', 
                    bbox_to_anchor=(0., -0.09, 1, 1), 
                    bbox_transform=ax.transAxes, borderpad=0) 

    cbr = plt.colorbar(barbs, cax=axins, orientation='horizontal', extendrect=None)
    cbr.set_label('Altura geopotencial [hPa]')
    im=gdf.plot(ax=ax, edgecolor='white', color='none',linewidth=.1)

    plt.title(title_dom_1+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S")+"UTC",fontsize = 9,x=0.43, y=0.95)
    plt.savefig(out_file_dom_1,bbox_inches = 'tight', pad_inches = 0)



    latitudes=D4_dataset.y.values
    longitudes=D4_dataset.x.values

    max_lat= get_like_index(latitudes,21.95 )
    min_lat=get_like_index(latitudes,15.3468056)
    min_lon=get_like_index(longitudes,-102.96770000)
    max_lon=get_like_index(longitudes,-93.19)

    bounds2 = (longitudes[min_lon],longitudes[max_lon],latitudes[min_lat],latitudes[max_lat]) 
    croped2=dmw.where((dmw.lon >= -102.96770000) & (dmw.lon <= -93.19) & (dmw.lat >=15.3468056 ) & (dmw.lat <= 21.95))

    croped2['wind_direction']-=180
    u = croped2.wind_speed * np.sin(croped2.wind_direction * np.pi/180)
    v = croped2.wind_speed * np.cos(croped2.wind_direction * np.pi/180)

    gdf2 = gpd.read_file('archivos_extra/destdv1gw/destdv1gw.shp')

    fig2 = plt.figure(figsize=(10.24,9.1), dpi=250)
    ax2= fig2.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=0.3,color='white',linestyle='--')
    gl2.top_labels = False
    gl2.right_labels = False
    gl2.ylabel_style = {'size': 6, 'color': 'white'}
    gl2.xlabel_style = {'size': 6, 'color': 'white'}
    gl2.ypadding = -10 
    gl2.xpadding = -10 

    ax2.set_xlim(longitudes[min_lon], longitudes[max_lon])
    ax2.set_ylim(latitudes[min_lat], latitudes[max_lat])
    ax2.coastlines(color='white', linewidth=0.5)


    ax2.imshow(D4T[max_lat:min_lat,min_lon:max_lon], extent=bounds2, cmap="Greys")


    barbs2 = ax2.barbs(croped2['lon'], croped2['lat'], u, v, croped2['pressure'], linewidth=0.2, length=3.5,cmap='rainbow')

    axins = inset_axes(ax2, 
                    width = '100%', 
                    height = '3%', 
                    loc = 'lower center', 
                    bbox_to_anchor=(0., -0.09, 1, 1), 
                    bbox_transform=ax2.transAxes, borderpad=0)

    cbr2 = plt.colorbar(barbs2, cax=axins, orientation='horizontal', extendrect=None)
    cbr2.set_label('Altura geopotencial [hPa]')
    im2=gdf2.plot(ax=ax2, edgecolor='white', color='none',linewidth=.1)
    plt.title(title_dom_2+ date_time_obj.strftime("%d/%m/%Y, %H:%M:%S") +"UTC",fontsize = 9,x=0.43, y=0.95)

    plt.savefig(out_file_dom_2,bbox_inches = 'tight', pad_inches = 0 ) 
