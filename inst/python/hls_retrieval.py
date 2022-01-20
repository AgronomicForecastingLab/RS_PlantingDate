import requests
import json
import time
import os
from pandas.io.json import json_normalize
from collections import defaultdict
import pandas as pd
import zipfile, io

values={"aoi": '{"type":"Feature","geometry":{"type":"Polygon","coordinates":[[[-92.8005803, 43.2758463], [-92.8005803, 43.2748463], [-92.7995803, 43.2748463], [-92.7995803, 43.2758463], [-92.8005803, 43.2758463]]]}}',
"Band": "['NDVI']",
"Enddate": "3/8/2019",
"Startdate": "3/2/2019",
"legendtype": "Relative",
"qafilter": "1",
"satellite": "Landsat,Sentinel",
"resolution":"0.00001"}
    
headers = {'Content-Type':'application/x-www-form-urlencoded'}
           
def harmonized_landsat_data(values, headers):
    try:
        url = 'https://ag-analytics.azure-api.net/hls-service/'
        response = requests.post(url, data = values, headers = headers).json()
        print(response)
        return response
    
    except Exception as e:
        raise e

hlsresponse = harmonized_landsat_data(values, headers)
