:::::::::::::::::
:: SHERPA TEST ::
:::::::::::::::::

SET city=London

:: Source receptor models
SET pm25model=input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc
SET pm10model=input/20151116_SR_no2_pm10_pm25/SR_PM10_Y.nc
SET no2eqmodel=input/20151116_SR_no2_pm10_pm25/SR_NO2eq_Y.nc
:: base case emissions
SET pm25emis=input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc
SET pm10emis=input/20151116_SR_no2_pm10_pm25/BC_emi_PM10_Y.nc
SET no2eqemis=input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc
:: base case concentrations
SET pm25conc=input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc
SET pm10conc=input/20151116_SR_no2_pm10_pm25/BC_conc_PM10_Y.nc
SET no2eqconc=input/20151116_SR_no2_pm10_pm25/BC_conc_NO2_NO2eq_Y_mgm3.nc

:: run module 1
echo "running module 1 for NO2eq"
sherpa.exe 1 %no2eqemis% input/London_region.nc input/user_reduction_snap7.txt %no2eqconc% %no2eqmodel% output/NO2/London/
echo "module 1 finished"

echo "running module 1 for PM25"
sherpa.exe 1 %pm25eqemis% input/London_region.nc input/user_reduction_snap7.txt %pm25conc% %pm25model% output/PM25/London/
echo "module 1 finished"

echo "running module 1 for PM10"
sherpa.exe 1 %pm10emis% input/London_region.nc input/user_reduction_snap7.txt %pm10conc% %pm25model% output/PM10/London/
echo "module 1 finished"

echo "running module 6 for NO2"
sherpa.exe 6 %no2eqemis% ../input/EMI_RED_ATLAS_NUTS2.nc "51.51" "-0.13" ../input/user_reduction_all50.txt %no2eqconc% %no2eqmodel% output/NO2/%city%/
echo "module 6 finished"

echo "running module 6 for PM25"
sherpa.exe 6 %pm25emis% ../input/EMI_RED_ATLAS_NUTS2.nc "51.51" "-0.13" ../input/user_reduction_all50.txt %pm25conc% %pm25model% output/PM25/London/
echo "module 6 finished"

echo "running module 6 for PM10"
sherpa.exe 6 %pm10emis% ../input/EMI_RED_ATLAS_NUTS2.nc "51.51" "-0.13" ../input/user_reduction_all50.txt %pm10conc% %pm10model% output/PM10/London/
echo "module 6 finished"

echo "running module 3a with 1 precursor"
sherpa.exe 31 %pm25emis% input/London_region.nc ../input/potency_reduction_module3a1P.txt %pm25conc% %pm25model% output/PM25/London/
echo "module 3a with 1 precursor finished"

echo "running module 3a with 2 precursors"
sherpa.exe 31 %pm25emis% input/London_region.nc ../input/potency_reduction_module3a2P.txt %pm25conc% %pm25model% output/PM25/London/
echo "module 3a with 2 precursors finished"

echo "running module 3b"
sherpa.exe 32 %pm25emis% input/London_region.nc ../input/potency_reduction_module3b.txt %pm25conc% %pm25model% output/PM25/London/
echo "module 3b finished"

echo "running module 2"
sherpa.exe 2 %pm25emis% ../input/EMI_RED_ATLAS_NUTS0.nc input/user_reduction_snap7.txt %modelpath% output/
echo "module 2 finished"

REM pause

REM echo "running module 4"
REM sherpa.exe 4 %pm25emis% input/London_region.nc ../input/user_reduction_snap7.txt %pm25conc% %pm25model% output/PM25/London/
REM echo "module 4 finished"
REM pause

REM echo "running module 5"
REM sherpa.exe 5 %pm25emis% ../input/EMI_RED_ATLAS_NUTS0.nc ../input/user_reduction_snap7.txt %pm25conc% %pm25model% output/PM25/London/
REM echo "module 5 finished"
REM pause

