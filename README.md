# SHERPA-simulation

## Working with SRR at both 0.1x0.1 and 0.1x0.05 deg resolution
Now the code has been adapted creating a new 'sherpa_globals.py', to allow simulations for

- SRR at 0.1x0.1 degrees (using 'sherpa_globals_cams42_01_01.py')
- SRR at 0.1x0.05 degrees (using 'sherpa_globals_cams42_01_005.py')

To use one of the 2 options, please copy 'sherpa_globals_cams42_01_**.py' to 'sherpa_globals.py'

Note that the downscaling feature (presented in the next subsection) works both for 0.1x0.1 and 0.1x0.05 SRR.

## Preparing results for the downscaling module
In this branch I adapt the 0.1x0.1 deg resolutionSRR to work with a new downscaling module. This module requires a different structure of emissions and concentration files.

For emissions:

- a unique file for all air quality indexes, containing all precursors and (for PM) PPM2.5, PPMco, PPM10
- same content should be available for basecase (initial) emissions and delta (resulting from the SRR simulation) emissions 
- the code will work with 12 GNFR sectors
- note that PPM10 = PPM2.5 + PMco

For concentrations

- you also need basecase files for PPM2.5 and PPM10 concentrations
- the PPM2.5 basecase is directly extracted from the EMEP results, the PPM10 is on the contrary computed as the sum of PPM2.5 + PPMco
- not that the SRR for PM downscaling will be run twice: once to simulate total PM2.5 (or PM10) and once to simulate the PPM2.5 (or PPM10) reductions. This second run is performed only reducing PPM2.5 or PPM10 emissions.

Now the code uses different tagas for emissions definitions (PPM25 or PPM10, depending on emissions to be reduced). 
Also now a new variable has been added when running module1 (downscale_request). 
If the variable downscale_request=0, a run with normal reductions is done. If downscale_request=1 then a run is done only reducing PPM emissions. 
This second run creates concentrations change for the PPM component. Finally, now filename are automatically created, depending on the concentrations to be modelled.
But also depending on the downscale_request. When downscale_request=1 indeed a 'primary' tag is added to the module1 output file.

# Module available for simulation
With this code you can run simulations with the SHERPA SRR (built i.e. with the git repository SHERPA-training), to evaluate the impact of emission reduction scenarios on concentrations.
This version of the SHERPA-simulation code has been tested withe the EDGAR2015 emission inventory, and with SRR built from EMEP air quality model simulations.

The modules that can be used in the code are:

-  Module 1 (scenario assessment): to simulate the impact on air quality of a specific emission reduction scenario (defined also through the previous two steps)
-  Module 3 (Source allocation): to understand how the air quality in a given area is influenced by different sources. This module runs in two modes: precursor-based source allocation, and sector-based one. Note that this module works running module 4 and module 1, in sequence.
-  Module 6 (Governance): to analyze how one should coordinate with the surrounding regions to optimally improve air quality;
-  Module 8 (health impact): to evaluate PM2.5 health-related impact, when running module 1
-  Module 9 (aggregation): to aggregate emissions and concentrations, at NUTS or FUAs level.

# Publications

- Degraeuwe, B., Pisoni, E., Thunis, P.
Prioritising the sources of pollution in European cities: Do air quality modelling applications provide consistent responses?
(2020) Geoscientific Model Development, 13 (11), pp. 5725-5736. 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85096928556&doi=10.5194%2fgmd-13-5725-2020&partnerID=40&md5=d4cb7a3413af75830b782af529db3727

- Sartini, L., Antonelli, M., Pisoni, E., Thunis, P.
From emissions to source allocation: Synergies and trade-offs between top-down and bottom-up information
(2020) Atmospheric Environment: X, 7, art. no. 100088, . 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85090038106&doi=10.1016%2fj.aeaoa.2020.100088&partnerID=40&md5=c41e93233b2f100e7dcda69e064a90ea

- Belis, C.A., Pisoni, E., Degraeuwe, B., Peduzzi, E., Thunis, P., Monforti-Ferrario, F., Guizzardi, D.
Urban pollution in the Danube and Western Balkans regions: The impact of major PM2.5 sources
(2019) Environment International, 133, art. no. 105158, . 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85073224266&doi=10.1016%2fj.envint.2019.105158&partnerID=40&md5=81f54f1ee6fc53184059b655797f3ffa

- Pisoni, E., Thunis, P., Clappier, A.
Application of the SHERPA source-receptor relationships, based on the EMEP MSC-W model, for the assessment of air quality policy scenarios
(2019) Atmospheric Environment: X, 4, art. no. 100047, . 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85072582264&doi=10.1016%2fj.aeaoa.2019.100047&partnerID=40&md5=34146ae8e90b2bc98ed4babfea0991a3

- Peduzzi, E., Pisoni, E., Clappier, A., Thunis, P.
Multi-level policies for air quality: implications of national and sub-national emission reductions on population exposure
(2018) Air Quality, Atmosphere and Health, 11 (9), pp. 1121-1135. 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85053388835&doi=10.1007%2fs11869-018-0613-1&partnerID=40&md5=ec09ca92720060078aac98a3408f85a5

- Monforti-Ferrario, F., Kona, A., Peduzzi, E., Pernigotti, D., Pisoni, E.
The impact on air quality of energy saving measures in the major cities signatories of the Covenant of Mayors initiative
(2018) Environment International, 118, pp. 222-234. 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85048115990&doi=10.1016%2fj.envint.2018.06.001&partnerID=40&md5=39ee982bcf359ed03505d2e3d33ec995

- Pisoni, E., Albrecht, D., Mara, T.A., Rosati, R., Tarantola, S., Thunis, P.
Application of uncertainty and sensitivity analysis to the air quality SHERPA modelling tool
(2018) Atmospheric Environment, 183, pp. 84-93. 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85045688029&doi=10.1016%2fj.atmosenv.2018.04.006&partnerID=40&md5=9432074fa33ac072995b6076f64cce3c

- Thunis, P.
On the validity of the incremental approach to estimate the impact of cities on air quality
(2018) Atmospheric Environment, 173, pp. 210-222. 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85034083743&doi=10.1016%2fj.atmosenv.2017.11.012&partnerID=40&md5=f877863da1d86c9bb5d844870ce349aa

- Thunis, P., Degraeuwe, B., Pisoni, E., Ferrari, F., Clappier, A.
On the design and assessment of regional air quality plans: The SHERPA approach
(2016) Journal of Environmental Management, 183, pp. 952-958. 
https://www.scopus.com/inward/record.uri?eid=2-s2.0-84994012263&doi=10.1016%2fj.jenvman.2016.09.049&partnerID=40&md5=1561546680304fcf57e914bdf441d452
