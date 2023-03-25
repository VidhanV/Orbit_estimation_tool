#!/usr/bin/env python
# coding: utf-8

# In[35]:


import pandas as pd
from datetime import datetime
from decimal import Decimal
columns= ["epoch","x","y","z","vx","vy","vz"]
pos_df= pd.read_csv(r'C:\Users\VIVO\Documents\posvec_obs.csv', usecols=columns)
date=[]
for i in pos_df["epoch"]:
    d = datetime.utcfromtimestamp(i)
    date.append(d)


pos_df["Date"]= date
pos_df.dropna(inplace=True)
display(pos_df)


# In[36]:


import datetime
orbit_epoch = datetime.datetime(2021, 12, 19)


# In[21]:


import orekit
orekit.initVM()


# In[4]:


from orekit.pyhelpers import download_orekit_data_curdir, setup_orekit_curdir
download_orekit_data_curdir()
setup_orekit_curdir()


# In[22]:


from org.orekit.frames import FramesFactory
from org.orekit.utils import IERSConventions
gcrf = FramesFactory.getGCRF()
itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
eme2000 = FramesFactory.getEME2000()
from org.orekit.time import TimeScalesFactory
utc = TimeScalesFactory.getUTC()
from org.orekit.models.earth import ReferenceEllipsoid
wgs84Ellipsoid = ReferenceEllipsoid.getWgs84(itrf)
from org.orekit.bodies import CelestialBodyFactory
moon = CelestialBodyFactory.getMoon()
sun = CelestialBodyFactory.getSun()


# In[37]:


from org.hipparchus.geometry.euclidean.threed import Vector3D
from orekit.pyhelpers import datetime_to_absolutedate
from org.orekit.utils import PVCoordinates
import math
i = math.ceil(len(pos_df.index) / 2)
points_for_iod = pos_df.iloc[[0, i, -1]]

from org.orekit.estimation.measurements import Position
from org.hipparchus.geometry.euclidean.threed import Vector3D
from orekit.pyhelpers import datetime_to_absolutedate
from org.orekit.utils import PVCoordinates
from org.orekit.estimation.measurements import PV
pos_1 = points_for_iod.iloc[0].to_list()
p_vector_1 = Vector3D(float(pos_1[0]),float(pos_1[1]),float(pos_1[2]))
v_vector_1 = Vector3D(float(pos_1[3]),float(pos_1[4]),float(pos_1[5]))
date_1 = datetime_to_absolutedate(pos_1[7])
tmp_pv1 = PVCoordinates(p_vector_1, v_vector_1)
pv1 = itrf.getTransformTo(eme2000, date_1).transformPVCoordinates(tmp_pv1)

pos_2 = points_for_iod.iloc[1].to_list()
p_vector_2 = Vector3D(float(pos_2[0]),float(pos_2[1]),float(pos_2[2]))
v_vector_2 = Vector3D(float(pos_2[3]),float(pos_2[4]),float(pos_2[5]))
date_2 = datetime_to_absolutedate(pos_2[7])
tmp_pv2 = PVCoordinates(p_vector_2, v_vector_2)
pv2 = itrf.getTransformTo(eme2000, date_2).transformPVCoordinates(tmp_pv2)

pos_3 = points_for_iod.iloc[2].to_list()
p_vector_3 = Vector3D(float(pos_3[0]),float(pos_3[1]),float(pos_3[2]))
v_vector_3 = Vector3D(float(pos_3[3]),float(pos_3[4]),float(pos_3[5]))
date_3 = datetime_to_absolutedate(pos_3[7])
tmp_pv3 = PVCoordinates(p_vector_3, v_vector_3)
pv3 = itrf.getTransformTo(eme2000, date_3).transformPVCoordinates(tmp_pv3)


# In[53]:


from org.orekit.estimation.iod import IodGibbs
from org.orekit.utils import Constants as orekit_constants
iod_gibbs = IodGibbs(orekit_constants.EIGEN5C_EARTH_MU)
initialOrbit = iod_gibbs.estimate(eme2000,
                                      pv1.getPosition(), date_1,
                                      pv2.getPosition(), date_2,
                                      pv3.getPosition(), date_3)
display(initialOrbit)
initialCartOrbit = OrbitType.CARTESIAN.convertType(initialOrbit)
display(initialCartOrbit)
display(initialCartOrbit.getDate())


# In[39]:


# Estimator parameters
estimator_position_scale = 1.0 # m
estimator_convergence_thres = 1e-2
estimator_max_iterations = 50
estimator_max_evaluations = 60

# Orbit propagator parameters
prop_min_step = 0.001 # s
prop_max_step = 300.0 # s
prop_position_error = 10.0 # m

from org.orekit.propagation.conversion import DormandPrince853IntegratorBuilder
integratorBuilder = DormandPrince853IntegratorBuilder(prop_min_step, prop_max_step, prop_position_error)

from org.orekit.propagation.conversion import NumericalPropagatorBuilder
from org.orekit.orbits import PositionAngle, OrbitType
initialCartOrbit = OrbitType.CARTESIAN.convertType(initialOrbit)
propagatorBuilder = NumericalPropagatorBuilder(initialCartOrbit,
                                               integratorBuilder,
                                               PositionAngle.TRUE,
                                               estimator_position_scale)

from org.orekit.forces.gravity.potential import GravityFieldFactory
gravityProvider = GravityFieldFactory.getConstantNormalizedProvider(7, 7)
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
gravityAttractionModel = HolmesFeatherstoneAttractionModel(itrf, gravityProvider)
propagatorBuilder.addForceModel(gravityAttractionModel)

from org.hipparchus.linear import QRDecomposer
matrix_decomposer = QRDecomposer(1e-11)
from org.hipparchus.optim.nonlinear.vector.leastsquares import GaussNewtonOptimizer
optimizer = GaussNewtonOptimizer(matrix_decomposer, False)

from org.orekit.estimation.leastsquares import BatchLSEstimator
estimator = BatchLSEstimator(optimizer, propagatorBuilder)
estimator.setParametersConvergenceThreshold(estimator_convergence_thres)
estimator.setMaxIterations(estimator_max_iterations)
estimator.setMaxEvaluations(estimator_max_evaluations)


# In[40]:


from org.orekit.estimation.measurements import Position, ObservableSatellite

observableSatellite = ObservableSatellite(0) # Propagator index = 0

for index, pv_gcrf in pos_df.iterrows():
    tmp_p = Vector3D(pv_gcrf[['x','y','z']].to_list())
    tmp_v = Vector3D(pv_gcrf[['vx','vy','vz']].to_list())
    tmp_pvc = PVCoordinates(tmp_p, tmp_v)
    pvc = itrf.getTransformTo(eme2000, datetime_to_absolutedate(pv_gcrf[7])).transformPVCoordinates(tmp_pvc)
    orekit_position = Position(
        datetime_to_absolutedate(pv_gcrf[7]),
        pvc.getPosition(),
        10.0,
        1.0,  # Base weight
        observableSatellite
    )
    estimator.addMeasurement(orekit_position)

estimatedPropagatorArray = estimator.estimate()


# In[41]:


from org.orekit.orbits import KeplerianOrbit
estimatedPropagator = estimatedPropagatorArray[0]
estimatedInitialState = estimatedPropagator.getInitialState()
tmp_orbit = OrbitType.KEPLERIAN.convertType(estimatedInitialState.getOrbit())
orbit = KeplerianOrbit.cast_(tmp_orbit)
display(orbit)

