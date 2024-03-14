
import random
from typing import Union, Iterable

import numpy as np
from numpy import ndarray

from models import TransitionModel,ObservationModel,StateModel

#
# Add your Robot Simulator here
#
class RobotSim:
    def __init__(self, trueState, sm, tm, om):
        self.__currentState = trueState
        self.__sm = sm
        self.__tm = tm
        self.__om = om

    def new_state(self) -> int:
        #matrix = self.__tm.get_T()
        #max_index = np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)
        #new_state = max_index
        #return new_state
        prob = self.__tm.get_T()[self.__currentState]
        self.__currentState = np.random.choice(self.__sm.get_num_of_states(), p=prob)
        return self.__currentState

    def new_sensor(self) -> int:
        prob = []
        for num in range(self.__om.get_nr_of_readings()):
            prob = np.append(prob, self.__om.get_o_reading_state(num, self.__currentState))
        reading = np.random.choice(self.__sm.get_num_of_readings(), p=prob)
        #prob = self.__om.get_o_reading(self.__sm.state_to_reading(self.__currentState))
        #prob = [prob[i][i] for i in range(0, len(prob), 4)]
        #prob = prob + [1 - sum(prob)]

        if reading >= len(prob) - 1:
            new_reading = None
        else:
            new_reading = reading
        #matrix = self.__om.get_o_reading_state(self.__currentState)
        #max_index = np.unravel_index(np.argmax(matrix, axis=None), matrix.shape)
        #new_reading = self.__sm.position_to_reading(max_index[0], max_index[1])
        #return new_reading
        #new_reading = self.__sm.state_to_reading(self.__currentState)
        return new_reading
        
#
# Add your Filtering approach here (or within the Localiser, that is your choice!)
#
class HMMFilter:
    def __init__(self, sm, tm, om):
        self.__sm = sm
        self.__tm = tm
        self.__om = om

    def update_f(self, sensor_reading, fVec):
        #print(self.__om.get_o_reading(sensor_reading))
        #print(self.__tm.get_T_transp())
        #print(fVec)
        fVec_new = np.matmul(self.__om.get_o_reading(sensor_reading), np.matmul(self.__tm.get_T_transp(), fVec))
        fVec_new = fVec_new / np.sum(fVec_new)
        self.__currentF = fVec_new
        return self.__currentF

        
        
        
