# MultiRoomAcousticsSimulator

The code in this repository contains the simulater that is described in this article: https://israelcohen.com/wp-content/uploads/2021/08/3059_1.pdf

The simulator produces a RIR (Room Impulse Response) that describes the transition of a sound source between adjacent rooms so that the source is in one room and the receiver / microphone in the other.
The parameters are documented in the body of the 'test.m' file and they generate one example. 

Note: If you change the parameters, please make sure that the locations of the source and the receiver are within the boundaries of the room. 
It is assumed that room of the source is in (0,0,0) to (L_xs, L_ys, L_zs). 
Aמג אhe room of the microphone (receiver) is adjacent to it and is between (L_xs, 0,0) and (L_xs + L_xr, L_yr, L_zr).

## How to activate
1. In Matlab, you need to first compile the cpp function into a 'mex' file by using: "mex -g StIM_rir_generator.cpp"
2. Run 'test.m' in order to produce the desired RIR.
