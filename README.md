# SoundComm

Transfer Text Message over sound signal by means of digital communication techniques

Doesn't support realtime processing.

inside the Transmitter folder use TransmitterApp.mlapp to create a sound signal which contains your text message. user can specify frequency,SymbolRate and FEC there, also live spectrum button can be used to see 0-20kHz live spectrum for choosing proper band

in Receiver folder use RecieverApp.mlapp for precessing recieved sound signal. the message file needs to be saved on the receiver device and its path is one of the inputs of Receiver. This receiver operates in three modes:

Rates about 1kBits/sec over 15KHz carrier have been Tested Successfully


DFE-RLS

DFE-LMS

viterbi demodulator without equalizer
