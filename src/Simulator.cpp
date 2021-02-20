#include "Simulator.hh"
#include "xdrfile.h"
Simulator::Simulator(const string _name, const int _nd, const int _np, const int _ts, const float _ts_delta, const Device _device) :
  name(_name),
  nd(_nd),
  np(_np),
  ts(_ts),
  ts_delta(_ts_delta),
  device(_device)
{}