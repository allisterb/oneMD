#include "Exceptions.hh"

NotImplementedException::NotImplementedException() : std::logic_error("Not yet implemented.") { };
