#include <stdexcept>

class NotImplementedException : public std::logic_error
{
  public: 
    NotImplementedException();
};