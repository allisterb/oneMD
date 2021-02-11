#include "Figlet.hh"
#include <sstream>
using namespace std;

int
main() {
  Figlet::standard.print("Fractions");
  for ( int i = 2; i <= 4; ++i ) {
    ostringstream ss;
    ss << "5/" << i << " = " << 5.0/i;
    Figlet::small.print(ss.str().c_str());
  }
  cout << "ALL DONE!\n";
  return 0;
}