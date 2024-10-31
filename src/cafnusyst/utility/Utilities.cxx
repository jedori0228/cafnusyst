#include "cafnusyst/utility/Utilities.h"

namespace cafnusyst{

int UniqueName(){
  static int N = 0;
  return N++;
}

} // END namespace cafnusyst
