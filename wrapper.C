#include "lvis_mol_lens.C"
#include "read_inputs/read_inputs.h"
#include "read_inputs/read_inputs.cxx"

#ifdef __MAKECINT__
// add here any needed #pragma
#endif

int wrapper(int arg = 0, string region = "", string par_type = "", int ite = 0, double step = 0)
{
  return lvis_mol_lens(arg, region, par_type, ite, step);
}
