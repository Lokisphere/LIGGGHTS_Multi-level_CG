#ifdef FIX_CLASS
FixStyle(mcg,FixMCG)
#else

#ifndef LMP_FIX_MCG_H
#define LMP_FIX_MCG_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixMCG : public Fix {

 public:
  FixMCG(class LAMMPS *, int, char **);
  ~FixMCG();

  int setmask();
  void end_of_step();
  void pre_exchange();

 private:

  // control
  bigint mcg_every;
  int refine_flag;

  // core logic
  void grow_particles();
  int  mcg_refine_and_coarsen();
  int refine();
  int coarsen();
  void refine_particle(int);
  void delete_marked_atoms(const std::vector<int>&);

};

}

#endif
#endif

