/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */


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

// ========== Modification by Tarun (Part 1/1) =====================
// ========== Additional member functions included by Tarun ======
  virtual void mcg_refine();
  virtual void mcg_coarsen();
  virtual void refine_particle(int);
  virtual void delete_marked_atoms(const std::vector<int> &);
  int setmask();
  virtual void end_of_step();
  virtual void grow_particles();

private:
  //FixPropertyAtom *fix_template_;
  int *dlist;
// ========== Modification by Tarun (Part 1/1) ends ==============

// protected:
//  int extra;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Angle coeffs are not set

No angle coefficients have been assigned in the data file or via the
angle_coeff command.

E: All angle coeffs are not set

All angle coefficients must be set in the data file or by the
angle_coeff command before running a simulation.

*/
