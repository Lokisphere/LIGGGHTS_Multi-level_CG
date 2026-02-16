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

#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h"

// ========== Modification by Tarun (Part-1)  ================
// ========== Additional header files included by Tarun ======
#include "group.h"
#include "region.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "atom_vec_sphere.h"
#include "fix_property_atom.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "modify.h"
#include "fix_mcg.h"

// ========== Modification by Tarun (Part-1) ends =============

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMCG::FixMCG(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"MCG") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix MCG command");
}

/* ---------------------------------------------------------------------- */

FixMCG::~FixMCG()
{
//Destructor
}
/* ---------------------------------------------------------------------- */

int FixMCG::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
}

void FixMCG::end_of_step()
{
    grow_particles();
    if (update->ntimestep % 100 == 0) {
        mcg_refine();
        mcg_coarsen();
    
    atom->nghost = 0;
    neighbor->ago = 0;
    
    }
}

void FixMCG::grow_particles()
{
    for (int i = 0; i < atom->nlocal; i++) {

        double r = atom->radius[i];

        // ---- fine particle growth ----
        if (r >= 0.0018 && r < 0.0025) {

            atom->radius[i] += 0.000001;

            if (atom->radius[i] > 0.0025)
                atom->radius[i] = 0.0025;
        }

        // ---- coarse particle growth ----
        double parent_r = 0.005;
        double start_r  = 0.732 * parent_r;

        if (r >= start_r && r < parent_r) {

            atom->radius[i] += 0.000002;   // faster growth optional

            if (atom->radius[i] > parent_r)
                atom->radius[i] = parent_r;
        }

        // ---- update mass consistently ----
        double rho = atom->density[i];
        double rr  = atom->radius[i];
        atom->rmass[i] = (4.0/3.0) * M_PI * rho * rr*rr*rr;
    }
}


void FixMCG::mcg_refine()
{
    int nlocal_old = atom->nlocal;
    std::vector<int> delete_list;

    for (int i = 0; i < nlocal_old; i++) {

        if (atom->x[i][1] >= 0.2 &&
            atom->x[i][1] <= 0.25 &&
            fabs(atom->radius[i] - 0.005) < 1e-9) {

            refine_particle(i);
            delete_list.push_back(i);
        }
    }

    if (!delete_list.empty())
        delete_marked_atoms(delete_list);
}


/* ---------------------------------------------------------------------- */

void FixMCG::refine_particle(int ith)
{
    const double scale = 0.366;
    const double r_parent = atom->radius[ith];
    const double delta = r_parent * scale;

    int itype = atom->type[ith];

    for (int ix = -1; ix <= 1; ix += 2) {
        for (int iy = -1; iy <= 1; iy += 2) {
            for (int iz = -1; iz <= 1; iz += 2) {

                double xnew[3];
                xnew[0] = atom->x[ith][0] + ix * delta;
                xnew[1] = atom->x[ith][1] + iy * delta;
                xnew[2] = atom->x[ith][2] + iz * delta;

                // ---- OFFICIAL creation path ----
                atom->avec->create_atom(itype, xnew);
                int m = atom->nlocal - 1;

                // kinematics
                atom->v[m][0] = atom->v[ith][0];
                atom->v[m][1] = atom->v[ith][1];
                atom->v[m][2] = atom->v[ith][2];

                atom->omega[m][0] = atom->omega[ith][0];
                atom->omega[m][1] = atom->omega[ith][1];
                atom->omega[m][2] = atom->omega[ith][2];

                // physical properties
                atom->radius[m]  = r_parent * scale;
                atom->density[m] = atom->density[ith];
                atom->rmass[m]   = atom->rmass[ith] / 8.0;

                // initialize fix-owned arrays
                for (int j = 0; j < modify->nfix; j++) {
                    if (modify->fix[j]->create_attribute) {
                        modify->fix[j]->pre_set_arrays();
                        modify->fix[j]->set_arrays(m);
                    }
                }
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixMCG::delete_marked_atoms(const std::vector<int> &delete_list)
{
    int nlocal = atom->nlocal;

    std::vector<int> mark(nlocal,0);
    for (int idx : delete_list)
        if (idx < nlocal)
            mark[idx] = 1;

    int i = 0;
    while (i < atom->nlocal) {

        if (mark[i]) {

            int last = atom->nlocal - 1;

            atom->avec->copy(last, i, 1);

            mark[i] = mark[last];

            atom->nlocal--;

        } else {
            i++;
        }
    }
}

void FixMCG::mcg_coarsen()
{
    const double child_r  = 0.005 * 0.5;
    const double parent_r = 0.005;		//0.732
    const double tol = 1e-9;

    bool did_coarsen = false;

    while (true)
    {
        int nlocal = atom->nlocal;
        std::vector<int> cluster;

        // --------------------------------------------------
        // Collect first 8 fine particles below plane
        // --------------------------------------------------
        for (int i = 0; i < nlocal; i++)
        {
            if (fabs(atom->radius[i] - child_r) < tol &&
                atom->x[i][1] < 0.15 && atom->x[i][1] > 0)
            {
                cluster.push_back(i);
                if (cluster.size() == 8)
                    break;
            }
        }

        if (cluster.size() < 8)
            break;  // no more full groups available

        did_coarsen = true;

        // --------------------------------------------------
        // Compute mass + momentum
        // --------------------------------------------------
        double total_mass = 0.0;
        double vnew[3] = {0,0,0};
        double omeganew[3] = {0,0,0};

        for (int idx : cluster)
        {
            double m = atom->rmass[idx];
            total_mass += m;

            vnew[0] += m * atom->v[idx][0];
            vnew[1] += m * atom->v[idx][1];
            vnew[2] += m * atom->v[idx][2];

            omeganew[0] += m * atom->omega[idx][0];
            omeganew[1] += m * atom->omega[idx][1];
            omeganew[2] += m * atom->omega[idx][2];
        }

        vnew[0] /= total_mass;
        vnew[1] /= total_mass;
        vnew[2] /= total_mass;

        omeganew[0] /= total_mass;
        omeganew[1] /= total_mass;
        omeganew[2] /= total_mass;

        // --------------------------------------------------
        // Place coarse particle at first child location
        // --------------------------------------------------
        int trigger = cluster[0];

        double xnew[3] = {
            atom->x[trigger][0],
            atom->x[trigger][1],
            atom->x[trigger][2]
        };

        int itype = atom->type[trigger];

        atom->avec->create_atom(itype, xnew);
        int m = atom->nlocal - 1;

        atom->v[m][0] = vnew[0];
        atom->v[m][1] = vnew[1];
        atom->v[m][2] = vnew[2];

        atom->omega[m][0] = omeganew[0];
        atom->omega[m][1] = omeganew[1];
        atom->omega[m][2] = omeganew[2];

	double start_r = 0.732 * parent_r;
	atom->radius[m]  = start_r;
        atom->density[m] = atom->density[trigger];
        double rho = atom->density[m];
	atom->rmass[m] = (4.0/3.0) * M_PI * rho * start_r*start_r*start_r;

        // --------------------------------------------------
        // Delete 8 fine particles
        // --------------------------------------------------
        delete_marked_atoms(cluster);
    }

    if (!did_coarsen) return;
}


