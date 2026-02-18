#include <cmath>
#include <algorithm>
#include <string.h>

#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "error.h"
#include "neighbor.h"
#include "modify.h"
#include "fix_mcg.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMCG::FixMCG(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp, narg, arg)
{
  if (strcmp(style,"mcg") != 0 && narg < 3)
    error->all(FLERR,"Illegal fix MCG command");

  mcg_every = 100;     // refinement interval
  refine_flag = 0;
}

/* ---------------------------------------------------------------------- */

FixMCG::~FixMCG() {}

/* ---------------------------------------------------------------------- */

int FixMCG::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    mask |= PRE_EXCHANGE;
    return mask;
}

void FixMCG::end_of_step()
{
    grow_particles();

    if (update->ntimestep % mcg_every == 0)
        refine_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixMCG::pre_exchange()
{
    if (!refine_flag)
        return;

    refine_flag = 0;

    mcg_refine_and_coarsen();

    // ----------- SAFE GLOBAL NATOMS RECOMPUTE -----------
    bigint nlocal = atom->nlocal;

    MPI_Allreduce(&nlocal,
                  &atom->natoms,
                  1,
                  MPI_LMP_BIGINT,
                  MPI_SUM,
                  world);

    // ----------- TAG EXTENSION -----------
    if (atom->tag_enable)
    {
        atom->tag_extend();

        if (atom->map_style)
        {
            atom->nghost = 0;
            atom->map_init();
            atom->map_set();
        }
    }

    neighbor->ago = 0;   // force rebuild
}


/* ---------------------------------------------------------------------- */

void FixMCG::grow_particles()
{
  for (int i = 0; i < atom->nlocal; i++) {

      double r = atom->radius[i];

      // fine growth
      if (r >= 0.0018 && r < 0.0025)
          atom->radius[i] = std::min(r + 0.2e-6, 0.0025);

      // coarse growth
      double parent_r = 0.005;
      double start_r  = 0.51 * parent_r;

      if (r >= start_r && r < parent_r)
          atom->radius[i] = std::min(r + 0.4e-6, parent_r);

      double rho = atom->density[i];
      double rr  = atom->radius[i];
      atom->rmass[i] = (4.0/3.0)*M_PI*rho*rr*rr*rr;
  }
}

/* ---------------------------------------------------------------------- */

int FixMCG::mcg_refine_and_coarsen()
{
  int created = 0;
  created += refine();
  created += coarsen();
  return created;
}

/* ---------------------------------------------------------------------- */

int FixMCG::refine()
{
  int nlocal_old = atom->nlocal;
  std::vector<int> del;
  int created = 0;

  for (int i = 0; i < nlocal_old; i++) {

      if (atom->x[i][1] >= 0.2 &&
          atom->x[i][1] <= 0.25 &&
          fabs(atom->radius[i] - 0.005) < 1e-9)
      {
          refine_particle(i);
          del.push_back(i);
          created += 8;
      }
  }

  if (!del.empty())
      delete_marked_atoms(del);

  return created;
}

/* ---------------------------------------------------------------------- */

void FixMCG::refine_particle(int ith)
{
  const double scale = 0.366;
  const double r_parent = atom->radius[ith];
  const double delta = r_parent * scale;

  int itype = atom->type[ith];

  for (int ix=-1; ix<=1; ix+=2)
  for (int iy=-1; iy<=1; iy+=2)
  for (int iz=-1; iz<=1; iz+=2) {

      double xnew[3] = {
          atom->x[ith][0] + ix*delta,
          atom->x[ith][1] + iy*delta,
          atom->x[ith][2] + iz*delta
      };

      atom->avec->create_atom(itype,xnew);
      int m = atom->nlocal-1;

      for(int k=0;k<3;k++){
          atom->v[m][k] = atom->v[ith][k];
          atom->omega[m][k] = atom->omega[ith][k];
      }

      atom->radius[m] = r_parent*scale;
      atom->density[m] = atom->density[ith];
      atom->rmass[m] = atom->rmass[ith]/8.0;
            
      // ---- CRITICAL FIX INITIALIZATION ----
      for(int j = 0; j < modify->nfix; j++)
      {
          if(modify->fix[j]->create_attribute)
          {
              modify->fix[j]->pre_set_arrays();
              modify->fix[j]->set_arrays(m);
          }
      }
  }
}

/* ---------------------------------------------------------------------- */

int FixMCG::coarsen()
{
  const double child_r = 0.0025;
  const double parent_r = 0.005;
  const double tol = 1e-9;

  int created = 0;

  while(true)
  {
      std::vector<int> cluster;

      for(int i=0;i<atom->nlocal;i++)
      {
          if(fabs(atom->radius[i]-child_r)<tol &&
             atom->x[i][1] < 0.15)
          {
              cluster.push_back(i);
              if(cluster.size()==8) break;
          }
      }

      if(cluster.size()<8) break;

      created += 1;

      double mass=0,vnew[3]={0,0,0},omeganew[3]={0,0,0};

      for(int idx:cluster){
          double m=atom->rmass[idx];
          mass+=m;
          for(int k=0;k<3;k++){
              vnew[k]+=m*atom->v[idx][k];
              omeganew[k]+=m*atom->omega[idx][k];
          }
      }

      for(int k=0;k<3;k++){
          vnew[k]/=mass;
          omeganew[k]/=mass;
      }

      int trigger=cluster[0];
      double xnew[3]={
          atom->x[trigger][0],
          atom->x[trigger][1],
          atom->x[trigger][2]
      };

      atom->avec->create_atom(atom->type[trigger],xnew);
      int m=atom->nlocal-1;
      
      for (int j = 0; j < modify->nfix; j++)
	{
	    if (modify->fix[j]->create_attribute)
	    {
	        modify->fix[j]->pre_set_arrays();
		 modify->fix[j]->set_arrays(m);
	    }
	}

      for(int k=0;k<3;k++){
          atom->v[m][k]=vnew[k];
          atom->omega[m][k]=omeganew[k];
      }

      double start_r=0.52*parent_r;
      atom->radius[m]=start_r;
      atom->density[m]=atom->density[trigger];

      double rho=atom->density[m];
      atom->rmass[m]=(4.0/3.0)*M_PI*rho*start_r*start_r*start_r;

      delete_marked_atoms(cluster);
  }

  return created;
}

/* ---------------------------------------------------------------------- */

void FixMCG::delete_marked_atoms(const std::vector<int>& del)
{
    int nlocal = atom->nlocal;

    std::vector<int> mark(nlocal,0);
    for(int idx : del)
        if(idx < nlocal)
            mark[idx] = 1;

    int i = 0;

    while(i < atom->nlocal)
    {
        if(mark[i])
        {
            int last = atom->nlocal - 1;

            // copy core atom data
            atom->avec->copy(last, i, 1);

            // copy ALL fix-owned arrays
            for(int j = 0; j < modify->nfix; j++)
                modify->fix[j]->copy_arrays(last, i, 1);

            mark[i] = mark[last];

            atom->nlocal--;
        }
        else
        {
            i++;
        }
    }
}



