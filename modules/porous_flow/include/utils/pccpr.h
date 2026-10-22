#pragma once
/*
  PCCPR - Constrained Pressure Residual preconditioner with quasi-IMPES
  decoupling, as a registered PETSc PC type ("cpr").

  Usage from any PETSc code (or MOOSE petsc_options):
    -pc_type cpr
    -pc_cpr_block_size 2          dofs per node (default: matrix block size)
    -pc_cpr_pressure_field 0      which dof within the node is pressure
    -pc_cpr_weights qimpes        qimpes | sum | raw
    -pc_cpr_type multiplicative   multiplicative | additive
    -pc_cpr_pressure_ksp_type / -pc_cpr_pressure_pc_type   (default preonly + hypre/boomeramg)
    -pc_cpr_global_ksp_type   / -pc_cpr_global_pc_type     (default preonly + bjacobi/ilu)

  Without recompiling the application, the PC can be loaded at run time:
    PETSC_OPTIONS="-dll_append /path/to/libpccpr.so" ./app ...
*/
#include <petscpc.h>

PETSC_EXTERN PetscErrorCode PCCPRRegister(void);
PETSC_EXTERN PetscErrorCode PCCPRSetBlockSize(PC, PetscInt);
PETSC_EXTERN PetscErrorCode PCCPRSetPressureField(PC, PetscInt);
/* Optional general layout: dofs[] holds, for each locally owned node, the bs
   global dof indices of that node (node-major).  May be called before the PC
   type is set to "cpr"; the layout is attached to the PC object. */
PETSC_EXTERN PetscErrorCode PCCPRSetNodeDofs(PC, PetscInt bs, PetscInt nnodes, const PetscInt dofs[]);
PETSC_EXTERN PetscErrorCode PCCPRGetPressureKSP(PC, KSP *);
PETSC_EXTERN PetscErrorCode PCCPRGetGlobalKSP(PC, KSP *);
